"""Command-line interface for GTDB downloader"""

import argparse
import concurrent.futures
import json
import re
import sys
import tarfile
import time
from pathlib import Path
from typing import Optional, List, Tuple, Dict

import requests

from gtdb_downloader.config import get_base_dir, get_marker_genes_url, GTDB_VERSIONS
from gtdb_downloader.downloader import (
    check_aria2c_available,
    download_file,
    download_files_aria2,
    download_metadata,
    generate_download_links,
    resolve_download_link,
)
from gtdb_downloader.metadata import MetadataParser, GENOME_TYPE_CATEGORIES


MGNIFY_API_BASE = "https://www.ebi.ac.uk/metagenomics/api/v1"

# Mapping files share a common column layout. ``is_representative`` is "1" when the
# genome is its species-cluster representative, "0" otherwise; ``genome_length`` is
# always the full genome size from metadata, even in the marker-gene mapping.
MAPPING_HEADER = "accession\tgenome_path\tgtdb_taxonomy\tis_representative\tgenome_length\n"
MG_MAPPING_HEADER = (
    "accession\tconcatenated_fna_path\tgtdb_taxonomy\tis_representative\tgenome_length\n"
)
MG_DIR_MAPPING_HEADER = (
    "accession\tmarker_genes_dir\tgtdb_taxonomy\tis_representative\tgenome_length\n"
)


def _sanitize_name(name: str) -> str:
    """Sanitize folder name: replace spaces with underscores and strip unwanted characters"""
    if not name:
        return name
    # Replace spaces with underscores
    out = name.replace(" ", "_")
    # Remove leading/trailing slashes
    out = out.strip("/\\")
    return out


def _get_symlink_name(genome_filename: str, is_species_rep: bool, flag_rep: bool) -> str:
    """Return symlink filename, optionally tagging species representatives."""
    if not (flag_rep and is_species_rep):
        return genome_filename
    if genome_filename.endswith(".fna.gz"):
        return genome_filename[:-7] + ".speciesrep.fna.gz"
    return genome_filename + ".speciesrep.fna.gz"


_ACCESSION_FILE_SUFFIXES = (".txt", ".tsv", ".csv", ".list", ".lst")


def _looks_like_accession_file(value: str) -> bool:
    """Heuristic: does this --accessions value name a file rather than an accession?"""
    if "/" in value or "\\" in value or value.startswith("~"):
        return True
    return Path(value).suffix.lower() in _ACCESSION_FILE_SUFFIXES


def _read_accession_file(path: Path) -> List[str]:
    """Read accessions from a text/TSV file (first column, '#' comments ignored)."""
    accessions: List[str] = []
    with open(path, "r", encoding="utf-8") as handle:
        for line in handle:
            content = line.split("#", 1)[0].strip()
            if not content:
                continue
            accessions.append(content.split()[0])
    return accessions


def _load_accession_queries(values: List[str]) -> List[str]:
    """
    Expand --accessions values into a de-duplicated list of accessions.

    Each value is either a path to a file listing accessions (one per line) or
    a literal accession, optionally comma-separated.

    Raises:
        FileNotFoundError: If a value looks like a file path but does not exist
    """
    queries: List[str] = []
    seen: set = set()

    for value in values:
        candidate = Path(value).expanduser()
        if candidate.is_file():
            raw_items = _read_accession_file(candidate)
        elif _looks_like_accession_file(value):
            raise FileNotFoundError(f"Accession list file not found: {value}")
        else:
            raw_items = [value]

        for item in raw_items:
            for accession in item.split(","):
                accession = accession.strip()
                if accession and accession not in seen:
                    seen.add(accession)
                    queries.append(accession)

    return queries


def _get_default_mapping_path(base_dir: Path, version: str) -> Path:
    """Return the default accession-to-path mapping file location."""
    return base_dir / version / "accession_path_map.tsv"


def _get_shared_raw_genomes_dir(base_dir: Path) -> Path:
    """Return the shared raw genomes directory used by all GTDB versions."""
    return base_dir / "raw"


def _iter_legacy_raw_genomes_dirs(base_dir: Path) -> List[Path]:
    """Return legacy per-version raw genome directories that still exist."""
    legacy_dirs: List[Path] = []
    for child in sorted(base_dir.iterdir()) if base_dir.exists() else []:
        if not child.is_dir():
            continue
        legacy_raw = child / "genomes" / "raw"
        if legacy_raw.is_dir():
            legacy_dirs.append(legacy_raw)
    return legacy_dirs


def _populate_shared_raw_from_legacy(base_dir: Path, shared_raw_dir: Path, verbose: bool = False) -> int:
    """
    Backfill shared raw storage from legacy per-version raw directories.

    Creates symlinks for missing files so existing downloads are reused without
    duplicating data on disk.
    """
    linked_count = 0
    for legacy_raw in _iter_legacy_raw_genomes_dirs(base_dir):
        for legacy_genome in sorted(legacy_raw.glob("*.fna.gz")):
            shared_path = shared_raw_dir / legacy_genome.name
            if shared_path.exists():
                continue
            try:
                shared_path.symlink_to(legacy_genome.resolve())
                linked_count += 1
            except Exception as exc:
                if verbose:
                    print(
                        f"Warning: Could not link legacy genome {legacy_genome} -> {shared_path}: {exc}",
                        file=sys.stderr,
                    )
    return linked_count


def _resolve_mapping_path(base_dir: Path, version: str, mapping_file: Optional[Path]) -> Path:
    """Resolve mapping file path, using CWD for relative custom paths."""
    default_path = _get_default_mapping_path(base_dir, version)
    if mapping_file is None:
        return default_path
    if mapping_file.is_absolute():
        return mapping_file
    return Path.cwd() / mapping_file


def _resolve_failed_path(base_dir: Path, version: str, failed_file: Optional[Path]) -> Optional[Path]:
    """Resolve failed-genome output path, using version directory for relative paths."""
    if failed_file is None:
        return None
    default_dir = base_dir / version
    if failed_file.is_absolute():
        return failed_file
    return default_dir / failed_file


def _normalize_ncbi_accession(accession: str) -> Optional[str]:
    """Normalize accession to NCBI GCA_/GCF_ form."""
    if accession.startswith(("RS_", "GB_")):
        accession = accession[3:]
    if accession.startswith(("GCA_", "GCF_")):
        return accession
    return None


def _extract_ncbi_accession(genome_id: str, genome_metadata: Optional[dict]) -> Optional[str]:
    """Extract a normalized NCBI accession from metadata or genome id."""
    if genome_metadata:
        metadata_accession = genome_metadata.get("accession")
        if isinstance(metadata_accession, str):
            normalized = _normalize_ncbi_accession(metadata_accession)
            if normalized:
                return normalized
    return _normalize_ncbi_accession(genome_id)


def _fetch_ncbi_datasets_status(accession: str, timeout_seconds: int = 12) -> Tuple[str, Optional[str]]:
    """
    Fetch NCBI Datasets page and extract the "Status:" line if present.

    Returns:
        Tuple of (datasets_url, extracted_status_text_or_none)
    """
    datasets_url = f"https://www.ncbi.nlm.nih.gov/datasets/genome/{accession}/"
    query_urls = [
        datasets_url,
        f"{datasets_url}?report=assembly",
        f"https://api.ncbi.nlm.nih.gov/datasets/v2/genome/accession/{accession}/dataset_report",
    ]
    headers = {
        "User-Agent": "gtdb-downloader/0.1 (+https://github.com/)"
    }

    def _extract_status_from_text(text: str) -> Optional[str]:
        compact = re.sub(r"\s+", " ", text)
        status_line = re.search(
            r"Status:\s*(.{1,260}?)\s{1,}(?:This record|Actions|Download|datasets|API|FTP|$)",
            compact,
            re.IGNORECASE,
        )
        if status_line:
            return status_line.group(1).strip()

        # Some pages expose "RefSeq <accession> is suppressed" without a clean Status block.
        suppressed_phrase = re.search(
            rf"((?:RefSeq|GenBank)?\s*{re.escape(accession)}\s+is\s+suppressed)",
            compact,
            re.IGNORECASE,
        )
        if suppressed_phrase:
            return f"Status: {suppressed_phrase.group(1)}"

        if accession.lower() in compact.lower() and "suppressed" in compact.lower():
            return f"Status: {accession} appears suppressed"
        return None

    for url in query_urls:
        try:
            response = requests.get(url, timeout=timeout_seconds, headers=headers)
        except Exception:
            continue

        body = response.text or ""
        direct_status = _extract_status_from_text(body)
        if direct_status:
            return datasets_url, direct_status

        # Strip tags and re-check plain text for server-side rendered fragments.
        plain_text = re.sub(r"<[^>]+>", " ", body)
        plain_text = re.sub(r"\s+", " ", plain_text)
        plain_status = _extract_status_from_text(plain_text)
        if plain_status:
            return datasets_url, plain_status

        # Parse JSON payloads when present.
        try:
            payload = response.json()
        except Exception:
            payload = None
        if payload is not None:
            payload_text = json.dumps(payload, ensure_ascii=False)
            json_status = _extract_status_from_text(payload_text)
            if json_status:
                return datasets_url, json_status

            if "suppressed" in payload_text.lower() and accession.lower() in payload_text.lower():
                return datasets_url, f"Status: {accession} appears suppressed"

    return datasets_url, None


def _fetch_ncbi_status_batch(
    accessions: List[str],
    *,
    timeout_seconds: int = 4,
    max_workers: int = 24,
) -> Dict[str, Optional[str]]:
    """Fetch NCBI Datasets status lines for many accessions in parallel."""
    if not accessions:
        return {}

    unique_accessions = sorted(set(accessions))
    workers = max(1, min(max_workers, len(unique_accessions)))
    results: Dict[str, Optional[str]] = {}
    started = time.monotonic()
    print(
        f"Checking NCBI status for {len(unique_accessions)} accessions "
        f"(workers={workers}, timeout={timeout_seconds}s)...",
        flush=True,
    )

    with concurrent.futures.ThreadPoolExecutor(max_workers=workers) as executor:
        future_to_accession = {
            executor.submit(_fetch_ncbi_datasets_status, accession, timeout_seconds): accession
            for accession in unique_accessions
        }
        done = 0
        for future in concurrent.futures.as_completed(future_to_accession):
            accession = future_to_accession[future]
            try:
                _, status_text = future.result()
            except Exception:
                status_text = None
            results[accession] = status_text
            done += 1
            if done == len(unique_accessions) or done % 25 == 0:
                print(f"  NCBI status progress: {done}/{len(unique_accessions)}", flush=True)

    elapsed = time.monotonic() - started
    print(f"NCBI status checks completed in {elapsed:.1f}s", flush=True)
    return results


def _collect_existing_genome_mappings(genomes_dir: Path) -> Dict[str, Path]:
    """Scan the raw genome directory and build accession-to-path mappings."""
    mappings: Dict[str, Path] = {}
    pattern = re.compile(r"^(GC[AF]_\d+\.\d+)")

    if not genomes_dir.exists():
        return mappings

    for genome_path in sorted(genomes_dir.glob("*.fna.gz")):
        match = pattern.match(genome_path.name)
        if not match:
            continue
        mappings[match.group(1)] = genome_path

    return mappings


def _normalize_mapping_accession(accession: str) -> Optional[str]:
    """Normalize metadata accession to mapping accession key."""
    if accession.startswith(("RS_", "GB_")):
        accession = accession[3:]
    if accession.startswith(("GCA_", "GCF_")):
        return accession
    return None


def _collect_target_accessions_for_mapping(
    *,
    version: str,
    datasets: List[str],
    taxon: Optional[str],
    only_rep: bool,
    base_dir: Path,
    mirror: str,
    verbose: bool,
    accessions: Optional[List[str]] = None,
    ignore_prefix: bool = False,
) -> Optional[set]:
    """
    Collect normalized accession keys for custom mapping subsets.

    Returns None when no filtering is needed (full map).
    """
    if taxon is None and not only_rep and not accessions:
        return None

    target_accessions: set = set()
    version_dir = setup_version_dir(version, base_dir)

    for dataset in datasets:
        metadata_file = download_metadata(version, dataset, version_dir, mirror=mirror, verbose=verbose)
        if metadata_file is None:
            continue
        parser = MetadataParser(metadata_file)
        if taxon is None and not accessions:
            genome_ids = list(parser.data.keys())
        else:
            genome_ids = parser.get_genomes_by_taxon(taxon) if taxon else []
            if accessions:
                matched_ids, _ = parser.get_genomes_by_accessions(
                    accessions, ignore_prefix=ignore_prefix
                )
                genome_ids = list(dict.fromkeys([*genome_ids, *matched_ids]))
        if only_rep:
            genome_ids = [gid for gid in genome_ids if parser.is_species_cluster_representative(gid)]
        for genome_id in genome_ids:
            genome_metadata = parser.get_genome_metadata(genome_id)
            if not genome_metadata:
                continue
            accession = genome_metadata.get("accession")
            if not isinstance(accession, str):
                continue
            normalized = _normalize_mapping_accession(accession)
            if normalized:
                target_accessions.add(normalized)

    return target_accessions


def _collect_taxonomy_lookup_for_mapping(
    *,
    version: str,
    datasets: List[str],
    base_dir: Path,
    mirror: str,
    verbose: bool,
    ensure_metadata: bool = False,
) -> Tuple[Dict[str, Tuple[str, str]], Dict[str, str]]:
    """Build accession -> (taxonomy, is_representative) and accession -> genome_length
    mappings from GTDB metadata.

    ``is_representative`` is "1" if the genome is its species-cluster
    representative, "0" otherwise.
    """
    lookup: Dict[str, Tuple[str, str]] = {}
    length_lookup: Dict[str, str] = {}
    version_dir = setup_version_dir(version, base_dir)

    for dataset in datasets:
        metadata_pattern = f"{dataset}_metadata_{version}.tsv.gz"
        metadata_file = version_dir / metadata_pattern
        if not metadata_file.exists():
            if not ensure_metadata:
                continue
            downloaded = download_metadata(version, dataset, version_dir, mirror=mirror, verbose=verbose)
            if downloaded is None:
                continue
            metadata_file = downloaded

        parser = MetadataParser(metadata_file)
        for genome_metadata in parser.data.values():
            accession = genome_metadata.get("accession")
            taxonomy = genome_metadata.get("gtdb_taxonomy")
            if not isinstance(accession, str) or not isinstance(taxonomy, str):
                continue
            normalized = _normalize_mapping_accession(accession)
            if not normalized:
                continue
            is_rep = "1" if parser.is_species_cluster_representative(accession) else "0"
            # Keep first seen values for a stable map.
            lookup.setdefault(normalized, (taxonomy, is_rep))
            length_lookup.setdefault(normalized, str(genome_metadata.get("genome_size", "")))

    return lookup, length_lookup


def build_mapping_file(
    version: str,
    base_dir: Optional[Path] = None,
    mapping_file: Optional[Path] = None,
    include_accessions: Optional[set] = None,
    taxonomy_lookup: Optional[Dict[str, Tuple[str, str]]] = None,
    genome_length_lookup: Optional[Dict[str, str]] = None,
    show_progress: bool = False,
) -> Path:
    """Create or refresh the accession-to-path mapping file from existing genomes."""
    if base_dir is None:
        base_dir = get_base_dir()

    genomes_dir = _get_shared_raw_genomes_dir(base_dir)
    genomes_dir.mkdir(parents=True, exist_ok=True)
    _populate_shared_raw_from_legacy(base_dir, genomes_dir, verbose=show_progress)
    resolved_mapping_path = _resolve_mapping_path(base_dir, version, mapping_file)
    mappings = _collect_existing_genome_mappings(genomes_dir)
    if include_accessions is not None:
        mappings = {acc: p for acc, p in mappings.items() if acc in include_accessions}
    resolved_mapping_path.parent.mkdir(parents=True, exist_ok=True)
    tmp_path = resolved_mapping_path.with_suffix(resolved_mapping_path.suffix + ".tmp")

    if show_progress:
        print(f"Building mapping file from: {genomes_dir}")

    with open(tmp_path, "w", encoding="utf-8") as handle:
        handle.write(MAPPING_HEADER)
        for count, (accession, genome_path) in enumerate(sorted(mappings.items()), start=1):
            taxonomy = ""
            is_representative = "0"
            if taxonomy_lookup is not None:
                taxonomy, is_representative = taxonomy_lookup.get(accession, ("", "0"))
            genome_length = ""
            if genome_length_lookup is not None:
                genome_length = genome_length_lookup.get(accession, "")
            handle.write(
                f"{accession}\t{genome_path}\t{taxonomy}\t{is_representative}\t{genome_length}\n"
            )
            if show_progress and count % 1000 == 0:
                print(f"  Mapped {count} genomes...")

    tmp_path.replace(resolved_mapping_path)

    if show_progress:
        print(f"Finished mapping {len(mappings)} genomes")
    return resolved_mapping_path


def _chunked(items: List[Dict[str, object]], chunk_size: int) -> List[List[Dict[str, object]]]:
    """Split a list into fixed-size chunks."""
    return [items[i:i + chunk_size] for i in range(0, len(items), chunk_size)]


def _render_progress(
    current: int,
    total: int,
    *,
    prefix: str,
    width: int = 28,
    done: bool = False,
) -> None:
    """Render progress to stdout (TTY bar, line-based fallback for logs)."""
    if total <= 0:
        return

    pct = int((current / total) * 100)
    if sys.stdout.isatty():
        filled = min(width, int((current / total) * width))
        bar = "#" * filled + "-" * (width - filled)
        end = "\n" if done else "\r"
        print(
            f"{prefix} [{bar}] {pct:3d}% ({current}/{total})",
            end=end,
            flush=True,
        )
        return

    # In non-interactive logs (e.g., batch jobs), emit a line every 5%.
    if done or current == total or current == 1 or current % max(1, total // 20) == 0:
        print(f"{prefix} {pct:3d}% ({current}/{total})", flush=True)


def _get_accession_keys(accession: str, ignore_prefix: bool) -> List[str]:
    """Return accession keys to use for local lookup."""
    if ignore_prefix and accession.startswith(("GCA_", "GCF_")):
        return [accession[4:]]
    return [accession]


def _index_existing_genomes(genomes_dir: Path, ignore_prefix: bool) -> Dict[str, Path]:
    """Index existing raw genomes by accession for fast presence checks."""
    indexed: Dict[str, Path] = {}
    for accession, genome_path in _collect_existing_genome_mappings(genomes_dir).items():
        for key in _get_accession_keys(accession, ignore_prefix):
            indexed.setdefault(key, genome_path)
    return indexed


def _find_existing_genome_path(
    existing_genomes: Dict[str, Path],
    genome_metadata: dict,
    ignore_prefix: bool,
) -> Optional[Path]:
    """Return an existing local genome path, optionally ignoring GCA/GCF prefix differences."""
    accession = genome_metadata.get("accession")
    if not accession:
        return None

    normalized = accession[3:] if accession.startswith(("RS_", "GB_")) else accession
    for key in _get_accession_keys(normalized, ignore_prefix):
        path = existing_genomes.get(key)
        if path is not None:
            return path
    return None


def setup_version_dir(version: str, base_dir: Path) -> Path:
    """Setup directory for a specific GTDB version"""
    version_dir = base_dir / version
    version_dir.mkdir(parents=True, exist_ok=True)
    return version_dir


def download_genomes(
    version: str,
    taxon: Optional[str] = None,
    accessions: Optional[List[str]] = None,
    dataset: str = "bac120",
    mirror: str = "europe",
    base_dir: Optional[Path] = None,
    output_dir: Optional[Path] = None,
    flat: Optional[str] = None,
    flag_rep: bool = False,
    only_rep: bool = False,
    ignore_prefix: bool = False,
    failed_file: Optional[Path] = None,
    verbose: bool = False,
    dry_run: bool = False,
    resolved_accessions: Optional[set] = None,
    genome_type: Optional[str] = None,
) -> bool:
    """
    Download genomes selected by taxon and/or by explicit accession

    Args:
        version: GTDB version
        taxon: Taxon to search for (optional if accessions are given)
        accessions: Explicit list of assembly accessions to download
        dataset: Dataset type (bac120 or ar53)
        mirror: Mirror to use for download
        base_dir: Base directory for GTDB data
        output_dir: Output directory for symlink taxonomy structure (not genomes)
        flat: Rank at which to build a flat symlink structure
        flag_rep: Tag species representatives in symlink names
        only_rep: Restrict the selection to species representatives
        ignore_prefix: Treat GCA_/GCF_ prefixes as interchangeable
        failed_file: Where to write the failed-genome report
        verbose: Verbose output
        dry_run: Don't actually download, just show what would be downloaded
        resolved_accessions: Optional set updated with the requested accessions
            that were found in this dataset (used to report unknown accessions
            once all datasets have been searched)
        genome_type: Restrict to a genome source category ("isolate", "mag",
            "sag", "env"); "all" or None keeps every genome

    Returns:
        True if successful, False otherwise
    """
    if not taxon and not accessions:
        print("Error: Either a taxon or a list of accessions is required", file=sys.stderr)
        return False

    if base_dir is None:
        base_dir = get_base_dir()
    
    version_dir = setup_version_dir(version, base_dir)
    
    # Genomes are shared across versions at base_dir / raw
    genomes_dir = _get_shared_raw_genomes_dir(base_dir)
    genomes_dir.mkdir(parents=True, exist_ok=True)
    linked_count = _populate_shared_raw_from_legacy(base_dir, genomes_dir, verbose=verbose)
    if verbose and linked_count:
        print(f"Linked {linked_count} legacy genomes into shared raw directory")
    
    # Symlink taxonomy structure goes to output_dir (or base_dir if not specified)
    if output_dir is None:
        taxonomy_dir = base_dir / version / "genomes" / "taxonomy"
    else:
        taxonomy_dir = output_dir
    
    taxonomy_dir.mkdir(parents=True, exist_ok=True)
    
    # Download metadata if needed
    print("Preparing metadata...")
    metadata_file = download_metadata(version, dataset, version_dir, mirror=mirror, verbose=verbose)
    if metadata_file is None:
        return False
    print(f"Metadata ready: {metadata_file}")
    
    # Parse metadata
    print("Loading metadata into memory...")
    try:
        parser = MetadataParser(metadata_file)
    except Exception as e:
        print(f"Error parsing metadata: {e}", file=sys.stderr)
        return False
    print("Metadata loaded")
    
    # Find matching genomes
    matching_genomes: List[str] = []

    if taxon:
        print(f"Filtering genomes for taxon query: {taxon}")
        matching_genomes.extend(parser.get_genomes_by_taxon(taxon))
        if not matching_genomes:
            message = f"No genomes found for taxon: {taxon} (dataset: {dataset})"
            if not accessions:
                print(message, file=sys.stderr)
                return False
            print(message)

    if accessions:
        print(f"Matching {len(accessions)} requested accessions against {dataset} metadata...")
        accession_genomes, matched_queries = parser.get_genomes_by_accessions(
            accessions, ignore_prefix=ignore_prefix
        )
        if resolved_accessions is not None:
            resolved_accessions.update(matched_queries)
        print(f"Matched {len(accession_genomes)} genomes by accession")
        matching_genomes.extend(accession_genomes)

    # De-duplicate while keeping selection order (taxon hits may overlap accessions).
    matching_genomes = list(dict.fromkeys(matching_genomes))

    label_parts = []
    if taxon:
        label_parts.append(f"taxon: {taxon}")
    if accessions:
        label_parts.append("requested accessions")
    selection_label = " + ".join(label_parts)

    if not matching_genomes:
        print(f"No genomes in dataset {dataset} matched the requested accessions")
        # Accession-only queries commonly match just one dataset; that is not a failure.
        return True

    if only_rep:
        matching_genomes = [
            genome_id
            for genome_id in matching_genomes
            if parser.is_species_cluster_representative(genome_id)
        ]
        if not matching_genomes:
            print(
                f"No species representative genomes found for {selection_label}",
                file=sys.stderr,
            )
            return False

    if genome_type and genome_type != "all":
        matching_genomes = parser.filter_by_genome_type(matching_genomes, genome_type)
        if not matching_genomes:
            print(
                f"No {genome_type} genomes found for {selection_label}",
                file=sys.stderr,
            )
            return False

    print(f"Found {len(matching_genomes)} genomes for {selection_label}")


    if verbose:
        print(f"Dataset: {dataset}")
        print(f"Version: {version}")
        print(f"Genomes directory: {genomes_dir}")
        print(f"Symlink directory: {taxonomy_dir}")
        print("\nGenomes to download:")
        for gid in matching_genomes[:10]:
            print(f"  - {gid}")
        if len(matching_genomes) > 10:
            print(f"  ... and {len(matching_genomes) - 10} more")
    
    if dry_run:
        print("\n[DRY RUN] Download would proceed for the above genomes")
        return True
    
    downloadable: List[Dict[str, object]] = []
    representative_by_cluster: Dict[str, List[str]] = {}
    failed_count = 0
    failed_genomes: List[str] = []
    failed_attempted_urls: Dict[str, List[str]] = {}
    failed_status_notes: Dict[str, str] = {}
    suppressed_genomes: List[str] = []
    status_cache: Dict[str, Optional[str]] = {}
    ncbi_checked_genomes = 0
    ncbi_missing_status_genomes = 0
    existing_genomes = _index_existing_genomes(genomes_dir, ignore_prefix)
    total_genomes = len(matching_genomes)
    show_prep_progress = not verbose and total_genomes > 200
    last_progress_update = 0.0

    if show_prep_progress:
        print(f"\nPreparing {total_genomes} genome entries before download...")
        _render_progress(0, total_genomes, prefix="Preparing", done=False)

    for i, genome_id in enumerate(matching_genomes, 1):
        if verbose:
            print(f"\n[{i}/{len(matching_genomes)}] Processing {genome_id}...")
        elif show_prep_progress:
            now = time.monotonic()
            if i == total_genomes or now - last_progress_update >= 0.2:
                _render_progress(i, total_genomes, prefix="Preparing", done=(i == total_genomes))
                last_progress_update = now

        genome_metadata = parser.get_genome_metadata(genome_id)
        if genome_metadata is None:
            if verbose:
                print("  Skipped: Could not find genome metadata")
            failed_count += 1
            failed_genomes.append(genome_id)
            failed_attempted_urls.setdefault(genome_id, [])
            continue

        tax_info = parser.get_taxon_path(genome_id)
        if tax_info is None:
            if verbose:
                print("  Skipped: Could not find taxonomy info")
            failed_count += 1
            failed_genomes.append(genome_id)
            failed_attempted_urls.setdefault(genome_id, [])
            continue

        taxonomy_str, _ = tax_info

        try:
            download_urls = generate_download_links(genome_metadata, ignore_prefix=ignore_prefix)
            download_url = download_urls[0]
            genome_filename = download_url.split("/")[-1]
        except Exception as e:
            if verbose:
                print(f"  Skipped: Could not generate download link: {e}")
            failed_count += 1
            failed_genomes.append(genome_id)
            failed_attempted_urls.setdefault(genome_id, [])
            continue

        genome_path = genomes_dir / genome_filename
        existing_genome_path = _find_existing_genome_path(existing_genomes, genome_metadata, ignore_prefix)
        local_present = False
        if existing_genome_path is not None:
            genome_path = existing_genome_path
            local_present = True
        elif genome_path.exists():
            # Fallback check for already-downloaded files not captured by accession index.
            local_present = True
        is_species_rep = parser.is_species_cluster_representative(genome_id)
        cluster_rep = parser.get_species_cluster_representative(genome_id) or "unknown_cluster"
        downloadable.append({
            "genome_id": genome_id,
            "taxonomy_str": taxonomy_str,
            "download_url": download_url,
            "genome_path": genome_path,
            "local_present": local_present,
            "is_species_rep": is_species_rep,
            "cluster_rep": cluster_rep,
            "genome_metadata": genome_metadata,
            "download_urls": download_urls,
            "attempted_urls": [],
        })

        if flag_rep and is_species_rep:
            representative_by_cluster.setdefault(cluster_rep, []).append(genome_id)

        if verbose and not local_present:
            print(f"  Queued download: {download_url}")

    if flag_rep:
        for cluster_rep, rep_genomes in representative_by_cluster.items():
            if len(rep_genomes) > 1:
                print(
                    (
                        f"Warning: More than one genome qualifies as cluster representative "
                        f"for {cluster_rep}: {', '.join(rep_genomes)}"
                    ),
                    file=sys.stderr,
                )

    batch_results: Dict[Path, bool] = {}
    pending_downloads = [
        item for item in downloadable if not item["local_present"]  # type: ignore[index]
    ]

    if pending_downloads:
        use_aria2 = check_aria2c_available()
        downloader_name = "aria2c" if use_aria2 else "wget"
        print(
            f"\nStarting download of {len(pending_downloads)} genomes with chunked fallback retries using {downloader_name}..."
        )
        chunk_size = 250 if use_aria2 else 25
        total_chunks = (len(pending_downloads) + chunk_size - 1) // chunk_size

        for chunk_index, chunk in enumerate(_chunked(pending_downloads, chunk_size), start=1):
            print(
                f"\nPrimary download chunk {chunk_index}/{total_chunks} "
                f"({len(chunk)} genomes)..."
            )

            primary_downloads: List[Tuple[str, Path]] = []
            for item in chunk:
                if item["local_present"]:  # type: ignore[index]
                    continue
                url = item["download_urls"][0]  # type: ignore[index]
                item["attempted_urls"].append(url)  # type: ignore[index]
                primary_downloads.append((url, item["genome_path"]))  # type: ignore[index]

            if use_aria2:
                tmp_dir = genomes_dir / ".tmp"
                chunk_results = download_files_aria2(
                    primary_downloads,
                    verbose=verbose,
                    tmp_dir=tmp_dir,
                )
            else:
                chunk_results = {}
                for url, genome_path in primary_downloads:
                    chunk_results[genome_path] = download_file(
                        url,
                        genome_path,
                        verbose=verbose,
                        use_aria2=False,
                    )

            batch_results.update(chunk_results)

            failed_chunk_items = [
                item
                for item in chunk
                if not item["local_present"] and not batch_results.get(item["genome_path"], False)  # type: ignore[index]
            ]

            if not failed_chunk_items:
                continue

            print(
                f"Resolving fallbacks for {len(failed_chunk_items)} failed genomes in chunk {chunk_index}/{total_chunks}..."
            , flush=True)

            # Run NCBI status checks in parallel before fallback resolution.
            accessions_to_query: List[str] = []
            for item in failed_chunk_items:
                accession = _extract_ncbi_accession(  # type: ignore[arg-type]
                    item["genome_id"], item["genome_metadata"]
                )
                if accession and accession not in status_cache:
                    accessions_to_query.append(accession)
            if accessions_to_query:
                status_cache.update(_fetch_ncbi_status_batch(accessions_to_query))

            fallback_downloads: List[Tuple[str, Path]] = []
            for idx, item in enumerate(failed_chunk_items, start=1):
                genome_id = item["genome_id"]
                genome_metadata = item["genome_metadata"]
                if verbose:
                    print(f"  Resolving fallback for {genome_id}...")
                elif idx == len(failed_chunk_items) or idx % 50 == 0:
                    print(f"  Fallback decision progress: {idx}/{len(failed_chunk_items)}", flush=True)

                accession = _extract_ncbi_accession(genome_id, genome_metadata)  # type: ignore[arg-type]
                if accession:
                    ncbi_checked_genomes += 1

                # The NCBI directory listing is authoritative: if a genome file
                # exists on the FTP we download it regardless of the status
                # heuristic, which has false positives (it reports plenty of
                # *current* genomes as "suppressed"). Only when resolution
                # genuinely finds no file do we record the NCBI status as the
                # failure reason.
                try:
                    fallback_url = resolve_download_link(
                        genome_metadata,
                        verbose=verbose,
                        ignore_prefix=ignore_prefix,
                    )  # type: ignore[arg-type]
                    fallback_filename = fallback_url.split("/")[-1]
                    fallback_path = genomes_dir / fallback_filename
                    item["download_url"] = fallback_url
                    item["genome_path"] = fallback_path
                    fallback_accession = fallback_filename.split("_genomic.fna.gz", 1)[0].split("_", 2)
                    if len(fallback_accession) >= 2:
                        normalized_accession = "_".join(fallback_accession[:2])
                        for key in _get_accession_keys(normalized_accession, ignore_prefix):
                            existing_genomes[key] = fallback_path
                    print(f"  Fallback resolved for {genome_id}: {fallback_filename}")
                except Exception as e:
                    # Unresolvable -- now (and only now) consult the NCBI status
                    # to explain why, and flag genuinely-suppressed genomes.
                    if accession:
                        status_text = status_cache.get(accession)
                        if status_text:
                            failed_status_notes[genome_id] = f"status:{status_text}"
                            if "suppressed" in status_text.lower():
                                datasets_url = f"https://www.ncbi.nlm.nih.gov/datasets/genome/{accession}/"
                                print(
                                    f"Status: {status_text} ({datasets_url}); no downloadable file found",
                                    file=sys.stderr,
                                )
                                suppressed_genomes.append(genome_id)  # type: ignore[arg-type]
                        else:
                            failed_status_notes[genome_id] = "status:unavailable_or_not_found"
                            ncbi_missing_status_genomes += 1
                    else:
                        failed_status_notes[genome_id] = "status:no_accession"
                    if verbose:
                        print(f"  Fallback resolution failed for {genome_id}: {e}")
                    continue

                if fallback_path.exists():
                    batch_results[fallback_path] = True
                    item["local_present"] = True
                    continue

                item["attempted_urls"].append(fallback_url)  # type: ignore[index]
                fallback_downloads.append((fallback_url, fallback_path))

            if not fallback_downloads:
                continue

            if use_aria2:
                tmp_dir = genomes_dir / ".tmp"
                fallback_results = download_files_aria2(
                    fallback_downloads,
                    verbose=verbose,
                    tmp_dir=tmp_dir,
                )
            else:
                fallback_results = {}
                for url, genome_path in fallback_downloads:
                    fallback_results[genome_path] = download_file(
                        url,
                        genome_path,
                        verbose=verbose,
                        use_aria2=False,
                    )

            batch_results.update(fallback_results)

    print("\nFinalizing symlinks/results...", flush=True)
    downloaded_count = 0
    show_finalize_progress = not verbose and len(downloadable) > 200
    finalize_total = len(downloadable)
    finalize_last_update = 0.0
    if show_finalize_progress:
        _render_progress(0, finalize_total, prefix="Finalizing", done=False)

    for idx, item in enumerate(downloadable, start=1):
        if show_finalize_progress:
            now = time.monotonic()
            if idx == finalize_total or now - finalize_last_update >= 0.2:
                _render_progress(idx, finalize_total, prefix="Finalizing", done=(idx == finalize_total))
                finalize_last_update = now

        genome_id = item["genome_id"]
        taxonomy_str = item["taxonomy_str"]
        genome_path = item["genome_path"]
        local_present = bool(item.get("local_present", False))
        is_species_rep = item["is_species_rep"]

        if not local_present and not genome_path.exists() and not batch_results.get(genome_path, False):
            if verbose:
                print(f"  Failed to download {genome_id}")
            failed_count += 1
            failed_genomes.append(genome_id)
            failed_attempted_urls[genome_id] = list(item.get("attempted_urls", []))  # type: ignore[arg-type]
            continue

        # Create symlink in taxonomy structure
        try:
            if flat:
                # Create a flat structure at the requested rank (e.g., species -> s__...)
                comp = parser.get_taxon_component_at_rank(taxonomy_str, flat)
                if comp:
                    # use the component as a single folder name (keep prefix)
                    tax_folder = taxonomy_dir / _sanitize_name(comp)
                else:
                    # fallback to full path
                    tax_parts = parser.parse_taxonomy_to_path(taxonomy_str)
                    tax_parts = [_sanitize_name(p) for p in tax_parts]
                    tax_folder = taxonomy_dir.joinpath(*tax_parts) if tax_parts else taxonomy_dir
            else:
                tax_parts = parser.parse_taxonomy_to_path(taxonomy_str)
                tax_parts = [_sanitize_name(p) for p in tax_parts]
                tax_folder = taxonomy_dir.joinpath(*tax_parts) if tax_parts else taxonomy_dir

            tax_folder.mkdir(parents=True, exist_ok=True)

            link_name = _get_symlink_name(genome_path.name, is_species_rep, flag_rep)
            link_path = tax_folder / link_name
            if not link_path.exists():
                link_path.symlink_to(genome_path.resolve())
                if verbose:
                    print(f"  Symlinked: {link_path}")
        except Exception as e:
            if verbose:
                print(f"  Warning: Could not create symlink: {e}")
        
        downloaded_count += 1

    taxonomy_lookup, genome_length_lookup = _collect_taxonomy_lookup_for_mapping(
        version=version,
        datasets=["bac120", "ar53"],
        base_dir=base_dir,
        mirror=mirror,
        verbose=verbose,
        ensure_metadata=False,
    )
    mapping_path = build_mapping_file(
        version,
        base_dir=base_dir,
        taxonomy_lookup=taxonomy_lookup,
        genome_length_lookup=genome_length_lookup,
        show_progress=verbose,
    )
    print(f"Mapping file written to: {mapping_path}")

    resolved_failed_path = _resolve_failed_path(base_dir, version, failed_file)
    if resolved_failed_path is not None:
        resolved_failed_path.parent.mkdir(parents=True, exist_ok=True)
        with open(resolved_failed_path, "w", encoding="utf-8") as handle:
            for genome_id in sorted(set(failed_genomes)):
                attempted_urls = failed_attempted_urls.get(genome_id, [])
                if attempted_urls:
                    handle.write("\t".join([genome_id, *attempted_urls]) + "\n")
                else:
                    handle.write(f"{genome_id}\n")
        print(f"Failed genome list written to: {resolved_failed_path}")
        status_path = resolved_failed_path.with_name(resolved_failed_path.stem + ".status.tsv")
        with open(status_path, "w", encoding="utf-8") as handle:
            for genome_id in sorted(set(failed_genomes)):
                note = failed_status_notes.get(genome_id, "status:not_checked")
                handle.write(f"{genome_id}\t{note}\n")
        print(f"Failed status list written to: {status_path}")
    
    print(f"\n✓ Downloaded: {downloaded_count}")
    print(f"✗ Failed: {failed_count}")
    if failed_count:
        print(
            f"NCBI status checks: {ncbi_checked_genomes} "
            f"(no status line or unreachable: {ncbi_missing_status_genomes})"
        )
    if suppressed_genomes:
        print(f"Suppressed at NCBI Datasets: {len(set(suppressed_genomes))}")
    print(f"Genomes stored in: {genomes_dir}")
    print(f"Taxonomy structure in: {taxonomy_dir}")
    
    return failed_count == 0


def download_genomes_for_taxon(taxon: str, version: str, **kwargs) -> bool:
    """Backwards-compatible wrapper around :func:`download_genomes`."""
    return download_genomes(version, taxon=taxon, **kwargs)


def _concat_marker_genes_for_dataset(
    extracted_dir: Path,
    verbose: bool = False,
) -> Optional[Path]:
    """Create one FASTA file per genome containing one record per marker gene.

    Processes one marker file at a time — only that marker's sequences are kept
    in memory at once. Each output file gets 120 records (one per marker present),
    in consistent sorted-filename order, suitable for per-marker MSA preparation.
    RS_/GB_ prefixes are stripped from output filenames.
    Returns the output directory, or None if marker genes are not found.
    """
    fna_dir = extracted_dir / "fna"
    if not fna_dir.exists():
        return None

    fna_files = sorted(fna_dir.glob("*.fna"))
    if not fna_files:
        print(f"No .fna files found in {fna_dir}", file=sys.stderr)
        return None

    output_dir = extracted_dir / "concatenated"
    output_dir.mkdir(parents=True, exist_ok=True)

    print(f"Processing {len(fna_files)} marker files from {fna_dir}...")
    genome_count = 0

    for i, fna_file in enumerate(fna_files, 1):
        marker_name = fna_file.stem
        if verbose:
            print(f"  [{i}/{len(fna_files)}] {fna_file.name}")
        elif i == 1 or i == len(fna_files) or i % 10 == 0:
            print(f"  Marker {i}/{len(fna_files)}: {fna_file.name}", flush=True)

        # Read this one marker's sequences into memory (~200 MB peak, then discarded)
        marker_seqs: Dict[str, str] = {}
        current_acc: Optional[str] = None
        current_seq: List[str] = []
        with open(fna_file, "r", encoding="utf-8") as fh:
            for line in fh:
                line = line.rstrip()
                if line.startswith(">"):
                    if current_acc is not None:
                        marker_seqs[current_acc] = "".join(current_seq)
                    raw_acc = line[1:].split()[0]
                    current_acc = raw_acc[3:] if raw_acc.startswith(("RS_", "GB_")) else raw_acc
                    current_seq = []
                elif line:
                    current_seq.append(line)
            if current_acc is not None:
                marker_seqs[current_acc] = "".join(current_seq)

        if i == 1:
            genome_count = len(marker_seqs)
            print(f"  {genome_count} genomes found", flush=True)

        # Write/append each genome's record immediately — "w" on first marker
        # to initialise clean files, "a" on all subsequent markers.
        mode = "w" if i == 1 else "a"
        for acc, seq in marker_seqs.items():
            out_path = output_dir / f"{acc}.fna"
            with open(out_path, mode, encoding="utf-8") as fh:
                fh.write(f">{marker_name}\n{seq}\n")

    print(f"Done: {genome_count} genome files written to {output_dir}")
    return output_dir


def _build_mg_mapping_file(
    version: str,
    datasets: List[str],
    base_dir: Path,
    mapping_file: Path,
    mirror: str,
    verbose: bool,
) -> Optional[Path]:
    """Build a combined marker-gene path mapping file across datasets.

    Writes the same columns as the genome mapping file, but the path column points
    to the concatenated marker gene file for each accession.
    genome_length is the full genome size from metadata, not the marker gene file size.
    """
    version_dir = setup_version_dir(version, base_dir)
    marker_genes_dir = version_dir / "marker_genes"

    rows: List[Tuple[str, Path, str, str, str]] = []

    for dataset in datasets:
        extracted_dir = marker_genes_dir / f"{dataset}_marker_genes_all_{version}"
        concat_dir = extracted_dir / "concatenated"
        fna_dir = extracted_dir / "fna"

        if not fna_dir.exists():
            print(
                f"Marker gene directory not found for {dataset}: {fna_dir}\n"
                f"Run --download-marker-genes first.",
                file=sys.stderr,
            )
            continue

        if not concat_dir.exists():
            print(
                f"Concatenated marker gene files not found for {dataset}: {concat_dir}\n"
                f"Run --build-mg first.",
                file=sys.stderr,
            )
            continue

        raw_accessions = _collect_accessions_from_marker_fna(fna_dir)
        if not raw_accessions:
            print(f"No accessions found in {fna_dir}", file=sys.stderr)
            continue

        metadata_file = download_metadata(version, dataset, version_dir, mirror=mirror, verbose=verbose)
        if metadata_file is None:
            print(f"Warning: Could not get metadata for {dataset}", file=sys.stderr)
            continue

        parser = MetadataParser(metadata_file)
        for raw_acc in raw_accessions:
            genome_metadata = parser.get_genome_metadata(raw_acc)
            if genome_metadata is None:
                if verbose:
                    print(f"  No metadata for {raw_acc}")
                continue
            accession = _normalize_mapping_accession(
                genome_metadata.get("accession", raw_acc)
            ) or raw_acc
            # Use the normalized accession (RS_/GB_ stripped) to match the output filenames
            normalized_raw = raw_acc[3:] if raw_acc.startswith(("RS_", "GB_")) else raw_acc
            genome_file = concat_dir / f"{normalized_raw}.fna"
            if not genome_file.exists():
                if verbose:
                    print(f"  No concatenated file for {normalized_raw}, skipping")
                continue
            taxonomy = genome_metadata.get("gtdb_taxonomy", "")
            genome_size = str(genome_metadata.get("genome_size", ""))
            is_rep = "1" if parser.is_species_cluster_representative(raw_acc) else "0"
            rows.append((accession, genome_file, taxonomy, is_rep, genome_size))

    if not rows:
        print(
            "No concatenated marker gene files found. "
            "Run --build-mg first to create per-genome files.",
            file=sys.stderr,
        )
        return None

    resolved_path = mapping_file if mapping_file.is_absolute() else Path.cwd() / mapping_file
    resolved_path.parent.mkdir(parents=True, exist_ok=True)
    tmp_path = resolved_path.with_suffix(resolved_path.suffix + ".tmp")

    with open(tmp_path, "w", encoding="utf-8") as fh:
        fh.write(MG_MAPPING_HEADER)
        for accession, mg_path, taxonomy, is_rep, genome_size in sorted(rows):
            fh.write(f"{accession}\t{mg_path}\t{taxonomy}\t{is_rep}\t{genome_size}\n")

    tmp_path.replace(resolved_path)
    return resolved_path


def _collect_accessions_from_marker_fna(fna_dir: Path) -> List[str]:
    """
    Read FASTA headers from marker gene .fna files to collect all accessions present.

    The files are organized per-marker (one file per PF family, all genomes inside),
    so we only need to read one file to get the full accession list.
    """
    accessions: List[str] = []
    fna_files = sorted(fna_dir.glob("*.fna"))
    if not fna_files:
        return accessions
    with open(fna_files[0], "r", encoding="utf-8") as fh:
        for line in fh:
            if line.startswith(">"):
                accessions.append(line[1:].strip().split()[0])
    return accessions


def _build_marker_genes_mapping(
    version: str,
    dataset: str,
    extracted_dir: Path,
    version_dir: Path,
    base_dir: Path,
    mirror: str,
    verbose: bool,
) -> Optional[Path]:
    """Build accession→local_path mapping TSV for an extracted marker genes directory.

    Marker genes are stored per-marker (one file per PF family containing all genomes),
    so local_path points to the shared extracted directory for every accession.
    """
    metadata_file = download_metadata(version, dataset, version_dir, mirror=mirror, verbose=verbose)
    if metadata_file is None:
        print(f"Warning: Could not get metadata for {dataset}, skipping mapping file", file=sys.stderr)
        return None

    parser = MetadataParser(metadata_file)

    fna_dir = extracted_dir / "fna"
    if not fna_dir.exists():
        print(f"Warning: fna subdirectory not found in {extracted_dir}", file=sys.stderr)
        return None

    raw_accessions = _collect_accessions_from_marker_fna(fna_dir)
    if not raw_accessions:
        print(f"Warning: No accessions found in {fna_dir}", file=sys.stderr)
        return None

    rows: List[Tuple[str, Path, str, str, str]] = []
    for raw_acc in raw_accessions:
        genome_metadata = parser.get_genome_metadata(raw_acc)
        if genome_metadata is None:
            if verbose:
                print(f"  Warning: No metadata found for {raw_acc}")
            continue

        accession = _normalize_mapping_accession(
            genome_metadata.get("accession", raw_acc)
        ) or raw_acc
        genome_size = str(genome_metadata.get("genome_size", ""))
        taxonomy = genome_metadata.get("gtdb_taxonomy", "")
        is_rep = "1" if parser.is_species_cluster_representative(raw_acc) else "0"
        rows.append((accession, extracted_dir, taxonomy, is_rep, genome_size))

    mapping_path = version_dir / f"{dataset}_marker_genes_map.tsv"
    tmp_path = mapping_path.with_suffix(".tsv.tmp")
    with open(tmp_path, "w", encoding="utf-8") as fh:
        fh.write(MG_DIR_MAPPING_HEADER)
        for accession, marker_path, taxonomy, is_rep, genome_size in sorted(rows):
            fh.write(f"{accession}\t{marker_path}\t{taxonomy}\t{is_rep}\t{genome_size}\n")
    tmp_path.replace(mapping_path)

    print(f"Mapped {len(rows)} genomes to marker gene paths")
    return mapping_path


def _list_mgnify_catalogues() -> List[Dict]:
    """Return available MGnify genome catalogues from the API."""
    try:
        resp = requests.get(
            f"{MGNIFY_API_BASE}/genome-catalogues",
            params={"page_size": 200},
            timeout=30,
        )
        resp.raise_for_status()
        return resp.json().get("data", [])
    except Exception as e:
        print(f"Error fetching MGnify catalogues: {e}", file=sys.stderr)
        return []


def _fetch_mgnify_accessions(catalogue_id: str, verbose: bool = False) -> Optional[set]:
    """Fetch all NCBI GCA/GCF base accessions for a MGnify catalogue.

    Downloads the FTP genomes-all_metadata.tsv for the catalogue (single file,
    much faster than API pagination) and reads the Genome_accession column.
    Returns base accessions without version suffix (e.g. GCA_018713305) since
    MGnify omits version numbers — GTDB matching strips versions on that side.
    Returns None on error, or an empty set if no NCBI accessions are available
    for this catalogue.
    """
    try:
        resp = requests.get(
            f"{MGNIFY_API_BASE}/genome-catalogues/{catalogue_id}",
            timeout=30,
        )
        resp.raise_for_status()
        ftp_url = resp.json().get("data", {}).get("attributes", {}).get("ftp-url", "").rstrip("/")
    except Exception as e:
        print(f"Error fetching catalogue info for {catalogue_id}: {e}", file=sys.stderr)
        return None

    if not ftp_url:
        print(f"No FTP URL for catalogue {catalogue_id}", file=sys.stderr)
        return None

    metadata_url = f"{ftp_url}/genomes-all_metadata.tsv"
    if verbose:
        print(f"  Downloading: {metadata_url}")

    try:
        resp = requests.get(metadata_url, timeout=120, stream=True)
        resp.raise_for_status()
    except Exception as e:
        print(f"Error downloading MGnify metadata TSV: {e}", file=sys.stderr)
        return None

    accessions: set = set()
    genome_acc_col: Optional[int] = None

    for i, line in enumerate(resp.iter_lines(decode_unicode=True)):
        if i == 0:
            headers = line.split("\t")
            try:
                genome_acc_col = headers.index("Genome_accession")
            except ValueError:
                print(
                    f"Column 'Genome_accession' not found in MGnify metadata TSV.\n"
                    f"Available columns: {', '.join(headers[:15])}",
                    file=sys.stderr,
                )
                return None
            continue

        parts = line.split("\t")
        if genome_acc_col is not None and genome_acc_col < len(parts):
            acc = parts[genome_acc_col].strip()
            if acc.startswith(("GCA_", "GCF_")):
                accessions.add(acc.rsplit(".", 1)[0])

    print(f"  Found {len(accessions)} NCBI accessions in MGnify catalogue", flush=True)
    return accessions


def _filter_gtdb_metadata_by_mgnify(
    version: str,
    datasets: List[str],
    base_dir: Path,
    mgnify_accessions: set,
    output_path: Path,
    mirror: str,
    verbose: bool,
) -> int:
    """Intersect MGnify accessions against the full GTDB metadata (all genomes, not just
    downloaded ones). Genome paths are filled in from the local mapping file where available.

    Output columns: accession, genome_path, gtdb_taxonomy, is_representative, genome_length
    genome_path is empty when the genome has not been downloaded locally.
    Returns count of matching rows written.
    """
    version_dir = base_dir / version

    # Load local genome paths from the mapping file (best-effort, may not exist)
    path_lookup: Dict[str, str] = {}
    mapping_path = _get_default_mapping_path(base_dir, version)
    if mapping_path.exists():
        with open(mapping_path, "r", encoding="utf-8") as fh:
            for line in fh:
                parts = line.rstrip("\n").split("\t")
                if parts[0] == "accession" or len(parts) < 2:
                    continue
                path_lookup[parts[0]] = parts[1]

    output_path.parent.mkdir(parents=True, exist_ok=True)
    tmp_path = output_path.with_suffix(output_path.suffix + ".tmp")
    count = 0

    with open(tmp_path, "w", encoding="utf-8") as fout:
        fout.write(MAPPING_HEADER)
        for dataset in datasets:
            metadata_file = download_metadata(version, dataset, version_dir, mirror=mirror, verbose=verbose)
            if metadata_file is None:
                print(f"Warning: metadata not available for {dataset}, skipping", file=sys.stderr)
                continue

            print(f"Scanning {dataset} metadata...", flush=True)
            parser = MetadataParser(metadata_file)
            for genome_metadata in parser.data.values():
                raw_acc = genome_metadata.get("accession", "")
                normalized = _normalize_mapping_accession(raw_acc)
                if not normalized:
                    continue
                if normalized.rsplit(".", 1)[0] not in mgnify_accessions:
                    continue
                taxonomy = genome_metadata.get("gtdb_taxonomy", "")
                genome_size = str(genome_metadata.get("genome_size", ""))
                is_rep = "1" if parser.is_species_cluster_representative(raw_acc) else "0"
                genome_path = path_lookup.get(normalized, "")
                fout.write(
                    f"{normalized}\t{genome_path}\t{taxonomy}\t{is_rep}\t{genome_size}\n"
                )
                count += 1

    tmp_path.replace(output_path)
    return count


def main():
    """Main entry point for CLI"""
    parser = argparse.ArgumentParser(
        description="Download GTDB genomes by taxonomy",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Download metadata for r226
  gtdb-dl --gtdb r226 --download

  # Download all Bacillota genomes
  gtdb-dl --gtdb r226 --taxon "Bacillota"

  # Download with verbose output
  gtdb-dl --gtdb r226 --taxon "d__Bacteria;p__Bacillota" -v

  # Download a specific list of assemblies
  gtdb-dl --gtdb r226 --accessions GCF_000970205.1 GCA_023390935.1

  # Download the assemblies listed in a file (one accession per line)
  gtdb-dl --gtdb r226 --accessions my_accessions.txt

  # Set custom base directory
  GTDBDL_DATA=/data/gtdb gtdb-dl --gtdb r226 --taxon "Archaea"
        """
    )
    
    parser.add_argument(
        "--gtdb",
        required=True,
        choices=list(GTDB_VERSIONS.keys()),
        help="GTDB version to use (e.g., r207, r214, r220, r226)"
    )
    
    parser.add_argument(
        "--taxon",
        help=(
            "Taxon to search for (e.g., 'Bacillota', 'Bacteria,Archaea', "
            "or full GTDB taxonomy path)"
        )
    )
    
    parser.add_argument(
        "--accessions",
        nargs="+",
        metavar="ACCESSION|FILE",
        help=(
            "Download a specific list of assemblies instead of (or in addition to) a taxon. "
            "Each value is either an accession or a file listing one accession per line "
            "(first column used, '#' comments ignored). Accessions may be given as "
            "GCF_000970205.1, GCA_000970205, or RS_GCF_000970205.1."
        )
    )

    parser.add_argument(
        "--dataset",
        choices=["bac120", "ar53", "all"],
        default="all",
        help="Dataset type (default: all)"
    )
    
    parser.add_argument(
        "--mirror",
        choices=["europe", "asia-pacific1", "asia-pacific2"],
        default="europe",
        help="Mirror to download from (default: europe)"
    )
    
    parser.add_argument(
        "--flat",
        choices=["domain","phylum","class","order","family","genus","species","d","p","c","o","f","g","s"],
        help="Create a flat symlink structure at the given rank (e.g. --flat species)")

    parser.add_argument(
        "--flag-rep",
        action="store_true",
        help="Append .speciesrep.fna.gz to symlinks for species-cluster representatives"
    )

    parser.add_argument(
        "--only-rep",
        action="store_true",
        help="Only include genomes marked as species representatives in metadata (gtdb_representative=t)"
    )

    parser.add_argument(
        "--genome-type",
        choices=["all", "isolate", "mag", "sag", "env"],
        default="all",
        help=(
            "Filter by genome source: isolate (cultured), mag (metagenome-assembled), "
            "sag (single-cell), env (environmental sample). Default: all"
        )
    )

    parser.add_argument(
        "--ignore-prefix",
        action="store_true",
        help="Treat GCA_ and GCF_ accession prefixes as interchangeable for local file checks and URL resolution"
    )
    
    parser.add_argument(
        "--output",
        "-o",
        type=Path,
        help="Output directory for symlink taxonomy structure (default: ~/.gtdb_downloader/{version}/genomes/taxonomy)"
    )

    parser.add_argument(
        "--mapping-file",
        nargs="?",
        const=Path("accession_path_map.tsv"),
        type=Path,
        help=(
            "Write or refresh a TSV mapping file (col1 accession, col2 local raw genome path, col3 taxonomy, "
            "col4 is_representative [1/0]). "
            "With no file value, print the global mapping file path for the selected version. "
            "Relative custom paths are resolved from the current working directory."
        )
    )

    parser.add_argument(
        "--mapping-file-mg",
        nargs="?",
        const=Path("accession_mg_path_map.tsv"),
        type=Path,
        help=(
            "Write a TSV mapping file pointing to per-genome concatenated marker gene files "
            "(accession, concatenated_fna_path, gtdb_taxonomy, is_representative, genome_length). "
            "genome_length is the full genome size from metadata, not the marker gene file size. "
            "Requires --download-marker-genes and --build-mg to be run first."
        )
    )

    parser.add_argument(
        "--failed-file",
        nargs="?",
        const=Path("failed_genomes.txt"),
        type=Path,
        help=(
            "Write failed-genome TSV: column 1 genome ID, columns 2+ attempted URLs. "
            "Relative paths are written under ~/.gtdb_downloader/{version}/."
        )
    )
    
    parser.add_argument(
        "--download",
        action="store_true",
        help="Only download metadata, do not download genomes"
    )

    parser.add_argument(
        "--download-marker-genes",
        action="store_true",
        help="Download marker gene tarballs for the selected GTDB version into {base_dir}/{version}/marker_genes/"
    )

    parser.add_argument(
        "--build-mg",
        action="store_true",
        help=(
            "Concatenate per-marker FASTA files into one file per genome, written to "
            "{base_dir}/{version}/marker_genes/{dataset}_marker_genes_all_{version}/concatenated/. "
            "Marker genes must be downloaded first with --download-marker-genes."
        )
    )

    parser.add_argument(
        "--mgnify-catalogues",
        action="store_true",
        help="List available MGnify genome catalogues (biome-specific genome collections)."
    )

    parser.add_argument(
        "--mgnify-filter",
        metavar="CATALOGUE_ID",
        help=(
            "Filter the GTDB accession mapping file to genomes present in a MGnify catalogue "
            "(e.g. human-gut-v2-0, marine-v1-0). "
            "Use --mgnify-catalogues to list available catalogues. "
            "Reads the default mapping file for the selected GTDB version."
        )
    )

    parser.add_argument(
        "--mgnify-output",
        metavar="FILE",
        type=Path,
        help=(
            "Output path for the MGnify-filtered mapping TSV "
            "(default: {catalogue_id}_accession_path_map.tsv in current directory)."
        )
    )

    parser.add_argument(
        "--verbose",
        "-v",
        action="store_true",
        help="Verbose output"
    )
    
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Show what would be downloaded without actually downloading"
    )
    
    parser.add_argument(
        "--base-dir",
        type=Path,
        help="Base directory for GTDB data (can also set GTDBDL_DATA env var)"
    )
    
    args = parser.parse_args()
    argv = sys.argv[1:]
    mapping_flag_without_value = False
    for i, token in enumerate(argv):
        if token == "--mapping-file":
            if i == len(argv) - 1 or argv[i + 1].startswith("-"):
                mapping_flag_without_value = True
            break

    base_dir = args.base_dir or get_base_dir()

    if args.verbose:
        print(f"Using base directory: {base_dir}")

    accession_queries: List[str] = []
    if args.accessions:
        try:
            accession_queries = _load_accession_queries(args.accessions)
        except OSError as e:
            print(f"Error: {e}", file=sys.stderr)
            return 1
        if not accession_queries:
            print("Error: No accessions found in --accessions input", file=sys.stderr)
            return 1
        if args.verbose:
            print(f"Requested accessions: {len(accession_queries)}")

    if args.mapping_file is not None and not args.taxon and not accession_queries and not args.download:
        if mapping_flag_without_value:
            print(f"Global mapping file: {_get_default_mapping_path(base_dir, args.gtdb)}")
            return 0

        datasets = ["bac120", "ar53"] if args.dataset == "all" else [args.dataset]
        include_accessions = _collect_target_accessions_for_mapping(
            version=args.gtdb,
            datasets=datasets,
            taxon=None,
            only_rep=args.only_rep,
            base_dir=base_dir,
            mirror=args.mirror,
            verbose=args.verbose,
        )
        taxonomy_lookup, genome_length_lookup = _collect_taxonomy_lookup_for_mapping(
            version=args.gtdb,
            datasets=datasets,
            base_dir=base_dir,
            mirror=args.mirror,
            verbose=args.verbose,
            ensure_metadata=True,
        )
        mapping_path = build_mapping_file(
            args.gtdb,
            base_dir=base_dir,
            mapping_file=args.mapping_file,
            include_accessions=include_accessions,
            taxonomy_lookup=taxonomy_lookup,
            genome_length_lookup=genome_length_lookup,
            show_progress=True,
        )
        count = sum(1 for _ in open(mapping_path, "r", encoding="utf-8")) - 1  # exclude header
        print(f"Mapping file written to: {mapping_path}")
        print(f"Mapped genomes: {count}")
        return 0

    if args.mapping_file_mg is not None and not args.taxon and not args.download:
        datasets = ["bac120", "ar53"] if args.dataset == "all" else [args.dataset]
        mapping_path = _build_mg_mapping_file(
            version=args.gtdb,
            datasets=datasets,
            base_dir=base_dir,
            mapping_file=args.mapping_file_mg,
            mirror=args.mirror,
            verbose=args.verbose,
        )
        if mapping_path is None:
            return 1
        count = sum(1 for _ in open(mapping_path, "r", encoding="utf-8")) - 1  # exclude header
        print(f"Marker gene mapping file written to: {mapping_path}")
        print(f"Mapped genomes: {count}")
        return 0

    if args.mgnify_catalogues:
        catalogues = _list_mgnify_catalogues()
        if not catalogues:
            return 1
        print("Available MGnify genome catalogues:")
        for cat in catalogues:
            cat_id = cat.get("id", "")
            attrs = cat.get("attributes", {})
            name = attrs.get("name", "")
            biome = attrs.get("catalogue-biome-label", "")
            genome_count = attrs.get("genome-count", "?")
            print(f"  {cat_id:40s}  {genome_count:>8} genomes  {name or biome}")
        return 0

    if args.mgnify_filter:
        output_path = args.mgnify_output or Path(f"{args.mgnify_filter}_accession_path_map.tsv")
        if not output_path.is_absolute():
            output_path = Path.cwd() / output_path

        datasets = ["bac120", "ar53"] if args.dataset == "all" else [args.dataset]

        print(f"Fetching MGnify catalogue: {args.mgnify_filter}...")
        mgnify_accessions = _fetch_mgnify_accessions(args.mgnify_filter, verbose=args.verbose)
        if mgnify_accessions is None:
            return 1
        print(f"Fetched {len(mgnify_accessions)} accessions from MGnify.")

        count = _filter_gtdb_metadata_by_mgnify(
            version=args.gtdb,
            datasets=datasets,
            base_dir=base_dir,
            mgnify_accessions=mgnify_accessions,
            output_path=output_path,
            mirror=args.mirror,
            verbose=args.verbose,
        )
        print(f"Filtered mapping written to: {output_path}")
        print(f"Matching genomes: {count} (out of {len(mgnify_accessions)} MGnify accessions)")
        return 0

    # Handle --download flag (metadata only)
    if args.download:
        version_dir = setup_version_dir(args.gtdb, base_dir)
        datasets = ["bac120", "ar53"] if args.dataset == "all" else [args.dataset]
        ok = True
        for ds in datasets:
            print(f"Downloading metadata for {args.gtdb} ({ds}) from {args.mirror}...")
            metadata_file = download_metadata(
                args.gtdb,
                ds,
                version_dir,
                mirror=args.mirror,
                verbose=args.verbose
            )
            if metadata_file:
                print(f"✓ Metadata downloaded successfully for {ds}")
            else:
                print(f"✗ Failed to download metadata for {ds}", file=sys.stderr)
                ok = False
        return 0 if ok else 1

    # Handle --download-marker-genes flag
    if args.download_marker_genes:
        datasets = ["bac120", "ar53"] if args.dataset == "all" else [args.dataset]
        version_dir = setup_version_dir(args.gtdb, base_dir)
        marker_genes_dir = version_dir / "marker_genes"
        marker_genes_dir.mkdir(parents=True, exist_ok=True)
        ok = True
        for ds in datasets:
            try:
                url = get_marker_genes_url(args.gtdb, ds, mirror=args.mirror)
            except ValueError as e:
                print(f"Error: {e}", file=sys.stderr)
                ok = False
                continue

            filename = url.split("/")[-1]
            dest = marker_genes_dir / filename

            # Download only if tarball not already present (safe to skip after interrupted extraction)
            if dest.exists():
                print(f"Already downloaded: {dest}")
            else:
                print(f"Downloading marker genes ({ds}) from {url}...")
                if not download_file(url, dest, verbose=args.verbose):
                    print(f"✗ Failed to download marker genes for {ds}", file=sys.stderr)
                    ok = False
                    continue
                print(f"✓ Marker genes downloaded: {dest}")

            # Determine extracted directory name from tarball top-level entry
            with tarfile.open(dest, "r:gz") as tar:
                top_dirs = {m.name.split("/")[0] for m in tar.getmembers() if m.name.split("/")[0]}
            extracted_name = top_dirs.pop() if len(top_dirs) == 1 else filename.removesuffix(".tar.gz")
            extracted_dir = marker_genes_dir / extracted_name

            if not extracted_dir.exists():
                print(f"Extracting {filename}...")
                with tarfile.open(dest, "r:gz") as tar:
                    tar.extractall(marker_genes_dir)
                print(f"✓ Extracted to: {extracted_dir}")
            else:
                print(f"Already extracted: {extracted_dir}")

            mapping_path = _build_marker_genes_mapping(
                version=args.gtdb,
                dataset=ds,
                extracted_dir=extracted_dir,
                version_dir=version_dir,
                base_dir=base_dir,
                mirror=args.mirror,
                verbose=args.verbose,
            )
            if mapping_path:
                print(f"Marker genes mapping file: {mapping_path}")
            else:
                ok = False
        return 0 if ok else 1

    # Handle --build-mg flag
    if args.build_mg:
        datasets = ["bac120", "ar53"] if args.dataset == "all" else [args.dataset]
        version_dir = setup_version_dir(args.gtdb, base_dir)
        marker_genes_dir = version_dir / "marker_genes"
        ok = True
        for ds in datasets:
            extracted_dir = marker_genes_dir / f"{ds}_marker_genes_all_{args.gtdb}"
            if not (extracted_dir / "fna").exists():
                print(
                    f"Marker genes not found for {ds}: {extracted_dir / 'fna'}\n"
                    f"Download them first with: gtdb-dl --gtdb {args.gtdb} --download-marker-genes",
                    file=sys.stderr,
                )
                ok = False
                continue
            print(f"Building concatenated marker gene files for {ds}...")
            out_dir = _concat_marker_genes_for_dataset(extracted_dir, verbose=args.verbose)
            if out_dir is None:
                ok = False
            else:
                print(f"✓ Concatenated files written to: {out_dir}")
        return 0 if ok else 1

    # Handle taxon- and/or accession-based download
    if args.taxon or accession_queries:
        datasets = ["bac120", "ar53"] if args.dataset == "all" else [args.dataset]
        overall_success = True
        resolved_accessions: set = set()

        selection_parts = []
        if args.taxon:
            selection_parts.append(f"taxon: {args.taxon}")
        if accession_queries:
            selection_parts.append(f"{len(accession_queries)} requested accessions")
        selection_desc = ", ".join(selection_parts)

        for ds in datasets:
            print(f"Downloading genomes for {selection_desc} (dataset: {ds})")
            success = download_genomes(
                args.gtdb,
                taxon=args.taxon,
                accessions=accession_queries or None,
                dataset=ds,
                mirror=args.mirror,
                base_dir=base_dir,
                output_dir=args.output,
                flat=args.flat,
                flag_rep=args.flag_rep,
                only_rep=args.only_rep,
                ignore_prefix=args.ignore_prefix,
                failed_file=args.failed_file,
                verbose=args.verbose,
                dry_run=args.dry_run,
                resolved_accessions=resolved_accessions,
                genome_type=args.genome_type,
            )
            overall_success = overall_success and success

            # Pure accession queries need no further datasets once everything matched
            # (parsing the bac120 metadata is expensive).
            if (
                accession_queries
                and not args.taxon
                and len(resolved_accessions) == len(accession_queries)
            ):
                break

        if accession_queries:
            unknown = [acc for acc in accession_queries if acc not in resolved_accessions]
            if unknown:
                print(
                    f"\n✗ {len(unknown)} requested accessions not found in "
                    f"GTDB {args.gtdb} metadata:",
                    file=sys.stderr,
                )
                for accession in unknown[:20]:
                    print(f"  - {accession}", file=sys.stderr)
                if len(unknown) > 20:
                    print(f"  ... and {len(unknown) - 20} more", file=sys.stderr)
                if not args.ignore_prefix:
                    print(
                        "  (try --ignore-prefix to match GCA_/GCF_ counterparts)",
                        file=sys.stderr,
                    )
                overall_success = False

        # Optional project-specific mapping file after download runs.
        if args.mapping_file is not None and not mapping_flag_without_value:
            include_accessions = _collect_target_accessions_for_mapping(
                version=args.gtdb,
                datasets=datasets,
                taxon=args.taxon,
                only_rep=args.only_rep,
                base_dir=base_dir,
                mirror=args.mirror,
                verbose=args.verbose,
                accessions=accession_queries or None,
                ignore_prefix=args.ignore_prefix,
            )
            taxonomy_lookup, genome_length_lookup = _collect_taxonomy_lookup_for_mapping(
                version=args.gtdb,
                datasets=datasets,
                base_dir=base_dir,
                mirror=args.mirror,
                verbose=args.verbose,
                ensure_metadata=True,
            )
            mapping_path = build_mapping_file(
                args.gtdb,
                base_dir=base_dir,
                mapping_file=args.mapping_file,
                include_accessions=include_accessions,
                taxonomy_lookup=taxonomy_lookup,
                genome_length_lookup=genome_length_lookup,
                show_progress=True,
            )
            print(f"Requested mapping file written to: {mapping_path}")
        elif args.mapping_file is not None and mapping_flag_without_value:
            print(f"Global mapping file: {_get_default_mapping_path(base_dir, args.gtdb)}")

        return 0 if overall_success else 1

    # If neither --download nor --taxon, show help
    parser.print_help()
    return 1


if __name__ == "__main__":
    sys.exit(main())
