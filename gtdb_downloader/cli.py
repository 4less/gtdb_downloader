"""Command-line interface for GTDB downloader"""

import argparse
import re
import sys
import time
from pathlib import Path
from typing import Optional, List, Tuple, Dict

import requests

from gtdb_downloader.config import get_base_dir, GTDB_VERSIONS
from gtdb_downloader.downloader import (
    check_aria2c_available,
    download_file,
    download_files_aria2,
    download_metadata,
    generate_download_links,
    resolve_download_link,
)
from gtdb_downloader.metadata import MetadataParser


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
    """Resolve mapping file path, using the default location when omitted or relative."""
    default_path = _get_default_mapping_path(base_dir, version)
    if mapping_file is None:
        return default_path
    if mapping_file.is_absolute():
        return mapping_file
    return default_path.parent / mapping_file


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
    try:
        response = requests.get(datasets_url, timeout=timeout_seconds)
        response.raise_for_status()
    except Exception:
        return datasets_url, None

    # Strip tags for a resilient plain-text search.
    plain_text = re.sub(r"<[^>]+>", " ", response.text)
    plain_text = re.sub(r"\s+", " ", plain_text)
    match = re.search(r"Status:\s*(.{1,220}?)\s{1,}(?:This record|Actions|Download|datasets|API|FTP|$)", plain_text, re.IGNORECASE)
    if match:
        return datasets_url, match.group(1).strip()

    idx = plain_text.lower().find("status:")
    if idx >= 0:
        snippet = plain_text[idx:idx + 240].strip()
        return datasets_url, snippet

    return datasets_url, None


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


def build_mapping_file(
    version: str,
    base_dir: Optional[Path] = None,
    mapping_file: Optional[Path] = None,
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
    resolved_mapping_path.parent.mkdir(parents=True, exist_ok=True)
    tmp_path = resolved_mapping_path.with_suffix(resolved_mapping_path.suffix + ".tmp")

    if show_progress:
        print(f"Building mapping file from: {genomes_dir}")

    with open(tmp_path, "w", encoding="utf-8") as handle:
        for count, (accession, genome_path) in enumerate(sorted(mappings.items()), start=1):
            handle.write(f"{accession}\t{genome_path}\n")
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


def download_genomes_for_taxon(
    taxon: str,
    version: str,
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
    dry_run: bool = False
) -> bool:
    """
    Download all genomes for a specific taxon
    
    Args:
        taxon: Taxon to search for
        version: GTDB version
        dataset: Dataset type (bac120 or ar53)
        mirror: Mirror to use for download
        base_dir: Base directory for GTDB data
        output_dir: Output directory for symlink taxonomy structure (not genomes)
        verbose: Verbose output
        dry_run: Don't actually download, just show what would be downloaded
    
    Returns:
        True if successful, False otherwise
    """
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
    print(f"Filtering genomes for taxon query: {taxon}")
    matching_genomes = parser.get_genomes_by_taxon(taxon)
    
    if not matching_genomes:
        print(f"No genomes found for taxon: {taxon}", file=sys.stderr)
        return False

    if only_rep:
        matching_genomes = [
            genome_id
            for genome_id in matching_genomes
            if parser.is_species_cluster_representative(genome_id)
        ]
        if not matching_genomes:
            print(
                f"No species representative genomes found for taxon: {taxon}",
                file=sys.stderr,
            )
            return False
    
    print(f"Found {len(matching_genomes)} genomes for taxon: {taxon}")
    
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
            )

            fallback_downloads: List[Tuple[str, Path]] = []
            for item in failed_chunk_items:
                genome_id = item["genome_id"]
                genome_metadata = item["genome_metadata"]
                if verbose:
                    print(f"  Resolving fallback for {genome_id}...")

                accession = _extract_ncbi_accession(genome_id, genome_metadata)  # type: ignore[arg-type]
                if accession:
                    ncbi_checked_genomes += 1
                    if accession not in status_cache:
                        _, status_text = _fetch_ncbi_datasets_status(accession)
                        status_cache[accession] = status_text
                    status_text = status_cache.get(accession)
                    if status_text:
                        failed_status_notes[genome_id] = f"status:{status_text}"
                    else:
                        failed_status_notes[genome_id] = "status:unavailable_or_not_found"
                        ncbi_missing_status_genomes += 1
                    if status_text and "suppressed" in status_text.lower():
                        datasets_url = f"https://www.ncbi.nlm.nih.gov/datasets/genome/{accession}/"
                        suppressed_msg = (
                            f"Status: {status_text} ({datasets_url}); "
                            f"stopping further lookup for this genome"
                        )
                        print(suppressed_msg, file=sys.stderr)
                        suppressed_genomes.append(genome_id)  # type: ignore[arg-type]
                        continue
                else:
                    failed_status_notes[genome_id] = "status:no_accession"

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

    downloaded_count = 0
    for item in downloadable:
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

    mapping_path = build_mapping_file(version, base_dir=base_dir, show_progress=verbose)
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
        "--dataset",
        choices=["bac120", "ar53"],
        default="bac120",
        help="Dataset type (default: bac120)"
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
            "Write or refresh a TSV mapping file of accession to local raw genome path. "
            "Can be used without --taxon to build the mapping from existing downloaded genomes."
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
    
    base_dir = args.base_dir or get_base_dir()
    
    if args.verbose:
        print(f"Using base directory: {base_dir}")

    if args.mapping_file is not None and not args.taxon and not args.download:
        mapping_path = build_mapping_file(
            args.gtdb,
            base_dir=base_dir,
            mapping_file=args.mapping_file,
            show_progress=True,
        )
        count = len(_collect_existing_genome_mappings(_get_shared_raw_genomes_dir(base_dir)))
        print(f"Mapping file written to: {mapping_path}")
        print(f"Mapped genomes: {count}")
        return 0
    
    # Handle --download flag (metadata only)
    if args.download:
        print(f"Downloading metadata for {args.gtdb} ({args.dataset}) from {args.mirror}...")
        version_dir = setup_version_dir(args.gtdb, base_dir)
        metadata_file = download_metadata(
            args.gtdb,
            args.dataset,
            version_dir,
            mirror=args.mirror,
            verbose=args.verbose
        )
        if metadata_file:
            print(f"✓ Metadata downloaded successfully")
            return 0
        else:
            print(f"✗ Failed to download metadata", file=sys.stderr)
            return 1
    
    # Handle taxon-based download
    if args.taxon:
        print(f"Downloading genomes for taxon: {args.taxon}")
        success = download_genomes_for_taxon(
            args.taxon,
            args.gtdb,
            dataset=args.dataset,
            mirror=args.mirror,
            base_dir=base_dir,
            output_dir=args.output,
            flat=args.flat,
            flag_rep=args.flag_rep,
            only_rep=args.only_rep,
            ignore_prefix=args.ignore_prefix,
            failed_file=args.failed_file,
            verbose=args.verbose,
            dry_run=args.dry_run
        )
        if args.mapping_file is not None:
            mapping_path = build_mapping_file(
                args.gtdb,
                base_dir=base_dir,
                mapping_file=args.mapping_file,
                show_progress=True,
            )
            print(f"Requested mapping file written to: {mapping_path}")
        return 0 if success else 1

    # If neither --download nor --taxon, show help
    parser.print_help()
    return 1


if __name__ == "__main__":
    sys.exit(main())
