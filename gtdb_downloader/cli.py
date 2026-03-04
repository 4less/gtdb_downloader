"""Command-line interface for GTDB downloader"""

import argparse
import re
import sys
from pathlib import Path
from typing import Optional, List, Tuple, Dict

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


def _resolve_mapping_path(base_dir: Path, version: str, mapping_file: Optional[Path]) -> Path:
    """Resolve mapping file path, using the default location when omitted or relative."""
    default_path = _get_default_mapping_path(base_dir, version)
    if mapping_file is None:
        return default_path
    if mapping_file.is_absolute():
        return mapping_file
    return default_path.parent / mapping_file


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

    genomes_dir = base_dir / version / "genomes" / "raw"
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
        if path is not None and path.exists():
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
    ignore_prefix: bool = False,
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
    
    # Genomes always go to base_dir / version / genomes / raw
    genomes_dir = base_dir / version / "genomes" / "raw"
    genomes_dir.mkdir(parents=True, exist_ok=True)
    
    # Symlink taxonomy structure goes to output_dir (or base_dir if not specified)
    if output_dir is None:
        taxonomy_dir = base_dir / version / "genomes" / "taxonomy"
    else:
        taxonomy_dir = output_dir
    
    taxonomy_dir.mkdir(parents=True, exist_ok=True)
    
    # Download metadata if needed
    metadata_file = download_metadata(version, dataset, version_dir, mirror=mirror, verbose=verbose)
    if metadata_file is None:
        return False
    
    # Parse metadata
    try:
        parser = MetadataParser(metadata_file)
    except Exception as e:
        print(f"Error parsing metadata: {e}", file=sys.stderr)
        return False
    
    # Find matching genomes
    matching_genomes = parser.get_genomes_by_taxon(taxon)
    
    if not matching_genomes:
        print(f"No genomes found for taxon: {taxon}", file=sys.stderr)
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
    existing_genomes = _index_existing_genomes(genomes_dir, ignore_prefix)

    for i, genome_id in enumerate(matching_genomes, 1):
        if verbose:
            print(f"\n[{i}/{len(matching_genomes)}] Processing {genome_id}...")

        genome_metadata = parser.get_genome_metadata(genome_id)
        if genome_metadata is None:
            if verbose:
                print("  Skipped: Could not find genome metadata")
            failed_count += 1
            continue

        tax_info = parser.get_taxon_path(genome_id)
        if tax_info is None:
            if verbose:
                print("  Skipped: Could not find taxonomy info")
            failed_count += 1
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
            continue

        genome_path = genomes_dir / genome_filename
        existing_genome_path = _find_existing_genome_path(existing_genomes, genome_metadata, ignore_prefix)
        if existing_genome_path is not None:
            genome_path = existing_genome_path
        is_species_rep = parser.is_species_cluster_representative(genome_id)
        cluster_rep = parser.get_species_cluster_representative(genome_id) or "unknown_cluster"
        downloadable.append({
            "genome_id": genome_id,
            "taxonomy_str": taxonomy_str,
            "download_url": download_url,
            "genome_path": genome_path,
            "is_species_rep": is_species_rep,
            "cluster_rep": cluster_rep,
            "genome_metadata": genome_metadata,
            "download_urls": download_urls,
        })

        if flag_rep and is_species_rep:
            representative_by_cluster.setdefault(cluster_rep, []).append(genome_id)

        if verbose and not genome_path.exists():
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
        item for item in downloadable if not item["genome_path"].exists()  # type: ignore[union-attr]
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

            primary_downloads: List[Tuple[str, Path]] = [
                (
                    item["download_urls"][0],  # type: ignore[index]
                    item["genome_path"],
                )
                for item in chunk
                if not item["genome_path"].exists()  # type: ignore[union-attr]
            ]

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
                if not item["genome_path"].exists() and not batch_results.get(item["genome_path"], False)  # type: ignore[union-attr]
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
                    continue

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
        is_species_rep = item["is_species_rep"]

        if not genome_path.exists() and not batch_results.get(genome_path, False):
            if verbose:
                print(f"  Failed to download {genome_id}")
            failed_count += 1
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
    
    print(f"\n✓ Downloaded: {downloaded_count}")
    print(f"✗ Failed: {failed_count}")
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
        help="Taxon to search for (e.g., 'Bacillota' or full GTDB taxonomy path)"
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
        count = len(_collect_existing_genome_mappings(base_dir / args.gtdb / "genomes" / "raw"))
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
            ignore_prefix=args.ignore_prefix,
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
