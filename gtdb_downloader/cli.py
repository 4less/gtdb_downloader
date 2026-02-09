"""Command-line interface for GTDB downloader"""

import argparse
import sys
from pathlib import Path
from typing import Optional

from gtdb_downloader.config import get_base_dir, GTDB_VERSIONS
from gtdb_downloader.downloader import download_metadata, download_file, generate_download_link
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
    
    downloaded_count = 0
    failed_count = 0
    
    for i, genome_id in enumerate(matching_genomes, 1):
        if verbose:
            print(f"\n[{i}/{len(matching_genomes)}] Processing {genome_id}...")
        
        # Get metadata row for this genome
        genome_metadata = parser.get_genome_metadata(genome_id)
        if genome_metadata is None:
            if verbose:
                print(f"  Skipped: Could not find genome metadata")
            failed_count += 1
            continue
        
        # Get taxonomy info
        tax_info = parser.get_taxon_path(genome_id)
        if tax_info is None:
            if verbose:
                print(f"  Skipped: Could not find taxonomy info")
            failed_count += 1
            continue
        
        taxonomy_str, det_dataset = tax_info
        
        # Generate download link
        try:
            download_url = generate_download_link(genome_metadata)
            # Filename is the last part of the URL (after last /)
            genome_filename = download_url.split("/")[-1]
        except Exception as e:
            if verbose:
                print(f"  Skipped: Could not generate download link: {e}")
            failed_count += 1
            continue
        
        # Download genome
        genome_path = genomes_dir / genome_filename
        
        if genome_path.exists():
            if verbose:
                print(f"  Already exists: {genome_path}")
        else:
            if verbose:
                print(f"  Downloading from: {download_url}")
            
            if not download_file(download_url, genome_path, verbose=verbose):
                if verbose:
                    print(f"  Failed to download")
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

            link_path = tax_folder / genome_filename
            if not link_path.exists():
                link_path.symlink_to(genome_path.resolve())
                if verbose:
                    print(f"  Symlinked: {link_path}")
        except Exception as e:
            if verbose:
                print(f"  Warning: Could not create symlink: {e}")
        
        downloaded_count += 1
    
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

  # Download all Firmicutes genomes
  gtdb-dl --gtdb r226 --taxon "Firmicutes"

  # Download with verbose output
  gtdb-dl --gtdb r226 --taxon "d__Bacteria;p__Firmicutes" -v

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
        help="Taxon to search for (e.g., 'Firmicutes' or full GTDB taxonomy path)"
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
        "--output",
        "-o",
        type=Path,
        help="Output directory for symlink taxonomy structure (default: ~/.gtdb_downloader/{version}/genomes/taxonomy)"
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
            verbose=args.verbose,
            dry_run=args.dry_run
        )
        return 0 if success else 1
    
    # If neither --download nor --taxon, show help
    parser.print_help()
    return 1


if __name__ == "__main__":
    sys.exit(main())
