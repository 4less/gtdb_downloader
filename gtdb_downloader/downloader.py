"""Download utilities for GTDB"""

import subprocess
import sys
from pathlib import Path
from typing import Optional, List


def check_aria2c_available() -> bool:
    """Check if aria2c is available in PATH"""
    try:
        result = subprocess.run(
            ["which", "aria2c"],
            capture_output=True,
            text=True,
            timeout=5
        )
        return result.returncode == 0
    except Exception:
        return False


def check_wget_available() -> bool:
    """Check if wget is available in PATH"""
    try:
        result = subprocess.run(
            ["which", "wget"],
            capture_output=True,
            text=True,
            timeout=5
        )
        return result.returncode == 0
    except Exception:
        return False


def generate_download_link(genome_metadata: dict) -> str:
    """
    Generate a download link for a specific genome from NCBI
    
    Args:
        genome_metadata: Dictionary with genome metadata (from MetadataParser)
    
    Returns:
        Download URL for the genome
    
    Raises:
        ValueError: If required fields are missing
    """
    accession = genome_metadata.get("accession")
    assembly_name = genome_metadata.get("ncbi_assembly_name")
    
    if not accession or not assembly_name:
        raise ValueError(f"Missing required fields for genome: accession={accession}, assembly_name={assembly_name}")
    
    # Parse accession: RS_GCF_034719275.1 or GB_GCA_034719275.1
    if accession.startswith("RS_"):
        # RefSeq
        acc = accession[3:]  # Remove RS_ prefix
    elif accession.startswith("GB_"):
        # GenBank
        acc = accession[3:]  # Remove GB_ prefix
    else:
        # Assume it's already just the accession
        acc = accession
    
    # Extract parts: GCF_034719275 -> 034, 719, 275
    # Format: GC[AF]_XXXXXX...
    if not acc.startswith(("GCA_", "GCF_")):
        raise ValueError(f"Unknown accession format: {acc}")
    
    # Extract numeric part after GCA_ or GCF_
    parts = acc.split("_")
    if len(parts) < 2:
        raise ValueError(f"Invalid accession format: {acc}")
    
    # Get the numeric identifier (e.g., 034719275 from GCF_034719275.1)
    numeric_id = parts[1].split(".")[0]  # Remove version number if present
    
    # Split into groups of 3 digits from the left
    # 034719275 -> 034/719/275
    d1 = numeric_id[:3]           # 034
    d2 = numeric_id[3:6] if len(numeric_id) >= 6 else ""  # 719
    d3 = numeric_id[6:] if len(numeric_id) > 6 else ""    # 275
    
    if not d2 or not d3:
        raise ValueError(f"Accession numeric ID too short: {numeric_id}")
    
    # Construct full accession with version
    full_accession = acc
    
    # Build URL: https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/034/719/275/GCA_034719275.1_ASM3471927v1/GCA_034719275.1_ASM3471927v1_genomic.fna.gz
    base_url = "https://ftp.ncbi.nlm.nih.gov/genomes/all"
    
    # Include GCA/GCF in path
    prefix = acc.split("_")[0]  # GCA or GCF
    
    # Directory path
    dir_path = f"{base_url}/{prefix}/{d1}/{d2}/{d3}/{full_accession}_{assembly_name}"
    
    # Filename
    filename = f"{full_accession}_{assembly_name}_genomic.fna.gz"
    
    # Complete URL
    url = f"{dir_path}/{filename}"
    
    return url


def download_file_aria2(
    url: str,
    output_path: Path,
    max_connections: int = 4,
    verbose: bool = False
) -> bool:
    """
    Download a file using aria2c
    
    Args:
        url: URL to download
        output_path: Path where to save the file
        max_connections: Number of parallel connections
        verbose: Print verbose output
    
    Returns:
        True if successful, False otherwise
    """
    output_path.parent.mkdir(parents=True, exist_ok=True)
    
    cmd = [
        "aria2c",
        f"--max-connection-per-server={max_connections}",
        f"--split={max_connections}",
        "-x", str(max_connections),
        "-k", "1M",
        "-o", output_path.name,
        "-d", str(output_path.parent),
    ]
    
    if not verbose:
        cmd.append("--quiet")
    
    cmd.append(url)
    
    try:
        result = subprocess.run(cmd, capture_output=not verbose, text=True, timeout=3600)
        return result.returncode == 0
    except subprocess.TimeoutExpired:
        print(f"Download timeout for {url}", file=sys.stderr)
        return False
    except Exception as e:
        print(f"Error downloading {url}: {e}", file=sys.stderr)
        return False


def download_file_wget(url: str, output_path: Path, verbose: bool = False) -> bool:
    """
    Download a file using wget
    
    Args:
        url: URL to download
        output_path: Path where to save the file
        verbose: Print verbose output
    
    Returns:
        True if successful, False otherwise
    """
    output_path.parent.mkdir(parents=True, exist_ok=True)
    
    cmd = [
        "wget",
        "-O", str(output_path),
        url
    ]
    
    if not verbose:
        cmd.insert(1, "-q")
    
    try:
        result = subprocess.run(cmd, capture_output=not verbose, text=True, timeout=3600)
        return result.returncode == 0
    except subprocess.TimeoutExpired:
        print(f"Download timeout for {url}", file=sys.stderr)
        return False
    except Exception as e:
        print(f"Error downloading {url}: {e}", file=sys.stderr)
        return False


def download_file(
    url: str,
    output_path: Path,
    verbose: bool = False,
    use_aria2: bool = True
) -> bool:
    """
    Download a file using available tools (aria2c preferred, fallback to wget)
    
    Args:
        url: URL to download
        output_path: Path where to save the file
        verbose: Print verbose output
        use_aria2: Try to use aria2c if available
    
    Returns:
        True if successful, False otherwise
    """
    if use_aria2 and check_aria2c_available():
        return download_file_aria2(url, output_path, verbose=verbose)
    elif check_wget_available():
        return download_file_wget(url, output_path, verbose=verbose)
    else:
        print("Error: Neither aria2c nor wget found in PATH", file=sys.stderr)
        print("Please install aria2 or wget to download files", file=sys.stderr)
        return False


def download_metadata(
    version: str,
    dataset: str,
    output_dir: Path,
    mirror: str = "europe",
    verbose: bool = False
) -> Optional[Path]:
    """
    Download GTDB metadata file for a specific version
    
    Args:
        version: GTDB version (e.g., 'r226')
        dataset: Dataset type ('bac120' or 'ar53')
        output_dir: Directory to store the metadata
        mirror: Mirror to use ('europe', 'asia-pacific1', 'asia-pacific2')
        verbose: Print verbose output
    
    Returns:
        Path to the metadata file, or None if download failed
    """
    from gtdb_downloader.config import get_metadata_url, GTDB_VERSIONS
    
    if version not in GTDB_VERSIONS:
        print(f"Error: Unknown GTDB version {version}", file=sys.stderr)
        print(f"Available versions: {', '.join(GTDB_VERSIONS.keys())}", file=sys.stderr)
        return None
    
    output_dir.mkdir(parents=True, exist_ok=True)
    
    try:
        url = get_metadata_url(version, dataset, mirror)
    except ValueError as e:
        print(f"Error: {e}", file=sys.stderr)
        return None
    
    filename = url.split("/")[-1]
    metadata_file = output_dir / filename
    
    # Check if already downloaded
    if metadata_file.exists():
        if verbose:
            print(f"Metadata file already exists: {metadata_file}")
        return metadata_file
    
    if verbose:
        print(f"Downloading metadata from: {url}")
    
    if not download_file(url, metadata_file, verbose=verbose):
        print(f"Failed to download metadata", file=sys.stderr)
        return None
    
    if verbose:
        print(f"Successfully downloaded metadata to: {metadata_file}")
    
    return metadata_file
