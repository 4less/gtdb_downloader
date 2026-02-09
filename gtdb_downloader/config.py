"""Configuration for GTDB downloader"""

import os
from pathlib import Path

# Default base directory for GTDB data
DEFAULT_BASE_DIR = Path.home() / ".gtdb_downloader"

# Mirror URLs for GTDB
MIRRORS = {
    "europe": "https://data.gtdb.aau.ecogenomic.org/releases/",
    "asia-pacific1": "https://data.gtdb.ecogenomic.org/releases/",
    "asia-pacific2": "https://data.ace.uq.edu.au/public/gtdb/data/releases/",
}

# GTDB versions and their release paths
GTDB_VERSIONS = {
    "r207": "release207/207.0",
    "r214": "release214/214.1",
    "r220": "release220/220.0",
    "r226": "release226/226.0",
}

# GTDB datasets
DATASETS = {
    "bac120": "bac120_metadata_{}.tsv.gz",
    "ar53": "ar53_metadata_{}.tsv.gz",
}

# Get base directory from environment or use default
def get_base_dir():
    """Get the base directory for GTDB data storage"""
    base_dir = os.environ.get("GTDBDL_DATA", DEFAULT_BASE_DIR)
    return Path(base_dir)


def get_metadata_url(version: str, dataset: str, mirror: str = "europe") -> str:
    """
    Get the download URL for metadata
    
    Args:
        version: GTDB version (e.g., 'r226')
        dataset: Dataset type ('bac120' or 'ar53')
        mirror: Mirror to use ('europe', 'asia-pacific1', 'asia-pacific2')
    
    Returns:
        Download URL for the metadata file
    """
    if mirror not in MIRRORS:
        raise ValueError(f"Unknown mirror: {mirror}. Available: {', '.join(MIRRORS.keys())}")
    
    if version not in GTDB_VERSIONS:
        raise ValueError(f"Unknown version: {version}. Available: {', '.join(GTDB_VERSIONS.keys())}")
    
    if dataset not in DATASETS:
        raise ValueError(f"Unknown dataset: {dataset}. Available: {', '.join(DATASETS.keys())}")
    
    base_url = MIRRORS[mirror]
    release_path = GTDB_VERSIONS[version]
    filename = DATASETS[dataset].format(version)
    
    return f"{base_url}{release_path}/{filename}"

