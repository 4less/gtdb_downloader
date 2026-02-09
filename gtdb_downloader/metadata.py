"""Metadata parsing and handling"""

import csv
import gzip
import tarfile
from pathlib import Path
from typing import Dict, List, Tuple, Optional
import io


class MetadataParser:
    """Parse GTDB metadata files"""
    
    def __init__(self, metadata_path: Path):
        """
        Initialize metadata parser
        
        Args:
            metadata_path: Path to metadata file or tar.gz
        """
        self.metadata_path = Path(metadata_path)
        self.data = {}
        self._load_metadata()
    
    def _load_metadata(self):
        """Load metadata from file"""
        if self.metadata_path.suffix == ".gz":
            self._load_from_gzip()
        else:
            self._load_from_csv()
    
    def _load_from_gzip(self):
        """Load from gzipped TSV file"""
        with gzip.open(self.metadata_path, "rt") as f:
            reader = csv.DictReader(f, delimiter="\t")
            for row in reader:
                genome_id = row.get("accession") or row.get("Genome")
                if genome_id:
                    self.data[genome_id] = row
    
    def _load_from_csv(self):
        """Load from CSV file"""
        with open(self.metadata_path, "r") as f:
            reader = csv.DictReader(f, delimiter="\t")
            for row in reader:
                genome_id = row.get("accession") or row.get("Genome")
                if genome_id:
                    self.data[genome_id] = row
    
    def get_genomes_by_taxon(self, taxon: str, field: str = "gtdb_taxonomy") -> List[str]:
        """
        Get all genomes that match a specific taxon
        
        Args:
            taxon: Taxon name (e.g., 'd__Bacteria;p__Firmicutes')
            field: Taxonomy field name (usually 'gtdb_taxonomy')
        
        Returns:
            List of genome IDs
        """
        matching_genomes = []
        taxon_lower = taxon.lower()
        
        for genome_id, row in self.data.items():
            taxonomy = row.get(field, "").lower()
            if taxon_lower in taxonomy:
                matching_genomes.append(genome_id)
        
        return matching_genomes
    
    def get_genome_metadata(self, genome_id: str) -> Optional[dict]:
        """
        Get the full metadata row for a genome
        
        Args:
            genome_id: Genome accession
        
        Returns:
            Metadata dictionary or None if not found
        """
        return self.data.get(genome_id)
    
    def get_taxon_path(self, genome_id: str, field: str = "gtdb_taxonomy") -> Optional[Tuple[str, str]]:
        """
        Get taxonomy path for a genome
        
        Args:
            genome_id: Genome accession
            field: Taxonomy field name
        
        Returns:
            Tuple of (full_taxonomy_string, dataset_type) or None
        """
        if genome_id not in self.data:
            return None
        
        row = self.data[genome_id]
        taxonomy = row.get(field, "")
        
        # Determine dataset type from genome ID or field
        dataset = "bac120" if not genome_id.startswith("GCF") else "ar53"
        if "Archaea" in taxonomy:
            dataset = "ar53"
        
        return taxonomy, dataset
    
    def parse_taxonomy_to_path(self, taxonomy_str: str) -> List[str]:
        """
        Parse GTDB taxonomy string into folder hierarchy
        
        Args:
            taxonomy_str: GTDB taxonomy string (e.g., 'd__Bacteria;p__Firmicutes;c__Clostridia')
        
        Returns:
            List of folder names
        """
        if not taxonomy_str:
            return []
        
        parts = []
        for rank in taxonomy_str.split(";"):
            if rank and rank != "":
                # Remove rank prefix (d__, p__, c__, etc.)
                name = rank.split("__")[-1].strip()
                if name:
                    parts.append(name)
        
        return parts

    def get_taxonomy_components(self, taxonomy_str: str) -> List[str]:
        """
        Return taxonomy components including rank prefixes (e.g., 'd__Bacteria')

        Args:
            taxonomy_str: GTDB taxonomy string

        Returns:
            List of components with prefixes
        """
        if not taxonomy_str:
            return []
        return [part.strip() for part in taxonomy_str.split(";") if part.strip()]

    def get_taxon_component_at_rank(self, taxonomy_str: str, rank: str) -> Optional[str]:
        """
        Return the taxonomy component at the requested rank (with prefix).

        Args:
            taxonomy_str: GTDB taxonomy string
            rank: requested rank (e.g., 's', 'species', 'g', 'genus')

        Returns:
            The matching component string (e.g., 's__Escherichia coli') or None
        """
        if not taxonomy_str:
            return None

        # Normalize rank to prefix (single letter)
        rank_map = {
            "domain": "d",
            "d": "d",
            "phylum": "p",
            "p": "p",
            "class": "c",
            "c": "c",
            "order": "o",
            "o": "o",
            "family": "f",
            "f": "f",
            "genus": "g",
            "g": "g",
            "species": "s",
            "s": "s",
        }

        r = rank.lower()
        if r not in rank_map:
            return None
        prefix = rank_map[r] + "__"

        for comp in self.get_taxonomy_components(taxonomy_str):
            if comp.startswith(prefix):
                return comp
        return None
