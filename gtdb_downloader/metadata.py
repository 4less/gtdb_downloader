"""Metadata parsing and handling"""

import csv
import gzip
import tarfile
from pathlib import Path
from typing import Dict, Iterable, List, Set, Tuple, Optional
import io


GENOME_TYPE_CATEGORIES = {
    "isolate": "none",
    "mag": "derived from metagenome",
    "sag": "derived from single cell",
    "env": "derived from environmental sample",
}


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
        self._accession_index: Optional[Dict[str, str]] = None
        self._accession_suffix_index: Optional[Dict[str, str]] = None
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
        Get all genomes that match a specific taxon exactly.
        
        Args:
            taxon: Taxon name or taxonomy path (e.g., 'Bacillota' or
                'd__Bacteria;p__Bacillota;c__Clostridia')
            field: Taxonomy field name (usually 'gtdb_taxonomy')
        
        Returns:
            List of genome IDs
        """
        matching_genomes = []
        query = taxon.strip()
        if not query:
            return matching_genomes

        # Support comma-separated OR queries, e.g. "Bacteria,Archaea".
        comma_queries = [part.strip() for part in query.split(",") if part.strip()]
        if len(comma_queries) > 1:
            seen = set()
            for subquery in comma_queries:
                for genome_id in self.get_genomes_by_taxon(subquery, field=field):
                    if genome_id not in seen:
                        seen.add(genome_id)
                        matching_genomes.append(genome_id)
            return matching_genomes

        # Normalise underscores to spaces in the part after the rank prefix so
        # that file-name-style queries like s__Bacteroides_ovatus match the
        # taxonomy string s__Bacteroides ovatus.
        def _norm_underscores(component: str) -> str:
            if "__" in component:
                prefix, name = component.split("__", 1)
                return f"{prefix}__{name.replace('_', ' ')}"
            return component.replace("_", " ")

        query_components = [_norm_underscores(part.strip()) for part in query.split(";") if part.strip()]
        query_components_lower = [part.lower() for part in query_components]



        def _strip_rank_prefix(component: str) -> str:
            if "__" in component:
                return component.split("__", 1)[1]
            return component

        query_names_lower = [_strip_rank_prefix(part).lower() for part in query_components]

        for genome_id, row in self.data.items():
            taxonomy = row.get(field, "")
            if not taxonomy:
                continue

            taxonomy_components = self.get_taxonomy_components(taxonomy)
            taxonomy_components_lower = [part.lower() for part in taxonomy_components]
            taxonomy_names_lower = [_strip_rank_prefix(part).lower() for part in taxonomy_components]

            # Path query: require exact component-by-component match in order,
            # with optional rank prefixes in the query (prefixes must match when supplied).
            if len(query_components_lower) > 1:
                if len(query_components_lower) > len(taxonomy_components_lower):
                    continue

                is_match = True
                for i, query_part in enumerate(query_components_lower):
                    tax_part = taxonomy_components_lower[i]
                    tax_name = taxonomy_names_lower[i]

                    if "__" in query_part:
                        if query_part != tax_part:
                            is_match = False
                            break
                    else:
                        if query_names_lower[i] != tax_name:
                            is_match = False
                            break

                if is_match:
                    matching_genomes.append(genome_id)
                continue

            # Single taxon query: exact match to a taxonomy component only (not substring).
            query_part = query_components_lower[0]
            query_name = query_names_lower[0]
            if "__" in query_part:
                if query_part in taxonomy_components_lower:
                    matching_genomes.append(genome_id)
            else:
                if query_name in taxonomy_names_lower:
                    matching_genomes.append(genome_id)
        
        return matching_genomes
    
    @staticmethod
    def _strip_accession_version(accession: str) -> str:
        """Drop the trailing assembly version (GCF_000970205.1 -> GCF_000970205)."""
        return accession.split(".", 1)[0]

    def _build_accession_index(self) -> None:
        """
        Build lookup tables from accession spellings to genome IDs.

        Two tables are built: one keyed by the normalized accession (with and
        without the assembly version) and one keyed by the numeric part only,
        so GCA_/GCF_ prefixes can be treated as interchangeable on request.
        """
        exact: Dict[str, str] = {}
        by_suffix: Dict[str, str] = {}

        for genome_id, row in self.data.items():
            raw_accession = row.get("accession") or row.get("Genome") or genome_id
            normalized = self._normalize_accession(str(raw_accession).upper())
            if not normalized:
                continue

            exact.setdefault(normalized, genome_id)
            exact.setdefault(self._strip_accession_version(normalized), genome_id)

            if normalized.startswith(("GCA_", "GCF_")):
                suffix = normalized[4:]
                by_suffix.setdefault(suffix, genome_id)
                by_suffix.setdefault(self._strip_accession_version(suffix), genome_id)

        self._accession_index = exact
        self._accession_suffix_index = by_suffix

    def get_genomes_by_accessions(
        self,
        accessions: Iterable[str],
        ignore_prefix: bool = False,
    ) -> Tuple[List[str], Set[str]]:
        """
        Look up genomes by accession.

        Accepts GTDB (RS_/GB_ prefixed) and plain NCBI accessions, with or
        without the assembly version suffix, in any capitalization.

        Args:
            accessions: Requested accessions
            ignore_prefix: Treat GCA_ and GCF_ as interchangeable

        Returns:
            Tuple of (matching genome IDs, set of requested accessions that matched)
        """
        if self._accession_index is None or self._accession_suffix_index is None:
            self._build_accession_index()

        genome_ids: List[str] = []
        matched_queries: Set[str] = set()
        seen_genomes: Set[str] = set()

        for query in accessions:
            key = self._normalize_accession(str(query).upper())
            if not key:
                continue

            genome_id = None
            for candidate in (key, self._strip_accession_version(key)):
                genome_id = self._accession_index.get(candidate)
                if genome_id is not None:
                    break

            if genome_id is None and ignore_prefix and key.startswith(("GCA_", "GCF_")):
                suffix = key[4:]
                for candidate in (suffix, self._strip_accession_version(suffix)):
                    genome_id = self._accession_suffix_index.get(candidate)
                    if genome_id is not None:
                        break

            if genome_id is None:
                continue

            matched_queries.add(query)
            if genome_id not in seen_genomes:
                seen_genomes.add(genome_id)
                genome_ids.append(genome_id)

        return genome_ids, matched_queries

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

    @staticmethod
    def _normalize_accession(accession: Optional[str]) -> str:
        """
        Normalize GTDB/NCBI accessions for stable comparisons.

        Examples:
            RS_GCF_000001405.40 -> GCF_000001405.40
            GB_GCA_000001405.28 -> GCA_000001405.28
        """
        if not accession:
            return ""
        value = accession.strip()
        if value.startswith("RS_") or value.startswith("GB_"):
            value = value[3:]
        return value

    def filter_by_genome_type(self, genome_ids: List[str], genome_type: str) -> List[str]:
        """Filter genome IDs by ncbi_genome_category. genome_type must be a key of GENOME_TYPE_CATEGORIES."""
        category = GENOME_TYPE_CATEGORIES.get(genome_type)
        if category is None:
            return genome_ids
        return [
            gid for gid in genome_ids
            if self.data.get(gid, {}).get("ncbi_genome_category", "none") == category
        ]

    def get_species_cluster_representative(self, genome_id: str) -> Optional[str]:
        """
        Return the species-cluster representative accession for a genome, if available.
        """
        row = self.data.get(genome_id)
        if not row:
            return None

        for key in (
            "gtdb_genome_representative",
            "species_cluster_representative",
            "gtdb_species_representative",
        ):
            rep = self._normalize_accession(row.get(key))
            if rep:
                return rep
        return None

    def is_species_cluster_representative(self, genome_id: str) -> bool:
        """
        Determine whether the genome is marked as a species-cluster representative.
        """
        row = self.data.get(genome_id)
        if not row:
            return False

        # Prefer explicit boolean markers when available.
        for key in (
            "gtdb_representative",
            "species_representative",
            "is_representative",
        ):
            value = row.get(key, "")
            if value is None:
                continue
            if str(value).strip().lower() in {"t", "true", "1", "yes", "y"}:
                return True

        accession = self._normalize_accession(row.get("accession") or row.get("Genome"))
        representative = self.get_species_cluster_representative(genome_id)
        if accession and representative and accession == representative:
            return True

        return False
