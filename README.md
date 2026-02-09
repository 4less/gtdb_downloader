# GTDB Downloader

A Python script to download GTDB (Genome Taxonomy Database) genomes by taxonomy. Downloads genomes to a central directory with automatic symlink structure based on taxonomy classification.

## Disclaimer
This script was developed with the help of ChatGPT.

## Features

- Download GTDB metadata and genomes for multiple versions (r207, r214, r220, r226)
- Search for genomes by taxon name
- Automatic metadata download on first use
- Multiple mirror support (Europe, Asia-Pacific)
- Download genomes using aria2 (with wget fallback)
- Store genomes in a single folder with taxonomy-based symlink structure
- Custom base directory support via environment variable

## Installation

### Prerequisites

Install aria2 (recommended for faster downloads) and/or wget:

**Ubuntu/Debian:**
```bash
sudo apt-get install aria2 wget
```

**macOS:**
```bash
brew install aria2 wget
```

**Fedora/RHEL:**
```bash
sudo dnf install aria2 wget
```

### Install gtdb-dl

Install in development/editable mode so it's added to your PATH and updates reflect immediately:

```bash
cd /path/to/gtdb_downloader
pip install -e .
```

This will:
1. Install the `gtdb-dl` command to your PATH
2. Create a base directory at `~/.gtdb_downloader` for storing GTDB data
3. Allow for easy updates (reinstall with `pip install -e .` after pulling changes)

### Verify Installation

```bash
gtdb-dl --help
```

## Usage

### Download metadata only

```bash
# Download metadata for a specific version (default: Europe mirror)
gtdb-dl --gtdb r226 --download

# Download ar53 (archaeal) dataset
gtdb-dl --gtdb r226 --download --dataset ar53

# Use different mirror
gtdb-dl --gtdb r226 --download --mirror asia-pacific1
```

### Download all genomes for a taxon

```bash
# Search by simple taxon name
gtdb-dl --gtdb r226 --taxon "Firmicutes"

# Full GTDB taxonomy path
gtdb-dl --gtdb r226 --taxon "d__Bacteria;p__Firmicutes;c__Clostridia"

# With verbose output
gtdb-dl --gtdb r226 --taxon "Firmicutes" -v

# Dry run (show what would be downloaded)
gtdb-dl --gtdb r226 --taxon "Firmicutes" --dry-run

# Use different mirror
gtdb-dl --gtdb r226 --taxon "Firmicutes" --mirror asia-pacific2

# Create symlink structure in custom directory
gtdb-dl --gtdb r226 --taxon "Firmicutes" --output /path/to/my_project/genomes
```

### Custom base directory

By default, data is stored in `~/.gtdb_downloader`. Change this with an environment variable:

```bash
# Use custom directory for this session
GTDBDL_DATA=/mnt/large_storage/gtdb gtdb-dl --gtdb r226 --taxon "Firmicutes"

# Or set permanently in your shell config (.bashrc, .zshrc, etc.)
export GTDBDL_DATA=/mnt/large_storage/gtdb
gtdb-dl --gtdb r226 --taxon "Firmicutes"
```

### Mirror options

Available mirrors:
- **europe** (default): https://data.gtdb.aau.ecogenomic.org/releases/
- **asia-pacific1**: https://data.gtdb.ecogenomic.org/releases/
- **asia-pacific2**: https://data.ace.uq.edu.au/public/gtdb/data/releases/

### Output structure

When you download genomes for a taxon, genomes are stored centrally and symlinks are organized by taxonomy:

**Data directory (central storage):**
```
~/.gtdb_downloader/r226/genomes/raw/
├── GCA_023390935.1_ASM2339093v1_genomic.fna.gz
├── GCA_946488845.1_ERR4303201_bin.25_MetaWRAP_v1.1_MAG_genomic.fna.gz
└── ... (all genomes from all projects)
```

**Output directory (taxonomy structure with symlinks):**

Default (if `--output` not specified):
```
~/.gtdb_downloader/r226/genomes/taxonomy/
└── Bacteria/
    └── Bacillota/
        └── Clostridia/
            └── Clostridiales/
                └── Clostridiaceae/
                    └── Clostridium/
                        └── Clostridium cuniculi/
                            └── GCA_023390935.1_ASM2339093v1_genomic.fna.gz -> ~/.gtdb_downloader/r226/genomes/raw/...
```

With `--output /path/to/my_project`:
```
/path/to/my_project/
└── Bacteria/
    └── Bacillota/
        └── Clostridia/
            └── ... (same hierarchy)
                └── Clostridium cuniculi/
                    └── GCA_023390935.1_ASM2339093v1_genomic.fna.gz -> ~/.gtdb_downloader/r226/genomes/raw/...
```

**Benefits:**
- **Efficient storage**: Genomes stored once, shared across all projects
- **Multiple projects**: Each project can have its own taxonomy-organized folder
- **Easy browsing**: Navigate through taxonomy hierarchy without file duplication
- **Metadata file**: `~/.gtdb_downloader/r226/bac120_metadata_r226.tsv.gz`

## Command Reference

```
usage: gtdb-dl [-h] --gtdb {r207,r214,r220,r226} [--taxon TAXON]
               [--dataset {bac120,ar53}] [--mirror {europe,asia-pacific1,asia-pacific2}]
               [--output OUTPUT] [--download] [--verbose] [--dry-run]
               [--base-dir BASE_DIR]

optional arguments:
  -h, --help                              show this help message and exit
  --gtdb {r207,r214,r220,r226}            GTDB version to use (required)
  --taxon TAXON                           Taxon to search for
  --dataset {bac120,ar53}                 Dataset type (default: bac120)
  --mirror {europe,asia-pacific1,asia-pacific2}  Mirror to download from (default: europe)
  --output OUTPUT, -o OUTPUT              Output directory for symlink taxonomy structure
  --download                              Only download metadata, don't download genomes
  --verbose, -v                          Verbose output
  --dry-run                               Show what would be downloaded
  --base-dir BASE_DIR                     Custom base directory (or use GTDBDL_DATA env var)
```

## Troubleshooting

### "Neither aria2c nor wget found in PATH"

Install one of these download tools:

```bash
# aria2 (recommended)
sudo apt-get install aria2

# or wget
sudo apt-get install wget
```

### "gtdb-dl: command not found"

Make sure the package is installed in editable mode:

```bash
cd /path/to/gtdb_downloader
pip install -e .
```

Check that the command is in your PATH:

```bash
which gtdb-dl
```

### Large downloads timing out

Use aria2 for better parallel downloads. It's automatically used if available, or you can install it:

```bash
sudo apt-get install aria2
```

## Development

To work on the code:

```bash
cd /path/to/gtdb_downloader
pip install -e .
```

Make changes and test:

```bash
gtdb-dl --help
```

Changes take effect immediately since it's installed in editable mode.

## License

MIT

## References

- GTDB Website: https://gtdb.ecogenomic.org/
- GTDB Downloads: https://data.gtdb.ecogenomic.org/releases/
