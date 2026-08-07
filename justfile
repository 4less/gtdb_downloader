# gtdb_downloader tasks -- run `just` to list recipes

python := env_var_or_default("PYTHON", "python3")

# List available recipes
default:
    @just --list

# Install gtdb-dl locally in editable mode
install:
    {{python}} -m pip install -e .
    @echo "Installed: $(command -v gtdb-dl)"

# Remove the local installation
uninstall:
    {{python}} -m pip uninstall -y gtdb-downloader

# Show the installed command's help
check:
    gtdb-dl --help
