#!/bin/bash
# Regenerate every figure, then independently verify figure CSVs against the raw data.
set -e
for script in figure1/*.py figure2/*.py figure3/*.py figure4/[0-9]*.py figure5/[0-9]*.py figure6/*.py figure7/*.py supplementary/*.py statistics/permanova.py; do
    echo "Running script: $script"
    (cd "$(dirname "$script")" && python "$(basename "$script")" > /dev/null)
done
echo "Verifying figures against raw data:"
python verify_figures.py
