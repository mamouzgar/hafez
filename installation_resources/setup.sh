#!/bin/bash

# --- Step 1: Install system dependencies ---
brew install gcc pkg-config icu4c udunits abseil cmake cmake-docs

# --- Step 2: Export required path variables ---
export PKG_CONFIG_PATH="/opt/homebrew/opt/icu4c/lib/pkgconfig"

# --- Step 3: Create conda environment ---
conda env create --no-default-packages -f environment.yml

# --- Step 4: Activate environment (for interactive sessions) ---
echo "To activate the environment, run: conda activate hafez_env"
