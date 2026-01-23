#!/bin/bash
#
# DynamicBind - Master Overnight Structure Preparation Script
#
# Run this script on ALCF Polaris login node to prepare all 25 protein structures
# for GaMD simulation. Designed to run overnight (~8 hours) without hanging.
#
# Usage:
#   ./run_overnight.sh              # Normal run
#   ./run_overnight.sh --reset      # Reset and reprocess all
#   nohup ./run_overnight.sh &      # Run in background
#
# The script will:
# 1. Set up resource limits for login node
# 2. Activate the gamd_env conda environment
# 3. Download 25 protein structures from PDB, AlphaFold, and PED
# 4. Process structures (extract models/chains, clean, fix)
# 5. Save prepared structures for GaMD simulation
#

set -e  # Exit on error

# ============================================================================
# CONFIGURATION
# ============================================================================

# Detect script location
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
PROJECT_ROOT="$( cd "$SCRIPT_DIR/../.." && pwd )"

# Output directory
OUTPUT_DIR="${PROJECT_ROOT}/data/structures"

# Log file
LOG_FILE="${PROJECT_ROOT}/structure_prep_$(date +%Y%m%d_%H%M%S).log"

# PDBFixer timeout (seconds per structure)
PDBFIXER_TIMEOUT=120

# ============================================================================
# RESOURCE LIMITS FOR LOGIN NODE
# ============================================================================

# Prevent OpenBLAS/MKL from spawning too many threads
export OPENBLAS_NUM_THREADS=4
export OMP_NUM_THREADS=4
export MKL_NUM_THREADS=4
export NUMEXPR_MAX_THREADS=4
export VECLIB_MAXIMUM_THREADS=4

# Limit Python's memory usage if needed
# ulimit -v 16000000  # 16GB virtual memory limit

echo "=============================================="
echo "DynamicBind Structure Preparation"
echo "=============================================="
echo "Project root: $PROJECT_ROOT"
echo "Output dir: $OUTPUT_DIR"
echo "Log file: $LOG_FILE"
echo "Start time: $(date)"
echo ""

# ============================================================================
# CONDA ENVIRONMENT SETUP
# ============================================================================

# Try different conda activation methods (for different Polaris configs)
activate_conda() {
    # Method 1: Direct miniconda
    if [ -f "$HOME/miniconda3/etc/profile.d/conda.sh" ]; then
        source "$HOME/miniconda3/etc/profile.d/conda.sh"
        echo "Activated conda from miniconda3"
        return 0
    fi

    # Method 2: Anaconda3
    if [ -f "$HOME/anaconda3/etc/profile.d/conda.sh" ]; then
        source "$HOME/anaconda3/etc/profile.d/conda.sh"
        echo "Activated conda from anaconda3"
        return 0
    fi

    # Method 3: Module system
    if command -v module &> /dev/null; then
        module load conda/2023-10-04 2>/dev/null || true
        echo "Loaded conda module"
        return 0
    fi

    # Method 4: Already in path
    if command -v conda &> /dev/null; then
        echo "Using conda from PATH"
        return 0
    fi

    echo "ERROR: Could not find conda installation"
    return 1
}

echo "Setting up conda environment..."
activate_conda || exit 1

# Activate gamd_env
if conda activate gamd_env 2>/dev/null; then
    echo "Activated gamd_env"
elif conda activate base 2>/dev/null; then
    echo "WARNING: gamd_env not found, using base environment"
fi

# Verify Python
echo "Python: $(which python)"
python --version

# ============================================================================
# RUN STRUCTURE PREPARATION
# ============================================================================

echo ""
echo "Starting structure preparation..."
echo "This may take several hours for all 25 structures."
echo ""

# Parse arguments
RESET_FLAG=""
if [ "$1" == "--reset" ]; then
    RESET_FLAG="--reset"
    echo "Reset mode: will reprocess all structures"
fi

# Create output directory
mkdir -p "$OUTPUT_DIR"

# Run the preparation script
cd "$PROJECT_ROOT"

python -u src/phase1_md/simulation/prepare_structure.py \
    --target all \
    --output "$OUTPUT_DIR" \
    --timeout $PDBFIXER_TIMEOUT \
    $RESET_FLAG \
    2>&1 | tee "$LOG_FILE"

# ============================================================================
# SUMMARY
# ============================================================================

echo ""
echo "=============================================="
echo "STRUCTURE PREPARATION COMPLETE"
echo "=============================================="
echo "End time: $(date)"
echo ""

# Count prepared files
if [ -d "$OUTPUT_DIR/prepared" ]; then
    PREP_COUNT=$(find "$OUTPUT_DIR/prepared" -name "*.pdb" | wc -l)
    echo "Prepared structures: $PREP_COUNT"
    echo ""
    echo "Files in prepared directory:"
    ls -la "$OUTPUT_DIR/prepared/"
fi

echo ""
echo "Log saved to: $LOG_FILE"
echo ""

# Check for any failures
if [ -f "$OUTPUT_DIR/checkpoint.json" ]; then
    echo "Checkpoint status:"
    cat "$OUTPUT_DIR/checkpoint.json"
fi
