#!/bin/bash
#
# Master Structure Preparation Runner for Polaris
#
# Usage:
#   ./run_master.sh                    # Normal run
#   nohup ./run_master.sh > master.log 2>&1 &   # Background overnight
#

echo "=============================================="
echo "MASTER STRUCTURE PREPARATION"
echo "=============================================="
echo "Start: $(date)"
echo ""

# Resource limits for login node
export OPENBLAS_NUM_THREADS=2
export OMP_NUM_THREADS=2
export MKL_NUM_THREADS=2
export NUMEXPR_MAX_THREADS=2

# Activate conda
if [ -f "$HOME/miniconda3/etc/profile.d/conda.sh" ]; then
    source "$HOME/miniconda3/etc/profile.d/conda.sh"
elif [ -f "$HOME/anaconda3/etc/profile.d/conda.sh" ]; then
    source "$HOME/anaconda3/etc/profile.d/conda.sh"
fi

conda activate gamd_env 2>/dev/null || conda activate base

echo "Python: $(which python)"
echo "Conda env: $CONDA_DEFAULT_ENV"
echo ""

# Run the master script
cd ~/DynamicBind_GaMD
python ~/DynamicBind/scripts/preprocessing/master_prepare_structures.py

echo ""
echo "=============================================="
echo "COMPLETE: $(date)"
echo "=============================================="

# Show results
echo ""
echo "Prepared files:"
ls -la ~/DynamicBind_GaMD/prepared/*.pdb 2>/dev/null | head -30

echo ""
echo "Checkpoint:"
cat ~/DynamicBind_GaMD/checkpoint.json 2>/dev/null
