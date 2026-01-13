#!/bin/bash
#==============================================================================
# DynamicBind - Polaris Environment Setup
#
# This script sets up the conda environment on ALCF Polaris.
# Run this ONCE before submitting any jobs.
#
# Usage:
#   ./setup_environment.sh
#==============================================================================

set -e

echo "=================================================="
echo "DynamicBind - Polaris Environment Setup"
echo "=================================================="

# Configuration
ENV_NAME="dynamicbind_env"
ENV_PATH="/eagle/projects/dynamicbind/envs/${ENV_NAME}"
PROJECT_DIR="/eagle/projects/dynamicbind/DynamicBind"

# Load base conda
module use /soft/modulefiles
module load conda/2024-04-29

echo "Creating conda environment at: ${ENV_PATH}"

# Create environment from YAML if it exists
if [ -f "${PROJECT_DIR}/environments/dynamicbind_polaris.yml" ]; then
    conda env create -f ${PROJECT_DIR}/environments/dynamicbind_polaris.yml -p ${ENV_PATH}
else
    # Create from scratch
    conda create -p ${ENV_PATH} python=3.10 -y

    # Activate
    conda activate ${ENV_PATH}

    # Install core dependencies
    echo "Installing OpenMM..."
    conda install -c conda-forge openmm cudatoolkit=12.2 -y

    echo "Installing MDAnalysis..."
    conda install -c conda-forge mdanalysis -y

    echo "Installing PDBFixer..."
    conda install -c conda-forge pdbfixer -y

    echo "Installing BioPython..."
    conda install -c conda-forge biopython -y

    echo "Installing scientific Python stack..."
    conda install -c conda-forge numpy scipy pandas matplotlib scikit-learn -y

    echo "Installing PyYAML..."
    conda install -c conda-forge pyyaml -y

    echo "Installing requests..."
    conda install -c conda-forge requests -y

    # Install PyTorch with CUDA support (for Phase 2/3)
    echo "Installing PyTorch..."
    conda install pytorch pytorch-cuda=12.1 -c pytorch -c nvidia -y

    # Install PyTorch Geometric (for GNN pocket scorer)
    echo "Installing PyTorch Geometric..."
    conda install pyg -c pyg -y
fi

echo ""
echo "=================================================="
echo "Environment Setup Complete!"
echo "=================================================="
echo ""
echo "To activate:"
echo "  module load conda/2024-04-29"
echo "  conda activate ${ENV_PATH}"
echo ""
echo "To test OpenMM:"
echo "  python -c \"import openmm; print(openmm.__version__)\""
echo "  python -m openmm.testInstallation"
echo ""
