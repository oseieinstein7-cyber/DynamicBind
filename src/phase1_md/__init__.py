"""
Phase 1: Molecular Dynamics Simulation Module

This module provides tools for running Gaussian Accelerated Molecular Dynamics (GaMD)
simulations on Intrinsically Disordered Proteins (IDPs) to capture their conformational
dynamics and identify cryptic binding pockets.

Key Components:
- simulation/: Core simulation engines (OpenMM + GaMD)
- analysis/: Trajectory analysis tools
- utils/: Helper functions for structure preparation

Supported Targets:
- Amyloid-Beta 42 (Alzheimer's)
- Tau K18 (Alzheimer's)
- c-Myc (Cancer)
- Alpha-Synuclein (Parkinson's)
- p53 (Cancer)
"""

from .simulation import run_gamd
from .analysis import analyze_trajectory

__version__ = "0.1.0"
__all__ = ["run_gamd", "analyze_trajectory"]
