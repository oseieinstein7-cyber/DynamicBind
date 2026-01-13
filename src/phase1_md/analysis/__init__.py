"""
Analysis submodule for Phase 1 MD trajectories.

Contains:
- Trajectory analysis (RMSD, RMSF, Rg)
- Secondary structure analysis
- Conformational clustering
- GaMD reweighting
"""

from .trajectory_analysis import TrajectoryAnalyzer, analyze_trajectory

__all__ = ["TrajectoryAnalyzer", "analyze_trajectory"]
