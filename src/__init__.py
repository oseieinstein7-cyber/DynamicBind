"""
DynamicBind: SE(3)-Equivariant Diffusion Model for IDP Drug Design

A computational pipeline for de novo drug design targeting
Intrinsically Disordered Proteins (IDPs).

Phases:
1. phase1_md: Molecular Dynamics simulations (GaMD)
2. phase2_pockets: Cryptic pocket detection
3. phase3_diffusion: SE(3)-equivariant diffusion model
4. phase4_validation: Docking and ADMET validation
"""

__version__ = "0.1.0"
__author__ = "DynamicBind Team"
