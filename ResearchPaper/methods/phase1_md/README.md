# Phase 1: Molecular Dynamics Simulation Methods

## Overview

Phase 1 generates conformational ensembles of Intrinsically Disordered Proteins (IDPs) using Gaussian Accelerated Molecular Dynamics (GaMD). This captures the dynamic nature of IDPs that is missed by static crystal structures.

## Why GaMD?

Standard MD simulations of IDPs require microsecond-to-millisecond timescales to observe relevant conformational changes. GaMD accelerates this by:

1. Adding a harmonic boost potential that smooths the energy landscape
2. Enabling faster transitions between conformational states
3. Allowing recovery of original free energy profiles through reweighting

## Software Stack

| Software | Version | Purpose |
|----------|---------|---------|
| OpenMM | ≥8.0 | MD engine with GPU acceleration |
| PDBFixer | latest | Structure preparation |
| MDAnalysis | ≥2.4 | Trajectory analysis |
| BioPython | ≥1.80 | PDB file handling |

## Simulation Protocol

### 1. Structure Preparation

```
Input: PDB file from RCSB or PED database
       ↓
Step 1: Download raw structure
       ↓
Step 2: Extract target chain
       ↓
Step 3: Fix missing atoms/residues (PDBFixer)
       ↓
Step 4: Add hydrogens at pH 7.0
       ↓
Output: Cleaned PDB ready for simulation
```

### 2. System Setup

| Parameter | Value | Rationale |
|-----------|-------|-----------|
| Force Field | AMBER14 | Well-validated for proteins |
| Water Model | TIP3P | Standard explicit water |
| Box Padding | 1.2-1.5 nm | Sufficient for IDP extension |
| Ionic Strength | 0.15 M NaCl | Physiological conditions |
| Temperature | 310 K | Physiological |
| Pressure | 1 bar | Standard conditions |

### 3. Energy Minimization

- Algorithm: L-BFGS
- Max iterations: 5000
- Convergence: When energy change < 1 kJ/mol

### 4. Equilibration

| Phase | Ensemble | Duration | Purpose |
|-------|----------|----------|---------|
| Phase 1 | NVT | 100 ps | Temperature equilibration |
| Phase 2 | NPT | 1 ns | Pressure/density equilibration |

### 5. GaMD Equilibration

| Step | Duration | Purpose |
|------|----------|---------|
| Statistics collection | 5 ns | Calculate Vmax, Vmin, Vavg, σ |
| Boost parameter calculation | - | Determine k0 and E |
| GaMD equilibration | 5 ns | Verify stable boosted dynamics |

### 6. Production GaMD

| Target | Duration | Frames | Size (est.) |
|--------|----------|--------|-------------|
| c-Myc | 2 μs | 200,000 | ~50 GB |
| p53 | 1.5 μs | 150,000 | ~40 GB |
| α-Synuclein | 1 μs | 100,000 | ~25 GB |
| Aβ42 | 1 μs | 100,000 | ~15 GB |
| Tau K18 | 1 μs | 100,000 | ~30 GB |

## GaMD Parameters

### Boost Potential

The GaMD boost potential is:

```
ΔV(r) = ½k(E - V(r))²    when V(r) < E
ΔV(r) = 0                 when V(r) ≥ E
```

Where:
- V(r) = original potential energy
- E = threshold energy (set to Vmax)
- k = force constant

### Dual-Boost Parameters

| Parameter | Description | Typical Value |
|-----------|-------------|---------------|
| σ0_total | Upper limit of total boost SD | 6.0 kcal/mol |
| σ0_dihedral | Upper limit of dihedral boost SD | 6.0 kcal/mol |

## Output Files

```
output/
├── solvated.pdb          # Solvated system
├── minimized.pdb         # After energy minimization
├── equilibrated.pdb      # After equilibration
├── production_traj.dcd   # Production trajectory
├── production_log.csv    # Energy/temperature log
├── gamd_parameters.txt   # GaMD boost parameters
├── final.pdb             # Final structure
└── analysis/
    ├── rmsd.csv          # RMSD vs time
    ├── rmsf.csv          # Per-residue RMSF
    ├── rg.csv            # Radius of gyration
    ├── clusters.csv      # Cluster assignments
    └── cluster_*.pdb     # Representative structures
```

## Quality Control Metrics

### Convergence Checks
- [ ] RMSD plateau reached
- [ ] Rg stabilized
- [ ] GaMD boost statistics converged

### Expected Values for IDPs
| Metric | Expected Range | Concern if |
|--------|----------------|------------|
| RMSD | 5-15 Å | < 3 Å (too rigid) |
| Rg | 15-40 Å | Constant (not sampling) |
| Boost | < 10 kcal/mol | > 20 (over-boosted) |

## Computational Resources

### ALCF Polaris Specifications
- 4x NVIDIA A100 (40GB) per node
- ~1 μs/day throughput for typical IDP system
- Checkpoint every 100 ps for fault tolerance

### Resource Allocation (5000 node-hours)
| Target | Priority | Hours | Rationale |
|--------|----------|-------|-----------|
| c-Myc | Critical | 400 | 70% of cancers |
| p53 | High | 350 | Well-characterized |
| α-Syn | High | 250 | Major unmet need |
| Aβ42 | Medium | 200 | Complex aggregation |
| Tau K18 | Medium | 200 | IDP challenge |

## References

1. Miao Y, et al. (2015) J. Chem. Theory Comput. 11, 3584-3595.
   "Gaussian Accelerated Molecular Dynamics: Unconstrained Enhanced Sampling"

2. Pang YT, et al. (2017) J. Chem. Theory Comput. 13, 9-19.
   "Gaussian Accelerated Molecular Dynamics in NAMD"

3. Eastman P, et al. (2017) PLoS Comput. Biol. 13, e1005659.
   "OpenMM 7: Rapid development of high performance algorithms"
