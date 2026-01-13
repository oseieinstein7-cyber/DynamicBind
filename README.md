# DynamicBind

## An SE(3)-Equivariant Diffusion Model for De Novo Drug Design Targeting Intrinsically Disordered Proteins (IDPs)

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

---

## Executive Summary

Roughly 30-50% of eukaryotic proteins are **Intrinsically Disordered Proteins (IDPs)**. Unlike normal proteins that have a fixed "lock" shape for drug "keys," IDPs are shapeshifters—they constantly change structure. Because of this, they are considered **"undruggable,"** yet they drive major diseases:

- **Cancer**: c-Myc (~70% of cancers), p53
- **Alzheimer's**: Tau, Amyloid-Beta (Aβ42)
- **Parkinson's**: Alpha-synuclein

### Our Approach

Instead of docking drugs to static protein snapshots, DynamicBind:

1. **Simulates** the moving "video" of proteins via Gaussian Accelerated Molecular Dynamics (GaMD)
2. **Identifies** fleeting cryptic pockets that open for milliseconds
3. **Generates** novel drug-like molecules using SE(3)-Equivariant Diffusion Models
4. **Validates** candidates through computational docking and ADMET analysis

---

## Project Structure

```
DynamicBind/
├── data/                          # All input/output data
│   ├── structures/                # PDB files by target
│   │   ├── abeta42/
│   │   ├── tau_k18/
│   │   ├── cmyc/
│   │   ├── alpha_synuclein/
│   │   └── p53/
│   ├── trajectories/              # MD simulation outputs
│   ├── pockets/                   # Detected cryptic pockets
│   └── ligands/                   # Generated/validated ligands
│
├── src/                           # Source code
│   ├── phase1_md/                 # Molecular Dynamics simulations
│   ├── phase2_pockets/            # Cryptic pocket detection
│   ├── phase3_diffusion/          # SE(3) diffusion model
│   ├── phase4_validation/         # Docking & ADMET validation
│   └── utils/                     # Shared utilities
│
├── configs/                       # Configuration files
│   ├── targets/                   # Target protein configs
│   ├── simulation/                # MD simulation parameters
│   ├── model/                     # Model hyperparameters
│   └── hpc/                       # HPC job configurations
│
├── scripts/                       # Execution scripts
│   ├── polaris/                   # ALCF Polaris submission scripts
│   ├── local/                     # Local execution scripts
│   └── preprocessing/             # Data preprocessing
│
├── models/                        # Trained models
├── notebooks/                     # Jupyter analysis notebooks
├── tests/                         # Unit and integration tests
├── environments/                  # Conda environment files
└── logs/                          # Execution logs
```

---

## Target Proteins & PDB Sources

### 1. Amyloid-Beta 42 (Alzheimer's)
| PDB ID | Description | Rationale |
|--------|-------------|-----------|
| 2NAO | NMR fibril structure (S-bend) | Disease-relevant fibril architecture |
| 1Z0Q | Solution NMR, aqueous | Collapsed but soluble starting point |
| 1IYT | Solid-state NMR fibril | Fibril-competent monomer extraction |

### 2. Tau K18 (Alzheimer's)
| Source | Description | Rationale |
|--------|-------------|-----------|
| PED00192 | IDP ensemble | Conformational diversity |
| PED00443 | IDP ensemble | Additional states |
| 2MZ7 | NMR structure | Experimental reference |

### 3. c-Myc (Cancer - ~70% of cancers)
| PDB ID | Description | Rationale |
|--------|-------------|-----------|
| 6G6K | Apo c-Myc:Max complex | "Ready-to-bind" helical conformation |
| 1NKP | DNA-bound complex | Fully folded alpha-helical state |
| 2A93 | NMR zipper ensemble | Conformational flexibility |

### 4. Alpha-Synuclein (Parkinson's)
| PDB ID | Description | Rationale |
|--------|-------------|-----------|
| 8A9L | Cryo-EM PD/DLB brain fibrils | Actual Lewy pathology fold |
| 8FPT | Solid-state NMR LBD fibrils | Complementary polymorph |
| 2KKW | Micelle-bound, partially folded | Pre-fibrillar cryptic pocket state |

### 5. p53 (Cancer)
| PDB ID | Description | Rationale |
|--------|-------------|-----------|
| 2AHI | WT crystal structure (1.85 Å) | Clean starting point for mutations |
| 2FEJ | NMR ensemble | Natural pocket breathing dynamics |

---

## Methodology Overview

### Phase 1: Molecular Dynamics Simulation
- **Tool**: OpenMM + GaMD (Gaussian Accelerated MD)
- **Goal**: Generate conformational ensembles capturing IDP dynamics
- **Output**: Trajectory files showing protein shape changes over time

### Phase 2: Cryptic Pocket Detection
- **Tools**: CryptoSite, LIGSITE, MDAnalysis, Custom GNN scorer
- **Innovation**: Fast "frame triage" scorer ranking frames by druggability
- **Output**: 3D coordinates of transient binding pockets

### Phase 3: SE(3)-Equivariant Diffusion Model
- **Base**: DiffDock/TargetDiff architecture
- **Innovation**: Modified to accept ensemble of protein shapes
- **Output**: Novel drug-like molecules targeting dynamic pockets

### Phase 4: Validation
- **Docking**: AutoDock Vina / Glide
- **ADMET**: SwissADME for toxicity/absorption analysis
- **FEP**: Free Energy Perturbation for binding affinity

---

## Computational Resources

- **Platform**: ALCF Polaris Supercomputer
- **Allocation**: 5,000 node-hours
- **Hardware**: 4x NVIDIA A100 GPUs per node

---

## Installation

```bash
# Clone repository
git clone https://github.com/yourusername/DynamicBind.git
cd DynamicBind

# Create conda environment
conda env create -f environments/dynamicbind.yml
conda activate dynamicbind

# For Polaris-specific environment
conda env create -f environments/dynamicbind_polaris.yml
```

---

## Quick Start

```bash
# 1. Download and prepare structures
python scripts/preprocessing/download_structures.py

# 2. Run Phase 1 MD simulation (local test)
python src/phase1_md/run_simulation.py --config configs/targets/abeta42.yaml --test

# 3. Submit to Polaris
qsub scripts/polaris/submit_gamd.pbs
```

---

## Citation

If you use DynamicBind in your research, please cite:

```bibtex
@software{dynamicbind2025,
  title={DynamicBind: SE(3)-Equivariant Diffusion for IDP Drug Design},
  author={Your Name},
  year={2025},
  url={https://github.com/yourusername/DynamicBind}
}
```

---

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

---

## Acknowledgments

- ALCF for Polaris compute allocation
- OpenMM and GROMACS development teams
- DiffDock and TargetDiff authors for foundational architectures
