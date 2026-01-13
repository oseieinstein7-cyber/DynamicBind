# DynamicBind Research Documentation

## Purpose

This folder contains **complete documentation** for the DynamicBind project. Everything needed to write a full research paper or create ISEF presentation boards is organized here.

**Goal**: A person who was not involved in the project must be able to write a complete research paper from this folder alone.

---

## Folder Structure

```
ResearchPaper/
├── figures/                    # All visualizations
│   ├── phase1_md/             # MD simulation figures
│   ├── phase2_pockets/        # Pocket detection visuals
│   ├── phase3_diffusion/      # Model architecture & results
│   ├── phase4_validation/     # Docking & ADMET results
│   ├── overview/              # Project overview diagrams
│   └── supplementary/         # Additional figures
│
├── data_tables/               # All numerical data
│   ├── simulation_stats/      # MD statistics
│   ├── pocket_analysis/       # Pocket properties
│   ├── binding_scores/        # Docking results
│   ├── admet_results/         # Drug-likeness data
│   └── comparative/           # Cross-target comparisons
│
├── methods/                   # Detailed methodology
│   ├── phase1_md/             # MD protocols
│   ├── phase2_pockets/        # Pocket detection methods
│   ├── phase3_diffusion/      # Model architecture details
│   ├── phase4_validation/     # Validation protocols
│   └── computational_resources/# HPC details
│
├── results/                   # Results by target protein
│   ├── abeta42/               # Amyloid-Beta results
│   ├── tau_k18/               # Tau K18 results
│   ├── cmyc/                  # c-Myc results
│   ├── alpha_synuclein/       # Alpha-synuclein results
│   ├── p53/                   # p53 results
│   └── comparative/           # Cross-target analysis
│
├── protocols/                 # Step-by-step guides
│   ├── setup/                 # Environment setup
│   ├── simulation/            # Running simulations
│   ├── analysis/              # Data analysis
│   └── validation/            # Validation procedures
│
├── presentations/             # ISEF & presentation materials
│   ├── isef_board/            # ISEF tri-fold board content
│   ├── poster/                # Scientific poster
│   └── slides/                # Presentation slides
│
├── supplementary/             # Extended materials
│   ├── raw_data/              # Unprocessed data links
│   ├── extended_methods/      # Detailed technical notes
│   └── additional_figures/    # Extra visualizations
│
├── daily_logs/                # Progress tracking
│   └── 2025/                  # Logs by year
│
├── references/                # Bibliography & citations
└── drafts/                    # Paper drafts
```

---

## Documentation Standards

### Figures
- **Format**: PNG (300 DPI minimum), SVG for diagrams
- **Naming**: `[phase]_[target]_[description]_v[version].png`
- **Example**: `phase1_abeta42_rmsd_trajectory_v1.png`
- **Requirements**: Always include axis labels, legends, scale bars

### Data Tables
- **Format**: CSV (raw), Markdown (readable)
- **Include**: Units, sample sizes, error margins
- **Naming**: `[phase]_[target]_[metric].csv`

### Methods
- **Detail Level**: Reproducible by another researcher
- **Include**: Software versions, parameters, commands

### Daily Logs
- **Format**: Markdown
- **Include**: Date, tasks completed, observations, next steps
- **Template**: See `daily_logs/TEMPLATE.md`

---

## Key Metrics to Document

### Phase 1: MD Simulations
- [ ] Simulation length (ns)
- [ ] Number of frames saved
- [ ] RMSD over time
- [ ] Radius of gyration
- [ ] Secondary structure evolution
- [ ] GaMD boost statistics

### Phase 2: Pocket Detection
- [ ] Number of pockets identified
- [ ] Pocket volumes
- [ ] Druggability scores
- [ ] Pocket persistence across frames
- [ ] Hydrophobicity profiles

### Phase 3: Diffusion Model
- [ ] Training loss curves
- [ ] Validation metrics
- [ ] Generated molecule statistics
- [ ] Diversity metrics

### Phase 4: Validation
- [ ] Docking scores (kcal/mol)
- [ ] ADMET predictions
- [ ] Binding pose analysis
- [ ] Comparison to known drugs

---

## ISEF Board Layout

```
┌─────────────────────────────────────────────────────────────────────────┐
│                              TITLE BANNER                                │
│   DynamicBind: AI-Powered Drug Design for "Undruggable" Cancer Proteins │
└─────────────────────────────────────────────────────────────────────────┘

┌──────────────────┐ ┌──────────────────┐ ┌──────────────────┐
│                  │ │                  │ │                  │
│   INTRODUCTION   │ │    HYPOTHESIS    │ │    MATERIALS     │
│                  │ │                  │ │                  │
│  - IDP Problem   │ │                  │ │  - Target list   │
│  - Cancer link   │ │                  │ │  - Software      │
│  - Our approach  │ │                  │ │  - Polaris HPC   │
│                  │ │                  │ │                  │
├──────────────────┤ ├──────────────────┤ ├──────────────────┤
│                  │ │                  │ │                  │
│    PROCEDURE     │ │     RESULTS      │ │   CONCLUSIONS    │
│                  │ │                  │ │                  │
│  Phase 1: MD     │ │  - MD movies     │ │  - Key findings  │
│  Phase 2: Pockets│ │  - Pocket maps   │ │  - Drug leads    │
│  Phase 3: AI     │ │  - New molecules │ │  - Future work   │
│  Phase 4: Valid. │ │  - Binding data  │ │                  │
│                  │ │                  │ │                  │
└──────────────────┘ └──────────────────┘ └──────────────────┘
```

---

## Contact

For questions about this documentation, contact the project lead.
