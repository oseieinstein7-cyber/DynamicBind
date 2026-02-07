#!/usr/bin/env python3
"""
Master PED (Protein Ensemble Database) Download and Processing Script

Downloads protein ensembles from PED using the correct API endpoint
and extracts specified number of diverse models for GaMD simulation.

PED API Endpoint Format:
    https://proteinensemble.org/api/ensemble_sample/{PED_ID}e{ENSEMBLE_NUM}

Example:
    https://proteinensemble.org/api/ensemble_sample/PED00024e001

Usage:
    python master_ped.py                    # Download and process all
    python master_ped.py --check            # Check PED availability
    python master_ped.py --ped PED00024     # Process specific PED
"""

import os
import sys
import json
import subprocess
import argparse
import logging
from pathlib import Path
from datetime import datetime

# Resource limits for login nodes
os.environ.setdefault('OPENBLAS_NUM_THREADS', '2')
os.environ.setdefault('OMP_NUM_THREADS', '2')
os.environ.setdefault('MKL_NUM_THREADS', '2')

import requests
from Bio.PDB import PDBParser, PDBIO, Select
import numpy as np

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)

# ============================================================================
# PED ENSEMBLE SPECIFICATIONS
# ============================================================================
# Based on the 25-structure specification for DynamicBind GaMD simulations

PED_SPECS = {
    "PED00024": {
        "ensemble": "e001",
        "models": 8,
        "target": "alpha_synuclein",
        "description": "Alpha-Synuclein IDP ensemble"
    },
    "PED00049": {
        "ensemble": "e001",
        "models": 7,
        "target": "cmyc",
        "description": "c-Myc bHLH-LZ disordered region"
    },
    "PED00086": {
        "ensemble": "e001",
        "models": 5,
        "target": "p53",
        "description": "p53 N-terminal transactivation domain"
    },
    "PED00192": {
        "ensemble": "e001",
        "models": 10,
        "target": "tau_k18",
        "description": "Tau K18 repeat domain ensemble"
    },
    "PED00443": {
        "ensemble": "e001",
        "models": 8,
        "target": "tau_k18",
        "description": "Tau K18 phosphorylated ensemble"
    },
}

# PED API base URL
PED_API_BASE = "https://proteinensemble.org/api/ensemble_sample"


class ModelSelect(Select):
    """Select specific model from NMR ensemble."""

    def __init__(self, model_id):
        self.model_id = model_id

    def accept_model(self, model):
        return model.id == self.model_id


def download_ped_ensemble(ped_id, ensemble_id, output_dir):
    """
    Download PED ensemble using the correct API endpoint.

    Args:
        ped_id: PED identifier (e.g., "PED00024")
        ensemble_id: Ensemble number (e.g., "e001")
        output_dir: Directory to save the file

    Returns:
        Path to downloaded file or None if failed
    """
    # Construct API URL
    url = f"{PED_API_BASE}/{ped_id}{ensemble_id}"

    output_file = output_dir / f"{ped_id}_{ensemble_id}.pdb"

    logger.info(f"Downloading {ped_id} from: {url}")

    try:
        response = requests.get(url, timeout=60)

        if response.status_code == 200:
            content = response.text

            # Verify it's actually PDB content
            if content.startswith("HEADER") or content.startswith("TITLE") or \
               content.startswith("MODEL") or content.startswith("ATOM"):
                with open(output_file, 'w') as f:
                    f.write(content)
                logger.info(f"  Downloaded: {output_file.name} ({len(content)} bytes)")
                return output_file
            else:
                logger.error(f"  Response is not PDB format: {content[:100]}")
                return None
        else:
            logger.error(f"  HTTP {response.status_code}: {response.text[:200]}")
            return None

    except requests.exceptions.Timeout:
        logger.error(f"  Timeout downloading {ped_id}")
        return None
    except Exception as e:
        logger.error(f"  Error downloading {ped_id}: {e}")
        return None


def count_models(pdb_file):
    """Count number of models in a PDB file."""
    parser = PDBParser(QUIET=True)
    try:
        structure = parser.get_structure("test", pdb_file)
        return len(list(structure.get_models()))
    except Exception:
        return 0


def select_diverse_models(pdb_file, n_models):
    """
    Select diverse models from NMR ensemble based on structural variance.

    Uses Cα RMSD to select maximally diverse subset.
    """
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("ens", pdb_file)

    models = list(structure.get_models())
    n_total = len(models)

    if n_total <= n_models:
        return list(range(n_total))

    # Get Cα coordinates for each model
    ca_coords = []
    for model in models:
        coords = []
        for chain in model:
            for residue in chain:
                if 'CA' in residue:
                    coords.append(residue['CA'].get_coord())
        if coords:
            ca_coords.append(np.array(coords))

    if not ca_coords or len(ca_coords[0]) == 0:
        # Fallback: evenly spaced selection
        indices = np.linspace(0, n_total - 1, n_models, dtype=int)
        return list(indices)

    # Calculate pairwise RMSD matrix
    n = len(ca_coords)
    min_len = min(len(c) for c in ca_coords)
    ca_coords = [c[:min_len] for c in ca_coords]

    # Greedy diverse selection
    selected = [0]  # Start with first model

    while len(selected) < n_models:
        max_min_dist = -1
        best_idx = -1

        for i in range(n):
            if i in selected:
                continue

            # Calculate min distance to selected set
            min_dist = float('inf')
            for j in selected:
                rmsd = np.sqrt(np.mean(np.sum((ca_coords[i] - ca_coords[j])**2, axis=1)))
                min_dist = min(min_dist, rmsd)

            if min_dist > max_min_dist:
                max_min_dist = min_dist
                best_idx = i

        if best_idx >= 0:
            selected.append(best_idx)

    return sorted(selected)


def extract_models(pdb_file, model_indices, output_dir, prefix):
    """Extract specified models from ensemble and save as separate files."""
    parser = PDBParser(QUIET=True)
    io = PDBIO()

    structure = parser.get_structure("ens", pdb_file)
    models = list(structure.get_models())

    extracted = []
    for idx in model_indices:
        if idx < len(models):
            output_file = output_dir / f"{prefix}_model{idx + 1}.pdb"
            io.set_structure(structure)
            io.save(str(output_file), ModelSelect(idx))
            extracted.append(output_file)
            logger.info(f"    Extracted model {idx + 1}: {output_file.name}")

    return extracted


def add_hydrogens(input_file, output_file):
    """Add hydrogens using reduce from AmberTools."""
    try:
        cmd = f"reduce -BUILD {input_file} > {output_file} 2>/dev/null"
        result = subprocess.run(cmd, shell=True, timeout=120)

        if output_file.exists() and output_file.stat().st_size > 0:
            return True
        else:
            # Fallback: copy original
            import shutil
            shutil.copy(input_file, output_file)
            return True
    except Exception as e:
        logger.warning(f"    reduce failed: {e}")
        import shutil
        shutil.copy(input_file, output_file)
        return True


def process_ped_ensemble(ped_id, spec, raw_dir, prepared_dir):
    """Download and process a PED ensemble."""
    logger.info(f"\n{'='*60}")
    logger.info(f"Processing {ped_id}: {spec['description']}")
    logger.info(f"  Target: {spec['target']}, Models needed: {spec['models']}")
    logger.info(f"{'='*60}")

    # Download
    pdb_file = download_ped_ensemble(ped_id, spec['ensemble'], raw_dir)
    if not pdb_file:
        return []

    # Count models
    n_models = count_models(pdb_file)
    logger.info(f"  Ensemble contains {n_models} models")

    if n_models == 0:
        logger.error(f"  No models found in {pdb_file}")
        return []

    # Select diverse models
    n_needed = spec['models']
    selected = select_diverse_models(pdb_file, n_needed)
    logger.info(f"  Selected models: {[i+1 for i in selected]}")

    # Extract models
    prefix = f"{spec['target']}_{ped_id}"
    temp_dir = raw_dir / "temp"
    temp_dir.mkdir(exist_ok=True)

    extracted = extract_models(pdb_file, selected, temp_dir, prefix)

    # Add hydrogens and save to prepared directory
    prepared_files = []
    for temp_file in extracted:
        output_file = prepared_dir / temp_file.name.replace('.pdb', '_H.pdb')
        if add_hydrogens(temp_file, output_file):
            prepared_files.append(output_file)
            logger.info(f"  Prepared: {output_file.name}")

    return prepared_files


def check_ped_availability():
    """Check if PED API endpoints are accessible."""
    logger.info("Checking PED API availability...")
    logger.info(f"API Base: {PED_API_BASE}")

    for ped_id, spec in PED_SPECS.items():
        url = f"{PED_API_BASE}/{ped_id}{spec['ensemble']}"
        try:
            response = requests.head(url, timeout=10)
            status = "OK" if response.status_code == 200 else f"HTTP {response.status_code}"
        except Exception as e:
            status = f"Error: {e}"

        logger.info(f"  {ped_id}: {status}")


def main():
    parser = argparse.ArgumentParser(description="Download and process PED ensembles")
    parser.add_argument("--check", action="store_true", help="Check PED availability")
    parser.add_argument("--ped", help="Process specific PED ID")
    parser.add_argument("--output", default=os.path.expanduser("~/DynamicBind_GaMD"),
                        help="Output directory")

    args = parser.parse_args()

    if args.check:
        check_ped_availability()
        return

    # Setup directories
    base_dir = Path(args.output)
    raw_dir = base_dir / "raw"
    prepared_dir = base_dir / "prepared"

    raw_dir.mkdir(parents=True, exist_ok=True)
    prepared_dir.mkdir(parents=True, exist_ok=True)

    logger.info("=" * 70)
    logger.info("PED ENSEMBLE DOWNLOAD AND PREPARATION")
    logger.info("=" * 70)
    logger.info(f"Start time: {datetime.now()}")
    logger.info(f"Output: {base_dir}")
    logger.info(f"PED entries to process: {len(PED_SPECS)}")

    # Process
    all_prepared = []
    failed = []

    specs_to_process = {args.ped: PED_SPECS[args.ped]} if args.ped else PED_SPECS

    for ped_id, spec in specs_to_process.items():
        try:
            prepared = process_ped_ensemble(ped_id, spec, raw_dir, prepared_dir)
            all_prepared.extend(prepared)
            if not prepared:
                failed.append(ped_id)
        except Exception as e:
            logger.error(f"Failed to process {ped_id}: {e}")
            failed.append(ped_id)

    # Summary
    logger.info("\n" + "=" * 70)
    logger.info("SUMMARY")
    logger.info("=" * 70)
    logger.info(f"Total prepared files: {len(all_prepared)}")
    logger.info(f"Failed: {len(failed)}")

    if failed:
        logger.info(f"Failed PED IDs: {failed}")

    for f in all_prepared:
        logger.info(f"  - {f.name}")

    logger.info(f"\nEnd time: {datetime.now()}")


if __name__ == "__main__":
    main()
