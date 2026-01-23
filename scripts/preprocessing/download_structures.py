#!/usr/bin/env python3
"""
DynamicBind - Structure Download and Preparation Script

Downloads and prepares all 25 protein structures for GaMD simulation:
- 5 Amyloid-Beta 42 structures (Alzheimer's)
- 4 Tau K18 structures (Alzheimer's)
- 4 c-Myc structures (Cancer)
- 5 Alpha-Synuclein structures (Parkinson's)
- 7 p53 structures including mutants (Cancer)

Usage:
    python download_structures.py                    # Prepare all targets
    python download_structures.py --target cmyc     # Prepare specific target
    python download_structures.py --target p53 --pdb 2AHI  # Specific PDB
    python download_structures.py --reset           # Reset and reprocess all
"""

import os
import sys
import argparse
import logging
from pathlib import Path

# Set resource limits BEFORE importing heavy libraries
os.environ.setdefault('OPENBLAS_NUM_THREADS', '4')
os.environ.setdefault('OMP_NUM_THREADS', '4')
os.environ.setdefault('MKL_NUM_THREADS', '4')

# Add src to path
PROJECT_ROOT = Path(__file__).parent.parent.parent
sys.path.insert(0, str(PROJECT_ROOT / "src"))

from phase1_md.simulation.prepare_structure import (
    StructurePreparer,
    TARGETS,
    STRUCTURE_SPECS
)

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    handlers=[
        logging.StreamHandler(),
        logging.FileHandler(PROJECT_ROOT / 'structure_prep.log')
    ]
)
logger = logging.getLogger(__name__)


def print_structure_summary():
    """Print summary of all structures to be prepared."""
    print("\n" + "="*70)
    print("STRUCTURE SUMMARY")
    print("="*70)

    for target_name, target in TARGETS.items():
        print(f"\n{target.name} ({target.disease}):")
        print(f"  {target.description}")
        print("  Structures:")

        for spec_id, spec in STRUCTURE_SPECS.items():
            if spec.target == target_name:
                features = []
                if spec.ensemble_models:
                    features.append(f"NMR={spec.ensemble_models}")
                if spec.fibril_chains:
                    features.append(f"fibril={spec.fibril_chains}")
                if spec.chains:
                    features.append(f"chains={spec.chains}")
                if spec.remove_dna:
                    features.append("remove_DNA")
                if spec.keep_ligand:
                    features.append(f"keep_lig={spec.keep_ligand}")
                if spec.generate_mutants:
                    features.append("mutants")
                if spec.remove_zinc:
                    features.append("no_Zn")

                feat_str = f" [{', '.join(features)}]" if features else ""
                print(f"    - {spec_id}: {spec.description}{feat_str}")

    print("\n" + "="*70 + "\n")


def main():
    parser = argparse.ArgumentParser(
        description="Download and prepare protein structures for DynamicBind",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
    python download_structures.py                    # All targets
    python download_structures.py --target cmyc     # Only c-Myc structures
    python download_structures.py --target p53 --pdb 2AHI  # Only p53 2AHI
    python download_structures.py --reset           # Reset and redo all
    python download_structures.py --list            # List all structures
        """
    )
    parser.add_argument(
        "-t", "--target",
        choices=list(TARGETS.keys()) + ["all"],
        default="all",
        help="Target protein to prepare (default: all)"
    )
    parser.add_argument(
        "-p", "--pdb",
        help="Specific PDB ID to prepare"
    )
    parser.add_argument(
        "-o", "--output",
        default=str(PROJECT_ROOT / "data" / "structures"),
        help="Output directory"
    )
    parser.add_argument(
        "--timeout",
        type=int,
        default=120,
        help="PDBFixer timeout in seconds (default: 120)"
    )
    parser.add_argument(
        "--reset",
        action="store_true",
        help="Reset checkpoint and reprocess all structures"
    )
    parser.add_argument(
        "--list",
        action="store_true",
        help="List all structures and exit"
    )
    parser.add_argument(
        "--skip-fix",
        action="store_true",
        help="Skip structure fixing (just download)"
    )

    args = parser.parse_args()

    # Print header
    logger.info("=" * 60)
    logger.info("DynamicBind Structure Preparation")
    logger.info("=" * 60)

    if args.list:
        print_structure_summary()
        return

    # Create preparer
    logger.info(f"Output directory: {args.output}")
    logger.info(f"PDBFixer timeout: {args.timeout}s")

    preparer = StructurePreparer(args.output, pdbfixer_timeout=args.timeout)

    # Reset if requested
    if args.reset:
        logger.info("Resetting checkpoint - will reprocess all structures")
        preparer.checkpoint = {"completed": [], "failed": [], "skipped": []}
        preparer._save_checkpoint()

    # Run preparation
    if args.target == "all":
        logger.info("Preparing all target proteins...")
        logger.info(f"Targets: {list(TARGETS.keys())}")
        logger.info(f"Total structures: {len(STRUCTURE_SPECS)}")

        print_structure_summary()

        results = preparer.prepare_all_targets()

        # Summary
        logger.info("\n" + "=" * 60)
        logger.info("SUMMARY")
        logger.info("=" * 60)

        total_files = 0
        for target, files in results.items():
            logger.info(f"\n{TARGETS[target].name}: {len(files)} structures")
            for f in files:
                logger.info(f"  - {Path(f).name}")
            total_files += len(files)

        logger.info(f"\nTotal prepared files: {total_files}")

    else:
        logger.info(f"Preparing target: {args.target}")
        if args.pdb:
            logger.info(f"PDB ID: {args.pdb}")

        prepared = preparer.prepare_target(args.target, args.pdb)

        logger.info("\n" + "=" * 60)
        logger.info(f"Prepared {len(prepared)} structures:")
        for f in prepared:
            logger.info(f"  - {f}")

    # Final status
    logger.info("\n" + "=" * 60)
    logger.info("Structure preparation complete!")
    logger.info(f"Completed: {len(preparer.checkpoint['completed'])}")
    logger.info(f"Failed: {len(preparer.checkpoint['failed'])}")
    logger.info("=" * 60)


if __name__ == "__main__":
    main()
