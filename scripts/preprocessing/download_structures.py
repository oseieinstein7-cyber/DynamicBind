#!/usr/bin/env python3
"""
DynamicBind - Structure Download and Preparation Script

Downloads all PDB structures for target proteins and prepares them for MD simulation.

Usage:
    python download_structures.py                    # Download all targets
    python download_structures.py --target cmyc     # Download specific target
    python download_structures.py --target cmyc --pdb 6G6K  # Download specific PDB
"""

import os
import sys
import argparse
import logging
from pathlib import Path

# Add src to path
PROJECT_ROOT = Path(__file__).parent.parent.parent
sys.path.insert(0, str(PROJECT_ROOT / "src"))

from phase1_md.simulation.prepare_structure import (
    StructurePreparer,
    TARGETS
)

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)


def main():
    parser = argparse.ArgumentParser(
        description="Download and prepare protein structures for DynamicBind"
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
        "--skip-fix",
        action="store_true",
        help="Skip structure fixing (just download)"
    )

    args = parser.parse_args()

    logger.info("=" * 60)
    logger.info("DynamicBind Structure Preparation")
    logger.info("=" * 60)

    # Create preparer
    preparer = StructurePreparer(args.output)

    if args.target == "all":
        logger.info("Preparing all target proteins...")
        logger.info(f"Targets: {list(TARGETS.keys())}")
        results = preparer.prepare_all_targets()

        # Summary
        logger.info("\n" + "=" * 60)
        logger.info("SUMMARY")
        logger.info("=" * 60)
        for target, files in results.items():
            logger.info(f"{target}: {len(files)} structures prepared")
            for f in files:
                logger.info(f"  - {Path(f).name}")
    else:
        logger.info(f"Preparing target: {args.target}")
        if args.pdb:
            logger.info(f"PDB ID: {args.pdb}")

        prepared = preparer.prepare_target(args.target, args.pdb)

        logger.info("\n" + "=" * 60)
        logger.info(f"Prepared {len(prepared)} structures:")
        for f in prepared:
            logger.info(f"  - {f}")

    logger.info("\n" + "=" * 60)
    logger.info("Structure preparation complete!")
    logger.info("=" * 60)


if __name__ == "__main__":
    main()
