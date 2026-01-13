"""
Structure Preparation Module

Handles downloading, cleaning, and preparing protein structures for MD simulation.
Supports:
- PDB download from RCSB
- PED ensemble download (for IDPs)
- Chain extraction
- Missing residue modeling
- Protonation state assignment
"""

import os
import logging
import requests
from pathlib import Path
from typing import Optional, List, Dict, Tuple
from dataclasses import dataclass

# Configure logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

# Try importing structure handling libraries
try:
    from Bio.PDB import PDBParser, PDBIO, Select
    from Bio.PDB.PDBExceptions import PDBConstructionWarning
    import warnings
    warnings.filterwarnings('ignore', category=PDBConstructionWarning)
    BIOPYTHON_AVAILABLE = True
except ImportError:
    BIOPYTHON_AVAILABLE = False
    logger.warning("BioPython not installed. Install with: pip install biopython")

try:
    from pdbfixer import PDBFixer
    from openmm.app import PDBFile
    PDBFIXER_AVAILABLE = True
except ImportError:
    PDBFIXER_AVAILABLE = False
    logger.warning("PDBFixer not installed. Install with: conda install -c conda-forge pdbfixer")


@dataclass
class TargetProtein:
    """Configuration for a target protein."""
    name: str
    disease: str
    pdb_ids: List[str]
    chains: Dict[str, str]  # pdb_id -> chain to extract
    description: str
    ped_ids: Optional[List[str]] = None  # For IDP ensembles


# Define all target proteins
TARGETS = {
    "abeta42": TargetProtein(
        name="Amyloid-Beta 42",
        disease="Alzheimer's",
        pdb_ids=["2NAO", "1Z0Q", "1IYT", "2BEG"],
        chains={"2NAO": "A", "1Z0Q": "A", "1IYT": "A", "2BEG": "A"},
        description="Amyloid-beta peptide implicated in Alzheimer's disease pathology"
    ),
    "tau_k18": TargetProtein(
        name="Tau K18",
        disease="Alzheimer's",
        pdb_ids=["2MZ7"],
        chains={"2MZ7": "A"},
        description="Tau protein microtubule-binding repeat region",
        ped_ids=["PED00192", "PED00443"]
    ),
    "cmyc": TargetProtein(
        name="c-Myc",
        disease="Cancer",
        pdb_ids=["6G6K", "1NKP", "2A93"],
        chains={"6G6K": "A", "1NKP": "A", "2A93": "A"},
        description="c-Myc oncogene transcription factor involved in ~70% of cancers"
    ),
    "alpha_synuclein": TargetProtein(
        name="Alpha-Synuclein",
        disease="Parkinson's",
        pdb_ids=["8A9L", "8FPT", "2KKW"],
        chains={"8A9L": "A", "8FPT": "A", "2KKW": "A"},
        description="Alpha-synuclein protein forming Lewy bodies in Parkinson's disease"
    ),
    "p53": TargetProtein(
        name="p53",
        disease="Cancer",
        pdb_ids=["2AHI", "2FEJ"],
        chains={"2AHI": "A", "2FEJ": "A"},
        description="p53 tumor suppressor DNA-binding domain"
    )
}


class ChainSelector(Select):
    """BioPython Select class for extracting specific chains."""

    def __init__(self, chain_ids: List[str]):
        self.chain_ids = chain_ids

    def accept_chain(self, chain):
        return chain.get_id() in self.chain_ids


class StructurePreparer:
    """
    Prepares protein structures for MD simulation.

    Handles:
    - Downloading from PDB/PED databases
    - Chain extraction
    - Missing residue/atom fixing
    - Protonation state assignment
    """

    def __init__(self, output_dir: str = "data/structures"):
        """
        Initialize structure preparer.

        Args:
            output_dir: Base directory for structure files
        """
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(parents=True, exist_ok=True)

        if BIOPYTHON_AVAILABLE:
            self.parser = PDBParser(QUIET=True)
            self.io = PDBIO()
        else:
            self.parser = None
            self.io = None

    def download_pdb(self, pdb_id: str, output_path: Optional[str] = None) -> str:
        """
        Download PDB file from RCSB.

        Args:
            pdb_id: 4-character PDB ID
            output_path: Optional output path

        Returns:
            Path to downloaded file
        """
        pdb_id = pdb_id.upper()

        if output_path is None:
            output_path = self.output_dir / f"{pdb_id}.pdb"
        else:
            output_path = Path(output_path)

        if output_path.exists():
            logger.info(f"PDB {pdb_id} already exists at {output_path}")
            return str(output_path)

        # Try RCSB PDB
        url = f"https://files.rcsb.org/download/{pdb_id}.pdb"
        logger.info(f"Downloading {pdb_id} from RCSB...")

        try:
            response = requests.get(url, timeout=30)
            response.raise_for_status()

            output_path.parent.mkdir(parents=True, exist_ok=True)
            with open(output_path, 'w') as f:
                f.write(response.text)

            logger.info(f"Saved to {output_path}")
            return str(output_path)

        except requests.RequestException as e:
            logger.error(f"Failed to download {pdb_id}: {e}")
            raise

    def download_ped_ensemble(
        self,
        ped_id: str,
        output_dir: Optional[str] = None
    ) -> List[str]:
        """
        Download IDP ensemble from Protein Ensemble Database (PED).

        Args:
            ped_id: PED entry ID (e.g., PED00192)
            output_dir: Output directory for ensemble files

        Returns:
            List of paths to downloaded conformer files
        """
        if output_dir is None:
            output_dir = self.output_dir / ped_id
        else:
            output_dir = Path(output_dir)

        output_dir.mkdir(parents=True, exist_ok=True)

        # PED ensemble download URL
        base_url = f"https://proteinensemble.org/api/entries/{ped_id}/ensemble"

        logger.info(f"Downloading PED ensemble {ped_id}...")

        try:
            response = requests.get(base_url, timeout=60)
            response.raise_for_status()

            # Save full ensemble
            ensemble_path = output_dir / f"{ped_id}_ensemble.pdb"
            with open(ensemble_path, 'w') as f:
                f.write(response.text)

            logger.info(f"Saved ensemble to {ensemble_path}")
            return [str(ensemble_path)]

        except requests.RequestException as e:
            logger.warning(f"Failed to download from PED API: {e}")
            logger.info("PED entries may require manual download from proteinensemble.org")
            return []

    def extract_chain(
        self,
        pdb_path: str,
        chain_id: str,
        output_path: Optional[str] = None
    ) -> str:
        """
        Extract specific chain from PDB file.

        Args:
            pdb_path: Path to input PDB
            chain_id: Chain ID to extract
            output_path: Output path (optional)

        Returns:
            Path to output file with extracted chain
        """
        if not BIOPYTHON_AVAILABLE:
            raise ImportError("BioPython required for chain extraction")

        pdb_path = Path(pdb_path)
        if output_path is None:
            output_path = pdb_path.parent / f"{pdb_path.stem}_chain{chain_id}.pdb"

        logger.info(f"Extracting chain {chain_id} from {pdb_path.name}...")

        structure = self.parser.get_structure("protein", str(pdb_path))

        self.io.set_structure(structure)
        self.io.save(str(output_path), ChainSelector([chain_id]))

        logger.info(f"Saved chain {chain_id} to {output_path}")
        return str(output_path)

    def fix_structure(
        self,
        pdb_path: str,
        output_path: Optional[str] = None,
        add_hydrogens: bool = True,
        add_missing_residues: bool = True,
        ph: float = 7.0
    ) -> str:
        """
        Fix PDB structure using PDBFixer.

        Fixes:
        - Missing heavy atoms
        - Missing residues (optional)
        - Missing hydrogens (optional)
        - Non-standard residues

        Args:
            pdb_path: Input PDB path
            output_path: Output path
            add_hydrogens: Add missing hydrogens
            add_missing_residues: Add missing residues
            ph: pH for protonation states

        Returns:
            Path to fixed structure
        """
        if not PDBFIXER_AVAILABLE:
            raise ImportError("PDBFixer required for structure fixing")

        pdb_path = Path(pdb_path)
        if output_path is None:
            output_path = pdb_path.parent / f"{pdb_path.stem}_fixed.pdb"

        logger.info(f"Fixing structure {pdb_path.name}...")

        fixer = PDBFixer(filename=str(pdb_path))

        # Find and replace non-standard residues
        fixer.findNonstandardResidues()
        fixer.replaceNonstandardResidues()

        # Find and add missing residues
        if add_missing_residues:
            fixer.findMissingResidues()
            # Only add terminal missing residues, not internal gaps
            # (internal gaps might be intentional for IDPs)
            keys_to_remove = []
            for key in fixer.missingResidues:
                # Keep only if it's a short gap (< 5 residues)
                if len(fixer.missingResidues[key]) > 5:
                    keys_to_remove.append(key)
            for key in keys_to_remove:
                del fixer.missingResidues[key]

        # Find and add missing atoms (including hydrogens)
        fixer.findMissingAtoms()
        fixer.addMissingAtoms()

        if add_hydrogens:
            fixer.addMissingHydrogens(ph)

        # Remove water (will re-solvate later)
        fixer.removeHeterogens(keepWater=False)

        # Save fixed structure
        with open(output_path, 'w') as f:
            PDBFile.writeFile(fixer.topology, fixer.positions, f)

        logger.info(f"Saved fixed structure to {output_path}")
        return str(output_path)

    def prepare_target(
        self,
        target_name: str,
        pdb_id: Optional[str] = None
    ) -> List[str]:
        """
        Prepare all structures for a target protein.

        Args:
            target_name: Target name (e.g., 'cmyc', 'p53')
            pdb_id: Specific PDB ID (optional, prepares all if None)

        Returns:
            List of paths to prepared structures
        """
        if target_name not in TARGETS:
            raise ValueError(f"Unknown target: {target_name}. "
                           f"Available: {list(TARGETS.keys())}")

        target = TARGETS[target_name]
        prepared_files = []

        # Create target output directory
        target_dir = self.output_dir / target_name
        pdb_dir = target_dir / "pdb"
        processed_dir = target_dir / "processed"
        pdb_dir.mkdir(parents=True, exist_ok=True)
        processed_dir.mkdir(parents=True, exist_ok=True)

        # Determine which PDB IDs to process
        pdb_ids = [pdb_id] if pdb_id else target.pdb_ids

        for pid in pdb_ids:
            try:
                # Download
                raw_path = self.download_pdb(pid, pdb_dir / f"{pid}.pdb")

                # Extract chain
                chain = target.chains.get(pid, "A")
                chain_path = self.extract_chain(
                    raw_path,
                    chain,
                    processed_dir / f"{pid}_chain{chain}.pdb"
                )

                # Fix structure
                fixed_path = self.fix_structure(
                    chain_path,
                    processed_dir / f"{pid}_chain{chain}_fixed.pdb"
                )

                prepared_files.append(fixed_path)
                logger.info(f"Successfully prepared {pid}")

            except Exception as e:
                logger.error(f"Failed to prepare {pid}: {e}")
                continue

        # Download PED ensembles if available
        if target.ped_ids:
            for ped_id in target.ped_ids:
                try:
                    ensemble_files = self.download_ped_ensemble(
                        ped_id,
                        target_dir / "ped_ensembles" / ped_id
                    )
                    prepared_files.extend(ensemble_files)
                except Exception as e:
                    logger.warning(f"Failed to download PED {ped_id}: {e}")

        return prepared_files

    def prepare_all_targets(self) -> Dict[str, List[str]]:
        """
        Prepare structures for all target proteins.

        Returns:
            Dictionary mapping target names to list of prepared file paths
        """
        results = {}

        for target_name in TARGETS:
            logger.info(f"\n{'='*60}")
            logger.info(f"Preparing {TARGETS[target_name].name}")
            logger.info(f"Disease: {TARGETS[target_name].disease}")
            logger.info(f"{'='*60}")

            prepared = self.prepare_target(target_name)
            results[target_name] = prepared

        return results


def main():
    """Main function for structure preparation."""
    import argparse

    parser = argparse.ArgumentParser(description="Prepare protein structures for MD")
    parser.add_argument(
        "-t", "--target",
        choices=list(TARGETS.keys()) + ["all"],
        default="all",
        help="Target protein to prepare"
    )
    parser.add_argument(
        "-p", "--pdb",
        help="Specific PDB ID to prepare"
    )
    parser.add_argument(
        "-o", "--output",
        default="data/structures",
        help="Output directory"
    )

    args = parser.parse_args()

    preparer = StructurePreparer(args.output)

    if args.target == "all":
        preparer.prepare_all_targets()
    else:
        preparer.prepare_target(args.target, args.pdb)


if __name__ == "__main__":
    main()
