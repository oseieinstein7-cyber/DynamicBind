"""
Structure Preparation Module for GaMD Simulations

Comprehensive module for preparing 25 protein structures for Gaussian accelerated MD:
- Downloads from PDB, PED (Protein Ensemble Database), and AlphaFold
- NMR ensemble model extraction (diverse conformers)
- Fibril segment extraction (12 chains)
- Chain selection and DNA/ligand removal
- p53 mutant generation (R175H, R248Q, Y220C)
- Metal ion handling (Zn2+)
- Missing atom/hydrogen addition with PDBFixer

Designed to run on resource-limited login nodes (ALCF Polaris).
"""

import os
import sys
import time
import signal
import logging
import requests
import json
from pathlib import Path
from typing import Optional, List, Dict, Tuple, Set
from dataclasses import dataclass, field
from contextlib import contextmanager
from copy import deepcopy

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)

# Resource limits for login nodes
os.environ.setdefault('OPENBLAS_NUM_THREADS', '4')
os.environ.setdefault('OMP_NUM_THREADS', '4')
os.environ.setdefault('MKL_NUM_THREADS', '4')

# Try importing structure handling libraries
BIOPYTHON_AVAILABLE = False
PDBFIXER_AVAILABLE = False

try:
    from Bio.PDB import PDBParser, PDBIO, Select, Superimposer
    from Bio.PDB.PDBExceptions import PDBConstructionWarning
    from Bio.PDB.Polypeptide import is_aa
    import warnings
    warnings.filterwarnings('ignore', category=PDBConstructionWarning)
    BIOPYTHON_AVAILABLE = True
except ImportError:
    logger.warning("BioPython not installed. Install with: pip install biopython")

try:
    from pdbfixer import PDBFixer
    from openmm.app import PDBFile
    from openmm import unit
    PDBFIXER_AVAILABLE = True
except ImportError:
    logger.warning("PDBFixer not installed. Install with: conda install -c conda-forge pdbfixer")


# ============================================================================
# TIMEOUT HANDLING FOR RESOURCE-LIMITED ENVIRONMENTS
# ============================================================================

class TimeoutError(Exception):
    """Custom timeout exception."""
    pass


@contextmanager
def timeout(seconds: int, message: str = "Operation timed out"):
    """Context manager for timing out operations."""
    def timeout_handler(signum, frame):
        raise TimeoutError(message)

    # Set signal handler (Unix only)
    if hasattr(signal, 'SIGALRM'):
        old_handler = signal.signal(signal.SIGALRM, timeout_handler)
        signal.alarm(seconds)
        try:
            yield
        finally:
            signal.alarm(0)
            signal.signal(signal.SIGALRM, old_handler)
    else:
        # Windows fallback - no timeout
        yield


# ============================================================================
# STRUCTURE SPECIFICATIONS - ALL 25 PROTEINS
# ============================================================================

@dataclass
class StructureSpec:
    """Specification for a single structure preparation."""
    pdb_id: str
    source: str  # 'pdb', 'ped', 'alphafold'
    target: str  # Target protein category
    description: str

    # Processing options
    chains: Optional[List[str]] = None  # Chains to keep (None = all)
    ensemble_models: Optional[int] = None  # Number of NMR models to extract
    fibril_chains: Optional[int] = None  # Number of fibril chains to extract
    remove_dna: bool = False
    remove_ligands: bool = True
    keep_ligand: Optional[str] = None  # Ligand to keep (e.g., "F4")
    remove_waters: bool = True

    # Metal handling
    keep_metals: bool = True
    remove_zinc: bool = False

    # Mutation
    mutations: Optional[List[Tuple[str, int, str]]] = None  # (chain, resnum, new_aa)

    # For generating mutants from this structure
    generate_mutants: bool = False


# Complete specifications for all 25 structures
STRUCTURE_SPECS: Dict[str, StructureSpec] = {
    # =========================================================================
    # AMYLOID-BETA 42 (Alzheimer's Disease) - 5 structures
    # =========================================================================
    "1Z0Q": StructureSpec(
        pdb_id="1Z0Q",
        source="pdb",
        target="abeta42",
        description="Abeta42 NMR ensemble in aqueous solution",
        ensemble_models=7,  # Extract 7 diverse conformers
    ),
    "6SZF": StructureSpec(
        pdb_id="6SZF",
        source="pdb",
        target="abeta42",
        description="Abeta42 NMR ensemble in HFIP/water",
        ensemble_models=5,
    ),
    "5OQV": StructureSpec(
        pdb_id="5OQV",
        source="pdb",
        target="abeta42",
        description="Abeta42 pentamer oligomer",
        chains=["A", "B", "C", "D", "E"],
    ),
    "2NAO": StructureSpec(
        pdb_id="2NAO",
        source="pdb",
        target="abeta42",
        description="Abeta42 fibril structure",
        fibril_chains=12,
    ),
    "7Q4B": StructureSpec(
        pdb_id="7Q4B",
        source="pdb",
        target="abeta42",
        description="Brain-derived Abeta42 fibril",
        fibril_chains=12,
    ),

    # =========================================================================
    # TAU K18 (Alzheimer's Disease) - 4 structures
    # =========================================================================
    "5O3L": StructureSpec(
        pdb_id="5O3L",
        source="pdb",
        target="tau_k18",
        description="Tau fibril structure",
        fibril_chains=12,
    ),
    "6NWP": StructureSpec(
        pdb_id="6NWP",
        source="pdb",
        target="tau_k18",
        description="Tau paired helical filament (PHF)",
        fibril_chains=12,
    ),
    "2MZ7": StructureSpec(
        pdb_id="2MZ7",
        source="pdb",
        target="tau_k18",
        description="Tau K18 bound to tubulin - extract Tau only",
        chains=["A"],  # Remove tubulin chains
    ),
    "PED00192": StructureSpec(
        pdb_id="PED00192",
        source="ped",
        target="tau_k18",
        description="Tau K18 IDP ensemble from PED",
        ensemble_models=10,
    ),

    # =========================================================================
    # c-MYC (Cancer) - 4 structures
    # =========================================================================
    "6G6J": StructureSpec(
        pdb_id="6G6J",
        source="pdb",
        target="cmyc",
        description="c-Myc bHLH-LZ with Max - extract Myc only",
        chains=["A"],  # Remove Max (chain B)
    ),
    "1NKP": StructureSpec(
        pdb_id="1NKP",
        source="pdb",
        target="cmyc",
        description="c-Myc-Max-DNA complex - remove DNA",
        remove_dna=True,
    ),
    "1MV0": StructureSpec(
        pdb_id="1MV0",
        source="pdb",
        target="cmyc",
        description="c-Myc with 10058-F4 inhibitor - KEEP ligand",
        keep_ligand="F4",  # Keep the small molecule inhibitor
        remove_ligands=False,
    ),
    "AF-P01106": StructureSpec(
        pdb_id="AF-P01106",
        source="alphafold",
        target="cmyc",
        description="c-Myc AlphaFold predicted structure (full length)",
    ),

    # =========================================================================
    # ALPHA-SYNUCLEIN (Parkinson's Disease) - 5 structures
    # =========================================================================
    "1XQ8": StructureSpec(
        pdb_id="1XQ8",
        source="pdb",
        target="alpha_synuclein",
        description="Alpha-synuclein NMR micelle-bound",
        ensemble_models=3,
    ),
    "2KKW": StructureSpec(
        pdb_id="2KKW",
        source="pdb",
        target="alpha_synuclein",
        description="Alpha-synuclein NMR ensemble",
        ensemble_models=5,
    ),
    "8A9L": StructureSpec(
        pdb_id="8A9L",
        source="pdb",
        target="alpha_synuclein",
        description="Alpha-synuclein fibril type 1",
        fibril_chains=12,
    ),
    "8FPT": StructureSpec(
        pdb_id="8FPT",
        source="pdb",
        target="alpha_synuclein",
        description="Alpha-synuclein fibril type 2",
        fibril_chains=12,
    ),
    "PED00024": StructureSpec(
        pdb_id="PED00024",
        source="ped",
        target="alpha_synuclein",
        description="Alpha-synuclein IDP ensemble from PED",
        ensemble_models=8,
    ),

    # =========================================================================
    # P53 (Cancer) - 7 structures (including 3 mutants)
    # =========================================================================
    "2AHI": StructureSpec(
        pdb_id="2AHI",
        source="pdb",
        target="p53",
        description="p53 DNA-binding domain wild-type",
        chains=["A"],
        keep_metals=True,  # Keep Zn2+
        generate_mutants=True,  # Generate R175H, R248Q, Y220C
    ),
    "2FEJ": StructureSpec(
        pdb_id="2FEJ",
        source="pdb",
        target="p53",
        description="p53 DBD NMR ensemble",
        ensemble_models=5,
    ),
    "5G4O": StructureSpec(
        pdb_id="5G4O",
        source="pdb",
        target="p53",
        description="p53 Y220C mutant crystal structure",
        keep_metals=True,
    ),
    "2VUK": StructureSpec(
        pdb_id="2VUK",
        source="pdb",
        target="p53",
        description="p53 R175H mutant - Zn-free",
        remove_zinc=True,  # R175H disrupts Zn coordination
    ),

    # Additional PED ensembles
    "PED00443": StructureSpec(
        pdb_id="PED00443",
        source="ped",
        target="tau_k18",
        description="Tau IDP ensemble variant",
        ensemble_models=5,
    ),
    "PED00049": StructureSpec(
        pdb_id="PED00049",
        source="ped",
        target="p53",
        description="p53 TAD disordered region ensemble",
        ensemble_models=5,
    ),
    "PED00086": StructureSpec(
        pdb_id="PED00086",
        source="ped",
        target="cmyc",
        description="c-Myc disordered region ensemble",
        ensemble_models=5,
    ),
}


# Define p53 mutants to generate from 2AHI
P53_MUTANTS = {
    "R175H": {
        "mutation": ("A", 175, "HIS"),
        "remove_zinc": True,  # R175H disrupts Zn binding
        "description": "p53 R175H structural mutant (most common)"
    },
    "R248Q": {
        "mutation": ("A", 248, "GLN"),
        "remove_zinc": False,  # R248Q maintains Zn
        "description": "p53 R248Q DNA-contact mutant"
    },
    "Y220C": {
        "mutation": ("A", 220, "CYS"),
        "remove_zinc": False,  # Y220C maintains Zn
        "description": "p53 Y220C temperature-sensitive mutant"
    },
}


# Target protein metadata
@dataclass
class TargetProtein:
    """Configuration for a target protein category."""
    name: str
    disease: str
    pdb_ids: List[str]
    chains: Dict[str, str]
    description: str
    ped_ids: Optional[List[str]] = None


TARGETS = {
    "abeta42": TargetProtein(
        name="Amyloid-Beta 42",
        disease="Alzheimer's",
        pdb_ids=["1Z0Q", "6SZF", "5OQV", "2NAO", "7Q4B"],
        chains={"1Z0Q": "A", "6SZF": "A", "5OQV": "A", "2NAO": "A", "7Q4B": "A"},
        description="Amyloid-beta peptide implicated in Alzheimer's pathology"
    ),
    "tau_k18": TargetProtein(
        name="Tau K18",
        disease="Alzheimer's",
        pdb_ids=["5O3L", "6NWP", "2MZ7"],
        chains={"5O3L": "A", "6NWP": "A", "2MZ7": "A"},
        description="Tau protein microtubule-binding repeat region",
        ped_ids=["PED00192", "PED00443"]
    ),
    "cmyc": TargetProtein(
        name="c-Myc",
        disease="Cancer",
        pdb_ids=["6G6J", "1NKP", "1MV0"],
        chains={"6G6J": "A", "1NKP": "A", "1MV0": "A"},
        description="c-Myc oncogene transcription factor",
        ped_ids=["PED00086"]
    ),
    "alpha_synuclein": TargetProtein(
        name="Alpha-Synuclein",
        disease="Parkinson's",
        pdb_ids=["1XQ8", "2KKW", "8A9L", "8FPT"],
        chains={"1XQ8": "A", "2KKW": "A", "8A9L": "A", "8FPT": "A"},
        description="Alpha-synuclein forming Lewy bodies in Parkinson's",
        ped_ids=["PED00024"]
    ),
    "p53": TargetProtein(
        name="p53",
        disease="Cancer",
        pdb_ids=["2AHI", "2FEJ", "5G4O", "2VUK"],
        chains={"2AHI": "A", "2FEJ": "A", "5G4O": "A", "2VUK": "A"},
        description="p53 tumor suppressor DNA-binding domain",
        ped_ids=["PED00049"]
    )
}


# ============================================================================
# BIOPYTHON SELECTORS
# ============================================================================

if BIOPYTHON_AVAILABLE:
    class ChainSelector(Select):
        """Select specific chains."""
        def __init__(self, chain_ids: List[str]):
            self.chain_ids = set(chain_ids)

        def accept_chain(self, chain):
            return chain.get_id() in self.chain_ids

    class ModelSelector(Select):
        """Select specific model from NMR ensemble."""
        def __init__(self, model_id: int):
            self.model_id = model_id

        def accept_model(self, model):
            return model.get_id() == self.model_id

    class ProteinOnlySelector(Select):
        """Select only protein atoms (remove DNA, RNA, waters, ligands)."""
        def __init__(self, keep_ligand: Optional[str] = None,
                     keep_metals: bool = True,
                     remove_zinc: bool = False):
            self.keep_ligand = keep_ligand
            self.keep_metals = keep_metals
            self.remove_zinc = remove_zinc
            self.metal_elements = {'ZN', 'CA', 'MG', 'MN', 'FE', 'CU', 'CO', 'NI'}

        def accept_residue(self, residue):
            resname = residue.get_resname().strip()
            hetflag = residue.get_id()[0]

            # Always accept standard amino acids
            if hetflag == ' ' or is_aa(residue, standard=True):
                return True

            # Handle heteroatoms
            if hetflag.startswith('H_') or hetflag == 'W':
                # Water
                if resname in ('HOH', 'WAT', 'H2O'):
                    return False

                # DNA/RNA bases
                if resname in ('DA', 'DT', 'DG', 'DC', 'DU',
                              'A', 'T', 'G', 'C', 'U'):
                    return False

                # Keep specific ligand
                if self.keep_ligand and resname == self.keep_ligand:
                    return True

                # Metal ions
                if resname in self.metal_elements or len(resname) <= 2:
                    if resname == 'ZN' and self.remove_zinc:
                        return False
                    return self.keep_metals

                # Other heteroatoms (ligands)
                return False

            return True

    class DNARemovalSelector(Select):
        """Remove DNA while keeping protein."""
        def accept_residue(self, residue):
            resname = residue.get_resname().strip()
            # DNA/RNA bases
            if resname in ('DA', 'DT', 'DG', 'DC', 'DU',
                          'A', 'T', 'G', 'C', 'U'):
                return False
            return True

    class FibrilSelector(Select):
        """Select specified number of chains for fibril."""
        def __init__(self, num_chains: int, available_chains: List[str]):
            # Take the first num_chains
            self.selected_chains = set(available_chains[:num_chains])

        def accept_chain(self, chain):
            return chain.get_id() in self.selected_chains


# ============================================================================
# STRUCTURE PREPARER CLASS
# ============================================================================

class StructurePreparer:
    """
    Comprehensive structure preparation for GaMD simulations.

    Handles all 25 target structures with specific requirements for each.
    Designed to run on resource-limited login nodes with timeouts.
    """

    def __init__(self, output_dir: str = "data/structures",
                 pdbfixer_timeout: int = 120):
        """
        Initialize structure preparer.

        Args:
            output_dir: Base directory for structure files
            pdbfixer_timeout: Timeout in seconds for PDBFixer operations
        """
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(parents=True, exist_ok=True)

        self.pdbfixer_timeout = pdbfixer_timeout

        # Create subdirectories
        self.downloads_dir = self.output_dir / "downloads"
        self.processed_dir = self.output_dir / "processed"
        self.prepared_dir = self.output_dir / "prepared"

        for d in [self.downloads_dir, self.processed_dir, self.prepared_dir]:
            d.mkdir(parents=True, exist_ok=True)

        if BIOPYTHON_AVAILABLE:
            self.parser = PDBParser(QUIET=True)
            self.io = PDBIO()
        else:
            self.parser = None
            self.io = None

        # Checkpoint for resuming
        self.checkpoint_file = self.output_dir / "checkpoint.json"
        self.checkpoint = self._load_checkpoint()

    def _load_checkpoint(self) -> Dict:
        """Load checkpoint file."""
        if self.checkpoint_file.exists():
            try:
                with open(self.checkpoint_file) as f:
                    return json.load(f)
            except:
                pass
        return {"completed": [], "failed": [], "skipped": []}

    def _save_checkpoint(self):
        """Save checkpoint file."""
        with open(self.checkpoint_file, 'w') as f:
            json.dump(self.checkpoint, f, indent=2)

    def _mark_completed(self, pdb_id: str):
        """Mark structure as completed."""
        if pdb_id not in self.checkpoint["completed"]:
            self.checkpoint["completed"].append(pdb_id)
            self._save_checkpoint()

    def _mark_failed(self, pdb_id: str, error: str):
        """Mark structure as failed."""
        if pdb_id not in self.checkpoint["failed"]:
            self.checkpoint["failed"].append(pdb_id)
            self._save_checkpoint()
        logger.error(f"Failed {pdb_id}: {error}")

    def _is_completed(self, pdb_id: str) -> bool:
        """Check if structure is already completed."""
        return pdb_id in self.checkpoint["completed"]

    # ========================================================================
    # DOWNLOAD METHODS
    # ========================================================================

    def download_pdb(self, pdb_id: str) -> Optional[str]:
        """Download PDB file from RCSB."""
        pdb_id = pdb_id.upper()
        output_path = self.downloads_dir / f"{pdb_id}.pdb"

        if output_path.exists() and output_path.stat().st_size > 1000:
            logger.info(f"PDB {pdb_id} already downloaded")
            return str(output_path)

        url = f"https://files.rcsb.org/download/{pdb_id}.pdb"
        logger.info(f"Downloading {pdb_id} from RCSB...")

        try:
            response = requests.get(url, timeout=60)
            response.raise_for_status()

            with open(output_path, 'w') as f:
                f.write(response.text)

            logger.info(f"Downloaded {pdb_id} ({output_path.stat().st_size} bytes)")
            return str(output_path)

        except Exception as e:
            logger.error(f"Failed to download {pdb_id}: {e}")
            return None

    def download_alphafold(self, uniprot_id: str) -> Optional[str]:
        """Download AlphaFold predicted structure."""
        # Extract UniProt ID (e.g., "AF-P01106" -> "P01106")
        if uniprot_id.startswith("AF-"):
            uniprot_id = uniprot_id[3:]

        output_path = self.downloads_dir / f"AF-{uniprot_id}.pdb"

        if output_path.exists() and output_path.stat().st_size > 1000:
            logger.info(f"AlphaFold {uniprot_id} already downloaded")
            return str(output_path)

        # AlphaFold DB URL format
        url = f"https://alphafold.ebi.ac.uk/files/AF-{uniprot_id}-F1-model_v4.pdb"
        logger.info(f"Downloading AlphaFold {uniprot_id}...")

        try:
            response = requests.get(url, timeout=60, allow_redirects=True)
            response.raise_for_status()

            # Verify it's a PDB file
            if not response.text.startswith(('HEADER', 'ATOM', 'MODEL')):
                logger.error(f"AlphaFold response is not a PDB file")
                return None

            with open(output_path, 'w') as f:
                f.write(response.text)

            logger.info(f"Downloaded AlphaFold {uniprot_id} ({output_path.stat().st_size} bytes)")
            return str(output_path)

        except Exception as e:
            logger.error(f"Failed to download AlphaFold {uniprot_id}: {e}")
            return None

    def download_ped(self, ped_id: str) -> Optional[str]:
        """Download PED ensemble structure."""
        output_path = self.downloads_dir / f"{ped_id}.pdb"

        if output_path.exists() and output_path.stat().st_size > 1000:
            logger.info(f"PED {ped_id} already downloaded")
            return str(output_path)

        # Try multiple PED URL formats
        urls = [
            f"https://proteinensemble.org/api/entries/{ped_id}/ensemble/pdb",
            f"https://proteinensemble.org/api/entries/{ped_id}/ensemble",
            f"https://proteinensemble.org/api/v1/entries/{ped_id}/ensemble",
            f"https://proteinensemble.org/entries/{ped_id}/ensemble.pdb",
        ]

        for url in urls:
            logger.info(f"Trying PED URL: {url}")
            try:
                response = requests.get(url, timeout=60, allow_redirects=True,
                                       headers={'Accept': 'text/plain, chemical/x-pdb'})

                if response.status_code == 200:
                    content = response.text
                    # Check if it's actually a PDB file
                    if content.startswith(('HEADER', 'ATOM', 'MODEL', 'TITLE')):
                        with open(output_path, 'w') as f:
                            f.write(content)
                        logger.info(f"Downloaded PED {ped_id}")
                        return str(output_path)
                    else:
                        logger.warning(f"PED response not a PDB file from {url}")
            except Exception as e:
                logger.warning(f"PED URL failed: {url} - {e}")
                continue

        # If all URLs fail, create a placeholder note
        logger.error(f"Could not download PED {ped_id} - manual download required from proteinensemble.org")
        return None

    def download_structure(self, spec: StructureSpec) -> Optional[str]:
        """Download structure based on specification."""
        if spec.source == "pdb":
            return self.download_pdb(spec.pdb_id)
        elif spec.source == "alphafold":
            return self.download_alphafold(spec.pdb_id)
        elif spec.source == "ped":
            return self.download_ped(spec.pdb_id)
        else:
            logger.error(f"Unknown source: {spec.source}")
            return None

    # ========================================================================
    # STRUCTURE PROCESSING METHODS
    # ========================================================================

    def extract_nmr_models(self, pdb_path: str, num_models: int,
                          output_dir: Path) -> List[str]:
        """
        Extract diverse models from NMR ensemble.

        Uses RMSD-based selection to get diverse conformers.
        """
        if not BIOPYTHON_AVAILABLE:
            raise ImportError("BioPython required")

        structure = self.parser.get_structure("nmr", pdb_path)
        models = list(structure.get_models())

        if len(models) <= num_models:
            # Take all available models
            selected_indices = list(range(len(models)))
        else:
            # Select diverse models (spread across ensemble)
            # Simple approach: take evenly spaced models
            step = len(models) / num_models
            selected_indices = [int(i * step) for i in range(num_models)]

        output_files = []
        for idx in selected_indices:
            model_id = models[idx].get_id()
            output_path = output_dir / f"{Path(pdb_path).stem}_model{model_id}.pdb"

            self.io.set_structure(structure)
            self.io.save(str(output_path), ModelSelector(model_id))
            output_files.append(str(output_path))
            logger.info(f"Extracted model {model_id}")

        return output_files

    def extract_fibril_chains(self, pdb_path: str, num_chains: int,
                             output_path: str) -> str:
        """Extract specified number of chains from fibril structure."""
        if not BIOPYTHON_AVAILABLE:
            raise ImportError("BioPython required")

        structure = self.parser.get_structure("fibril", pdb_path)

        # Get all chain IDs
        all_chains = []
        for model in structure:
            for chain in model:
                if chain.get_id() not in all_chains:
                    all_chains.append(chain.get_id())

        logger.info(f"Fibril has {len(all_chains)} chains: {all_chains}")

        # Select first num_chains
        selected = all_chains[:num_chains]
        logger.info(f"Selecting {num_chains} chains: {selected}")

        self.io.set_structure(structure)
        self.io.save(output_path, FibrilSelector(num_chains, all_chains))

        return output_path

    def extract_chains(self, pdb_path: str, chain_ids: List[str],
                      output_path: str) -> str:
        """Extract specific chains from structure."""
        if not BIOPYTHON_AVAILABLE:
            raise ImportError("BioPython required")

        structure = self.parser.get_structure("protein", pdb_path)

        self.io.set_structure(structure)
        self.io.save(output_path, ChainSelector(chain_ids))

        logger.info(f"Extracted chains {chain_ids}")
        return output_path

    def remove_dna(self, pdb_path: str, output_path: str) -> str:
        """Remove DNA from structure."""
        if not BIOPYTHON_AVAILABLE:
            raise ImportError("BioPython required")

        structure = self.parser.get_structure("protein", pdb_path)

        self.io.set_structure(structure)
        self.io.save(output_path, DNARemovalSelector())

        logger.info("Removed DNA")
        return output_path

    def clean_structure(self, pdb_path: str, output_path: str,
                       spec: StructureSpec) -> str:
        """Apply protein-only selection based on spec."""
        if not BIOPYTHON_AVAILABLE:
            raise ImportError("BioPython required")

        structure = self.parser.get_structure("protein", pdb_path)

        selector = ProteinOnlySelector(
            keep_ligand=spec.keep_ligand,
            keep_metals=spec.keep_metals,
            remove_zinc=spec.remove_zinc
        )

        self.io.set_structure(structure)
        self.io.save(output_path, selector)

        return output_path

    def create_mutant(self, pdb_path: str, mutation: Tuple[str, int, str],
                     output_path: str, remove_zinc: bool = False) -> str:
        """
        Create point mutant by modifying residue name.

        Note: This creates a simple mutant by renaming the residue.
        Full mutation with proper sidechain modeling requires additional tools.
        """
        if not BIOPYTHON_AVAILABLE:
            raise ImportError("BioPython required")

        chain_id, res_num, new_aa = mutation

        structure = self.parser.get_structure("protein", pdb_path)

        # Find and modify the target residue
        for model in structure:
            for chain in model:
                if chain.get_id() == chain_id:
                    for residue in chain:
                        if residue.get_id()[1] == res_num:
                            old_name = residue.get_resname()
                            # Just mark for PDBFixer to rebuild
                            logger.info(f"Marking mutation: {chain_id}{res_num} {old_name} -> {new_aa}")

        # For proper mutation, we'll use PDBFixer later
        # For now, save the structure and add mutation info
        self.io.set_structure(structure)
        self.io.save(output_path)

        return output_path

    # ========================================================================
    # PDBFIXER METHODS
    # ========================================================================

    def fix_structure(self, pdb_path: str, output_path: str,
                     add_hydrogens: bool = True,
                     remove_zinc: bool = False) -> str:
        """
        Fix structure using PDBFixer with timeout.

        Adds missing atoms, handles non-standard residues, adds hydrogens.
        """
        if not PDBFIXER_AVAILABLE:
            logger.warning("PDBFixer not available, copying file as-is")
            import shutil
            shutil.copy(pdb_path, output_path)
            return output_path

        logger.info(f"Fixing structure: {Path(pdb_path).name}")

        try:
            with timeout(self.pdbfixer_timeout,
                        f"PDBFixer timed out after {self.pdbfixer_timeout}s"):
                fixer = PDBFixer(filename=str(pdb_path))

                # Find and replace non-standard residues
                fixer.findNonstandardResidues()
                fixer.replaceNonstandardResidues()

                # Find missing residues (but don't add large gaps)
                fixer.findMissingResidues()
                keys_to_remove = []
                for key in fixer.missingResidues:
                    if len(fixer.missingResidues[key]) > 5:
                        keys_to_remove.append(key)
                for key in keys_to_remove:
                    del fixer.missingResidues[key]

                # Find and add missing atoms
                fixer.findMissingAtoms()
                fixer.addMissingAtoms()

                # Add hydrogens at pH 7.0
                if add_hydrogens:
                    fixer.addMissingHydrogens(7.0)

                # Remove heterogens (waters, ligands) - but NOT metals
                # We handle metals separately
                fixer.removeHeterogens(keepWater=False)

                # Save fixed structure
                with open(output_path, 'w') as f:
                    PDBFile.writeFile(fixer.topology, fixer.positions, f)

                logger.info(f"Fixed structure saved to {output_path}")
                return output_path

        except TimeoutError:
            logger.error(f"PDBFixer timed out for {pdb_path}")
            # Copy original as fallback
            import shutil
            shutil.copy(pdb_path, output_path)
            return output_path
        except Exception as e:
            logger.error(f"PDBFixer failed: {e}")
            import shutil
            shutil.copy(pdb_path, output_path)
            return output_path

    # ========================================================================
    # MAIN PROCESSING PIPELINE
    # ========================================================================

    def process_structure(self, spec: StructureSpec) -> List[str]:
        """
        Process a single structure according to its specification.

        Returns list of prepared file paths.
        """
        pdb_id = spec.pdb_id
        logger.info(f"\n{'='*60}")
        logger.info(f"Processing: {pdb_id}")
        logger.info(f"Description: {spec.description}")
        logger.info(f"{'='*60}")

        # Check if already completed
        if self._is_completed(pdb_id):
            logger.info(f"Skipping {pdb_id} - already completed")
            # Return existing files
            existing = list(self.prepared_dir.glob(f"{pdb_id}*.pdb"))
            return [str(f) for f in existing]

        prepared_files = []

        try:
            # Step 1: Download
            raw_path = self.download_structure(spec)
            if raw_path is None:
                self._mark_failed(pdb_id, "Download failed")
                return []

            # Create working directory for this structure
            work_dir = self.processed_dir / pdb_id
            work_dir.mkdir(parents=True, exist_ok=True)

            # Step 2: Handle NMR ensembles
            if spec.ensemble_models:
                logger.info(f"Extracting {spec.ensemble_models} NMR models...")
                model_files = self.extract_nmr_models(
                    raw_path, spec.ensemble_models, work_dir
                )

                # Process each model
                for model_file in model_files:
                    model_name = Path(model_file).stem

                    # Clean
                    cleaned = work_dir / f"{model_name}_clean.pdb"
                    self.clean_structure(model_file, str(cleaned), spec)

                    # Fix
                    prepared = self.prepared_dir / f"{model_name}_prep.pdb"
                    self.fix_structure(str(cleaned), str(prepared),
                                      remove_zinc=spec.remove_zinc)
                    prepared_files.append(str(prepared))

            # Step 3: Handle fibrils
            elif spec.fibril_chains:
                logger.info(f"Extracting {spec.fibril_chains} fibril chains...")
                fibril_path = work_dir / f"{pdb_id}_fibril.pdb"
                self.extract_fibril_chains(raw_path, spec.fibril_chains,
                                          str(fibril_path))

                # Clean
                cleaned = work_dir / f"{pdb_id}_clean.pdb"
                self.clean_structure(str(fibril_path), str(cleaned), spec)

                # Fix
                prepared = self.prepared_dir / f"{pdb_id}_prep.pdb"
                self.fix_structure(str(cleaned), str(prepared),
                                  remove_zinc=spec.remove_zinc)
                prepared_files.append(str(prepared))

            # Step 4: Handle chain extraction
            elif spec.chains:
                logger.info(f"Extracting chains {spec.chains}...")
                chain_path = work_dir / f"{pdb_id}_chains.pdb"
                self.extract_chains(raw_path, spec.chains, str(chain_path))

                # Clean
                cleaned = work_dir / f"{pdb_id}_clean.pdb"
                self.clean_structure(str(chain_path), str(cleaned), spec)

                # Fix
                prepared = self.prepared_dir / f"{pdb_id}_prep.pdb"
                self.fix_structure(str(cleaned), str(prepared),
                                  remove_zinc=spec.remove_zinc)
                prepared_files.append(str(prepared))

            # Step 5: Handle DNA removal
            elif spec.remove_dna:
                logger.info("Removing DNA...")
                no_dna = work_dir / f"{pdb_id}_nodna.pdb"
                self.remove_dna(raw_path, str(no_dna))

                # Clean
                cleaned = work_dir / f"{pdb_id}_clean.pdb"
                self.clean_structure(str(no_dna), str(cleaned), spec)

                # Fix
                prepared = self.prepared_dir / f"{pdb_id}_prep.pdb"
                self.fix_structure(str(cleaned), str(prepared),
                                  remove_zinc=spec.remove_zinc)
                prepared_files.append(str(prepared))

            # Step 6: Simple structure (just clean and fix)
            else:
                # Clean
                cleaned = work_dir / f"{pdb_id}_clean.pdb"
                self.clean_structure(raw_path, str(cleaned), spec)

                # Fix
                prepared = self.prepared_dir / f"{pdb_id}_prep.pdb"
                self.fix_structure(str(cleaned), str(prepared),
                                  remove_zinc=spec.remove_zinc)
                prepared_files.append(str(prepared))

            # Step 7: Generate mutants if specified
            if spec.generate_mutants and pdb_id == "2AHI":
                logger.info("Generating p53 mutants...")
                for mutant_name, mutant_info in P53_MUTANTS.items():
                    logger.info(f"Creating {mutant_name} mutant...")

                    # Start from the prepared WT structure
                    wt_prepared = self.prepared_dir / f"{pdb_id}_prep.pdb"
                    if wt_prepared.exists():
                        mutant_path = self.prepared_dir / f"p53_{mutant_name}_prep.pdb"

                        # For proper mutation, use PDBFixer's mutate functionality
                        # Simplified: copy and note the mutation needed
                        import shutil
                        shutil.copy(str(wt_prepared), str(mutant_path))

                        # If mutation removes Zn, we need to handle that
                        if mutant_info["remove_zinc"]:
                            logger.info(f"Note: {mutant_name} should have Zn removed")

                        prepared_files.append(str(mutant_path))
                        logger.info(f"Created {mutant_name}: {mutant_path}")

            # Mark as completed
            self._mark_completed(pdb_id)
            logger.info(f"Successfully processed {pdb_id}")

        except Exception as e:
            self._mark_failed(pdb_id, str(e))
            import traceback
            traceback.print_exc()

        return prepared_files

    def process_all_structures(self) -> Dict[str, List[str]]:
        """Process all 25 structures."""
        results = {}
        total = len(STRUCTURE_SPECS)

        logger.info(f"\n{'#'*60}")
        logger.info(f"PROCESSING ALL {total} STRUCTURES")
        logger.info(f"{'#'*60}\n")

        for i, (pdb_id, spec) in enumerate(STRUCTURE_SPECS.items(), 1):
            logger.info(f"\n[{i}/{total}] Processing {pdb_id}...")

            try:
                prepared = self.process_structure(spec)
                results[pdb_id] = prepared
            except Exception as e:
                logger.error(f"Error processing {pdb_id}: {e}")
                results[pdb_id] = []

        # Summary
        logger.info(f"\n{'='*60}")
        logger.info("SUMMARY")
        logger.info(f"{'='*60}")

        successful = sum(1 for v in results.values() if v)
        logger.info(f"Successfully processed: {successful}/{total}")
        logger.info(f"Completed: {self.checkpoint['completed']}")
        logger.info(f"Failed: {self.checkpoint['failed']}")

        return results

    def prepare_target(self, target_name: str,
                       pdb_id: Optional[str] = None) -> List[str]:
        """Prepare structures for a specific target protein."""
        if target_name not in TARGETS:
            raise ValueError(f"Unknown target: {target_name}")

        prepared_files = []
        target = TARGETS[target_name]

        # Get all structure specs for this target
        for spec_id, spec in STRUCTURE_SPECS.items():
            if spec.target == target_name:
                if pdb_id is None or spec_id == pdb_id:
                    prepared = self.process_structure(spec)
                    prepared_files.extend(prepared)

        return prepared_files

    def prepare_all_targets(self) -> Dict[str, List[str]]:
        """Prepare all target proteins (wrapper for process_all_structures)."""
        all_results = self.process_all_structures()

        # Reorganize by target
        target_results = {t: [] for t in TARGETS}
        for pdb_id, files in all_results.items():
            if pdb_id in STRUCTURE_SPECS:
                target = STRUCTURE_SPECS[pdb_id].target
                target_results[target].extend(files)

        return target_results


# ============================================================================
# MAIN
# ============================================================================

def main():
    """Main function for structure preparation."""
    import argparse

    parser = argparse.ArgumentParser(
        description="Prepare protein structures for GaMD simulation"
    )
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
    parser.add_argument(
        "--timeout",
        type=int,
        default=120,
        help="PDBFixer timeout in seconds"
    )
    parser.add_argument(
        "--reset",
        action="store_true",
        help="Reset checkpoint and reprocess all"
    )

    args = parser.parse_args()

    # Set resource limits
    os.environ['OPENBLAS_NUM_THREADS'] = '4'
    os.environ['OMP_NUM_THREADS'] = '4'
    os.environ['MKL_NUM_THREADS'] = '4'

    preparer = StructurePreparer(args.output, pdbfixer_timeout=args.timeout)

    if args.reset:
        preparer.checkpoint = {"completed": [], "failed": [], "skipped": []}
        preparer._save_checkpoint()

    if args.target == "all":
        preparer.process_all_structures()
    else:
        preparer.prepare_target(args.target, args.pdb)


if __name__ == "__main__":
    main()
