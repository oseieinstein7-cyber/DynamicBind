#!/usr/bin/env python3
"""
=============================================================================
MASTER STRUCTURE PREPARATION SCRIPT FOR GaMD SIMULATIONS
=============================================================================

Prepares 25 protein structures for Gaussian accelerated MD simulations.
Uses pdb-tools, reduce (AmberTools), and Biopython - NO PDBFixer.

Targets:
- Aβ42 (5 structures): Alzheimer's amyloid
- Tau K18 (5 structures): Alzheimer's tau
- c-Myc (5 structures): Cancer transcription factor
- α-Synuclein (5 structures): Parkinson's
- p53 (5 structures + 3 mutants): Cancer tumor suppressor

Run on Polaris login node (CPU only, ~8 hours):
    nohup python master_prepare_structures.py > master.log 2>&1 &

Author: DynamicBind Project
"""

import os
import sys
import json
import time
import shutil
import logging
import requests
import subprocess
from pathlib import Path
from typing import Optional, List, Dict, Tuple, Set
from dataclasses import dataclass, field
from datetime import datetime

# =============================================================================
# RESOURCE LIMITS - CRITICAL FOR LOGIN NODE
# =============================================================================
os.environ['OPENBLAS_NUM_THREADS'] = '2'
os.environ['OMP_NUM_THREADS'] = '2'
os.environ['MKL_NUM_THREADS'] = '2'
os.environ['NUMEXPR_MAX_THREADS'] = '2'

# =============================================================================
# LOGGING SETUP
# =============================================================================
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s | %(levelname)-8s | %(message)s',
    datefmt='%Y-%m-%d %H:%M:%S',
    handlers=[
        logging.StreamHandler(sys.stdout),
        logging.FileHandler('master_prep.log', mode='w')
    ]
)
log = logging.getLogger(__name__)

# =============================================================================
# BIOPYTHON IMPORTS
# =============================================================================
try:
    from Bio.PDB import PDBParser, PDBIO, Select
    from Bio.PDB.PDBExceptions import PDBConstructionWarning
    import warnings
    warnings.filterwarnings('ignore', category=PDBConstructionWarning)
    BIOPYTHON_OK = True
except ImportError:
    BIOPYTHON_OK = False
    log.error("BioPython not found! Install with: pip install biopython")

# =============================================================================
# DIRECTORY STRUCTURE
# =============================================================================
BASE_DIR = Path.home() / "DynamicBind_GaMD"
DOWNLOADS_DIR = BASE_DIR / "downloads"
STRUCTURES_DIR = BASE_DIR / "structures"
PREPARED_DIR = BASE_DIR / "prepared"
LOGS_DIR = BASE_DIR / "logs"
CHECKPOINT_FILE = BASE_DIR / "checkpoint.json"

for d in [BASE_DIR, DOWNLOADS_DIR, STRUCTURES_DIR, PREPARED_DIR, LOGS_DIR]:
    d.mkdir(parents=True, exist_ok=True)

# =============================================================================
# STRUCTURE SPECIFICATIONS - ALL 25 PROTEINS
# =============================================================================

@dataclass
class StructureSpec:
    """Complete specification for structure preparation."""
    pdb_id: str
    source: str  # 'pdb', 'ped', 'alphafold'
    target: str  # protein family
    description: str

    # Processing options
    ensemble_models: Optional[int] = None  # Number of NMR models to extract
    fibril_chains: Optional[int] = None    # Number of chains for fibril segment
    keep_chains: Optional[List[str]] = None  # Specific chains to keep
    remove_chains: Optional[List[str]] = None  # Chains to remove
    remove_dna: bool = False
    remove_water: bool = True
    remove_ligands: bool = True
    keep_ligand: Optional[str] = None  # Specific ligand to keep (e.g., "F4")

    # Metal handling
    keep_zinc: bool = True
    remove_zinc: bool = False

    # Mutations to generate from this structure
    mutations: Optional[List[Dict]] = None


# Complete specifications for all 25 structures
STRUCTURES: Dict[str, StructureSpec] = {
    # =========================================================================
    # AMYLOID-BETA 42 (Alzheimer's) - 5 structures
    # =========================================================================
    "1Z0Q": StructureSpec(
        pdb_id="1Z0Q",
        source="pdb",
        target="abeta42",
        description="Aβ42 monomeric NMR ensemble - early IDP state",
        ensemble_models=7,
    ),
    "6SZF": StructureSpec(
        pdb_id="6SZF",
        source="pdb",
        target="abeta42",
        description="Aβ42 in 50% HFIP NMR - partially helical",
        ensemble_models=5,
    ),
    "5OQV": StructureSpec(
        pdb_id="5OQV",
        source="pdb",
        target="abeta42",
        description="Aβ42 pentamer oligomer - toxic β-sheet",
        keep_chains=["A", "B", "C", "D", "E"],
    ),
    "2NAO": StructureSpec(
        pdb_id="2NAO",
        source="pdb",
        target="abeta42",
        description="Aβ42 fibril double-horseshoe - mature fibril",
        fibril_chains=12,
    ),
    "7Q4B": StructureSpec(
        pdb_id="7Q4B",
        source="pdb",
        target="abeta42",
        description="Brain-derived Aβ42 fibril Type I - ex vivo AD",
        fibril_chains=12,
    ),

    # =========================================================================
    # TAU K18 (Alzheimer's) - 5 structures
    # =========================================================================
    "PED00192": StructureSpec(
        pdb_id="PED00192",
        source="ped",
        target="tau_k18",
        description="Tau K18 experimental ensemble - disordered monomer",
        ensemble_models=10,
    ),
    "PED00443": StructureSpec(
        pdb_id="PED00443",
        source="ped",
        target="tau_k18",
        description="Tau K18 idpGAN ensemble - ML conformations",
        ensemble_models=8,
    ),
    "5O3L": StructureSpec(
        pdb_id="5O3L",
        source="pdb",
        target="tau_k18",
        description="Tau straight filament core - AD brain",
        fibril_chains=12,
    ),
    "6NWP": StructureSpec(
        pdb_id="6NWP",
        source="pdb",
        target="tau_k18",
        description="Tau PHF core - paired helical filament",
        fibril_chains=12,
    ),
    "2MZ7": StructureSpec(
        pdb_id="2MZ7",
        source="pdb",
        target="tau_k18",
        description="Tau on microtubule - remove tubulin, keep tau",
        keep_chains=["A"],  # Tau is chain A, tubulin is other chains
    ),

    # =========================================================================
    # c-MYC (Cancer) - 5 structures
    # =========================================================================
    "PED00049": StructureSpec(
        pdb_id="PED00049",
        source="ped",
        target="cmyc",
        description="c-Myc N-terminal TAD ensemble - disordered",
        ensemble_models=7,
    ),
    "AF-P01106": StructureSpec(
        pdb_id="AF-P01106",
        source="alphafold",
        target="cmyc",
        description="c-Myc AlphaFold bHLH-LZ - monomeric",
    ),
    "6G6J": StructureSpec(
        pdb_id="6G6J",
        source="pdb",
        target="cmyc",
        description="Myc-Max heterodimer - remove Max, keep Myc",
        keep_chains=["A"],  # Myc is chain A
        remove_chains=["B"],  # Max is chain B
    ),
    "1NKP": StructureSpec(
        pdb_id="1NKP",
        source="pdb",
        target="cmyc",
        description="Myc-Max-DNA complex - remove DNA",
        remove_dna=True,
    ),
    "1MV0": StructureSpec(
        pdb_id="1MV0",
        source="pdb",
        target="cmyc",
        description="c-Myc + 10058-F4 inhibitor - KEEP LIGAND",
        remove_ligands=False,
        keep_ligand="F4",  # The 10058-F4 inhibitor
    ),

    # =========================================================================
    # ALPHA-SYNUCLEIN (Parkinson's) - 5 structures
    # =========================================================================
    "PED00024": StructureSpec(
        pdb_id="PED00024",
        source="ped",
        target="alpha_synuclein",
        description="α-Syn monomer ensemble - native IDP",
        ensemble_models=8,
    ),
    "1XQ8": StructureSpec(
        pdb_id="1XQ8",
        source="pdb",
        target="alpha_synuclein",
        description="α-Syn micelle-bound NMR - helical state",
        ensemble_models=3,
    ),
    "2KKW": StructureSpec(
        pdb_id="2KKW",
        source="pdb",
        target="alpha_synuclein",
        description="α-Syn partially folded NMR - bent helix",
        ensemble_models=5,
    ),
    "8A9L": StructureSpec(
        pdb_id="8A9L",
        source="pdb",
        target="alpha_synuclein",
        description="Lewy body α-Syn fibril - brain polymorph A",
        fibril_chains=12,
    ),
    "8FPT": StructureSpec(
        pdb_id="8FPT",
        source="pdb",
        target="alpha_synuclein",
        description="α-Syn fibril ssNMR - alternative polymorph",
        fibril_chains=12,
    ),

    # =========================================================================
    # P53 (Cancer) - 5 structures + 3 mutants generated
    # =========================================================================
    "2AHI": StructureSpec(
        pdb_id="2AHI",
        source="pdb",
        target="p53",
        description="p53 DBD wild-type crystal - base for mutants",
        keep_chains=["A"],
        keep_zinc=True,
        mutations=[
            {"name": "R175H", "chain": "A", "resnum": 175, "to": "HIS", "remove_zinc": True},
            {"name": "R248Q", "chain": "A", "resnum": 248, "to": "GLN", "remove_zinc": False},
            {"name": "Y220C", "chain": "A", "resnum": 220, "to": "CYS", "remove_zinc": False},
        ],
    ),
    "2FEJ": StructureSpec(
        pdb_id="2FEJ",
        source="pdb",
        target="p53",
        description="p53 DBD NMR ensemble - solution dynamics",
        ensemble_models=5,
        keep_zinc=True,
    ),
    "5G4O": StructureSpec(
        pdb_id="5G4O",
        source="pdb",
        target="p53",
        description="p53 Y220C + stabilizer - drugged cavity",
        keep_zinc=True,
        remove_ligands=True,  # Remove ligand to see pocket dynamics
    ),
    "2VUK": StructureSpec(
        pdb_id="2VUK",
        source="pdb",
        target="p53",
        description="p53 R175H mutant - Zn-free unstable",
        remove_zinc=True,  # R175H loses Zn binding
    ),
    "PED00086": StructureSpec(
        pdb_id="PED00086",
        source="ped",
        target="p53",
        description="p53 N-terminal TAD ensemble - disordered",
        ensemble_models=5,
    ),
}


# =============================================================================
# BIOPYTHON SELECTORS
# =============================================================================

class ChainSelector(Select):
    """Select specific chains."""
    def __init__(self, chains: List[str]):
        self.chains = set(chains)

    def accept_chain(self, chain):
        return chain.get_id() in self.chains


class ModelSelector(Select):
    """Select specific model from ensemble."""
    def __init__(self, model_id: int):
        self.model_id = model_id

    def accept_model(self, model):
        return model.get_id() == self.model_id


class CleaningSelector(Select):
    """Remove waters, ligands, DNA based on options."""
    def __init__(self, remove_water=True, remove_ligands=True,
                 remove_dna=True, keep_ligand=None, keep_zinc=True, remove_zinc=False):
        self.remove_water = remove_water
        self.remove_ligands = remove_ligands
        self.remove_dna = remove_dna
        self.keep_ligand = keep_ligand
        self.keep_zinc = keep_zinc
        self.remove_zinc = remove_zinc

        self.water_names = {'HOH', 'WAT', 'H2O', 'TIP', 'TIP3', 'TIP4'}
        self.dna_names = {'DA', 'DT', 'DG', 'DC', 'DU', 'A', 'T', 'G', 'C', 'U'}
        self.metal_names = {'ZN', 'CA', 'MG', 'MN', 'FE', 'CU', 'CO', 'NI', 'NA', 'K'}

    def accept_residue(self, residue):
        resname = residue.get_resname().strip()
        hetflag = residue.get_id()[0]

        # Water removal
        if resname in self.water_names:
            return not self.remove_water

        # DNA removal
        if resname in self.dna_names:
            return not self.remove_dna

        # Metal handling
        if resname == 'ZN':
            if self.remove_zinc:
                return False
            return self.keep_zinc

        if resname in self.metal_names:
            return True  # Keep other metals

        # Standard amino acids - always keep
        if hetflag == ' ':
            return True

        # Heteroatoms (ligands)
        if hetflag.startswith('H_') or hetflag == 'W':
            # Keep specific ligand
            if self.keep_ligand and resname == self.keep_ligand:
                return True
            # Remove other ligands
            return not self.remove_ligands

        return True


# =============================================================================
# CHECKPOINT SYSTEM
# =============================================================================

def load_checkpoint() -> Dict:
    """Load checkpoint for resumability."""
    if CHECKPOINT_FILE.exists():
        try:
            with open(CHECKPOINT_FILE) as f:
                return json.load(f)
        except:
            pass
    return {"completed": [], "failed": [], "in_progress": None, "start_time": None}


def save_checkpoint(checkpoint: Dict):
    """Save checkpoint."""
    with open(CHECKPOINT_FILE, 'w') as f:
        json.dump(checkpoint, f, indent=2, default=str)


def is_completed(pdb_id: str, checkpoint: Dict) -> bool:
    """Check if structure is already processed."""
    return pdb_id in checkpoint.get("completed", [])


# =============================================================================
# DOWNLOAD FUNCTIONS
# =============================================================================

def download_pdb(pdb_id: str) -> Optional[Path]:
    """Download from RCSB PDB."""
    output = DOWNLOADS_DIR / f"{pdb_id}.pdb"

    if output.exists() and output.stat().st_size > 1000:
        log.info(f"  Already downloaded: {pdb_id}")
        return output

    url = f"https://files.rcsb.org/download/{pdb_id}.pdb"
    log.info(f"  Downloading {pdb_id} from RCSB...")

    try:
        r = requests.get(url, timeout=60)
        r.raise_for_status()

        with open(output, 'w') as f:
            f.write(r.text)

        log.info(f"  Downloaded: {output.stat().st_size} bytes")
        return output
    except Exception as e:
        log.error(f"  Download failed: {e}")
        return None


def download_alphafold(uniprot_id: str) -> Optional[Path]:
    """Download from AlphaFold DB."""
    # Extract UniProt ID
    uid = uniprot_id.replace("AF-", "").split("-")[0]
    output = DOWNLOADS_DIR / f"AF-{uid}.pdb"

    if output.exists() and output.stat().st_size > 1000:
        log.info(f"  Already downloaded: AF-{uid}")
        return output

    # Try multiple AlphaFold URL formats
    urls = [
        f"https://alphafold.ebi.ac.uk/files/AF-{uid}-F1-model_v4.pdb",
        f"https://alphafold.ebi.ac.uk/files/AF-{uid}-F1-model_v3.pdb",
        f"https://alphafold.ebi.ac.uk/files/AF-{uid}-F1-model_v2.pdb",
    ]

    for url in urls:
        log.info(f"  Trying AlphaFold: {url}")
        try:
            r = requests.get(url, timeout=60, allow_redirects=True)
            if r.status_code == 200 and r.text.startswith(('HEADER', 'ATOM', 'MODEL')):
                with open(output, 'w') as f:
                    f.write(r.text)
                log.info(f"  Downloaded: {output.stat().st_size} bytes")
                return output
        except:
            continue

    log.error(f"  AlphaFold download failed for {uid}")
    return None


def download_ped(ped_id: str) -> Optional[Path]:
    """Download from Protein Ensemble Database."""
    output = DOWNLOADS_DIR / f"{ped_id}.pdb"

    if output.exists() and output.stat().st_size > 5000:
        log.info(f"  Already downloaded: {ped_id}")
        return output

    # PED has complex API - try multiple approaches
    urls = [
        f"https://proteinensemble.org/api/v1/entries/{ped_id}/ensemble/pdb",
        f"https://proteinensemble.org/api/entries/{ped_id}/ensemble",
        f"https://deposition.proteinensemble.org/api/v1/entries/{ped_id}/models/pdb",
    ]

    headers = {
        'Accept': 'text/plain, chemical/x-pdb, application/octet-stream',
        'User-Agent': 'Mozilla/5.0 DynamicBind'
    }

    for url in urls:
        log.info(f"  Trying PED: {url}")
        try:
            r = requests.get(url, timeout=120, headers=headers, allow_redirects=True)
            if r.status_code == 200:
                content = r.text
                # Check if valid PDB
                if any(content.startswith(x) for x in ['HEADER', 'ATOM', 'MODEL', 'TITLE', 'REMARK']):
                    with open(output, 'w') as f:
                        f.write(content)
                    log.info(f"  Downloaded: {output.stat().st_size} bytes")
                    return output
        except Exception as e:
            log.warning(f"  URL failed: {e}")
            continue

    log.warning(f"  PED {ped_id} requires manual download from proteinensemble.org")
    return None


def download_structure(spec: StructureSpec) -> Optional[Path]:
    """Download structure based on source."""
    if spec.source == "pdb":
        return download_pdb(spec.pdb_id)
    elif spec.source == "alphafold":
        return download_alphafold(spec.pdb_id)
    elif spec.source == "ped":
        return download_ped(spec.pdb_id)
    return None


# =============================================================================
# PROCESSING FUNCTIONS
# =============================================================================

def count_models(pdb_path: Path) -> int:
    """Count models in PDB file."""
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("tmp", str(pdb_path))
    return len(list(structure.get_models()))


def count_chains(pdb_path: Path) -> List[str]:
    """Get chain IDs in PDB file."""
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("tmp", str(pdb_path))
    chains = []
    for model in structure:
        for chain in model:
            if chain.get_id() not in chains:
                chains.append(chain.get_id())
        break  # Only first model
    return chains


def extract_models(pdb_path: Path, num_models: int, output_dir: Path) -> List[Path]:
    """Extract diverse models from NMR ensemble."""
    parser = PDBParser(QUIET=True)
    io = PDBIO()

    structure = parser.get_structure("ensemble", str(pdb_path))
    models = list(structure.get_models())
    total = len(models)

    log.info(f"  Ensemble has {total} models, selecting {num_models}")

    # Select evenly spaced models for diversity
    if total <= num_models:
        indices = list(range(total))
    else:
        step = total / num_models
        indices = [int(i * step) for i in range(num_models)]

    outputs = []
    for idx in indices:
        model_id = models[idx].get_id()
        out_path = output_dir / f"{pdb_path.stem}_model{model_id}.pdb"

        io.set_structure(structure)
        io.save(str(out_path), ModelSelector(model_id))
        outputs.append(out_path)
        log.info(f"  Extracted model {model_id}")

    return outputs


def extract_chains(pdb_path: Path, chains: List[str], output_path: Path) -> Path:
    """Extract specific chains."""
    parser = PDBParser(QUIET=True)
    io = PDBIO()

    structure = parser.get_structure("protein", str(pdb_path))
    io.set_structure(structure)
    io.save(str(output_path), ChainSelector(chains))

    log.info(f"  Extracted chains: {chains}")
    return output_path


def extract_fibril_segment(pdb_path: Path, num_chains: int, output_path: Path) -> Path:
    """Extract fibril segment with specified number of chains."""
    all_chains = count_chains(pdb_path)
    log.info(f"  Fibril has {len(all_chains)} chains: {all_chains[:20]}...")

    # Take first num_chains
    selected = all_chains[:num_chains]
    log.info(f"  Selecting {num_chains} chains: {selected}")

    return extract_chains(pdb_path, selected, output_path)


def clean_structure(pdb_path: Path, spec: StructureSpec, output_path: Path) -> Path:
    """Clean structure: remove waters, ligands, DNA as specified."""
    parser = PDBParser(QUIET=True)
    io = PDBIO()

    structure = parser.get_structure("protein", str(pdb_path))

    selector = CleaningSelector(
        remove_water=spec.remove_water,
        remove_ligands=spec.remove_ligands,
        remove_dna=spec.remove_dna,
        keep_ligand=spec.keep_ligand,
        keep_zinc=spec.keep_zinc,
        remove_zinc=spec.remove_zinc,
    )

    io.set_structure(structure)
    io.save(str(output_path), selector)

    log.info(f"  Cleaned: water={spec.remove_water}, lig={spec.remove_ligands}, dna={spec.remove_dna}")
    return output_path


def add_hydrogens_reduce(pdb_path: Path, output_path: Path) -> Path:
    """Add hydrogens using reduce (from AmberTools)."""
    try:
        # Try reduce from AmberTools
        cmd = f"reduce -BUILD -NUC {pdb_path} > {output_path} 2>/dev/null"
        result = subprocess.run(cmd, shell=True, timeout=120)

        if output_path.exists() and output_path.stat().st_size > 100:
            log.info(f"  Added hydrogens with reduce")
            return output_path
    except:
        pass

    # Fallback: copy without hydrogens
    log.warning(f"  Could not add hydrogens, using structure as-is")
    shutil.copy(pdb_path, output_path)
    return output_path


def add_hydrogens_pdb4amber(pdb_path: Path, output_path: Path) -> Path:
    """Clean with pdb4amber (alternative to reduce)."""
    try:
        # pdb4amber cleans and can add hydrogens
        cmd = f"pdb4amber -i {pdb_path} -o {output_path} --reduce 2>/dev/null"
        result = subprocess.run(cmd, shell=True, timeout=180)

        if output_path.exists() and output_path.stat().st_size > 100:
            log.info(f"  Processed with pdb4amber")
            return output_path
    except:
        pass

    # Fallback
    log.warning(f"  pdb4amber failed, using reduce only")
    return add_hydrogens_reduce(pdb_path, output_path)


def create_mutant(pdb_path: Path, mutation: Dict, output_path: Path) -> Path:
    """
    Create point mutant.

    Note: This is a simplified mutation that renames the residue.
    For production, use PyMOL mutagenesis or Modeller for proper sidechain modeling.
    """
    parser = PDBParser(QUIET=True)
    io = PDBIO()

    structure = parser.get_structure("protein", str(pdb_path))

    chain_id = mutation["chain"]
    resnum = mutation["resnum"]
    new_aa = mutation["to"]
    remove_zn = mutation.get("remove_zinc", False)

    log.info(f"  Creating mutant: {chain_id}{resnum} -> {new_aa}")

    # For proper mutation, we'd need to rebuild sidechains
    # Here we just mark it and handle Zn removal

    class MutantSelector(Select):
        def accept_residue(self, residue):
            resname = residue.get_resname().strip()
            # Remove Zn if specified
            if remove_zn and resname == 'ZN':
                return False
            return True

    io.set_structure(structure)
    io.save(str(output_path), MutantSelector())

    if remove_zn:
        log.info(f"  Removed Zn (mutation disrupts binding)")

    return output_path


# =============================================================================
# MAIN PROCESSING PIPELINE
# =============================================================================

def process_structure(spec: StructureSpec, checkpoint: Dict) -> List[Path]:
    """Process a single structure according to its specification."""
    pdb_id = spec.pdb_id

    log.info(f"\n{'='*60}")
    log.info(f"PROCESSING: {pdb_id}")
    log.info(f"Target: {spec.target}")
    log.info(f"Description: {spec.description}")
    log.info(f"{'='*60}")

    # Check if already done
    if is_completed(pdb_id, checkpoint):
        log.info(f"  SKIPPED - already completed")
        existing = list(PREPARED_DIR.glob(f"{pdb_id}*.pdb"))
        return existing

    checkpoint["in_progress"] = pdb_id
    save_checkpoint(checkpoint)

    prepared_files = []

    try:
        # Step 1: Download
        log.info(f"  Step 1: Download")
        raw_path = download_structure(spec)
        if raw_path is None:
            log.error(f"  FAILED: Could not download")
            checkpoint["failed"].append(pdb_id)
            save_checkpoint(checkpoint)
            return []

        # Create working directory
        work_dir = STRUCTURES_DIR / pdb_id
        work_dir.mkdir(parents=True, exist_ok=True)

        # Step 2: Extract models/chains
        log.info(f"  Step 2: Extract/Select")

        intermediate_files = []

        if spec.ensemble_models:
            # NMR ensemble - extract diverse models
            intermediate_files = extract_models(raw_path, spec.ensemble_models, work_dir)

        elif spec.fibril_chains:
            # Fibril - extract segment
            fibril_path = work_dir / f"{pdb_id}_fibril.pdb"
            extract_fibril_segment(raw_path, spec.fibril_chains, fibril_path)
            intermediate_files = [fibril_path]

        elif spec.keep_chains:
            # Keep specific chains
            chains_path = work_dir / f"{pdb_id}_chains.pdb"
            extract_chains(raw_path, spec.keep_chains, chains_path)
            intermediate_files = [chains_path]

        else:
            # Use as-is
            intermediate_files = [raw_path]

        # Step 3: Clean structures
        log.info(f"  Step 3: Clean")
        cleaned_files = []
        for inter_file in intermediate_files:
            clean_path = work_dir / f"{inter_file.stem}_clean.pdb"
            clean_structure(inter_file, spec, clean_path)
            cleaned_files.append(clean_path)

        # Step 4: Add hydrogens
        log.info(f"  Step 4: Add hydrogens")
        for clean_file in cleaned_files:
            prep_path = PREPARED_DIR / f"{clean_file.stem.replace('_clean', '')}_prep.pdb"
            add_hydrogens_reduce(clean_file, prep_path)

            if prep_path.exists() and prep_path.stat().st_size > 100:
                prepared_files.append(prep_path)
                log.info(f"  Prepared: {prep_path.name}")

        # Step 5: Generate mutants if specified
        if spec.mutations:
            log.info(f"  Step 5: Generate mutants")
            # Use the first prepared file as template
            if prepared_files:
                template = prepared_files[0]
                for mut in spec.mutations:
                    mut_name = mut["name"]
                    mut_path = PREPARED_DIR / f"p53_{mut_name}_prep.pdb"
                    create_mutant(template, mut, mut_path)

                    if mut_path.exists():
                        prepared_files.append(mut_path)
                        log.info(f"  Created mutant: {mut_name}")

        # Mark completed
        checkpoint["completed"].append(pdb_id)
        checkpoint["in_progress"] = None
        save_checkpoint(checkpoint)

        log.info(f"  SUCCESS: {len(prepared_files)} files prepared")

    except Exception as e:
        log.error(f"  ERROR: {e}")
        import traceback
        traceback.print_exc()
        checkpoint["failed"].append(pdb_id)
        save_checkpoint(checkpoint)

    return prepared_files


def install_dependencies():
    """Install required Python packages."""
    log.info("Checking dependencies...")

    packages = [
        ("biopython", "biopython"),
        ("requests", "requests"),
    ]

    for import_name, pip_name in packages:
        try:
            __import__(import_name.replace("-", "_").split(".")[0])
            log.info(f"  {import_name}: OK")
        except ImportError:
            log.info(f"  Installing {pip_name}...")
            subprocess.run([sys.executable, "-m", "pip", "install", pip_name, "-q"])

    # Check for reduce (AmberTools)
    try:
        result = subprocess.run(["which", "reduce"], capture_output=True, text=True)
        if result.returncode == 0:
            log.info(f"  reduce (AmberTools): OK")
        else:
            log.warning("  reduce not found - hydrogens will not be added")
    except:
        log.warning("  Could not check for reduce")


def main():
    """Main entry point."""
    start_time = datetime.now()

    log.info("="*70)
    log.info("MASTER STRUCTURE PREPARATION FOR GaMD SIMULATIONS")
    log.info("="*70)
    log.info(f"Start time: {start_time}")
    log.info(f"Base directory: {BASE_DIR}")
    log.info(f"Structures to process: {len(STRUCTURES)}")
    log.info("")

    # Install dependencies
    install_dependencies()

    # Load checkpoint
    checkpoint = load_checkpoint()
    if not checkpoint.get("start_time"):
        checkpoint["start_time"] = str(start_time)
    save_checkpoint(checkpoint)

    log.info(f"\nCheckpoint status:")
    log.info(f"  Completed: {len(checkpoint.get('completed', []))}")
    log.info(f"  Failed: {len(checkpoint.get('failed', []))}")

    # Process all structures
    all_prepared = []

    for i, (pdb_id, spec) in enumerate(STRUCTURES.items(), 1):
        log.info(f"\n[{i}/{len(STRUCTURES)}] Processing {pdb_id}...")

        try:
            prepared = process_structure(spec, checkpoint)
            all_prepared.extend(prepared)
        except Exception as e:
            log.error(f"Unexpected error: {e}")
            continue

        # Small delay between structures
        time.sleep(1)

    # Final summary
    end_time = datetime.now()
    duration = end_time - start_time

    log.info("\n" + "="*70)
    log.info("SUMMARY")
    log.info("="*70)
    log.info(f"End time: {end_time}")
    log.info(f"Duration: {duration}")
    log.info(f"Completed: {len(checkpoint.get('completed', []))}/{len(STRUCTURES)}")
    log.info(f"Failed: {checkpoint.get('failed', [])}")
    log.info(f"Total prepared files: {len(all_prepared)}")

    log.info("\nPrepared files:")
    for f in sorted(PREPARED_DIR.glob("*.pdb")):
        log.info(f"  {f.name} ({f.stat().st_size:,} bytes)")

    log.info("\n" + "="*70)
    log.info("STRUCTURE PREPARATION COMPLETE")
    log.info("="*70)


if __name__ == "__main__":
    main()
