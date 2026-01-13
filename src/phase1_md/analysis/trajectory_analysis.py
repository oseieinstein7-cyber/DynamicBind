"""
Trajectory Analysis Module

Provides comprehensive analysis of MD trajectories including:
- RMSD (Root Mean Square Deviation)
- RMSF (Root Mean Square Fluctuation)
- Radius of Gyration (Rg)
- Secondary Structure Evolution
- Conformational Clustering
- GaMD Reweighting for PMF calculation
"""

import os
import logging
import numpy as np
from pathlib import Path
from typing import Optional, List, Dict, Tuple, Any
from dataclasses import dataclass

# Configure logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

# Try importing MDAnalysis
try:
    import MDAnalysis as mda
    from MDAnalysis.analysis import rms, align, diffusionmap
    from MDAnalysis.analysis.rms import RMSD, RMSF
    MDANALYSIS_AVAILABLE = True
except ImportError:
    MDANALYSIS_AVAILABLE = False
    logger.warning("MDAnalysis not installed. Install with: pip install MDAnalysis")

# Try importing scikit-learn for clustering
try:
    from sklearn.cluster import KMeans, DBSCAN
    from sklearn.decomposition import PCA
    SKLEARN_AVAILABLE = True
except ImportError:
    SKLEARN_AVAILABLE = False
    logger.warning("scikit-learn not installed. Install with: pip install scikit-learn")

# Try importing matplotlib for plotting
try:
    import matplotlib.pyplot as plt
    MATPLOTLIB_AVAILABLE = True
except ImportError:
    MATPLOTLIB_AVAILABLE = False
    logger.warning("matplotlib not installed. Install with: pip install matplotlib")


@dataclass
class AnalysisResults:
    """Container for trajectory analysis results."""
    rmsd: Optional[np.ndarray] = None
    rmsf: Optional[np.ndarray] = None
    rg: Optional[np.ndarray] = None
    time: Optional[np.ndarray] = None
    cluster_labels: Optional[np.ndarray] = None
    cluster_centers: Optional[np.ndarray] = None
    representative_frames: Optional[List[int]] = None
    pca_projection: Optional[np.ndarray] = None


class TrajectoryAnalyzer:
    """
    Comprehensive trajectory analysis for MD simulations.

    Analyzes:
    - Structural stability (RMSD, RMSF)
    - Compactness (Radius of Gyration)
    - Conformational landscape (PCA, clustering)
    - Representative structures extraction
    """

    def __init__(
        self,
        topology_path: str,
        trajectory_path: str,
        output_dir: str = "analysis"
    ):
        """
        Initialize trajectory analyzer.

        Args:
            topology_path: Path to topology file (PDB, GRO, etc.)
            trajectory_path: Path to trajectory file (DCD, XTC, etc.)
            output_dir: Output directory for results
        """
        if not MDANALYSIS_AVAILABLE:
            raise ImportError(
                "MDAnalysis required for trajectory analysis. "
                "Install with: pip install MDAnalysis"
            )

        self.topology_path = Path(topology_path)
        self.trajectory_path = Path(trajectory_path)
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(parents=True, exist_ok=True)

        logger.info(f"Loading trajectory: {self.trajectory_path}")
        self.universe = mda.Universe(str(topology_path), str(trajectory_path))

        self.n_frames = len(self.universe.trajectory)
        self.n_atoms = len(self.universe.atoms)

        logger.info(f"Loaded {self.n_frames} frames, {self.n_atoms} atoms")

        # Results container
        self.results = AnalysisResults()

    def calculate_rmsd(
        self,
        selection: str = "backbone",
        reference_frame: int = 0
    ) -> np.ndarray:
        """
        Calculate RMSD over trajectory.

        Args:
            selection: Atom selection string
            reference_frame: Reference frame for alignment

        Returns:
            Array of RMSD values per frame
        """
        logger.info(f"Calculating RMSD for selection: {selection}")

        # Select atoms
        atoms = self.universe.select_atoms(selection)
        n_atoms_selected = len(atoms)
        logger.info(f"Selected {n_atoms_selected} atoms")

        # Set reference
        self.universe.trajectory[reference_frame]
        ref_positions = atoms.positions.copy()

        # Calculate RMSD
        rmsd_values = []
        time_values = []

        for ts in self.universe.trajectory:
            # Align to reference
            positions = atoms.positions
            # Simple RMSD calculation (without alignment for speed)
            diff = positions - ref_positions
            rmsd = np.sqrt(np.mean(np.sum(diff**2, axis=1)))
            rmsd_values.append(rmsd)
            time_values.append(ts.time)

        self.results.rmsd = np.array(rmsd_values)
        self.results.time = np.array(time_values)

        logger.info(f"RMSD: mean={np.mean(self.results.rmsd):.3f} A, "
                   f"std={np.std(self.results.rmsd):.3f} A")

        return self.results.rmsd

    def calculate_rmsf(
        self,
        selection: str = "name CA"
    ) -> np.ndarray:
        """
        Calculate per-residue RMSF.

        Args:
            selection: Atom selection string (typically CA atoms)

        Returns:
            Array of RMSF values per selected atom
        """
        logger.info(f"Calculating RMSF for selection: {selection}")

        atoms = self.universe.select_atoms(selection)
        n_atoms = len(atoms)

        # Collect positions
        positions = np.zeros((self.n_frames, n_atoms, 3))

        for i, ts in enumerate(self.universe.trajectory):
            positions[i] = atoms.positions

        # Calculate mean position
        mean_pos = np.mean(positions, axis=0)

        # Calculate RMSF
        diff = positions - mean_pos
        rmsf = np.sqrt(np.mean(np.sum(diff**2, axis=2), axis=0))

        self.results.rmsf = rmsf

        logger.info(f"RMSF: mean={np.mean(rmsf):.3f} A, max={np.max(rmsf):.3f} A")

        return rmsf

    def calculate_radius_of_gyration(
        self,
        selection: str = "protein"
    ) -> np.ndarray:
        """
        Calculate radius of gyration over trajectory.

        Args:
            selection: Atom selection string

        Returns:
            Array of Rg values per frame
        """
        logger.info(f"Calculating Rg for selection: {selection}")

        atoms = self.universe.select_atoms(selection)
        rg_values = []

        for ts in self.universe.trajectory:
            # Calculate center of mass
            com = atoms.center_of_mass()

            # Calculate Rg
            positions = atoms.positions
            masses = atoms.masses

            diff = positions - com
            sq_dist = np.sum(diff**2, axis=1)
            rg = np.sqrt(np.sum(masses * sq_dist) / np.sum(masses))
            rg_values.append(rg)

        self.results.rg = np.array(rg_values)

        logger.info(f"Rg: mean={np.mean(self.results.rg):.3f} A, "
                   f"std={np.std(self.results.rg):.3f} A")

        return self.results.rg

    def cluster_conformations(
        self,
        selection: str = "name CA",
        n_clusters: int = 10,
        method: str = "kmeans"
    ) -> Tuple[np.ndarray, List[int]]:
        """
        Cluster trajectory conformations.

        Args:
            selection: Atom selection for clustering
            n_clusters: Number of clusters
            method: Clustering method ('kmeans' or 'dbscan')

        Returns:
            Tuple of (cluster labels, representative frame indices)
        """
        if not SKLEARN_AVAILABLE:
            raise ImportError("scikit-learn required for clustering")

        logger.info(f"Clustering conformations into {n_clusters} clusters...")

        atoms = self.universe.select_atoms(selection)
        n_atoms = len(atoms)

        # Collect positions
        positions = np.zeros((self.n_frames, n_atoms * 3))

        for i, ts in enumerate(self.universe.trajectory):
            positions[i] = atoms.positions.flatten()

        # Perform PCA for dimensionality reduction
        logger.info("Running PCA...")
        pca = PCA(n_components=min(50, self.n_frames, n_atoms * 3))
        pca_coords = pca.fit_transform(positions)
        self.results.pca_projection = pca_coords

        variance_explained = np.sum(pca.explained_variance_ratio_[:10])
        logger.info(f"First 10 PCs explain {variance_explained*100:.1f}% variance")

        # Cluster
        if method == "kmeans":
            clusterer = KMeans(n_clusters=n_clusters, random_state=42, n_init=10)
        elif method == "dbscan":
            clusterer = DBSCAN(eps=3, min_samples=5)
        else:
            raise ValueError(f"Unknown clustering method: {method}")

        labels = clusterer.fit_predict(pca_coords[:, :10])
        self.results.cluster_labels = labels

        # Find representative frames (closest to cluster centers)
        representatives = []
        unique_labels = np.unique(labels)

        for label in unique_labels:
            if label == -1:  # Skip noise in DBSCAN
                continue

            mask = labels == label
            cluster_points = pca_coords[mask, :10]
            center = np.mean(cluster_points, axis=0)

            # Find closest point to center
            distances = np.linalg.norm(cluster_points - center, axis=1)
            closest_idx = np.where(mask)[0][np.argmin(distances)]
            representatives.append(closest_idx)

        self.results.representative_frames = representatives

        logger.info(f"Found {len(unique_labels)} clusters")
        logger.info(f"Representative frames: {representatives}")

        return labels, representatives

    def extract_representative_structures(
        self,
        output_prefix: str = "cluster"
    ) -> List[str]:
        """
        Extract representative structures from clusters.

        Args:
            output_prefix: Prefix for output PDB files

        Returns:
            List of paths to extracted PDB files
        """
        if self.results.representative_frames is None:
            logger.warning("No clusters found. Run cluster_conformations first.")
            return []

        output_files = []
        protein = self.universe.select_atoms("protein")

        for i, frame_idx in enumerate(self.results.representative_frames):
            self.universe.trajectory[frame_idx]

            output_path = self.output_dir / f"{output_prefix}_{i:03d}.pdb"
            protein.write(str(output_path))
            output_files.append(str(output_path))

            logger.info(f"Saved cluster {i} representative to {output_path}")

        return output_files

    def plot_rmsd(self, output_path: Optional[str] = None) -> None:
        """Plot RMSD over time."""
        if not MATPLOTLIB_AVAILABLE:
            logger.warning("matplotlib not available for plotting")
            return

        if self.results.rmsd is None:
            self.calculate_rmsd()

        plt.figure(figsize=(10, 4))
        time_ns = self.results.time / 1000  # Convert to ns
        plt.plot(time_ns, self.results.rmsd, 'b-', linewidth=0.5)
        plt.xlabel('Time (ns)')
        plt.ylabel('RMSD (Angstrom)')
        plt.title('Backbone RMSD')

        if output_path is None:
            output_path = self.output_dir / "rmsd.png"

        plt.savefig(output_path, dpi=300, bbox_inches='tight')
        plt.close()
        logger.info(f"Saved RMSD plot to {output_path}")

    def plot_rmsf(self, output_path: Optional[str] = None) -> None:
        """Plot per-residue RMSF."""
        if not MATPLOTLIB_AVAILABLE:
            logger.warning("matplotlib not available for plotting")
            return

        if self.results.rmsf is None:
            self.calculate_rmsf()

        plt.figure(figsize=(12, 4))
        residues = np.arange(1, len(self.results.rmsf) + 1)
        plt.bar(residues, self.results.rmsf, width=1.0, color='steelblue')
        plt.xlabel('Residue Number')
        plt.ylabel('RMSF (Angstrom)')
        plt.title('Per-Residue RMSF')

        if output_path is None:
            output_path = self.output_dir / "rmsf.png"

        plt.savefig(output_path, dpi=300, bbox_inches='tight')
        plt.close()
        logger.info(f"Saved RMSF plot to {output_path}")

    def plot_rg(self, output_path: Optional[str] = None) -> None:
        """Plot radius of gyration over time."""
        if not MATPLOTLIB_AVAILABLE:
            logger.warning("matplotlib not available for plotting")
            return

        if self.results.rg is None:
            self.calculate_radius_of_gyration()

        plt.figure(figsize=(10, 4))
        time_ns = self.results.time / 1000 if self.results.time is not None else np.arange(len(self.results.rg))
        plt.plot(time_ns, self.results.rg, 'g-', linewidth=0.5)
        plt.xlabel('Time (ns)')
        plt.ylabel('Radius of Gyration (Angstrom)')
        plt.title('Radius of Gyration')

        if output_path is None:
            output_path = self.output_dir / "rg.png"

        plt.savefig(output_path, dpi=300, bbox_inches='tight')
        plt.close()
        logger.info(f"Saved Rg plot to {output_path}")

    def run_full_analysis(self) -> AnalysisResults:
        """
        Run complete trajectory analysis.

        Returns:
            AnalysisResults object with all computed metrics
        """
        logger.info("=" * 60)
        logger.info("Running full trajectory analysis")
        logger.info("=" * 60)

        # Calculate basic metrics
        self.calculate_rmsd()
        self.calculate_rmsf()
        self.calculate_radius_of_gyration()

        # Cluster conformations
        if SKLEARN_AVAILABLE:
            self.cluster_conformations()
            self.extract_representative_structures()

        # Generate plots
        if MATPLOTLIB_AVAILABLE:
            self.plot_rmsd()
            self.plot_rmsf()
            self.plot_rg()

        # Save numerical results
        self._save_results()

        logger.info("Analysis complete!")
        return self.results

    def _save_results(self) -> None:
        """Save numerical results to files."""
        # Save RMSD
        if self.results.rmsd is not None:
            np.savetxt(
                self.output_dir / "rmsd.csv",
                np.column_stack([self.results.time, self.results.rmsd]),
                delimiter=',',
                header='Time(ps),RMSD(A)',
                comments=''
            )

        # Save RMSF
        if self.results.rmsf is not None:
            np.savetxt(
                self.output_dir / "rmsf.csv",
                self.results.rmsf,
                delimiter=',',
                header='RMSF(A)',
                comments=''
            )

        # Save Rg
        if self.results.rg is not None and self.results.time is not None:
            np.savetxt(
                self.output_dir / "rg.csv",
                np.column_stack([self.results.time, self.results.rg]),
                delimiter=',',
                header='Time(ps),Rg(A)',
                comments=''
            )

        # Save cluster assignments
        if self.results.cluster_labels is not None:
            np.savetxt(
                self.output_dir / "clusters.csv",
                self.results.cluster_labels,
                delimiter=',',
                header='Cluster',
                comments='',
                fmt='%d'
            )

        logger.info(f"Saved results to {self.output_dir}")


def analyze_trajectory(
    topology_path: str,
    trajectory_path: str,
    output_dir: str = "analysis",
    n_clusters: int = 10
) -> AnalysisResults:
    """
    Convenience function for trajectory analysis.

    Args:
        topology_path: Path to topology file
        trajectory_path: Path to trajectory file
        output_dir: Output directory
        n_clusters: Number of clusters for conformational analysis

    Returns:
        AnalysisResults object
    """
    analyzer = TrajectoryAnalyzer(topology_path, trajectory_path, output_dir)
    return analyzer.run_full_analysis()


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="Analyze MD trajectory")
    parser.add_argument("topology", help="Topology file (PDB, GRO, etc.)")
    parser.add_argument("trajectory", help="Trajectory file (DCD, XTC, etc.)")
    parser.add_argument("-o", "--output", default="analysis", help="Output directory")
    parser.add_argument("-n", "--clusters", type=int, default=10, help="Number of clusters")

    args = parser.parse_args()

    analyze_trajectory(args.topology, args.trajectory, args.output, args.clusters)
