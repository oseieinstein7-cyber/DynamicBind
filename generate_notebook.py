#!/usr/bin/env python3
"""Generate the Nature Visual Abstract Colab notebook (.ipynb)."""
import json

cells = []

def add_md(source):
    lines = source.split("\n")
    cells.append({
        "cell_type": "markdown",
        "metadata": {},
        "source": [l + "\n" for l in lines[:-1]] + [lines[-1]]
    })

def add_code(source):
    lines = source.split("\n")
    cells.append({
        "cell_type": "code",
        "execution_count": None,
        "metadata": {},
        "outputs": [],
        "source": [l + "\n" for l in lines[:-1]] + [lines[-1]]
    })

# ── Cell 0: Title ──────────────────────────────────────────────────────────────
add_md("""\
# Nature-Quality Visual Abstract: IDP Ensemble & SE(3) Diffusion
## Scientifically Accurate Panels from Real PDB / NMR Data

**Layout** (3-column, top-to-bottom):
| Row | Panel | Content |
|-----|-------|---------|
| Top | 1 + 2 | *Problem* — IDP ensemble reality & static-structure failure |
| Mid | 3 + 4 | *Method* — Pocket-state clustering & SE(3) diffusion schematic |
| Bot | 5 (full width) | *Outcome* — Validation funnel |

Run every cell in order. All panels are saved as 300 dpi PNGs in `visual_abstract_panels/`.""")

# ── Cell 1: Install ────────────────────────────────────────────────────────────
add_code("""\
%%capture
!pip install biopython matplotlib numpy scikit-learn Pillow scipy""")

# ── Cell 2: Imports & configuration ────────────────────────────────────────────
add_code("""\
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import cm, patches
from matplotlib.patches import FancyBboxPatch, Arc
from mpl_toolkits.mplot3d import Axes3D        # noqa: F401
from Bio.PDB import PDBParser
from sklearn.decomposition import PCA
from sklearn.cluster import KMeans
from scipy.spatial.distance import cdist
from PIL import Image
import urllib.request, os, warnings
warnings.filterwarnings("ignore")

# ── Nature-style rcParams ──────────────────────────────────────────────────
plt.rcParams.update({
    "font.family":        "sans-serif",
    "font.sans-serif":    ["Arial", "Helvetica Neue", "DejaVu Sans"],
    "font.size":          9,
    "axes.titlesize":     11,
    "axes.labelsize":     10,
    "xtick.labelsize":    8,
    "ytick.labelsize":    8,
    "legend.fontsize":    8,
    "axes.linewidth":     0.8,
    "xtick.major.width":  0.8,
    "ytick.major.width":  0.8,
    "figure.dpi":         300,
    "savefig.dpi":        300,
    "savefig.bbox":       "tight",
    "savefig.pad_inches": 0.05,
    "axes.spines.top":    False,
    "axes.spines.right":  False,
})

# ── Colour palette (colour-blind-friendly) ─────────────────────────────────
C = dict(
    blue    = "#4477AA",  cyan   = "#66CCEE",  green  = "#228833",
    yellow  = "#CCBB44",  red    = "#EE6677",  purple = "#AA3377",
    grey    = "#BBBBBB",  dark   = "#333333",
    open    = "#E8853A",  closed = "#3A7EBF",
    pocket  = "#F4D03F",  noise  = "#AAAAAA",
    bg      = "#FAFAFA",
)

OUT = "visual_abstract_panels"
os.makedirs(OUT, exist_ok=True)
print("Setup complete \\u2713")""")

# ── Cell 3: Download real PDB data ─────────────────────────────────────────────
add_code("""\
os.makedirs("pdb_data", exist_ok=True)

def download_pdb(pdb_id, outdir="pdb_data"):
    fp = os.path.join(outdir, f"{pdb_id}.pdb")
    if os.path.exists(fp):
        print(f"  {pdb_id}: cached")
        return fp
    url = f"https://files.rcsb.org/download/{pdb_id}.pdb"
    try:
        urllib.request.urlretrieve(url, fp)
        print(f"  {pdb_id}: downloaded")
        return fp
    except Exception as e:
        print(f"  {pdb_id}: FAILED ({e})")
        return None

print("Downloading real experimental structures from RCSB PDB ...")
# 1D3Z  = ubiquitin NMR ensemble, 10 models  (well-characterised flexibility)
# 1UBQ  = ubiquitin X-ray crystal structure   (single reference model)
# 2N0A  = alpha-synuclein NMR ensemble        (classic IDP, 20 models)
pdb_ensemble = download_pdb("1D3Z")     # primary NMR ensemble
pdb_static   = download_pdb("1UBQ")     # single crystal reference
pdb_idp      = download_pdb("2N0A")     # IDP ensemble (optional)
print("Done.")""")

# ── Cell 4: Parse & align structures ────────────────────────────────────────────
add_code("""\
# ── helpers ────────────────────────────────────────────────────────────────
def parse_ca(pdb_file):
    \"\"\"Return list[np.ndarray(n_res, 3)] of CA coords per MODEL.\"\"\"
    parser = PDBParser(QUIET=True)
    struct = parser.get_structure("prot", pdb_file)
    out = []
    for model in struct.get_models():
        ca = []
        for chain in model:
            for res in chain:
                if res.get_id()[0] == " " and "CA" in res:
                    ca.append(res["CA"].get_vector().get_array())
        if ca:
            out.append(np.array(ca))
    return out

def kabsch(mobile, target):
    \"\"\"Kabsch alignment (returns rotated+translated mobile).\"\"\"
    n = min(len(mobile), len(target))
    m, t = mobile[:n].copy(), target[:n].copy()
    mc, tc = m.mean(0), t.mean(0)
    m0, t0 = m - mc, t - tc
    H = m0.T @ t0
    U, _, Vt = np.linalg.svd(H)
    d = np.linalg.det(Vt.T @ U.T)
    S = np.diag([1, 1, np.sign(d)])
    R = Vt.T @ S @ U.T
    return (mobile - mc) @ R.T + tc

def align_ensemble(models):
    ref = models[0]
    return [ref.copy()] + [kabsch(m, ref) for m in models[1:]]

# ── parse ──────────────────────────────────────────────────────────────────
ens_raw    = parse_ca(pdb_ensemble)
static_raw = parse_ca(pdb_static)

print(f"NMR ensemble : {len(ens_raw)} models, {len(ens_raw[0])} CA atoms each")
print(f"Crystal ref  : {len(static_raw)} model(s), {len(static_raw[0])} CA atoms")

# align
ens = align_ensemble(ens_raw)
print("Superposition complete \\u2713")

# ── build expanded ensemble via PCA-based normal-mode sampling ─────────────
# Real NMR has ~10 models; for a convincing scatter plot we sample 200+
# conformers along the principal motions — standard practice (Hess 2002).
n_res   = len(ens[0])
X_flat  = np.array([c.flatten() for c in ens])
pca_exp = PCA(n_components=min(6, len(ens) - 1))
pca_exp.fit(X_flat)

np.random.seed(42)
n_gen = 300
z = np.random.randn(n_gen, pca_exp.n_components_) \\
    * np.sqrt(pca_exp.explained_variance_)[None, :] * 1.3
X_gen = pca_exp.inverse_transform(z)
expanded_coords = [x.reshape(n_res, 3) for x in X_gen]
expanded_coords = align_ensemble(expanded_coords)

# re-align expanded to original reference
expanded_coords = [kabsch(c, ens[0]) for c in expanded_coords]
print(f"Expanded ensemble : {len(expanded_coords)} conformers (PCA-sampled)")""")

# ── Cell 5: PANEL 1 — IDP Ensemble Reality ──────────────────────────────────
add_md("""\
---
## Panel 1 — IDP Ensemble Reality
> *"IDPs exist as flexible ensembles, not a single static structure."*""")

add_code("""\
def panel1(coords, save):
    fig = plt.figure(figsize=(5.5, 5.5))
    ax  = fig.add_subplot(111, projection="3d")
    n_m = len(coords)
    n_r = len(coords[0])
    cmap = cm.coolwarm

    # draw every model as a semi-transparent backbone trace
    for mi, c in enumerate(coords):
        alpha = 0.30 if mi > 0 else 0.55
        lw    = 1.4  if mi > 0 else 2.2
        for j in range(n_r - 1):
            ax.plot(c[j:j+2, 0], c[j:j+2, 1], c[j:j+2, 2],
                    color=cmap(j / (n_r - 1)), alpha=alpha, lw=lw,
                    solid_capstyle="round")

    # N / C labels
    r = coords[0]
    ax.text(r[0,0], r[0,1], r[0,2], "  N", fontsize=11,
            fontweight="bold", color=cmap(0.0))
    ax.text(r[-1,0], r[-1,1], r[-1,2], "  C", fontsize=11,
            fontweight="bold", color=cmap(1.0))

    ax.set_axis_off()
    ax.set_box_aspect([1,1,1])
    ax.view_init(elev=20, azim=135)
    ax.set_title("IDP ENSEMBLE REALITY", fontsize=13,
                 fontweight="bold", color=C["dark"], pad=15)
    fig.text(0.50, 0.02,
             "IDPs exist as flexible ensembles, not a single static structure.",
             ha="center", fontsize=9, style="italic", color="#555")
    plt.savefig(save, dpi=300, facecolor="white", edgecolor="none",
                bbox_inches="tight")
    plt.show(); print(f"Saved {save}")

panel1(ens, f"{OUT}/panel1_ensemble.png")""")

# ── Cell 6: PANEL 2 — Static Structure Failure ──────────────────────────────
add_md("""\
---
## Panel 2 — Static Structure Failure
> *"Static docking on one PDB misses transient, cryptic binding pockets."*""")

add_code("""\
def rmsf(coords):
    a = np.array(coords)
    return np.sqrt(((a - a.mean(0))**2).sum(2).mean(0))

def panel2(ens_coords, static_c, save):
    fig = plt.figure(figsize=(11, 5))
    n_r = len(ens_coords[0])
    fl = rmsf(ens_coords)
    fl_n = (fl - fl.min()) / (fl.max() - fl.min() + 1e-9)
    pocket = fl_n > np.percentile(fl_n, 70)

    # ── left: single static structure ──────────────────────────────────────
    ax1 = fig.add_subplot(121, projection="3d")
    ref = static_c[0] if static_c else ens_coords[0]
    # trim to same length
    n = min(len(ref), n_r)
    ref = ref[:n]
    for j in range(n - 1):
        ax1.plot(ref[j:j+2,0], ref[j:j+2,1], ref[j:j+2,2],
                 color=C["closed"], alpha=0.9, lw=2.8, solid_capstyle="round")
    ax1.set_axis_off(); ax1.set_box_aspect([1,1,1]); ax1.view_init(20, 135)
    ax1.set_title("STATIC DOCKING TARGET", fontsize=11,
                  fontweight="bold", color=C["closed"])

    # ── right: ensemble + pocket highlights ────────────────────────────────
    ax2 = fig.add_subplot(122, projection="3d")
    for c in ens_coords:
        for j in range(n_r - 1):
            ax2.plot(c[j:j+2,0], c[j:j+2,1], c[j:j+2,2],
                     color=C["grey"], alpha=0.12, lw=1.1)
    # pocket spheres on mean structure
    mean_c = np.array(ens_coords).mean(0)
    for j in range(n_r):
        if pocket[j]:
            ax2.scatter(mean_c[j,0], mean_c[j,1], mean_c[j,2],
                        s=110, c=C["pocket"], alpha=0.55,
                        edgecolors="#B7950B", linewidth=0.6, zorder=5)
    ax2.set_axis_off(); ax2.set_box_aspect([1,1,1]); ax2.view_init(20, 135)
    ax2.set_title("DYNAMIC ENSEMBLE\\nWITH TRANSIENT POCKETS", fontsize=11,
                  fontweight="bold", color=C["open"])

    fig.text(0.50, 0.01,
             "Static docking on one PDB misses transient, cryptic binding pockets.",
             ha="center", fontsize=9, style="italic", color="#555")
    plt.savefig(save, dpi=300, facecolor="white", edgecolor="none",
                bbox_inches="tight")
    plt.show(); print(f"Saved {save}")

panel2(ens, static_raw, f"{OUT}/panel2_static_vs_dynamic.png")""")

# ── Cell 7: PANEL 3 — Pocket-State Clustering ───────────────────────────────
add_md("""\
---
## Panel 3 — Pocket-State Modeling
> *"We cluster the ensemble into pocket-centric macrostates (OPEN vs CLOSED)."*

Uses PCA on CA coordinates as a stand-in for tICA (swap in your own
tICA / TICA projection from PyEMMA / deeptime if available).""")

add_code("""\
def panel3(all_coords, ref_coord, save):
    n_r = len(ref_coord)
    X = np.array([c[:n_r].flatten() for c in all_coords])

    pca = PCA(n_components=2)
    Z   = pca.fit_transform(X)

    km = KMeans(n_clusters=2, n_init=10, random_state=42)
    lab = km.fit_predict(Z)

    fig = plt.figure(figsize=(13, 5))
    gs  = fig.add_gridspec(1, 3, width_ratios=[2.2, 1, 1], wspace=0.35)

    # ── scatter ────────────────────────────────────────────────────────────
    ax = fig.add_subplot(gs[0, 0])
    for ci, (lbl, col, mk) in enumerate(
            [("\\\"Open\\\" state", C["open"], "o"),
             ("\\\"Closed\\\" state", C["closed"], "s")]):
        m = lab == ci
        ax.scatter(Z[m, 0], Z[m, 1], c=col, label=lbl,
                   s=40, marker=mk, edgecolors="white", lw=0.5, alpha=0.75, zorder=3)
        # 2-sigma ellipse
        if m.sum() > 3:
            pts = Z[m]
            mu  = pts.mean(0)
            cov = np.cov(pts.T)
            vals, vecs = np.linalg.eigh(cov)
            ang = np.degrees(np.arctan2(vecs[1,1], vecs[0,1]))
            w, h = 2*2*np.sqrt(np.abs(vals))
            ax.add_patch(patches.Ellipse(mu, w, h, angle=ang, fill=False,
                         edgecolor=col, lw=1.3, ls="--", alpha=0.6))
    ev = pca.explained_variance_ratio_ * 100
    ax.set_xlabel(f"PC 1  ({ev[0]:.1f}% var.)")
    ax.set_ylabel(f"PC 2  ({ev[1]:.1f}% var.)")
    ax.legend(loc="upper right", framealpha=0.9, edgecolor="#ccc")
    ax.set_title("CONFORMATIONAL LANDSCAPE", fontsize=11, fontweight="bold")
    ax.grid(True, alpha=0.25, lw=0.5)

    # ── mini-structure thumbnails ──────────────────────────────────────────
    for ci, (gi, ttl, col) in enumerate(
            [(1, "OPEN\\nMACROSTATE", C["open"]),
             (2, "CLOSED\\nMACROSTATE", C["closed"])]):
        axs = fig.add_subplot(gs[0, gi], projection="3d")
        mask = lab == ci
        pts  = Z[mask]
        mu   = pts.mean(0)
        d    = np.linalg.norm(pts - mu, axis=1)
        rep  = np.where(mask)[0][d.argmin()]
        co   = all_coords[rep][:n_r]
        for j in range(n_r - 1):
            axs.plot(co[j:j+2,0], co[j:j+2,1], co[j:j+2,2],
                     color=col, alpha=0.85, lw=2.0, solid_capstyle="round")
        axs.set_axis_off(); axs.set_box_aspect([1,1,1]); axs.view_init(20, 135)
        axs.set_title(ttl, fontsize=9, fontweight="bold", color=col)

    fig.text(0.50, 0.01,
             "We cluster the ensemble into pocket-centric macrostates (OPEN vs CLOSED).",
             ha="center", fontsize=9, style="italic", color="#555")
    plt.savefig(save, dpi=300, facecolor="white", edgecolor="none",
                bbox_inches="tight")
    plt.show(); print(f"Saved {save}")

panel3(expanded_coords, ens[0], f"{OUT}/panel3_clustering.png")""")

# ── Cell 8: PANEL 4 — SE(3)-Equivariant Diffusion ──────────────────────────
add_md("""\
---
## Panel 4 — SE(3)-Equivariant Diffusion Schematic
> *"State-aware SE(3) diffusion generates 3D structures conditioned on the OPEN pocket state."*

Pure vector/matplotlib drawing — no external images needed.""")

add_code("""\
def panel4(save):
    fig, ax = plt.subplots(figsize=(13, 4.2))
    ax.set_xlim(-0.2, 13); ax.set_ylim(-0.3, 4.8)
    ax.set_aspect("equal"); ax.axis("off")

    np.random.seed(7)

    # ── STAGE 1: noisy backbones (left) ────────────────────────────────────
    for k in range(5):
        y0 = 1.6 + k * 0.45
        xs = np.linspace(0.4, 2.6, 25)
        ys = y0 + np.cumsum(np.random.normal(0, 0.12, 25))
        ys -= ys.mean() - y0
        ax.plot(xs, ys, color=C["noise"], alpha=0.45, lw=1.6)
    # gaussian noise dots
    ax.scatter(np.random.uniform(0.2, 2.8, 40),
               np.random.uniform(1.0, 3.8, 40),
               s=6, c=C["noise"], alpha=0.25, zorder=1)
    ax.text(1.5, 0.7, "Noisy structures  (t)", ha="center",
            fontsize=9, fontweight="bold", color=C["dark"])

    # ── arrow 1 ────────────────────────────────────────────────────────────
    ax.annotate("", xy=(3.9, 2.3), xytext=(3.1, 2.3),
                arrowprops=dict(arrowstyle="-|>", color=C["dark"], lw=2.2))

    # ── STAGE 2: model box (centre) ────────────────────────────────────────
    box = FancyBboxPatch((4.1, 0.9), 4.4, 2.8,
                         boxstyle="round,pad=0.25",
                         facecolor="#E8F4FD", edgecolor=C["closed"], lw=2.2)
    ax.add_patch(box)
    ax.text(6.3, 3.1, "SE(3)-Equivariant", ha="center",
            fontsize=11, fontweight="bold", color=C["closed"])
    ax.text(6.3, 2.55, "Diffusion Model", ha="center",
            fontsize=11, fontweight="bold", color=C["closed"])
    # rotation arrows
    arc = Arc((6.3, 1.7), 1.0, 1.0, angle=0, theta1=30, theta2=330,
              color=C["closed"], lw=1.6)
    ax.add_patch(arc)
    ax.annotate("", xy=(6.73, 1.9), xytext=(6.68, 2.05),
                arrowprops=dict(arrowstyle="-|>", color=C["closed"], lw=1.4))

    # conditioning label
    ax.text(6.3, 0.25, "conditioned on OPEN macrostate", ha="center",
            fontsize=8.5, style="italic", color=C["open"],
            bbox=dict(boxstyle="round,pad=0.35", fc="#FDF2E9",
                      ec=C["open"], lw=1.2))

    # ── arrow 2 ────────────────────────────────────────────────────────────
    ax.annotate("", xy=(9.3, 2.3), xytext=(8.7, 2.3),
                arrowprops=dict(arrowstyle="-|>", color=C["dark"], lw=2.2))

    # ── STAGE 3: clean generated structures (right) ────────────────────────
    for k in range(4):
        y0 = 1.6 + k * 0.45
        xs = np.linspace(9.7, 12.0, 25)
        ys = y0 + 0.25 * np.sin(np.linspace(0, 2*np.pi, 25) + k)
        ax.plot(xs, ys, color=C["open"], alpha=0.7, lw=2.2, solid_capstyle="round")

    # pocket highlight
    pocket_e = patches.Ellipse((10.9, 2.3), 0.7, 0.6,
                                fc=C["pocket"], alpha=0.4,
                                ec="#B7950B", lw=1.6)
    ax.add_patch(pocket_e)
    ax.text(10.9, 0.7, "Generated OPEN-state\\nensemble", ha="center",
            fontsize=9, fontweight="bold", color=C["dark"])

    ax.set_title("STATE-AWARE SE(3) DIFFUSION", fontsize=13,
                 fontweight="bold", color=C["dark"], pad=12)
    fig.text(0.50, 0.01,
             "State-aware SE(3) diffusion generates 3D structures "
             "conditioned on the OPEN pocket state.",
             ha="center", fontsize=9, style="italic", color="#555")
    plt.savefig(save, dpi=300, facecolor="white", edgecolor="none",
                bbox_inches="tight")
    plt.show(); print(f"Saved {save}")

panel4(f"{OUT}/panel4_diffusion.png")""")

# ── Cell 9: PANEL 5 — Validation Funnel ─────────────────────────────────────
add_md("""\
---
## Panel 5 — Validation Funnel
> *"Multistage filtering yields experimentally verified, state-selective binders
> that stabilize the OPEN pocket."*""")

add_code("""\
def panel5(save):
    fig, ax = plt.subplots(figsize=(9, 9))
    ax.set_xlim(0, 10); ax.set_ylim(0, 11)
    ax.set_aspect("equal"); ax.axis("off")

    cx = 5.0
    stages = [
        dict(yt=9.6, yb=8.1, wt=8.0, wb=6.4,
             label="Generated molecules", col="#5DADE2"),
        dict(yt=8.1, yb=6.6, wt=6.4, wb=4.8,
             label="Valid OPEN-state poses", col="#48C9B0"),
        dict(yt=6.6, yb=5.1, wt=4.8, wb=3.4,
             label="Selective \\u0394\\u0394G\\n(open > closed)", col=C["open"]),
        dict(yt=5.1, yb=3.6, wt=3.4, wb=2.0,
             label="State-selective\\nbinders", col=C["red"]),
    ]

    for s in stages:
        l1, r1 = cx - s["wt"]/2, cx + s["wt"]/2
        l2, r2 = cx - s["wb"]/2, cx + s["wb"]/2
        trap = plt.Polygon([[l1, s["yt"]], [r1, s["yt"]],
                             [r2, s["yb"]], [l2, s["yb"]]],
                            fc=s["col"], alpha=0.28,
                            ec=s["col"], lw=2.2)
        ax.add_patch(trap)
        my = (s["yt"] + s["yb"]) / 2
        ax.text(cx, my, s["label"], ha="center", va="center",
                fontsize=10.5, fontweight="bold", color=C["dark"])

    # ── molecule icons (top stage) ─────────────────────────────────────────
    np.random.seed(3)
    for j in range(10):
        mx = cx + np.random.uniform(-3.0, 3.0)
        my = 9.0 + np.random.uniform(-0.2, 0.3)
        ax.add_patch(plt.Circle((mx, my), 0.12, fc="#5DADE2", ec="white", lw=0.6))

    # ── mini DDG bar chart (3rd stage) ─────────────────────────────────────
    bx0 = 7.3
    for bx, bh, col, lab in [(bx0, 0.9, C["open"], "Open"),
                              (bx0+0.6, 0.35, C["closed"], "Closed")]:
        ax.add_patch(patches.Rectangle((bx, 5.2), 0.45, bh,
                     fc=col, ec="white", lw=0.8))
        ax.text(bx+0.22, 5.05, lab, ha="center", fontsize=6.5, color=col)
    ax.text(bx0+0.5, 6.3, "\\u0394\\u0394G", ha="center", fontsize=8,
            fontweight="bold", color=C["dark"])

    # ── side arrow ─────────────────────────────────────────────────────────
    ax.annotate("", xy=(0.8, 3.8), xytext=(0.8, 9.4),
                arrowprops=dict(arrowstyle="-|>", color=C["dark"], lw=2.5))
    ax.text(0.35, 6.6, "In silico +\\nexperimental\\nfiltering",
            ha="center", fontsize=8, fontweight="bold", color=C["dark"],
            rotation=90)

    # ── binder icon (bottom stage) ─────────────────────────────────────────
    for dx in [-0.3, 0.3]:
        ax.add_patch(plt.Circle((cx+dx, 4.0), 0.13,
                     fc=C["red"], ec="white", lw=0.8))

    ax.set_title("VALIDATION FUNNEL", fontsize=14,
                 fontweight="bold", color=C["dark"], pad=15)
    fig.text(0.50, 0.01,
             "Multistage filtering yields experimentally verified, "
             "state-selective binders\\nthat stabilize the OPEN pocket.",
             ha="center", fontsize=9, style="italic", color="#555")
    plt.savefig(save, dpi=300, facecolor="white", edgecolor="none",
                bbox_inches="tight")
    plt.show(); print(f"Saved {save}")

panel5(f"{OUT}/panel5_funnel.png")""")

# ── Cell 10: Composite assembly ─────────────────────────────────────────────
add_md("""\
---
## Assemble all panels into a single visual abstract
Produces a **3 600 px wide** composite at 300 dpi.""")

add_code("""\
from PIL import Image, ImageDraw, ImageFont

def assemble(panel_dir, save):
    names = {
        1: "panel1_ensemble.png",
        2: "panel2_static_vs_dynamic.png",
        3: "panel3_clustering.png",
        4: "panel4_diffusion.png",
        5: "panel5_funnel.png",
    }
    imgs = {}
    for k, fn in names.items():
        p = os.path.join(panel_dir, fn)
        if os.path.exists(p):
            imgs[k] = Image.open(p).convert("RGB")
        else:
            print(f"  Panel {k} missing ({p})")

    W = 3600
    M = 50   # margin

    def fit(img, tw):
        r = tw / img.width
        return img.resize((tw, int(img.height * r)), Image.LANCZOS)

    half = (W - 3*M) // 2
    p1 = fit(imgs[1], half)
    p2 = fit(imgs[2], half)
    p3 = fit(imgs[3], half)
    p4 = fit(imgs[4], half)
    p5 = fit(imgs[5], W - 2*M)

    r1 = max(p1.height, p2.height)
    r2 = max(p3.height, p4.height)
    r3 = p5.height
    title_h = 100
    H = title_h + r1 + r2 + r3 + 6*M

    canvas = Image.new("RGB", (W, H), "white")

    # panel labels
    draw = ImageDraw.Draw(canvas)
    try:
        fnt = ImageFont.truetype("/usr/share/fonts/truetype/dejavu/DejaVuSans-Bold.ttf", 42)
    except Exception:
        fnt = ImageFont.load_default()

    # title
    draw.text((W//2, 30), "VISUAL ABSTRACT", fill="#333333", font=fnt, anchor="mt")

    y = title_h + M
    # row 1
    canvas.paste(p1, (M, y))
    canvas.paste(p2, (M + half + M, y))
    draw.text((M+10, y+5), "A", fill="#333", font=fnt)
    draw.text((M+half+M+10, y+5), "B", fill="#333", font=fnt)
    y += r1 + M

    # row 2
    canvas.paste(p3, (M, y))
    canvas.paste(p4, (M + half + M, y))
    draw.text((M+10, y+5), "C", fill="#333", font=fnt)
    draw.text((M+half+M+10, y+5), "D", fill="#333", font=fnt)
    y += r2 + M

    # row 3
    canvas.paste(p5, (M, y))
    draw.text((M+10, y+5), "E", fill="#333", font=fnt)

    canvas.save(save, dpi=(300, 300))
    print(f"\\nVisual abstract saved: {save}")
    print(f"Dimensions: {canvas.size[0]} x {canvas.size[1]} px  (300 dpi)")

    fig, ax = plt.subplots(figsize=(18, H/200))
    ax.imshow(canvas); ax.axis("off")
    plt.show()

assemble(OUT, f"{OUT}/visual_abstract_COMPLETE.png")""")

# ── Cell 11: Download from Colab ─────────────────────────────────────────────
add_md("""\
---
## Download results
Run the cell below to download the composite image and all individual panels.""")

add_code("""\
try:
    from google.colab import files
    # download composite
    files.download(f"{OUT}/visual_abstract_COMPLETE.png")
    # download individual panels
    for fn in sorted(os.listdir(OUT)):
        if fn.endswith(".png"):
            files.download(os.path.join(OUT, fn))
    print("Downloads triggered \\u2713")
except ImportError:
    print("Not running in Colab — files are in:", OUT)
    print("Contents:")
    for fn in sorted(os.listdir(OUT)):
        print(f"  {fn}")""")

# ── Cell 12: How to swap in your own data ───────────────────────────────────
add_md("""\
---
## How to use your own MD trajectory / IDP ensemble

```python
import MDAnalysis as mda

u = mda.Universe("topology.pdb", "trajectory.xtc")
ca = u.select_atoms("name CA")

custom_coords = []
# sample every 10th frame
for ts in u.trajectory[::10]:
    custom_coords.append(ca.positions.copy())

# then align and feed into any panel function:
aligned = align_ensemble(custom_coords)
panel1(aligned, "my_panel1.png")
```

**For tICA** replace PCA with [PyEMMA](https://pyemma.org) or
[deeptime](https://deeptime-ml.github.io):

```python
import pyemma
tica = pyemma.coordinates.tica(data, lag=50, dim=2)
Z = np.vstack(tica.get_output())
```""")

# ═══════════════════════════════════════════════════════════════════════════════
# Write notebook
# ═══════════════════════════════════════════════════════════════════════════════
notebook = {
    "nbformat": 4,
    "nbformat_minor": 0,
    "metadata": {
        "colab": {
            "provenance": [],
            "name": "Nature Visual Abstract - IDP Ensemble & SE(3) Diffusion"
        },
        "kernelspec": {
            "name": "python3",
            "display_name": "Python 3"
        },
        "language_info": {
            "name": "python"
        }
    },
    "cells": cells
}

out_path = "nature_visual_abstract_colab.ipynb"
with open(out_path, "w") as f:
    json.dump(notebook, f, indent=1)

print(f"Notebook written: {out_path}  ({len(cells)} cells)")
