# CDBGC — Cell Density-Based Grid Clustering

A grid-based clustering algorithm with **automatic outlier detection**, designed for datasets where
clusters are non-spherical and where outliers would otherwise be absorbed into clusters by methods
such as K-means.

> Course project: **Individual Research Project (MATH40008)**, Department of Mathematics,
> Imperial College London — Year 1, Summer Term 2026. Sole author: Chen Yuze (Jack).
> Topic chosen from a list of five; the project was self-directed (no prescribed method).
> Status: coursework, not peer-reviewed. All reported metrics are on **synthetic data**.

---

## 1. Method

### 1.1 Grid construction
Take the min/max of each axis to form the bounding box of the dataset, then divide both axes into
`grid_n` intervals, giving a `grid_n × grid_n` mesh.

### 1.2 Density estimation (bilinear splatting)
For a point `(x, y)` inside its primary cell `[Lx, Rx] × [Ly, Ry]`, compute the relative position

```
tx = (x − Lx) / (Rx − Lx)        ty = (y − Ly) / (Ry − Ly)
```

Both are in `[0, 1]`. Define the weight vectors over
`[left/lower neighbour, primary cell, right/upper neighbour]`:

```
Wx = [ (1 − tx)/2 ,  0.5 ,  tx/2 ]
Wy = [ (1 − ty)/2 ,  0.5 ,  ty/2 ]
```

Each point distributes weight over its `3 × 3` neighbourhood via the outer product

```
Wab = Wx[a] × Wy[b],     a, b ∈ {1, 2, 3}
```

The nine weights sum to 1, so the total mass of the dataset is conserved. Spreading the mass instead
of counting points per cell prevents points lying near a cell boundary from flipping between cells
under small perturbations — the density field is smoothed.

### 1.3 Cluster formation and outlier detection
- Cells whose accumulated density is `≥ thr` are kept as **surviving cells**.
- Two surviving cells are **connected** if their grid distance is `≤ 1` in both directions
  (8-neighbourhood). Connected components of surviving cells = candidate clusters.
- A data point is assigned to a component if its **accumulated contribution to that component exceeds
  `1/4`**; otherwise it is labelled an outlier (`−1`).

The `1/4` rule follows from the weight structure: the largest weight a single point can place in any
one cell is `0.5 × 0.5 = 0.25` (its primary cell). Requiring more than `1/4` therefore means the
assignment cannot rest on a single cell — it must be supported by several cells of the same component.

### 1.4 Choosing `grid_n` and `thr`
`grid_n` and `thr` are selected by grid search. The objective is an **adjusted DBCV**:

```
Adjusted DBCV = Raw DBCV − β × r_outlier
r_outlier     = N_outlier / N_points,        β = 2.5
```

`DBCV` (Density-Based Clustering Validation) handles arbitrary cluster shapes but does **not** penalise
excessive outlier labelling — a configuration can inflate its raw score simply by declaring hard points
to be outliers. Subtracting a term proportional to the outlier ratio removes that incentive.

Configurations with fewer than 2 or more than 5 clusters, or containing any cluster of fewer than 5
points, are scored `−2.0` without computing DBCV.

### 1.5 Evaluation
External validation with the **Adjusted Rand Index (ARI)** against ground-truth labels, compared with
two baselines forced to `K = 2`: K-means (`nstart = 25`) and Ward.D2 hierarchical clustering.

---

## 2. Results

Synthetic benchmark: two interleaving half-circles ("two moons") of 200 points each, plus 15 uniformly
distributed outliers — 415 points total, `set.seed(42)`.

| Method | ARI |
|---|---|
| **Grid Density (CDBGC)** | **0.9709** |
| K-means (K = 2) | 0.2542 |
| Hierarchical, Ward.D2 (K = 2) | 0.238 |

Both baselines are forced to `K = 2` and have no outlier mechanism, so the outliers are absorbed into
the two clusters; CDBGC detects the number of clusters and the outliers itself.

---

## 3. Known limitations

- **Dimensionality.** The grid is partitioned strictly along the coordinate axes: a `10 × 10` split
  gives 100 cells in 2D, but `10³ = 1000` in 3D and `10^D` in `D` dimensions. The method is a 2D (at
  best low-dimensional) method and does not scale to high-dimensional data.
- The `1/4` assignment rule and `β = 2.5` were tuned on the benchmark above; both are dataset-dependent
  and would need re-validation on a new domain.
- `grid_n` and `thr` are found by exhaustive grid search, which becomes expensive for fine resolutions.
- Connectivity is defined on a fixed 8-neighbourhood, so two dense regions separated by a single
  sub-threshold cell are treated as different clusters.

---

## 4. Repository layout

```
R/01_grid_density_demo.R     Gaussian blobs + an outlier at the origin; walks through grid construction,
                             density accumulation, connected components and soft point assignment.
R/02_two_moons_tuning.R      Two-moons benchmark; hyper-parameter sweep over grid_n × thr with the
                             outlier-penalised DBCV objective, plus the tuning-surface and density plots.
R/03_ARI_vs_baselines.R      Two-moons benchmark; ARI comparison of CDBGC against K-means and
                             Ward.D2 hierarchical clustering.
poster/CDBGC_poster.pdf      The compiled A0 research poster.
latex/poster.tex             LaTeX source of the poster, plus its figures and bibliography.
                             Built with the Imperial College London poster template (XeLaTeX);
                             the class file is not redistributed here — see latex/README.md.
figures/                     Figures used in the poster: grid construction, density heat map, connected
                             components, adjusted-DBCV tuning surface, ground truth vs. output, K-means
                             comparison, and the end-to-end workflow diagram.
```

### Running the code

Requires R (4.x) and the packages `igraph`, `dbscan`, `ggplot2`, `patchwork`:

```r
install.packages(c("igraph", "dbscan", "ggplot2", "patchwork"))
```

Each script is self-contained (it generates its own data with a fixed seed) and can be run with
`Rscript R/03_ARI_vs_baselines.R`. `01_grid_density_demo.R` assumes the packages are already installed.
`03_ARI_vs_baselines.R` prints the ARI table quoted in section 2; `dbcv()` must be available from the
installed version of `dbscan`.

---

## 5. AI assistance (full disclosure)

- The algorithm idea, the density-weight formulation (`Wx`, `Wy`, `Wab`), the `1/4` assignment rule and
  the adjusted-DBCV objective are the author's own.
- ChatGPT and DeepSeek were used to **look up suitable validation indices** (which led to DBCV), to
  check whether an equivalent method already existed, to **scaffold parts of the R code** (e.g. which
  functions to use for binning and for connected components), to debug, and to draft the poster layout.
- This README was compiled with AI assistance from the code and poster.

## 6. Reference

Davoud Moulavi, Pablo A. Jaskowiak, Ricardo J. G. B. Campello, Arthur Zimek, Jörg Sander.
*Density-Based Clustering Validation*, pp. 839–847 (SDM 2014).
