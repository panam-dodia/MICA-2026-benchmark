# Benchmark protocol

## 1. Dataset

### Selection (scripts 01, 03, 06)
| Step | Criterion | Entries |
|---|---|---|
| 0 | RCSB search: electron microscopy, released on or after 2026-01-01, resolution ≤ 4.0 Å, no nucleic acid entity | 4,381 |
| 1 | Single particle reconstruction | 4,050 |
| 2 | Exactly one EMDB map | 4,046 |
| 3 | 1.5 ≤ resolution ≤ 4.0 Å | 4,042 |
| 4 | All polymer entities are protein | 4,041 |
| 5 | No non-standard residues | 3,697 |
| 6 | 300 ≤ deposited polymer residues ≤ 2,500 | 2,232 |
| 7 | ≤ 4 unique protein entities | 1,916 |
| 8 | Novelty: no protein entity (≥ 30 aa) with > 25% identity to any PDB entity released before 2026-01-01 | 109 |
| 8b | Exclude amyloid/fibril entries (whole word in title or PDB keywords) and entries whose longest chain has < 50 residues | 106 |
| 9 | One entry per 25%-identity cluster (union of 2026 hits; random representative, seed 2026) | 51 |
| 10 | Stratified random sample, 10 per size bin (300–700, 700–1,400, 1,400–2,500), seed 2026 | 30 |

* Sequence searches use the RCSB sequence service (MMseqs2), identity cutoff 0.25, E-value cutoff 0.1. Chains shorter
  than 30 residues are not searched (one 10-residue peptide in 9PZO).
* 1,807 of the 1,916 entries after step 7 (94.3%) have a > 25% identity homolog released before 2026.
* Step 8b removed three polymorphs of a designed Orb2 amyloid (9PX8, 9PZD, 9Q37).
* Redundancy search validation: every query must return its own entity (always a 2026 release); otherwise the search
  is retried. Validated hits are stored in `dataset/dedup_hits.json`, so reruns give identical clusters.
* The 300–700 bin had exactly 10 eligible entries, so all were taken.

### Size definition
The size filter and bins use `deposited_polymer_monomer_count`, i.e. **deposited polymer residues** (full sequence of
all chains, including residues not resolved in the map). All 30 entries also satisfy 300 ≤ modeled residues ≤ 2,500
(modeled 338 to 2,329). By modeled residues five entries change bin (21DI → 300–700; 10AW, 9T9R, 9R9A, 9OME →
700–1,400); results are reported for both definitions (`scripts/10_residue_counts.py`, `dataset/benchmark_30_sizes.csv`).

### Final set
30 entries, resolution 1.95–3.92 Å, 39 unique protein entities, 1–18 chains. For 9I8C (1 → 12 chains, C12),
9SZB (2 → 24, C12) and 12ZE (2 → 120, icosahedral) the deposited model is one asymmetric unit of a map showing the
full assembly; their biological assembly 1 is used for in-place metrics (Section 4).

### Files per entry (Zenodo, `dataset/entries/<EMDB>/`)
EMDB primary map (from the wwPDB EMDB mirror), RCSB FASTA, deposited mmCIF, and for the symmetric entries
`<pdb>_assembly1.cif`. Contour level and resolution come from the EMDB API (via MICA's EMDB helper); all EMDB
resolutions equal the PDB-reported resolution.

## 2. AlphaFold3 inputs
* 39 single-chain jobs (one per unique protein entity) on AlphaFold Server, as in the MICA paper.
* Fixed seed `modelSeeds: ["2026"]` (MICA's own job-file script leaves the seed empty, i.e. random).
* `useStructureTemplate: true`, `maxTemplateDate: 2025-02-03`, the latest date the server allows
  (AlphaFold Server JSON format documentation, google-deepmind/alphafold `server/README.md`). No 2026 structure can be
  a template.
* `model_0` of each job is used. The same files feed MICA, EModelX(+AF) and EMProt (each in its expected layout).
* Verified for all 39 jobs: model_0 present, seed and template date recorded by the server, model length equal to the
  submitted sequence.

## 3. Methods

All methods use default parameters. One method runs at a time on the GPU. A 12 h wall-clock limit per map applies to
every method (`scripts/watchdog.sh`); failures count as failures; maps are never cropped.

| Method | Version | Inputs | Final model scored |
|---|---|---|---|
| MICA | commit 726d1bb (2026-09-28) | map, FASTA, AF3 | Phenix-refined all-atom model (`MICA_pipeline.sh` steps) |
| EModelX(+AF) | commit 53713cf (2025-02-18) | map, FASTA, AF3 templates (protocol `temp_flex`) | Phenix-refined model (`--run_phenix`) |
| ModelAngelo | 1.0.18, commit 994945b | map, FASTA | `output/output.cif` |
| CryoAtom | CryoAtom2 2.1.1, commit 856e250 | map, FASTA | `output/output.cif` |
| EMProt | 1.2, commit 9489d52 | map, one FASTA record and one AF3 template per chain copy | `output/output.cif` |

Phenix refinement is applied only where it is part of the method's own pipeline (MICA, EModelX(+AF)).

### Installation notes
* **MICA**: conda part of `environment.yml` plus `cudatoolkit=11.1` (PyTorch 1.8.1 cu111 for an Ampere GPU). The pip
  part of `environment.yml` fails on Python 3.8 (unpinned `lmdb`) and would upgrade torch/numpy; instead pinned
  packages were installed with `--no-deps` (following MICA's own Dockerfile): mrcfile 1.3.0, superpose3d 1.1.1,
  tqdm 4.61.2, networkx 2.6.3, natsort 8.4.0, einops 0.8.1, open3d 0.13.0, rotary-embedding-torch 0.2.1, then
  matplotlib 3.3.4, requests, scikit-learn 1.0.2, addict 2.4.0, pyyaml 6.0.3 with numpy 1.19.1 and scipy 1.5.2 pinned.
  atom3d, wandb and scikit-image are not imported by MICA's inference code. Phenix 1.20.1 (the version in MICA's
  README; MICA's refinement parameter file is a 1.x parameter set).
* **EModelX(+AF)**: clone of the MICA environment (identical Python 3.8 / PyTorch 1.8.1 stack) plus scikit-image
  0.18.2, imageio 2.9.0, tifffile 2021.4.8, PyWavelets 1.1.1. Weights `models.zip` from the authors' Google Drive.
  `--download_afdb` is not used (it would replace the AF3 templates).
* **ModelAngelo**: repository `install_script.sh --download-weights` (Python 3.11, torch 2.9.1+cu128).
* **CryoAtom**: repository `install.sh` (run with `bash`); ESM2-650M and RNA-FM weights pre-cached so no download is
  needed at run time.
* **EMProt**: `conda env create -f environment.yml`, `pip install -e .`, weights from the authors' server.

### System-level adjustments (no change to algorithms or parameters)
1. **MICA and EModelX(+AF)**: `OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1`. Without it, the
   `multiprocessing.Pool` started after multithreaded numpy/torch work deadlocks (workers sleeping in `futex_wait`).
   The methods' own worker pools still provide the parallelism.
2. **MICA and EModelX(+AF)**: `MKL_THREADING_LAYER=GNU`. Merizo (MICA's domain splitting) imports torch before numpy and
   otherwise fails with "MKL_THREADING_LAYER=INTEL is incompatible with libgomp". `process_AF3_results.py` exits 0
   even when Merizo fails, so the run script checks that domains were produced.
3. **EModelX(+AF)**: `--mul_proc_num 14` instead of the default 30, because 30 forked workers exceed 50 GB RAM on large
   maps; 14 is the value MICA chooses automatically on the same machine (0.6 × 24 threads).
4. **CryoAtom and EMProt**: the WSL2 kernel has no vsyscall page. Precompiled helpers built on CentOS crash:
   EMProt's `dockx` and `domasgx` segfault at startup, and `getp` (both methods) crashed on real data. Following the
   authors' instructions, `getp` was rebuilt from the bundled source (CryoAtom `compile.sh`; EMProt `bin/srcs/getp`,
   `make`), and EMProt's `dockx`/`domasgx` were replaced with the recompiled versions shipped in `bin/recompiled/`
   (`dockx` needs `libfftw3.so.3`: `apt install libfftw3-double3`).
5. **Memory**: WSL2 memory raised from the default 31 GiB to 52 GB (+16 GB swap). The run scripts read the kernel
   `oom_kill` counter before and after each entry and mark the entry FAILED if it increased, because Phenix docking
   can lose worker processes to the OOM killer and still exit 0.
6. MICA and EModelX(+AF) start `phenix.real_space_refine` in the background; the run scripts wait for it before
   starting the next entry, so refinements never overlap.

### Failures (final)
| Method | Entry | Reason |
|---|---|---|
| MICA | 26FB | run.py OOM-killed after a clean 9.1 h docking (needs > ~66 GB RAM+swap) |
| MICA | 12ZE | Phenix docking of domain 1 of 4 unfinished after 11.7 h (12 h limit). MICA's `dock_in_map.py` has a timeout handler but `subprocess.run` is called without `timeout=`, so docking can run indefinitely. |
| EModelX(+AF) | 26FB | OOM-killed (~52 GB RSS) |
| EModelX(+AF) | 12ZE | `MemoryError`: 256 GiB for a 107,065 × 107,065 Cα distance matrix |
| EMProt | 9YFW, 12ZE | OOM-killed (~52 GB RSS) |
| ModelAngelo, CryoAtom | none | 12ZE took 4.5 h and 9.2 h respectively |

## 4. Evaluation (`scripts/score_raw.sh`, `scripts/parse_scores.py`)
| Metric | Computation | Reference |
|---|---|---|
| TM-score | `USalign model deposited -mm 1 -ter 0`, normalized by the deposited structure | deposited model |
| Aligned Cα length, sequence identity | same US-align run | deposited model |
| Cα match | `phenix.chain_comparison model reference` (model as target): % of model residues with a reference Cα within 3 Å | in-place reference |
| Sequence match | same run, SEQ MATCH (%) | in-place reference |
| Cα quality | Cα match × N_model / N_reference (Cryo2Struct definition) | in-place reference |
| Recall (extra) | `phenix.chain_comparison reference model` (reference as target): FOUND % | in-place reference |

* Cα match follows the MICA/Cryo2Struct definition ("percentage of Cα atoms in the model that have corresponding
  residues within 3 Å in the true structure"), i.e. precision; recall is reported separately.
* In-place reference = biological assembly 1 for 9I8C, 9SZB, 12ZE; otherwise the deposited model (assembly 1 equals
  the deposited model for the other 27 entries). Cα quality is not reported for the three symmetric entries (it
  mainly reflects how many copies a method builds).
* chain_comparison reports the better of `unique_part_of_target_only` True/False (rows `Unique_target` or
  `Entire_target`); both are parsed. Residue counts come from "Original residues in target / In query".
* 9V3Z's deposited model contains a ligand with a 5-character component ID (`A1BDB`) that Phenix 1.20.1 cannot read;
  `scripts/12_phenix_compatible_refs.py` writes a copy without those 19 ligand atoms for chain_comparison only.
* US-align cannot parse the mmCIF written by `phenix.real_space_refine` for EModelX(+AF); for US-align only, those
  models are rewritten as PDB with Phenix's own reader/writer (coordinates unchanged).
* Statistics (`scripts/analyze.py`): paired two-sided t test (as in the MICA paper) and Wilcoxon signed-rank test
  (normal approximation with tie correction) for MICA versus each method on the entries both built.

## 5. Runtime
Wall-clock time per successful map, final attempt only (`runs/timings/`). MICA's first 22 entries ran under the 31 GiB
WSL2 memory default, the rest under 52 GB. CPU scoring jobs occasionally overlapped with GPU runs; times are indicative.
