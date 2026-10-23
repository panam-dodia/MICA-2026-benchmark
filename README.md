# MICA 2026 benchmark

Benchmark of **MICA** (cryo-EM protein model building with AlphaFold3; Gyawali, Dhakal & Cheng, *Communications
Chemistry* 2025) against **EModelX(+AF)**, **ModelAngelo**, **CryoAtom** and **EMProt** on 30 single particle cryo-EM
maps released in EMDB/PDB in 2026, none of which has a protein chain with more than 25% sequence identity to any PDB
entry released before 2026.

* Code: this repository (MIT License)
* Data (maps, sequences, reference structures, AlphaFold3 predictions, raw scoring output): Zenodo, [doi.org/10.5281/zenodo.23584417](https://doi.org/10.5281/zenodo.23584417)
* Author: Panam Dodia ([github.com/panam-dodia](https://github.com/panam-dodia))

## Results (27 maps completed by all five methods)

| Method | Models built | TM-score mean (median) | Cα match % | Seq. match % | Cα quality | Recall % | Median minutes / map |
|---|---|---|---|---|---|---|---|
| MICA | 28 / 30 | 0.877 (0.968) | 92.4 | 94.5 | 91.1 | 89.4 | 42.2 |
| EModelX(+AF) | 28 / 30 | 0.890 (0.938) | 92.3 | 93.0 | 91.8 | 92.0 | 7.2 |
| ModelAngelo | 30 / 30 | 0.747 (0.879) | 94.3 | 98.2 | 81.7 | 83.7 | 5.8 |
| CryoAtom | 30 / 30 | 0.884 (0.984) | 88.6 | 97.2 | 95.2 | 96.2 | 6.8 |
| EMProt | 28 / 30 | 0.951 (0.983) | 82.1 | 97.0 | 95.6 | 95.9 | 6.4 |

Paired comparisons of MICA (TM-score, entries both methods built; paired t test / Wilcoxon):
better than ModelAngelo (0.880 vs 0.729, p = 0.014 / 0.016); not different from EModelX(+AF) (0.880 vs 0.879,
p = 0.97 / 0.91) or CryoAtom (0.880 vs 0.871, p = 0.85 / 0.65); lower than EMProt (0.877 vs 0.951, p = 0.019 / 0.017).
Full tables: [`eval/analysis.md`](eval/analysis.md); per entry scores: [`eval/scores.csv`](eval/scores.csv).
Protocol, decisions and every system level adjustment: [`docs/PROTOCOL.md`](docs/PROTOCOL.md).

## Repository layout

```
scripts/      pipeline (numbered dataset steps, one run script per method, scoring, analysis)
dataset/      candidate list, selection tables, novelty and redundancy search results, final 30 entries
af3/          AlphaFold Server job files (seed 2026, templates up to 2025-02-03) and job manifest
runs/         entry list used by the run scripts, EMProt per chain inputs, per method timing logs
eval/         parsed scores, summaries and the analysis tables
docs/         full protocol
```

After extracting the Zenodo archives into the repository root, `dataset/entries/`, `af3/results/` and `eval/raw/`
are added (they are ignored by Git).

## Requirements

* Linux (the benchmark ran on Ubuntu 26.04 under WSL2; AMD Ryzen 9 5900X, 64 GB RAM, NVIDIA RTX 3090 24 GB)
* conda; Python 3 for the dataset and analysis scripts (standard library only, plus `requests` for the download step,
  which uses MICA's EMDB helper to read contour levels)
* Phenix 1.20.1 (docking and refinement for MICA and EModelX(+AF); chain comparison metric)
* US-align (version 20261012 used)
* The five methods installed as described in `docs/PROTOCOL.md`

Locations are set in `scripts/config.sh` (benchmark root, method installs, run folders, Phenix, US-align, conda) and
can be overridden with environment variables.

## Reproducing the benchmark

```bash
# 1. Dataset (writes dataset/; uses the RCSB search and data APIs)
python3 scripts/01_fetch_candidates.py        # 2026 EM entries <= 4 A without nucleic acid, with metadata
python3 scripts/03_fetch_clusters.py          # sequences of candidate entries
python3 scripts/06_build_novel_set.py         # filters, 25% novelty, amyloid rule, redundancy, stratified sample (seed 2026)
MICA_REPO=/path/to/MICA python3 scripts/07_download_entries.py   # maps, FASTA, reference structures, contour levels
python3 scripts/11_download_assemblies.py     # biological assemblies for the 3 symmetric entries
python3 scripts/12_phenix_compatible_refs.py  # reference copies readable by Phenix 1.20.1
python3 scripts/10_residue_counts.py          # deposited and modeled residue counts

# 2. AlphaFold3
python3 scripts/08_make_af3_jobs.py           # job files (already provided in af3/)
#    upload af3/upload_day1.json and af3/upload_day2.json at alphafoldserver.com (daily quota 30 jobs), submit,
#    download the results and put each folder in af3/results/fold_<job>/ (or use the Zenodo copy)
python3 scripts/09_ready_entries.py           # entries whose AF3 results are present -> runs/ready_entries.tsv
python3 scripts/13_prepare_emprot_inputs.py   # one FASTA record and one template per chain copy for EMProt

# 3. Methods: one at a time; each batch is resumable and has a 12 hour per map watchdog
R=${RUNS_ROOT:-$HOME/bench/runs}; mkdir -p $R
bash scripts/run_mica_batch.sh > $R/mica_batch.log 2>&1 &
bash scripts/watchdog.sh $R/mica_batch.log $R/mica 12 > $R/mica_watchdog.log 2>&1 &
#    then the same for emodelx, modelangelo, cryoatom and emprot (run_<method>_batch.sh, <method>_batch.log)

# 4. Scoring and analysis
bash scripts/score_raw.sh mica       # also: emodelx_af, modelangelo, cryoatom, emprot
python3 scripts/parse_scores.py      # eval/scores.csv
python3 scripts/summarize_scores.py  # eval/summary.csv
python3 scripts/analyze.py           # eval/analysis.md, eval/per_entry_wide.csv
```

`bash scripts/status.sh` prints the state of all batches, failures and scoring at any time.

## Data licenses

* Code: MIT License.
* EMDB maps and PDB structures and sequences: CC0 (from the wwPDB archives).
* AlphaFold3 predictions: generated with AlphaFold Server and **subject to the AlphaFold Server Output Terms of Use**
  (https://alphafoldserver.com/output-terms); non commercial use only. They are not covered by the licenses above.

## Methods benchmarked

* MICA: Gyawali, Dhakal & Cheng, *Commun. Chem.* 8, 320 (2025). github.com/jianlin-cheng/MICA
* EModelX(+AF): Chen et al., *Nat. Commun.* 15, 8808 (2024). github.com/biomed-AI/EModelX
* ModelAngelo: Jamali et al., *Nature* 628, 450 (2024). github.com/3dem/model-angelo
* CryoAtom: Su et al., *Nat. Struct. Mol. Biol.* (2025). github.com/YangLab-SDU/CryoAtom
* EMProt: Li et al., *Nat. Struct. Mol. Biol.* (2025). github.com/huang-laboratory/EMProt
* Tools: Phenix (Liebschner et al. 2019), US-align (Zhang et al. 2022), AlphaFold3 (Abramson et al. 2024)
