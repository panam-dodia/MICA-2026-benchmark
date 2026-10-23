"""List benchmark entries whose AF3 jobs are all downloaded, for the method batch scripts.

Writes runs/ready_entries.tsv (tab-separated, LF line endings, sorted by residue count):
  emdb  pdb  contour  resolution  residues  jobs(comma-separated AF3 job names)
An entry is ready when every one of its AF3 jobs has af3/results/fold_<job>/fold_<job>_model_0.cif.
"""
import csv
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent


def main():
    details = {r["pdb_id"]: r for r in csv.DictReader(open(ROOT / "dataset" / "benchmark_30_details.csv", encoding="utf-8"))}
    jobs = {}
    for m in csv.DictReader(open(ROOT / "af3" / "job_manifest.csv", encoding="utf-8")):
        jobs.setdefault(m["pdb_id"], []).append(m["job_name"])
    ready, waiting = [], []
    for pdb, r in details.items():
        have = all((ROOT / "af3" / "results" / f"fold_{j}" / f"fold_{j}_model_0.cif").exists() for j in jobs[pdb])
        (ready if have else waiting).append(pdb)
    ready.sort(key=lambda p: int(details[p]["residue_count"]))
    out = ROOT / "runs" / "ready_entries.tsv"
    out.parent.mkdir(exist_ok=True)
    with open(out, "w", newline="\n", encoding="utf-8") as fh:
        fh.write("emdb\tpdb\tcontour\tresolution\tresidues\tjobs\n")
        for p in ready:
            r = details[p]
            fh.write(f"{r['emdb_id'].split('-')[1]}\t{p}\t{r['contour_level']}\t{r['emdb_resolution']}\t{r['residue_count']}\t{','.join(sorted(jobs[p]))}\n")
    print(f"ready: {len(ready)} -> {out}")
    print(f"waiting for AF3: {sorted(waiting)}")


if __name__ == "__main__":
    main()
