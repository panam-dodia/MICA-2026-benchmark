"""Create AlphaFold Server job files: one single-chain job per unique protein entity, as in the MICA paper.

Input:  dataset/benchmark_30_details.csv and each entry's RCSB FASTA.
Output: af3/jobs/<pdbid>_<n>.json      one job per file (AlphaFold Server upload format)
        af3/jobs_all.json              all jobs in one list
        af3/job_manifest.csv           job name -> entry, entity, length, expected MICA results folder

Job settings (recorded for the report): dialect "alphafoldserver", version 1, modelSeeds ["2026"],
useStructureTemplate true, maxTemplateDate "2025-02-03" (the server maximum; see docs/PROTOCOL.md).
Like MICA's utils/fasta_to_AF3_json.py, 'X' characters are removed from sequences.
"""
import csv
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
DATA = ROOT / "dataset"
OUT = ROOT / "af3"
SEED = "2026"
MAX_TEMPLATE_DATE = "2025-02-03"


def read_fasta(path):
    records, header, seq = [], None, []
    for line in path.read_text().splitlines():
        if line.startswith(">"):
            if header:
                records.append((header, "".join(seq)))
            header, seq = line[1:].strip(), []
        elif line.strip():
            seq.append(line.strip())
    if header:
        records.append((header, "".join(seq)))
    return records


def main():
    (OUT / "jobs").mkdir(parents=True, exist_ok=True)
    rows = list(csv.DictReader(open(DATA / "benchmark_30_details.csv", encoding="utf-8")))
    all_jobs, manifest = [], []
    for r in rows:
        entity_records = read_fasta(ROOT / r["fasta_file"])
        if len(entity_records) != int(r["unique_chain_count"]):
            print(f"WARNING {r['pdb_id']}: FASTA has {len(entity_records)} entities, expected {r['unique_chain_count']}")
        for header, seq in entity_records:
            entity = header.split("|")[0]                  # e.g. 9I8C_1
            name = entity.lower()                          # MICA folder name: AF3_results/9i8c_1/
            clean = seq.replace("X", "").replace("x", "")
            job = {
                "name": name,
                "modelSeeds": [SEED],
                "sequences": [{"proteinChain": {
                    "sequence": clean, "count": 1,
                    "useStructureTemplate": True, "maxTemplateDate": MAX_TEMPLATE_DATE}}],
                "dialect": "alphafoldserver",
                "version": 1,
            }
            (OUT / "jobs" / f"{name}.json").write_text(json.dumps([job], indent=2))
            all_jobs.append(job)
            manifest.append({
                "job_name": name, "pdb_id": r["pdb_id"], "emdb_id": r["emdb_id"], "entity": entity,
                "length": len(clean), "x_removed": len(seq) - len(clean), "chains": header.split("|")[1] if "|" in header else "",
                "mica_results_dir": f"dataset/entries/{r['emdb_id'].split('-')[1]}/AF3_results/{name}/",
                "submitted": "", "downloaded": "",
            })
    (OUT / "jobs_all.json").write_text(json.dumps(all_jobs, indent=2))
    with open(OUT / "job_manifest.csv", "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=list(manifest[0].keys()))
        w.writeheader()
        w.writerows(manifest)
    total = sum(m["length"] for m in manifest)
    print(f"{len(all_jobs)} jobs, {total} residues total, longest {max(m['length'] for m in manifest)}; files in {OUT}")


if __name__ == "__main__":
    main()
