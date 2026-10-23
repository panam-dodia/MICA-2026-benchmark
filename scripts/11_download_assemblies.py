"""Download biological assembly 1 for entries whose deposited model covers only part of the symmetric map.

Found 2026-10-19: for 9I8C (1 -> 12 chains, C12), 9SZB (2 -> 24, C12) and 12ZE (2 -> 120, icosahedral) the deposited
model is one asymmetric unit, while the map shows the whole assembly. In-place metrics (phenix.chain_comparison) use
assembly 1 as the reference for these entries, so a correct copy in any symmetric position counts. All other entries
have assembly 1 == deposited model. Writes dataset/entries/<EMDB>/<pdbid>_assembly1.cif.
"""
import csv
import gzip
import os
import shutil
import urllib.request
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
SYMMETRIC = {"9I8C", "9SZB", "12ZE"}
URL = "https://files.rcsb.org/download/{pdb}-assembly1.cif.gz"


def main():
    for r in csv.DictReader(open(ROOT / "dataset" / "benchmark_30_details.csv", encoding="utf-8")):
        if r["pdb_id"] not in SYMMETRIC:
            continue
        d = ROOT / "dataset" / "entries" / r["emdb_id"].split("-")[1]
        out = d / f"{r['pdb_id'].lower()}_assembly1.cif"
        gz = out.with_suffix(".cif.gz")
        with urllib.request.urlopen(URL.format(pdb=r["pdb_id"]), timeout=300) as resp, open(gz, "wb") as fh:
            shutil.copyfileobj(resp, fh)
        with gzip.open(gz, "rb") as fin, open(out, "wb") as fout:
            shutil.copyfileobj(fin, fout)
        gz.unlink()
        chains = {l.split()[6] for l in open(out) if l.startswith("ATOM")}
        print(f"{r['pdb_id']}: {out.relative_to(ROOT)}  {out.stat().st_size / 1e6:.1f} MB, {len(chains)} chains")


if __name__ == "__main__":
    main()
