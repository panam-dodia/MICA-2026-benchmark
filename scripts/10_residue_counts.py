"""Record both size measures for the 30 benchmark entries.

deposited_residues: RCSB rcsb_entry_info.deposited_polymer_monomer_count, the full deposited sequence of every chain,
                    including residues not resolved in the map. Used by 06_build_novel_set.py for the size filter and bins.
modeled_residues:   residues with a CA atom in the deposited mmCIF (what US-align uses as the reference length L).

Writes dataset/benchmark_30_sizes.csv with both counts and both bin assignments.
"""
import csv
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent


def size_bin(n):
    if n < 300:
        return "<300"
    return "300-700" if n < 700 else ("700-1400" if n < 1400 else "1400-2500")


def modeled_count(cif_path):
    """Count distinct (label_asym_id, label_seq_id, ins_code) residues with a CA atom (ATOM records only)."""
    seen = set()
    for line in open(cif_path):
        if line.startswith("ATOM"):
            f = line.split()
            if f[3] == "CA":
                seen.add((f[6], f[8], f[9]))
    return len(seen)


def main():
    rows = list(csv.DictReader(open(ROOT / "dataset" / "benchmark_30_details.csv", encoding="utf-8")))
    out = []
    for r in rows:
        dep = int(r["residue_count"])
        mod = modeled_count(ROOT / r["cif_file"])
        out.append({"pdb_id": r["pdb_id"], "emdb_id": r["emdb_id"],
                    "deposited_residues": dep, "modeled_residues": mod,
                    "modeled_fraction": round(mod / dep, 3),
                    "bin_deposited": r["size_bin"], "bin_modeled": size_bin(mod)})
    with open(ROOT / "dataset" / "benchmark_30_sizes.csv", "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=list(out[0].keys()))
        w.writeheader()
        w.writerows(out)
    changed = [o["pdb_id"] for o in out if o["bin_deposited"] != o["bin_modeled"]]
    counts = {b: sum(1 for o in out if o["bin_modeled"] == b) for b in ("300-700", "700-1400", "1400-2500")}
    print(f"modeled range {min(o['modeled_residues'] for o in out)}-{max(o['modeled_residues'] for o in out)}; "
          f"bins by modeled count {counts}; bin changes: {changed}")


if __name__ == "__main__":
    main()
