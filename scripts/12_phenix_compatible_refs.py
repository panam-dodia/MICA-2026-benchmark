"""Write Phenix-1.20-readable copies of reference structures that contain residue names longer than 3 characters.

Phenix 1.20 (2022) stores residue names in 3 characters. Structures released with the newer 5-character
chemical-component IDs (ligands) fail: "string is too long for target variable (maximum length is 3 characters,
5 given)", found for 9V3Z on 2026-10-19. For each reference mmCIF (deposited and assembly1), this writes
<name>_phenix.cif with every _atom_site row whose label_comp_id or auth_comp_id is longer than 3 characters removed.
The benchmark entries have no non-standard polymer residues (filter step 5), so only ligands are removed and the
protein coordinates are unchanged. Used for phenix.chain_comparison only; US-align reads the original files.
"""
import csv
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent


def filter_long_comp_ids(path):
    """Return (filtered_lines, removed_names) or (None, set()) if nothing needs removing."""
    lines = path.read_text().splitlines(keepends=True)
    out, removed = [], {}
    i = 0
    while i < len(lines):
        line = lines[i]
        if line.strip() == "loop_" and i + 1 < len(lines) and lines[i + 1].startswith("_atom_site."):
            out.append(line)
            i += 1
            cols = []
            while i < len(lines) and lines[i].startswith("_atom_site."):
                cols.append(lines[i].strip())
                out.append(lines[i])
                i += 1
            idx = [cols.index(c) for c in ("_atom_site.label_comp_id", "_atom_site.auth_comp_id") if c in cols]
            while i < len(lines) and not lines[i].startswith(("#", "loop_", "_")):
                f = lines[i].split()
                long_names = {f[k] for k in idx if k < len(f) and len(f[k]) > 3}
                if long_names:
                    for n in long_names:
                        removed[n] = removed.get(n, 0) + 1
                else:
                    out.append(lines[i])
                i += 1
            continue
        out.append(line)
        i += 1
    return (out, removed) if removed else (None, removed)


def main():
    rows = list(csv.DictReader(open(ROOT / "dataset" / "benchmark_30_details.csv", encoding="utf-8")))
    total = 0
    for r in rows:
        d = ROOT / "dataset" / "entries" / r["emdb_id"].split("-")[1]
        pdb = r["pdb_id"].lower()
        for src in (d / f"{pdb}.cif", d / f"{pdb}_assembly1.cif"):
            if not src.exists():
                continue
            filtered, removed = filter_long_comp_ids(src)
            dst = src.with_name(src.stem + "_phenix.cif")
            if filtered is None:
                if dst.exists():
                    dst.unlink()
                continue
            dst.write_text("".join(filtered))
            total += 1
            print(f"{r['pdb_id']}: {src.name} -> {dst.name}; removed atom_site rows by residue name: {removed}")
    print(f"{total} Phenix-compatible reference file(s) written")


if __name__ == "__main__":
    main()
