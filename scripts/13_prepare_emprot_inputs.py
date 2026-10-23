"""Prepare EMProt inputs: one FASTA record and one AF3 template per chain copy.

EMProt (huang-laboratory/EMProt v1.2) ignores FASTA headers, treats every record as one chain, and infers homo-oligomer
copies from repeated records (official 8AVV example: two identical records plus two identical template files). The RCSB
FASTA has one record per unique entity with its chain list in the header ("Chains A, B[auth C]"), so each entity is
expanded to one record per listed chain. The template for each copy is that entity's AF3 model_0.cif, the same file
MICA and EModelX(+AF) use. Chain counts therefore match the deposited model, as for every other method.

Writes runs/emprot_inputs/<EMDB>/seq.fa and templates.txt (one AF3 path per line, relative to the benchmark root).
"""
import csv
import re
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent


def read_fasta(path):
    recs, h, s = [], None, []
    for line in path.read_text().splitlines():
        if line.startswith(">"):
            if h:
                recs.append((h, "".join(s)))
            h, s = line[1:].strip(), []
        elif line.strip():
            s.append(line.strip())
    if h:
        recs.append((h, "".join(s)))
    return recs


def chain_ids(header):
    """'9I8C_1|Chains A, B[auth C]|...' -> ['A', 'B']  (label IDs; only the count matters to EMProt)."""
    field = header.split("|")[1]
    field = re.sub(r"^Chains?\s+", "", field.strip())
    return [re.sub(r"\[.*?\]", "", c).strip() for c in field.split(",") if c.strip()]


def main():
    entries = list(csv.DictReader(open(ROOT / "runs" / "ready_entries.tsv", encoding="utf-8"), delimiter="\t"))
    for e in entries:
        n, pdb = e["emdb"], e["pdb"].lower()
        recs = read_fasta(ROOT / "dataset" / "entries" / n / f"{pdb}.fasta")
        out = ROOT / "runs" / "emprot_inputs" / n
        out.mkdir(parents=True, exist_ok=True)
        fa, tmpl = [], []
        for header, seq in recs:
            job = header.split("|")[0].lower()                     # e.g. 9i8c_1
            af3 = Path("af3") / "results" / f"fold_{job}" / f"fold_{job}_model_0.cif"
            assert (ROOT / af3).exists(), af3
            for c in chain_ids(header):
                fa.append(f">{job}_chain{c}\n{seq}\n")
                tmpl.append(af3.as_posix())
        (out / "seq.fa").write_text("".join(fa), newline="\n")
        (out / "templates.txt").write_text("\n".join(tmpl) + "\n", newline="\n")
        print(f"{e['pdb']:5} entities={len(recs)} chain copies={len(fa)}")


if __name__ == "__main__":
    main()
