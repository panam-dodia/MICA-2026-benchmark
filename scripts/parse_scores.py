"""Parse raw scoring outputs (scripts/score_raw.sh) into one table: eval/scores.csv.

Per method and entry (rules: docs/PROTOCOL.md):
  tm_score        US-align TM-score normalized by the deposited structure ("Structure_2" line)
  aligned_len     US-align aligned length;  rmsd: US-align RMSD;  seq_identity: US-align Seq_ID
  ca_match        chain_comparison with the MODEL as target: FOUND % (model residues with a reference CA <= 3 A)
  seq_match       same run: SEQ MATCH(%)
  recall          chain_comparison with the REFERENCE as target: FOUND %
  n_model, n_ref  residue counts from "Original residues in target: A  In query: B" (precision run: target = model)
  ca_quality      ca_match * n_model / n_ref   (Cryo2Struct definition; n_ref = in-place reference, i.e. the
                  biological assembly for 9I8C, 9SZB, 12ZE). For those 3 entries it is left empty and written
                  to ca_quality_vs_assembly instead, because it mainly reflects how many copies were built.
Entries without a model are written with status "no model" so failures stay visible.
"""
import csv
import re
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
RAW = ROOT / "eval" / "raw"
METHODS = ["mica", "modelangelo", "cryoatom", "emprot", "emodelx_af"]
SYMMETRIC = {"9I8C", "9SZB", "12ZE"}


def parse_usalign(text):
    out = {}
    m = re.search(r"Aligned length=\s*(\d+), RMSD=\s*([\d.]+), Seq_ID=n_identical/n_aligned=\s*([\d.]+)", text)
    if m:
        out.update(aligned_len=int(m.group(1)), rmsd=float(m.group(2)), seq_identity=float(m.group(3)))
    m = re.search(r"TM-score=\s*([\d.]+) \(normalized by length of Structure_2: L=(\d+)", text)
    if m:
        out.update(tm_score=float(m.group(1)), usalign_ref_len=int(m.group(2)))
    return out


def parse_chaincomp(text):
    """Return FOUND %, SEQ MATCH %, target and query residue counts from a chain_comparison log."""
    out = {}
    m = re.search(r"Original residues in target:\s*(\d+)\s+In query:\s*(\d+)", text)
    if m:
        out.update(n_target=int(m.group(1)), n_query=int(m.group(2)))
    # result row is labelled Unique_target or Entire_target, whichever variant matched more (see evaluation.md)
    rows = [l for l in text.splitlines() if l.strip().startswith(("Unique_target", "Entire_target"))]
    if rows:
        f = rows[-1].split()
        # Unique_target RMSD Close_N Far_N Forward_N Reverse_N Mixed_N FOUND% CA_SCORE SEQ_MATCH% SEQ_SCORE MEAN_LEN FRAGS BAD
        out.update(found=float(f[7]), seq_match=float(f[9]), close_rmsd=float(f[1]), close_n=int(f[2]),
                   variant=f[0])
    return out


def main():
    entries = list(csv.DictReader(open(ROOT / "runs" / "ready_entries.tsv", encoding="utf-8"), delimiter="\t"))
    rows = []
    for method in METHODS:
        mdir = RAW / method
        if not mdir.exists():
            continue
        for e in entries:
            n, pdb = e["emdb"], e["pdb"]
            row = {"method": method, "emdb": n, "pdb": pdb, "symmetric_ref": pdb in SYMMETRIC}
            us, pr, rc = (mdir / f"{n}.{s}.txt" for s in ("usalign", "cc_precision", "cc_recall"))
            if not (us.exists() and pr.exists() and rc.exists()):
                row["status"] = "no model"
                rows.append(row)
                continue
            u, p, r = parse_usalign(us.read_text()), parse_chaincomp(pr.read_text()), parse_chaincomp(rc.read_text())
            row.update({k: u.get(k) for k in ("tm_score", "aligned_len", "rmsd", "seq_identity", "usalign_ref_len")})
            row.update(ca_match=p.get("found"), seq_match=p.get("seq_match"), recall=r.get("found"),
                       n_model=p.get("n_target"), n_ref=p.get("n_query"),
                       cc_variant_precision=p.get("variant"), cc_variant_recall=r.get("variant"))
            if row["ca_match"] is not None and row["n_model"] and row["n_ref"]:
                q = round(row["ca_match"] * row["n_model"] / row["n_ref"], 2)
                # For symmetric entries the value mostly reflects how many copies a method built (MICA follows the
                # FASTA chain count, ModelAngelo builds every copy), so it is n/a in the main column (evaluation.md).
                row["ca_quality_vs_assembly" if pdb in SYMMETRIC else "ca_quality"] = q
            missing = [k for k in ("tm_score", "ca_match", "seq_match", "recall") if row.get(k) is None]
            row["status"] = "ok" if not missing else "parse error: " + ",".join(missing)
            rows.append(row)
    fields = ["method", "emdb", "pdb", "status", "symmetric_ref", "tm_score", "aligned_len", "rmsd", "seq_identity",
              "ca_match", "seq_match", "ca_quality", "ca_quality_vs_assembly", "recall", "n_model", "n_ref", "usalign_ref_len", "cc_variant_precision", "cc_variant_recall"]
    out = ROOT / "eval" / "scores.csv"
    with open(out, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=fields)
        w.writeheader()
        w.writerows(rows)
    for r in rows:
        if r["status"] != "no model":
            print({k: r.get(k) for k in fields})
    print(f"wrote {len(rows)} rows -> {out}")


if __name__ == "__main__":
    main()
