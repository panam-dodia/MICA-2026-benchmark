"""Summarize eval/scores.csv per method.

Prints and writes eval/summary.csv with, per method and subset:
  n_scored, n_failed, and means of tm_score, ca_match, seq_match, ca_quality, aligned_len, seq_identity, recall.
Subsets:
  all_scored      every entry the method produced a model for
  non_symmetric   excluding 9I8C, 9SZB, 12ZE (deposited model = asymmetric unit of a symmetric map)
  common          entries every method listed in eval/scores.csv produced a model for (fair head-to-head)
Failures are counted, not imputed. ca_quality is n/a for symmetric entries (see docs/PROTOCOL.md).
"""
import csv
from pathlib import Path
from statistics import mean, median

ROOT = Path(__file__).resolve().parent.parent
METRICS = ["tm_score", "ca_match", "seq_match", "ca_quality", "aligned_len", "seq_identity", "recall"]


def fnum(x):
    try:
        return float(x)
    except (TypeError, ValueError):
        return None


def main():
    rows = list(csv.DictReader(open(ROOT / "eval" / "scores.csv", encoding="utf-8")))
    methods = list(dict.fromkeys(r["method"] for r in rows))
    ok = {m: {r["pdb"]: r for r in rows if r["method"] == m and r["status"] == "ok"} for m in methods}
    failed = {m: sorted(r["pdb"] for r in rows if r["method"] == m and r["status"] != "ok") for m in methods}
    common = set.intersection(*(set(v) for v in ok.values())) if ok else set()
    symmetric = {r["pdb"] for r in rows if r["symmetric_ref"] == "True"}

    out = []
    for m in methods:
        subsets = {"all_scored": set(ok[m]), "non_symmetric": set(ok[m]) - symmetric, "common": common}
        for name, ids in subsets.items():
            rec = {"method": m, "subset": name, "n_scored": len(ids), "n_failed": len(failed[m])}
            for k in METRICS:
                vals = [fnum(ok[m][p][k]) for p in ids if p in ok[m] and fnum(ok[m][p][k]) is not None]
                rec[f"mean_{k}"] = round(mean(vals), 3) if vals else None
                if k == "tm_score":
                    rec["median_tm_score"] = round(median(vals), 3) if vals else None
            out.append(rec)
    fields = ["method", "subset", "n_scored", "n_failed", "mean_tm_score", "median_tm_score", "mean_ca_match",
              "mean_seq_match", "mean_ca_quality", "mean_aligned_len", "mean_seq_identity", "mean_recall"]
    with open(ROOT / "eval" / "summary.csv", "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=fields)
        w.writeheader()
        w.writerows(out)
    for m in methods:
        print(f"{m}: scored {len(ok[m])}, failed/missing {len(failed[m])}: {failed[m]}")
    print(f"common subset ({len(common)}): {sorted(common)}")
    hdr = f"{'method':12} {'subset':14} {'n':>3} {'TM':>6} {'TMmed':>6} {'CAmatch':>8} {'SeqM':>6} {'CAqual':>7} {'AlnLen':>7} {'SeqID':>6} {'Recall':>7}"
    print(hdr)
    for r in out:
        g = lambda k: "" if r[k] is None else r[k]
        print(f"{r['method']:12} {r['subset']:14} {r['n_scored']:>3} {g('mean_tm_score'):>6} {g('median_tm_score'):>6} "
              f"{g('mean_ca_match'):>8} {g('mean_seq_match'):>6} {g('mean_ca_quality'):>7} {g('mean_aligned_len'):>7} "
              f"{g('mean_seq_identity'):>6} {g('mean_recall'):>7}")


if __name__ == "__main__":
    main()
