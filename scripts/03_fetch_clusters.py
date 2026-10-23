"""For entries passing the basic filters, fetch per-entity sequences and 30%-identity cluster IDs from RCSB.

Reads dataset/candidates_2026.csv, writes dataset/filtered_entities.json:
  {pdb_id: {"row": {...}, "entities": [{"entity_id", "chains", "seq", "len", "cluster30"}]}}
"""
import csv
import json
import time
from pathlib import Path

from importlib import util as _u
_spec = _u.spec_from_file_location("fetch", Path(__file__).with_name("01_fetch_candidates.py"))
fetch = _u.module_from_spec(_spec)
_spec.loader.exec_module(fetch)

DATA = Path(__file__).resolve().parent.parent / "dataset"

QUERY = """
query($ids: [String!]!) {
  entries(entry_ids: $ids) {
    rcsb_id
    polymer_entities {
      rcsb_id
      entity_poly { pdbx_seq_one_letter_code_can rcsb_entity_polymer_type }
      rcsb_polymer_entity_container_identifiers { auth_asym_ids }
      rcsb_cluster_membership { cluster_id identity }
    }
  }
}
"""


def passes(r):
    f = lambda k, c=float: c(r[k]) if r[k] not in ("", "None") else None
    return (f("n_emdb", int) == 1 and f("resolution") is not None and 1.5 <= f("resolution") <= 4.0
            and f("n_protein_entities", int) == f("n_polymer_entities", int)
            and f("n_nonstd_residues", int) == 0
            and f("n_residues_modeled", int) is not None and 300 <= f("n_residues_modeled", int) <= 5000
            and f("n_protein_entities", int) <= 4)


def main():
    rows = {r["pdb_id"]: r for r in csv.DictReader(open(DATA / "candidates_2026.csv", encoding="utf-8")) if passes(r)}
    ids = sorted(rows)
    print(f"entries passing filters: {len(ids)}")
    out = {}
    for i in range(0, len(ids), 200):
        res = fetch.post_json(fetch.GRAPHQL_URL, {"query": QUERY, "variables": {"ids": ids[i:i + 200]}})
        for e in res["data"]["entries"]:
            ents = []
            for p in e["polymer_entities"] or []:
                cl = {m["identity"]: m["cluster_id"] for m in (p.get("rcsb_cluster_membership") or [])}
                seq = (p["entity_poly"]["pdbx_seq_one_letter_code_can"] or "").replace("\n", "")
                ents.append({
                    "entity_id": p["rcsb_id"],
                    "chains": p["rcsb_polymer_entity_container_identifiers"]["auth_asym_ids"],
                    "seq": seq,
                    "len": len(seq),
                    "cluster30": cl.get(30),
                })
            out[e["rcsb_id"]] = {"row": rows[e["rcsb_id"]], "entities": ents}
        print(f"  {min(i + 200, len(ids))}/{len(ids)}")
        time.sleep(0.3)
    (DATA / "filtered_entities.json").write_text(json.dumps(out, indent=1))
    n_missing = sum(1 for v in out.values() for en in v["entities"] if en["cluster30"] is None)
    print(f"wrote {len(out)} entries; entities without a 30% cluster id: {n_missing}")


if __name__ == "__main__":
    main()
