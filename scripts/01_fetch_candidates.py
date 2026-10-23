"""Fetch 2026 cryo-EM protein-only PDB entries (<=4 A) and their metadata from RCSB.

Writes dataset/candidates_2026.csv with one row per PDB entry.
Standard library only.
"""
import csv
import json
import sys
import time
import urllib.request
from pathlib import Path

OUT_DIR = Path(__file__).resolve().parent.parent / "dataset"
SEARCH_URL = "https://search.rcsb.org/rcsbsearch/v2/query"
GRAPHQL_URL = "https://data.rcsb.org/graphql"
RELEASE_FROM = "2026-01-01"
MAX_RES = 4.0

SEARCH_QUERY = {
    "query": {"type": "group", "logical_operator": "and", "nodes": [
        {"type": "terminal", "service": "text", "parameters": {
            "attribute": "exptl.method", "operator": "exact_match", "value": "ELECTRON MICROSCOPY"}},
        {"type": "terminal", "service": "text", "parameters": {
            "attribute": "rcsb_accession_info.initial_release_date", "operator": "greater_or_equal", "value": RELEASE_FROM}},
        {"type": "terminal", "service": "text", "parameters": {
            "attribute": "rcsb_entry_info.resolution_combined", "operator": "less_or_equal", "value": MAX_RES}},
        {"type": "terminal", "service": "text", "parameters": {
            "attribute": "rcsb_entry_info.polymer_entity_count_nucleic_acid", "operator": "equals", "value": 0}},
    ]},
    "return_type": "entry",
    "request_options": {"return_all_hits": True},
}

GRAPHQL_QUERY = """
query($ids: [String!]!) {
  entries(entry_ids: $ids) {
    rcsb_id
    struct { title }
    rcsb_accession_info { initial_release_date }
    rcsb_entry_container_identifiers { emdb_ids }
    rcsb_entry_info {
      resolution_combined
      deposited_polymer_monomer_count
      polymer_entity_count_protein
      deposited_polymer_entity_instance_count
      deposited_nonpolymer_entity_instance_count
    }
    polymer_entities {
      entity_poly { rcsb_entity_polymer_type rcsb_sample_sequence_length rcsb_non_std_monomer_count }
      rcsb_polymer_entity_container_identifiers { auth_asym_ids }
    }
  }
}
"""


def post_json(url, payload, retries=4):
    data = json.dumps(payload).encode()
    for attempt in range(retries):
        try:
            req = urllib.request.Request(url, data=data, headers={"Content-Type": "application/json"})
            with urllib.request.urlopen(req, timeout=120) as resp:
                body = resp.read().decode()
                return json.loads(body) if body.strip() else None  # RCSB search: HTTP 204 = no hits
        except Exception as exc:  # network hiccups: back off and retry
            if attempt == retries - 1:
                raise
            print(f"  retry {attempt + 1} after error: {exc}", file=sys.stderr)
            time.sleep(2 ** attempt)


def main():
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    ids = [hit["identifier"] for hit in post_json(SEARCH_URL, SEARCH_QUERY)["result_set"]]
    print(f"search hits: {len(ids)}")

    rows = []
    batch = 200
    for i in range(0, len(ids), batch):
        chunk = ids[i:i + batch]
        result = post_json(GRAPHQL_URL, {"query": GRAPHQL_QUERY, "variables": {"ids": chunk}})
        for e in result["data"]["entries"]:
            info = e.get("rcsb_entry_info") or {}
            ents = e.get("polymer_entities") or []
            emdb = (e.get("rcsb_entry_container_identifiers") or {}).get("emdb_ids") or []
            res = info.get("resolution_combined") or [None]
            protein_ents = [p for p in ents if (p.get("entity_poly") or {}).get("rcsb_entity_polymer_type") == "Protein"]
            non_std = sum((p.get("entity_poly") or {}).get("rcsb_non_std_monomer_count") or 0 for p in protein_ents)
            rows.append({
                "pdb_id": e["rcsb_id"],
                "emdb_ids": ";".join(emdb),
                "n_emdb": len(emdb),
                "release_date": ((e.get("rcsb_accession_info") or {}).get("initial_release_date") or "")[:10],
                "resolution": res[0],
                "n_residues_modeled": info.get("deposited_polymer_monomer_count"),
                "n_chains": info.get("deposited_polymer_entity_instance_count"),
                "n_protein_entities": info.get("polymer_entity_count_protein"),
                "n_polymer_entities": len(ents),
                "n_ligand_instances": info.get("deposited_nonpolymer_entity_instance_count"),
                "max_entity_seq_len": max(((p.get("entity_poly") or {}).get("rcsb_sample_sequence_length") or 0) for p in protein_ents) if protein_ents else 0,
                "n_nonstd_residues": non_std,
                "title": ((e.get("struct") or {}).get("title") or "").replace("\n", " "),
            })
        print(f"  metadata {min(i + batch, len(ids))}/{len(ids)}")
        time.sleep(0.3)

    out = OUT_DIR / "candidates_2026.csv"
    with open(out, "w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)
    print(f"wrote {len(rows)} rows -> {out}")


if __name__ == "__main__":
    main()
