"""Build the novel 2026 benchmark set and the stratified random sample. No downloads of maps/models here.

Filters, in order (count printed after each):
  0. RCSB search: EM, released >= 2026-01-01, resolution <= 4.0 A, no nucleic acid entities
  1. single-particle reconstruction
  2. exactly one EMDB map
  3. 1.5 <= resolution <= 4.0 A
  4. all polymer entities are protein
  5. no non-standard residues
  6. 300 <= modeled residues <= 2500
  7. <= 4 unique protein entities (bounds AF3 jobs)
  8. novelty: no protein entity (>= MIN_LEN aa) has > 25% identity to any PDB entity released before 2026-01-01
 8b. exclude amyloid/fibril entries (title or PDB keywords) and entries whose longest chain is < 50 residues
  9. de-duplication: one entry per 25%-identity cluster within the novel set
 10. stratified random sample: 10 entries per size bin, fixed seed

Outputs (dataset/):
  filter_counts.txt, novelty_scan.jsonl (resumable cache), novel_set.csv, clusters.json, benchmark_30.csv, sample_seed.txt
"""
import csv
import json
import random
import re
import time
from pathlib import Path

from importlib import util as _u
_spec = _u.spec_from_file_location("fetch", Path(__file__).with_name("01_fetch_candidates.py"))
fetch = _u.module_from_spec(_spec)
_spec.loader.exec_module(fetch)

DATA = Path(__file__).resolve().parent.parent / "dataset"
SEED = 2026
IDENTITY = 0.25
NOVELTY_DATE = "2026-01-01"
MIN_LEN = 30  # shorter peptides are not searched (RCSB MMseqs2 search is unreliable for very short queries)
SIZE_BINS = [(300, 700), (700, 1400), (1400, 2500)]
PER_BIN = 10
MIN_LONGEST_CHAIN = 50
FIBRIL_RE = re.compile(r"\b(amyloids?|amyloidogenic|fibrils?|fibrillar)\b", re.IGNORECASE)  # whole words: not "fibrillin"

KEYWORDS_QUERY = """
query($ids: [String!]!) {
  entries(entry_ids: $ids) { rcsb_id struct_keywords { pdbx_keywords text } }
}
"""

EXTRA_QUERY = """
query($ids: [String!]!) {
  entries(entry_ids: $ids) {
    rcsb_id
    em_experiment { reconstruction_method aggregation_state }
  }
}
"""


def sequence_search(seq, before=None, on_or_after=None, return_ids=False):
    nodes = [{"type": "terminal", "service": "sequence", "parameters": {
        "evalue_cutoff": 0.1, "identity_cutoff": IDENTITY, "sequence_type": "protein", "value": seq}}]
    if before:
        nodes.append({"type": "terminal", "service": "text", "parameters": {
            "attribute": "rcsb_accession_info.initial_release_date", "operator": "less", "value": before}})
    if on_or_after:
        nodes.append({"type": "terminal", "service": "text", "parameters": {
            "attribute": "rcsb_accession_info.initial_release_date", "operator": "greater_or_equal", "value": on_or_after}})
    opts = {"return_all_hits": True} if return_ids else {"paginate": {"start": 0, "rows": 0}}
    q = {"query": {"type": "group", "logical_operator": "and", "nodes": nodes},
         "return_type": "polymer_entity", "request_options": opts}
    res = fetch.post_json(fetch.SEARCH_URL, q) or {}
    time.sleep(0.15)
    if return_ids:
        return [h["identifier"] for h in res.get("result_set", [])]
    return res.get("total_count", 0)


def as_num(v, cast=float):
    try:
        return cast(v)
    except (TypeError, ValueError):
        return None


def main():
    log = []

    def report(name, n):
        line = f"{name:70s} {n:5d}"
        print(line, flush=True)
        log.append(line)

    rows = list(csv.DictReader(open(DATA / "candidates_2026.csv", encoding="utf-8")))
    report("0. EM, released 2026, <=4.0 A, no nucleic acid (RCSB search)", len(rows))

    # 1. single particle: needs em_experiment, fetched for all candidates
    method = {}
    ids = [r["pdb_id"] for r in rows]
    for i in range(0, len(ids), 200):
        res = fetch.post_json(fetch.GRAPHQL_URL, {"query": EXTRA_QUERY, "variables": {"ids": ids[i:i + 200]}})
        for e in res["data"]["entries"]:
            em = e.get("em_experiment") or {}
            method[e["rcsb_id"]] = (em.get("reconstruction_method") or "", em.get("aggregation_state") or "")
        time.sleep(0.3)
    for r in rows:
        r["reconstruction_method"], r["aggregation_state"] = method.get(r["pdb_id"], ("", ""))
    rows = [r for r in rows if r["reconstruction_method"].upper() == "SINGLE PARTICLE"]
    report("1. single-particle reconstruction", len(rows))

    steps = [
        ("2. exactly one EMDB map", lambda r: as_num(r["n_emdb"], int) == 1),
        ("3. 1.5 <= resolution <= 4.0 A", lambda r: as_num(r["resolution"]) is not None and 1.5 <= as_num(r["resolution"]) <= 4.0),
        ("4. all polymer entities are protein", lambda r: r["n_protein_entities"] == r["n_polymer_entities"]),
        ("5. no non-standard residues", lambda r: as_num(r["n_nonstd_residues"], int) == 0),
        ("6. 300 <= deposited polymer residues <= 2500", lambda r: as_num(r["n_residues_modeled"], int) is not None and 300 <= as_num(r["n_residues_modeled"], int) <= 2500),
        ("7. <= 4 unique protein entities", lambda r: as_num(r["n_protein_entities"], int) <= 4),
    ]
    for name, fn in steps:
        rows = [r for r in rows if fn(r)]
        report(name, len(rows))

    # sequences come from the earlier fetch (covers 300-5000 residues, so a superset of this set)
    ents = json.loads((DATA / "filtered_entities.json").read_text())
    missing = [r["pdb_id"] for r in rows if r["pdb_id"] not in ents]
    assert not missing, f"no sequences cached for {missing[:5]}"

    # 8. novelty scan (resumable). Per entry: stop at the first entity with an older homolog.
    cache_path = DATA / "novelty_scan.jsonl"
    cache = {}
    if cache_path.exists():
        for line in cache_path.read_text().splitlines():
            if line.strip():
                d = json.loads(line)
                cache[d["pdb_id"]] = d
    with open(cache_path, "a") as fh:
        for i, r in enumerate(rows, 1):
            pdb = r["pdb_id"]
            if pdb in cache:
                continue
            checked, novel, short = [], True, []
            for e in ents[pdb]["entities"]:
                if e["len"] < MIN_LEN:
                    short.append(e["entity_id"])
                    continue
                n = sequence_search(e["seq"], before=NOVELTY_DATE)
                checked.append({"entity_id": e["entity_id"], "len": e["len"], "n_older_homologs": n})
                if n > 0:
                    novel = False
                    break
            cache[pdb] = {"pdb_id": pdb, "novel": novel, "entities": checked, "short_unsearched": short}
            fh.write(json.dumps(cache[pdb]) + "\n")
            fh.flush()
            if i % 100 == 0:
                print(f"   novelty scan {i}/{len(rows)}", flush=True)
    n_before = len(rows)
    rows = [r for r in rows if cache[r["pdb_id"]]["novel"]]
    report("8. novel: no entity >25% identity to PDB entries released before 2026", len(rows))
    log.append(f"   -> {n_before - len(rows)}/{n_before} ({100 * (n_before - len(rows)) / n_before:.1f}%) entries have an older >25% homolog")
    print(log[-1])

    # 8b. amyloid/fibril and short-chain exclusion
    kw = {}
    novel_list = [r["pdb_id"] for r in rows]
    for i in range(0, len(novel_list), 200):
        res = fetch.post_json(fetch.GRAPHQL_URL, {"query": KEYWORDS_QUERY, "variables": {"ids": novel_list[i:i + 200]}})
        for e in res["data"]["entries"]:
            k = e.get("struct_keywords") or {}
            kw[e["rcsb_id"]] = " ".join(x for x in (k.get("pdbx_keywords"), k.get("text")) if x)
    excluded = []
    for r in rows:
        pdb = r["pdb_id"]
        reasons = []
        hit = FIBRIL_RE.search(r["title"]) or FIBRIL_RE.search(kw.get(pdb, ""))
        if hit:
            reasons.append(f"amyloid/fibril ('{hit.group(0)}')")
        longest = max(e["len"] for e in ents[pdb]["entities"])
        if longest < MIN_LONGEST_CHAIN:
            reasons.append(f"longest chain {longest} < {MIN_LONGEST_CHAIN}")
        if reasons:
            excluded.append({"pdb_id": pdb, "reasons": "; ".join(reasons), "title": r["title"]})
    excluded_ids = {x["pdb_id"] for x in excluded}
    rows = [r for r in rows if r["pdb_id"] not in excluded_ids]
    report("8b. exclude amyloid/fibril and longest chain < 50 residues", len(rows))
    for x in excluded:
        line = f"   removed {x['pdb_id']}: {x['reasons']} | {x['title']}"
        print(line)
        log.append(line)
    (DATA / "excluded_8b.json").write_text(json.dumps(excluded, indent=1))

    # 9. de-duplicate: link novel entries that share a >25%-identity entity hit among 2026 releases
    novel_ids = {r["pdb_id"] for r in rows}
    parent = {p: p for p in novel_ids}

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    # Every query is itself a 2026 entity, so a valid search always returns its own entity.
    # A result without it is a failed/empty server reply: retry instead of trusting it.
    # Validated results are cached so reruns produce identical clusters.
    hits_path = DATA / "dedup_hits.json"
    hits_cache = json.loads(hits_path.read_text()) if hits_path.exists() else {}
    for pdb in sorted(novel_ids):
        for e in ents[pdb]["entities"]:
            if e["len"] < MIN_LEN:
                continue
            if e["entity_id"] not in hits_cache:
                for attempt in range(6):
                    hits = sequence_search(e["seq"], on_or_after=NOVELTY_DATE, return_ids=True)
                    if e["entity_id"] in hits:
                        break
                    print(f"   dedup search for {e['entity_id']} missing self-hit (attempt {attempt + 1}); retrying", flush=True)
                    time.sleep(3 * (attempt + 1))
                else:
                    raise RuntimeError(f"dedup search for {e['entity_id']} never returned its own entity")
                hits_cache[e["entity_id"]] = sorted(hits)
                hits_path.write_text(json.dumps(hits_cache, indent=1))
            for hit in hits_cache[e["entity_id"]]:
                other = hit.split("_")[0]
                if other in novel_ids and other != pdb:
                    parent[find(other)] = find(pdb)
    clusters = {}
    for p in novel_ids:
        clusters.setdefault(find(p), []).append(p)
    clusters = sorted(sorted(c) for c in clusters.values())
    (DATA / "clusters.json").write_text(json.dumps(clusters, indent=1))
    rng = random.Random(SEED)
    reps = {rng.choice(c) for c in clusters}  # one random representative per cluster (seeded)
    rows = [r for r in rows if r["pdb_id"] in reps]
    report("9. one entry per 25%-identity cluster", len(rows))

    with open(DATA / "novel_set.csv", "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(sorted(rows, key=lambda r: r["pdb_id"]))

    # 10. stratified random sample
    sample = []
    for lo, hi in SIZE_BINS:
        pool = sorted((r for r in rows if lo <= int(r["n_residues_modeled"]) < hi or (hi == 2500 and int(r["n_residues_modeled"]) == 2500)),
                      key=lambda r: r["pdb_id"])
        report(f"   size bin {lo}-{hi} residues: available", len(pool))
        pick = rng.sample(pool, min(PER_BIN, len(pool)))
        for r in pick:
            r["size_bin"] = f"{lo}-{hi}"
        sample += pick
    report("10. stratified random sample", len(sample))

    fields = ["pdb_id", "emdb_ids", "resolution", "n_residues_modeled", "n_protein_entities", "size_bin", "release_date", "title"]
    with open(DATA / "benchmark_30.csv", "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh)
        w.writerow(["pdb_id", "emdb_id", "resolution", "residue_count", "unique_chain_count", "size_bin", "release_date", "title"])
        for r in sorted(sample, key=lambda r: (int(r["n_residues_modeled"]), r["pdb_id"])):
            w.writerow([r[k] for k in fields])
    (DATA / "sample_seed.txt").write_text(
        f"seed={SEED}\nrng=random.Random(seed); used first for cluster representatives (sorted clusters), then per size bin (pools sorted by pdb_id) in order {SIZE_BINS}\n")
    (DATA / "filter_counts.txt").write_text("\n".join(log) + "\n")
    print("done", flush=True)


if __name__ == "__main__":
    main()
