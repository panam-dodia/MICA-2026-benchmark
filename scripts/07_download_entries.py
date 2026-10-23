"""Download inputs and ground truth for every entry in dataset/benchmark_30.csv.

Per entry, into dataset/entries/<EMDB number>/ (the layout MICA's README uses: input/<map id>/):
  emd_<N>.map        EMDB primary map (downloaded .map.gz, decompressed)
  <pdbid>.fasta      RCSB FASTA (official headers, e.g. >9I8C_1|Chain A|...)
  <pdbid>.cif        deposited model (ground truth for scoring)
Contour level and EMDB resolution come from MICA's utils/emdb_extractor.py (run as-is).

Resumable: existing, valid files are skipped. Writes dataset/benchmark_30_details.csv.
"""
import csv
import gzip
import os
import re
import shutil
import subprocess
import sys
import time
import urllib.request
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
DATA = ROOT / "dataset"
ENTRIES = DATA / "entries"
EMDB_EXTRACTOR = Path(os.environ.get("MICA_REPO", Path.home() / "bench" / "methods" / "MICA")) / "utils" / "emdb_extractor.py"

# wwPDB mirror of the EMDB archive (same files as EBI); EBI was ~0.3 MB/s and stalled on 2026-10-16
MAP_URLS = ["https://files.wwpdb.org/pub/emdb/structures/EMD-{n}/map/emd_{n}.map.gz",
            "https://ftp.ebi.ac.uk/pub/databases/emdb/structures/EMD-{n}/map/emd_{n}.map.gz"]
FASTA_URL = "https://www.rcsb.org/fasta/entry/{pdb}/download"
CIF_URL = "https://files.rcsb.org/download/{pdb}.cif.gz"


def download(urls, dest, retries=3):
    """Try each mirror in order; timeout=60 s applies per socket read, so a stalled connection fails fast."""
    urls = [urls] if isinstance(urls, str) else urls
    tmp = dest.with_suffix(dest.suffix + ".part")
    last = None
    for url in urls:
        for attempt in range(retries):
            try:
                with urllib.request.urlopen(url, timeout=60) as resp, open(tmp, "wb") as fh:
                    shutil.copyfileobj(resp, fh, length=8 * 1024 * 1024)
                tmp.replace(dest)
                return url
            except Exception as exc:
                last = exc
                print(f"    retry {attempt + 1} for {url}: {exc}", flush=True)
                time.sleep(5 * (attempt + 1))
    raise RuntimeError(f"all mirrors failed for {dest.name}: {last}")


def gunzip(src, dest):
    tmp = dest.with_suffix(dest.suffix + ".part")
    with gzip.open(src, "rb") as fin, open(tmp, "wb") as fout:
        shutil.copyfileobj(fin, fout, length=8 * 1024 * 1024)
    tmp.replace(dest)


def is_mrc(path):
    # MRC2014: bytes 208-211 are the machine stamp field label "MAP "
    try:
        with open(path, "rb") as fh:
            fh.seek(208)
            return fh.read(4) == b"MAP "
    except OSError:
        return False


def emdb_params(n):
    out = subprocess.run([sys.executable, "-I", str(EMDB_EXTRACTOR), "--emdb_id", n],
                         capture_output=True, text=True, timeout=120)
    contour = re.search(r"Contour Level:\s*([-\d.eE+]+)", out.stdout)
    res = re.search(r"Resolution:\s*([-\d.eE+]+)", out.stdout)
    return (contour.group(1) if contour else "", res.group(1) if res else "", out.stdout.strip().splitlines()[-1] if out.stdout.strip() else out.stderr.strip())


def fasta_header_issues(path):
    """MICA rule (README): header '>PDBID_n|...' with exactly one underscore in the whole header."""
    issues = []
    for line in path.read_text().splitlines():
        if line.startswith(">"):
            if not re.match(r"^>[0-9A-Za-z]+_\d+\|", line):
                issues.append(f"bad prefix: {line[:60]}")
            elif line.count("_") != 1:
                issues.append(f"{line.count('_')} underscores: {line[:80]}")
    return issues


def main():
    rows = list(csv.DictReader(open(DATA / "benchmark_30.csv", encoding="utf-8")))
    details = []
    for i, r in enumerate(rows, 1):
        pdb = r["pdb_id"]
        n = r["emdb_id"].split("-")[1]
        d = ENTRIES / n
        d.mkdir(parents=True, exist_ok=True)
        print(f"[{i}/{len(rows)}] {pdb} / EMD-{n}", flush=True)

        map_path = d / f"emd_{n}.map"
        if not is_mrc(map_path):
            gz = d / f"emd_{n}.map.gz"
            if not gz.exists():
                src = download([u.format(n=n) for u in MAP_URLS], gz)
                print(f"    map from {src.split('/')[2]}", flush=True)
            gunzip(gz, map_path)
            if not is_mrc(map_path):
                raise RuntimeError(f"{map_path} is not a valid MRC file")
            gz.unlink()

        fasta_path = d / f"{pdb.lower()}.fasta"
        if not fasta_path.exists() or fasta_path.stat().st_size == 0:
            download(FASTA_URL.format(pdb=pdb), fasta_path)

        cif_path = d / f"{pdb.lower()}.cif"
        if not cif_path.exists() or cif_path.stat().st_size == 0:
            gz = d / f"{pdb.lower()}.cif.gz"
            download(CIF_URL.format(pdb=pdb), gz)
            gunzip(gz, cif_path)
            gz.unlink()

        contour, emdb_res, extractor_msg = emdb_params(n)
        issues = fasta_header_issues(fasta_path)
        details.append({
            "pdb_id": pdb, "emdb_id": r["emdb_id"],
            "resolution": r["resolution"], "emdb_resolution": emdb_res, "contour_level": contour,
            "residue_count": r["residue_count"], "unique_chain_count": r["unique_chain_count"],
            "size_bin": r["size_bin"], "release_date": r["release_date"],
            "map_file": str(map_path.relative_to(ROOT)), "map_size_mb": round(map_path.stat().st_size / 1e6, 1),
            "fasta_file": str(fasta_path.relative_to(ROOT)), "cif_file": str(cif_path.relative_to(ROOT)),
            "fasta_header_issues": " | ".join(issues),
            "title": r["title"],
        })
        print(f"    map {details[-1]['map_size_mb']} MB, contour={contour or 'MISSING'}, emdb_res={emdb_res or 'MISSING'}"
              + (f", FASTA issues: {issues}" if issues else ""), flush=True)
        if not contour or not emdb_res:
            print(f"    emdb_extractor said: {extractor_msg}", flush=True)

    with open(DATA / "benchmark_30_details.csv", "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=list(details[0].keys()))
        w.writeheader()
        w.writerows(details)
    print("done", flush=True)


if __name__ == "__main__":
    main()
