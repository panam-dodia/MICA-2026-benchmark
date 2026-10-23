#!/bin/bash
# Score finished models and save the tools' raw output for later parsing (scripts/parse_scores.py).
#
#   US-align  <model> <deposited> -mm 1 -ter 0          -> <EMDB>.usalign.txt       TM-score, aligned length, Seq_ID
#   phenix.chain_comparison <model> <ref_inplace>       -> <EMDB>.cc_precision.txt  target = model: % of model residues
#                                                          with a reference CA within 3 A (the papers' "Calpha match")
#                                                          and sequence match
#   phenix.chain_comparison <ref_inplace> <model>       -> <EMDB>.cc_recall.txt     target = reference: % of reference
#                                                          residues matched (coverage)
# ref_inplace = biological assembly 1 (<pdbid>_assembly1.cif) when it exists, i.e. for 9I8C, 9SZB and 12ZE where the
# deposited model is one asymmetric unit of a symmetric map; otherwise the deposited model. If a <ref>_phenix.cif
# copy exists (ligands with 5-character names removed, which Phenix 1.20 cannot read), that copy is used.
# Rules: docs/PROTOCOL.md.
#
# Raw outputs: $BENCH/eval/raw/<method>/. Existing outputs are kept (resumable).
# Usage: bash score_raw.sh <method> [EMDB ...]     method: mica | modelangelo | cryoatom | ...
set -u
METHOD=${1:?usage: score_raw.sh <method> [EMDB ...]}; shift

source "$(dirname "${BASH_SOURCE[0]}")/config.sh"
RAW=$BENCH/eval/raw/$METHOD
mkdir -p "$RAW"

model_path() {   # final model chosen per method (see docs/PROTOCOL.md)
    local n=$1
    case "$METHOD" in
        mica)        ls "$RUNS_ROOT/mica/$n/output/"*_real_space_refined_000.pdb 2>/dev/null | head -1 ;;
        modelangelo) ls "$RUNS_ROOT/modelangelo/$n/output/output.cif" 2>/dev/null ;;
        cryoatom)    ls "$RUNS_ROOT/cryoatom/$n/output/output.cif" 2>/dev/null ;;
        emprot)      ls "$RUNS_ROOT/emprot/$n/output/output.cif" 2>/dev/null ;;
        emodelx_af)  ls "$RUNS_ROOT/emodelx/$n/output/"*real_space_refined_000.* 2>/dev/null | head -1 ;;
        *)           echo "" ;;
    esac
}

# phenix_env.sh references unset variables, so relax "set -u" while sourcing it
set +u; source "$PHENIX_ACT" > /dev/null 2>&1; set -u
command -v phenix.chain_comparison > /dev/null || { echo "phenix.chain_comparison not found after sourcing $PHENIX_ACT"; exit 1; }

header() { echo "# model: $MODEL"; echo "# reference: $1"; echo "# cmd: $2"; echo "# date: $(date '+%F %T')"; }

if [ $# -gt 0 ]; then ENTRIES="$*"; else ENTRIES=$(tail -n +2 "$BENCH/runs/ready_entries.tsv" | cut -f1); fi
for N in $ENTRIES; do
    PDB=$(awk -F'\t' -v n="$N" '$1 == n {print tolower($2)}' "$BENCH/runs/ready_entries.tsv")
    [ -n "$PDB" ] || { echo "$N: not in ready_entries.tsv"; continue; }
    E=$BENCH/dataset/entries/$N
    REF_DEP=$E/$PDB.cif
    REF_INPLACE=$REF_DEP
    [ -s "$E/${PDB}_assembly1.cif" ] && REF_INPLACE=$E/${PDB}_assembly1.cif
    # Phenix-1.20-compatible copy (5-character ligand names removed; scripts/12_phenix_compatible_refs.py)
    [ -s "${REF_INPLACE%.cif}_phenix.cif" ] && REF_INPLACE=${REF_INPLACE%.cif}_phenix.cif
    MODEL=$(model_path "$N")
    if [ -z "$MODEL" ] || [ ! -s "$MODEL" ]; then echo "$N: no $METHOD model"; continue; fi
    echo "$N ($PDB): $(basename "$MODEL")  [in-place ref: $(basename "$REF_INPLACE")]"

    # US-align cannot parse the mmCIF written by phenix.real_space_refine (EModelX(+AF): "Chain number 0").
    # For US-align only, such models are rewritten as PDB by Phenix itself (iotbx); coordinates are unchanged.
    US_MODEL=$MODEL
    if [ "$METHOD" = "emodelx_af" ] && [ "${MODEL##*.}" = "cif" ]; then
        US_MODEL=$RAW/$N.model_for_usalign.pdb
        [ -s "$US_MODEL" ] || phenix.python -c "import iotbx.pdb,sys; iotbx.pdb.input(file_name=sys.argv[1]).construct_hierarchy().write_pdb_file(sys.argv[2])" "$MODEL" "$US_MODEL"
    fi
    if [ ! -s "$RAW/$N.usalign.txt" ]; then
        { header "$REF_DEP" "USalign model deposited -mm 1 -ter 0 (model file: $(basename "$US_MODEL"))"; "$USALIGN" "$US_MODEL" "$REF_DEP" -mm 1 -ter 0; } > "$RAW/$N.usalign.txt" 2>&1
    fi
    if [ ! -s "$RAW/$N.cc_precision.txt" ]; then
        { header "$REF_INPLACE" "phenix.chain_comparison model reference (target=model)"; phenix.chain_comparison "$MODEL" "$REF_INPLACE"; } > "$RAW/$N.cc_precision.txt" 2>&1
    fi
    if [ ! -s "$RAW/$N.cc_recall.txt" ]; then
        { header "$REF_INPLACE" "phenix.chain_comparison reference model (target=reference)"; phenix.chain_comparison "$REF_INPLACE" "$MODEL"; } > "$RAW/$N.cc_recall.txt" 2>&1
    fi
    echo "    $(grep -E '^TM-score=.*Structure_2' "$RAW/$N.usalign.txt" | head -1)"
    echo "    precision: $(grep -E '^ *(Unique|Entire)_target' "$RAW/$N.cc_precision.txt" | tail -1)"
    echo "    recall:    $(grep -E '^ *(Unique|Entire)_target' "$RAW/$N.cc_recall.txt" | tail -1)"
done
