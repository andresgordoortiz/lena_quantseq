#!/bin/bash
# =============================================================================
# check_degs_annotation.sh
#
# Validate that none of the DEGs used in downstream analyses were affected
# by the GRCz11 → GRCz12tu annotation update.
#
# WHAT YOU NEED ON THE CLUSTER (lena_quantseq/):
#   check_degs_annotation.sh          ← this script
#   deg_symbols.txt                   ← gene list (generated locally, copy it over)
#   reference_GRCz11/*.gtf            ← already there
#   reference_GRCz12tu/*.gtf          ← already there
#
# USAGE:
#   scp check_degs_annotation.sh deg_symbols.txt clip:~/projects/lena_quantseq/
#   ssh clip
#   cd ~/projects/lena_quantseq
#   chmod +x check_degs_annotation.sh
#   ./check_degs_annotation.sh
#
# OUTPUT (in compare_annotations/):
#   deg_annotation_check.tsv   — per-gene TSV (copy back & import into R)
#   deg_annotation_summary.txt — human-readable summary
#
# Requirements: gawk (GNU awk with match() 3-arg form), standard coreutils.
# =============================================================================
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"

# ── Paths (relative to script location) ────────────────────────────────────
GRCz11_GTF="$SCRIPT_DIR/reference_GRCz11/GCF_000002035.6_GRCz11_genomic.gtf"
GRCz12tu_GTF="$SCRIPT_DIR/reference_GRCz12tu/GCF_049306965.1_GRCz12tu_genomic.gtf"
DEG_FILE="$SCRIPT_DIR/deg_symbols.txt"

OUTDIR="$SCRIPT_DIR/compare_annotations"
mkdir -p "$OUTDIR"

# ── Sanity checks ──────────────────────────────────────────────────────────
for f in "$GRCz11_GTF" "$GRCz12tu_GTF" "$DEG_FILE"; do
  if [ ! -f "$f" ]; then
    echo "ERROR: required file not found: $f" >&2
    exit 1
  fi
done

N_DEGS=$(wc -l < "$DEG_FILE" | tr -d ' ')
echo "Gene symbols to check: $N_DEGS  (from deg_symbols.txt)"
echo ""

# ── 1. Parse GTFs → gene_symbol <TAB> biotype ──────────────────────────────
echo "Step 1/4: Parsing GTF files..."

parse_gtf() {
  awk -F'\t' '$3 == "gene" {
    match($9, /gene_id "([^"]+)"/, gid)
    match($9, /gene "([^"]+)"/, gname)
    match($9, /gene_biotype "([^"]+)"/, bt)
    sym = (gname[1] != "") ? gname[1] : gid[1]
    bio = (bt[1] != "") ? bt[1] : "NA"
    print sym "\t" bio
  }' "$1" | sort -u -k1,1
}

GRCz11_GENES="$OUTDIR/_grch11_genes.tmp"
GRCz12_GENES="$OUTDIR/_grch12tu_genes.tmp"

parse_gtf "$GRCz11_GTF"  > "$GRCz11_GENES"
parse_gtf "$GRCz12tu_GTF" > "$GRCz12_GENES"

echo "   GRCz11 unique symbols:  $(wc -l < "$GRCz11_GENES" | tr -d ' ')"
echo "   GRCz12tu unique symbols: $(wc -l < "$GRCz12_GENES" | tr -d ' ')"

# ── 2. Build change-category files ─────────────────────────────────────────
echo "Step 2/4: Identifying annotation changes..."

REMOVED="$OUTDIR/_removed.tmp"
awk 'NR==FNR {a[$1]; next} !($1 in a) {print $1}' \
  "$GRCz12_GENES" "$GRCz11_GENES" | sort -u > "$REMOVED"

ADDED="$OUTDIR/_added.tmp"
awk 'NR==FNR {a[$1]; next} !($1 in a) {print $1}' \
  "$GRCz11_GENES" "$GRCz12_GENES" | sort -u > "$ADDED"

BIOTYPE_CHG="$OUTDIR/_biotype_changed.tmp"
awk -F'\t' 'NR==FNR {bt[$1]=$2; next}
  ($1 in bt) && bt[$1] != $2 {print $1 "\t" bt[$1] "\t" $2}' \
  "$GRCz11_GENES" "$GRCz12_GENES" | sort -u -k1,1 > "$BIOTYPE_CHG"

echo "   Removed in GRCz12tu:  $(wc -l < "$REMOVED" | tr -d ' ')"
echo "   Added in GRCz12tu:    $(wc -l < "$ADDED" | tr -d ' ')"
echo "   Biotype changed:      $(wc -l < "$BIOTYPE_CHG" | tr -d ' ')"

# ── 3. Cross-reference DEGs ────────────────────────────────────────────────
echo "Step 3/4: Cross-referencing gene symbols against annotation changes..."

OUTPUT="$OUTDIR/deg_annotation_check.tsv"

printf "gene\tin_GRCz11\tin_GRCz12tu\tremoved\tadded\tbiotype_changed\told_biotype\tnew_biotype\n" > "$OUTPUT"

awk -F'\t' '
  FILENAME == ARGV[1] { grch11[$1] = $2; next }
  FILENAME == ARGV[2] { grch12[$1] = $2; next }
  FILENAME == ARGV[3] { removed[$1] = 1; next }
  FILENAME == ARGV[4] { added[$1]   = 1; next }
  FILENAME == ARGV[5] { bt_old[$1]  = $2; bt_new[$1] = $3; next }
  {
    g = $1
    printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n", \
      g, \
      (g in grch11)  ? "TRUE" : "FALSE", \
      (g in grch12)  ? "TRUE" : "FALSE", \
      (g in removed) ? "TRUE" : "FALSE", \
      (g in added)   ? "TRUE" : "FALSE", \
      (g in bt_old)  ? "TRUE" : "FALSE", \
      (g in bt_old)  ? bt_old[g] : "NA", \
      (g in bt_new)  ? bt_new[g] : "NA"
  }
' "$GRCz11_GENES" "$GRCz12_GENES" "$REMOVED" "$ADDED" "$BIOTYPE_CHG" "$DEG_FILE" >> "$OUTPUT"

# ── 4. Summary ─────────────────────────────────────────────────────────────
echo "Step 4/4: Summary..."
echo ""

N_REM=$(awk -F'\t' 'NR>1 && $4=="TRUE"' "$OUTPUT" | wc -l | tr -d ' ')
N_ADD=$(awk -F'\t' 'NR>1 && $5=="TRUE"' "$OUTPUT" | wc -l | tr -d ' ')
N_BTC=$(awk -F'\t' 'NR>1 && $6=="TRUE"' "$OUTPUT" | wc -l | tr -d ' ')
N_MISS=$(awk -F'\t' 'NR>1 && $3=="FALSE"' "$OUTPUT" | wc -l | tr -d ' ')

SUMMARY="$OUTDIR/deg_annotation_summary.txt"
cat > "$SUMMARY" << EOF
=====================================================
DEG Annotation Validation
GRCz11 → GRCz12tu
$(date)
=====================================================

Total gene symbols checked:           $N_DEGS

DEGs REMOVED in GRCz12tu:             $N_REM
DEGs newly ADDED in GRCz12tu:         $N_ADD
DEGs with BIOTYPE change:             $N_BTC
DEGs not found in GRCz12tu at all:    $N_MISS

EOF

if [ "$N_REM" -gt 0 ]; then
  echo "--- DEGs removed in GRCz12tu ---" >> "$SUMMARY"
  awk -F'\t' 'NR>1 && $4=="TRUE" {print "  " $1}' "$OUTPUT" >> "$SUMMARY"
  echo "" >> "$SUMMARY"
fi
if [ "$N_ADD" -gt 0 ]; then
  echo "--- DEGs newly added in GRCz12tu ---" >> "$SUMMARY"
  awk -F'\t' 'NR>1 && $5=="TRUE" {print "  " $1}' "$OUTPUT" >> "$SUMMARY"
  echo "" >> "$SUMMARY"
fi
if [ "$N_BTC" -gt 0 ]; then
  echo "--- DEGs with biotype changes ---" >> "$SUMMARY"
  awk -F'\t' 'NR>1 && $6=="TRUE" {printf "  %s: %s -> %s\n", $1, $7, $8}' "$OUTPUT" >> "$SUMMARY"
  echo "" >> "$SUMMARY"
fi
if [ "$N_MISS" -gt 0 ]; then
  echo "--- DEGs not found in GRCz12tu (potential concern) ---" >> "$SUMMARY"
  awk -F'\t' 'NR>1 && $3=="FALSE" {print "  " $1}' "$OUTPUT" >> "$SUMMARY"
  echo "" >> "$SUMMARY"
fi

if [ "$N_REM" -eq 0 ] && [ "$N_ADD" -eq 0 ] && [ "$N_BTC" -eq 0 ] && [ "$N_MISS" -eq 0 ]; then
  echo "CONCLUSION: None of the $N_DEGS gene symbols were affected by" >> "$SUMMARY"
  echo "  the GRCz11 -> GRCz12tu annotation update." >> "$SUMMARY"
else
  echo "WARNING: Some DEGs were affected. See details above." >> "$SUMMARY"
fi

cat "$SUMMARY"

# Clean up temp files
rm -f "$OUTDIR"/_*.tmp

echo ""
echo "Output files:"
echo "  $OUTPUT"
echo "  $SUMMARY"
echo ""
echo "To import into R (copy deg_annotation_check.tsv back to your local machine):"
echo '  deg_check <- read.delim("compare_annotations/deg_annotation_check.tsv")'
echo '  table(deg_check$removed)'
echo '  subset(deg_check, removed | added | biotype_changed)'
