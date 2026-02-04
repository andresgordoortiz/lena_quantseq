#!/bin/bash

# Script to compare GRCz11 and GRCz12tu genome annotations
# Usage: ./compare_annotations.sh

set -e

# Define paths
GRCz11_GTF="/users/andres.ortiz/projects/lena_quantseq/reference_GRCz11/GCF_000002035.6_GRCz11_genomic.gtf"
GRCz12tu_GTF="/users/andres.ortiz/projects/lena_quantseq/reference_GRCz12tu/GCF_049306965.1_GRCz12tu_genomic.gtf"
OUTPUT_DIR="annotation_comparison_$(date +%Y%m%d_%H%M%S)"

# Create output directory
mkdir -p "$OUTPUT_DIR"
cd "$OUTPUT_DIR"

echo "=========================================="
echo "Genome Annotation Comparison"
echo "GRCz11 vs GRCz12tu"
echo "=========================================="
echo ""

# Extract gene information from GTF files
echo "Step 1: Extracting gene information from GTF files..."

# Extract genes from GRCz11 (gene_id and gene_name if available)
awk -F'\t' '$3 == "gene" {
    match($9, /gene_id "([^"]+)"/, gene_id)
    match($9, /gene "([^"]+)"/, gene_name)
    match($9, /gene_biotype "([^"]+)"/, biotype)

    gid = gene_id[1]
    gname = (gene_name[1] != "") ? gene_name[1] : "NA"
    gtype = (biotype[1] != "") ? biotype[1] : "NA"

    print gid "\t" gname "\t" gtype "\t" $1 "\t" $4 "\t" $5 "\t" $7
}' "$GRCz11_GTF" | sort -k1,1 > GRCz11_genes.txt

# Extract genes from GRCz12tu
awk -F'\t' '$3 == "gene" {
    match($9, /gene_id "([^"]+)"/, gene_id)
    match($9, /gene "([^"]+)"/, gene_name)
    match($9, /gene_biotype "([^"]+)"/, biotype)

    gid = gene_id[1]
    gname = (gene_name[1] != "") ? gene_name[1] : "NA"
    gtype = (biotype[1] != "") ? biotype[1] : "NA"

    print gid "\t" gname "\t" gtype "\t" $1 "\t" $4 "\t" $5 "\t" $7
}' "$GRCz12tu_GTF" | sort -k1,1 > GRCz12tu_genes.txt

echo "GRCz11 genes extracted: $(wc -l < GRCz11_genes.txt)"
echo "GRCz12tu genes extracted: $(wc -l < GRCz12tu_genes.txt)"
echo ""

# Step 2: Identify gene categories
echo "Step 2: Categorizing genes..."

# Genes only in GRCz11 (removed/deprecated in GRCz12tu)
awk 'NR==FNR {a[$1]; next} !($1 in a)' GRCz12tu_genes.txt GRCz11_genes.txt > genes_removed_in_GRCz12tu.txt

# Genes only in GRCz12tu (newly added)
awk 'NR==FNR {a[$1]; next} !($1 in a)' GRCz11_genes.txt GRCz12tu_genes.txt > genes_added_in_GRCz12tu.txt

# Genes in both (common genes)
awk 'NR==FNR {a[$1]=$0; next} ($1 in a) {print a[$1] "\t|\t" $0}' GRCz11_genes.txt GRCz12tu_genes.txt > genes_common.txt

echo "Genes removed in GRCz12tu: $(wc -l < genes_removed_in_GRCz12tu.txt)"
echo "Genes added in GRCz12tu: $(wc -l < genes_added_in_GRCz12tu.txt)"
echo "Genes in common: $(wc -l < genes_common.txt)"
echo ""

# Step 3: Analyze gene name changes (LOC/computational -> official gene names)
echo "Step 3: Analyzing gene nomenclature changes..."

# Genes that changed from LOC/computational to official names
awk -F'\t' '{
    split($0, parts, /\t\|\t/)
    if (length(parts) == 2) {
        split(parts[1], old, /\t/)
        split(parts[2], new, /\t/)

        old_name = old[2]
        new_name = new[2]

        # Check if old name was computational and new name is not
        if ((old_name ~ /^LOC/ || old_name ~ /^ENSDARG/ || old_name == "NA") &&
            new_name !~ /^LOC/ && new_name !~ /^ENSDARG/ && new_name != "NA") {
            print old[1] "\t" old_name "\t->\t" new_name "\t" new[4] ":" new[5] "-" new[6]
        }
    }
}' genes_common.txt > genes_with_new_names.txt

# Genes that went from having a name to LOC (unlikely but check)
awk -F'\t' '{
    split($0, parts, /\t\|\t/)
    if (length(parts) == 2) {
        split(parts[1], old, /\t/)
        split(parts[2], new, /\t/)

        old_name = old[2]
        new_name = new[2]

        # Check if old name was official and new name is computational
        if (old_name !~ /^LOC/ && old_name !~ /^ENSDARG/ && old_name != "NA" &&
            (new_name ~ /^LOC/ || new_name ~ /^ENSDARG/ || new_name == "NA")) {
            print old[1] "\t" old_name "\t->\t" new_name "\t" new[4] ":" new[5] "-" new[6]
        }
    }
}' genes_common.txt > genes_with_removed_names.txt

echo "Genes with new official names: $(wc -l < genes_with_new_names.txt)"
echo "Genes with names removed: $(wc -l < genes_with_removed_names.txt)"
echo ""

# Step 4: Analyze biotype changes
echo "Step 4: Analyzing gene biotype changes..."

awk -F'\t' '{
    split($0, parts, /\t\|\t/)
    if (length(parts) == 2) {
        split(parts[1], old, /\t/)
        split(parts[2], new, /\t/)

        if (old[3] != new[3]) {
            print old[1] "\t" old[2] "\t" old[3] "\t->\t" new[3]
        }
    }
}' genes_common.txt > genes_with_biotype_changes.txt

echo "Genes with biotype changes: $(wc -l < genes_with_biotype_changes.txt)"
echo ""

# Step 5: Analyze newly annotated genes (previously unannotated loci)
echo "Step 5: Identifying newly annotated loci..."

# New genes with official names (not LOC or ENSDARG)
awk -F'\t' '$2 !~ /^LOC/ && $2 !~ /^ENSDARG/ && $2 != "NA" {
    print $1 "\t" $2 "\t" $3 "\t" $4 ":" $5 "-" $6
}' genes_added_in_GRCz12tu.txt > newly_annotated_with_names.txt

# New genes with computational names
awk -F'\t' '($2 ~ /^LOC/ || $2 ~ /^ENSDARG/ || $2 == "NA") {
    print $1 "\t" $2 "\t" $3 "\t" $4 ":" $5 "-" $6
}' genes_added_in_GRCz12tu.txt > newly_annotated_computational.txt

echo "Newly annotated with official names: $(wc -l < newly_annotated_with_names.txt)"
echo "Newly annotated with computational IDs: $(wc -l < newly_annotated_computational.txt)"
echo ""

# Step 6: Create summary report
echo "Step 6: Generating summary report..."

cat > SUMMARY_REPORT.txt << EOF
========================================
Genome Annotation Comparison Summary
GRCz11 vs GRCz12tu
Generated: $(date)
========================================

OVERALL STATISTICS:
-------------------
Total genes in GRCz11:    $(wc -l < GRCz11_genes.txt)
Total genes in GRCz12tu:  $(wc -l < GRCz12tu_genes.txt)
Net change:               $(($(wc -l < GRCz12tu_genes.txt) - $(wc -l < GRCz11_genes.txt)))

GENE CHANGES:
-------------
Genes removed in GRCz12tu:           $(wc -l < genes_removed_in_GRCz12tu.txt)
Genes added in GRCz12tu:             $(wc -l < genes_added_in_GRCz12tu.txt)
Genes present in both versions:     $(wc -l < genes_common.txt)

ANNOTATION IMPROVEMENTS:
------------------------
Genes with new official names:       $(wc -l < genes_with_new_names.txt)
  (previously LOC/computational IDs)

Genes with names removed:            $(wc -l < genes_with_removed_names.txt)
  (reverted to computational IDs)

Genes with biotype changes:          $(wc -l < genes_with_biotype_changes.txt)

NEWLY ANNOTATED LOCI:
---------------------
New genes with official names:       $(wc -l < newly_annotated_with_names.txt)
New genes with computational IDs:    $(wc -l < newly_annotated_computational.txt)

========================================

OUTPUT FILES:
-------------
1. GRCz11_genes.txt                  - All genes in GRCz11
2. GRCz12tu_genes.txt                - All genes in GRCz12tu
3. genes_removed_in_GRCz12tu.txt     - Genes not in GRCz12tu
4. genes_added_in_GRCz12tu.txt       - New genes in GRCz12tu
5. genes_common.txt                  - Genes in both versions
6. genes_with_new_names.txt          - LOC/computational -> official names
7. genes_with_removed_names.txt      - Official names -> computational IDs
8. genes_with_biotype_changes.txt    - Genes with changed biotypes
9. newly_annotated_with_names.txt    - New genes with official names
10. newly_annotated_computational.txt - New genes with computational IDs

========================================
EOF

cat SUMMARY_REPORT.txt

echo ""
echo "Analysis complete!"
echo "Results saved in: $OUTPUT_DIR"
echo ""

# Step 7: Create biotype distribution comparison
echo "Step 7: Comparing biotype distributions..."

echo "GRCz11 Biotype Distribution:" > biotype_distribution.txt
awk -F'\t' '{print $3}' GRCz11_genes.txt | sort | uniq -c | sort -rn >> biotype_distribution.txt

echo "" >> biotype_distribution.txt
echo "GRCz12tu Biotype Distribution:" >> biotype_distribution.txt
awk -F'\t' '{print $3}' GRCz12tu_genes.txt | sort | uniq -c | sort -rn >> biotype_distribution.txt

echo "Biotype distribution saved to: biotype_distribution.txt"
echo ""
echo "=========================================="
echo "All done! Check the output directory for detailed results."
echo "=========================================="