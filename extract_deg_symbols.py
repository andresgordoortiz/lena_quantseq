import csv, glob

degs = set()

# Full DESeq2 outputs: filter padj < 0.05 and abs(LFC) >= 1.0
full_deseq = (
    glob.glob("results/q1_exp1_0vs15_activin_240min.csv")
    + glob.glob("results/q1_exp2_0vs15_dmso.csv")
    + glob.glob("results/q3_Activin_vs_Baseline.csv")
    + glob.glob("results/q3_SB50_*_vs_*.csv")
)

for f in full_deseq:
    with open(f) as fh:
        reader = csv.DictReader(fh)
        for row in reader:
            try:
                padj = float(row.get("padj", "1"))
                lfc = float(row.get("log2FoldChange", "0"))
                gene = row.get("gene", "").strip('"')
                if padj < 0.05 and abs(lfc) >= 1.0 and gene and gene != "NA":
                    degs.add(gene)
            except (ValueError, TypeError):
                pass

# Already-filtered gene lists: take gene column directly
filtered_files = (
    glob.glob("results/q1_shared_de_genes.csv")
    + glob.glob("results/q3_genelist_*.csv")
    + glob.glob("results/q3_always_not_blocked_genes.csv")
    + glob.glob("results/q2_nodal_divergence_*.csv")
    + glob.glob("results/q5_motility_gene_tracking.csv")
    + glob.glob("results/q5_motility_delayed_commitment_genes.csv")
    + glob.glob("results/q4_chipseq_candidates.csv")
    + glob.glob("results/q4_candidate_validation.csv")
    + glob.glob("results/q4_reversibility_classified.csv")
)
# NOTE: q3_gene_transfer_tracking.csv excluded — it contains ALL genes, not just DEGs

for f in filtered_files:
    with open(f) as fh:
        reader = csv.DictReader(fh)
        fields = reader.fieldnames or []
        for col in ["gene", "gene_id", "gene_symbol"]:
            if col in fields:
                fh.seek(0)
                next(fh)
                for row in csv.DictReader(fh):
                    g = row.get(col, "").strip('"')
                    if g and g != "NA":
                        degs.add(g)
                break

with open("compare_annotations/deg_symbols.txt", "w") as out:
    for g in sorted(degs):
        out.write(g + "\n")

print(f"Unique DEG symbols: {len(degs)}")
