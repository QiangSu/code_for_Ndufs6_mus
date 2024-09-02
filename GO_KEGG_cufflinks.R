library(clusterProfiler)
library(org.Mm.eg.db)
library(DOSE)
library(ggplot2)  # Load ggplot2 for theme function

# Step 1: Read the CSV file
file_path <- "C:/Users/Qiang/Desktop/SIAT/data/DEseq2_musKO_AAV_cufflinks_isoform_result_KO_downregulation.csv"
degs <- read.csv(file_path, header = TRUE, sep = ",")

# Verify the data is read correctly
print(head(degs))

# Step 2: Check and print column names
print(colnames(degs))

# Step 3: Identify the correct columns for gene_id and log2FoldChange
# Adjust the column names here based on the actual data
gene_id_col <- "tracking_id"  # Ensure this column exists in your data
log2fc_col <- "log2FoldChange"  # Ensure this column exists in your data

if(gene_id_col %in% colnames(degs) & log2fc_col %in% colnames(degs)) {
  gene_list <- setNames(degs[[log2fc_col]], degs[[gene_id_col]])
  # Verify the gene_list is created correctly
  print(head(gene_list))
} else {
  stop("The necessary columns (Transcript_id, log2FoldChange) are not present in the data frame.")
}

# Step 4: Convert gene identifiers to ENTREZID
gene_df <- bitr(names(gene_list), fromType = "ENSEMBLTRANS", toType = "ENTREZID", OrgDb = org.Mm.eg.db)
# Ensure that the conversion was successful
if (nrow(gene_df) == 0) {
  stop("Gene ID conversion failed. Please check your gene identifiers.")
}
gene_list <- gene_list[gene_df$ENSEMBLTRANS]
names(gene_list) <- gene_df$ENTREZID

# Step 5: Perform GO Enrichment Analysis
ego <- enrichGO(gene = names(gene_list),
                OrgDb = org.Mm.eg.db,
                keyType = "ENTREZID",  # Use ENTREZID after conversion
                pAdjustMethod = "BH",
                pvalueCutoff = 0.1,  # Less stringent p-value cutoff
                qvalueCutoff = 0.2,
                ont = "ALL")

# Check if any GO terms were enriched
if (is.null(ego) || length(ego) == 0) {
  cat("No significantly enriched GO terms found. Try less stringent thresholds or check your gene identifiers.\n")
} else {
  # Step 6: View Results
  print(summary(ego))
  
  # Adjust plot size and rotate labels to avoid overlapping
  barplot(ego, showCategory=20, title="GO Enrichment", font.size=8) + 
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  dotplot(ego, showCategory=20, title="GO Enrichment", font.size=8) + 
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  # Save the results
  #output_file <- "C:/Users/Qiang/Desktop/SIAT/data/GO_enrichment_results_DEseq2_musADULT_NB_cufflinks_isoform_result.csv"
  #write.csv(as.data.frame(ego), file = output_file)
}

# Step 7: Perform KEGG Enrichment Analysis
ekegg <- enrichKEGG(gene = names(gene_list),
                    organism = 'mmu',  # 'mmu' is the code for mouse
                    keyType = "kegg",  # Use ENTREZID for KEGG
                    pAdjustMethod = "BH",
                    pvalueCutoff = 0.1,  # Less stringent p-value cutoff
                    qvalueCutoff = 0.2)

# Check if any KEGG pathways were enriched
if (is.null(ekegg) || length(ekegg) == 0) {
  cat("No significantly enriched KEGG pathways found. Try less stringent thresholds or check your gene identifiers.\n")
} else {
  # Step 8: View Results
  print(summary(ekegg))
  
  # Adjust plot size and rotate labels to avoid overlapping
  barplot(ekegg, showCategory=20, title="KEGG Enrichment", font.size=6) + 
    theme(axis.text.x = element_text(angle = 0, hjust = 1))
  
  dotplot(ekegg, showCategory=20, title="KEGG Enrichment", font.size=6) + 
    theme(axis.text.x = element_text(angle = 0, hjust = 1))
  
  # Step 9: Save the Results
  #output_file_kegg <- "C:/Users/Qiang/Desktop/SIAT/data/KEGG_enrichment_results_DEseq2_musADULT_NB_cufflinks_isoform_result.csv"
  #write.csv(as.data.frame(ekegg), file = output_file_kegg)
}
