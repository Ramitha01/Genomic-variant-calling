# ----------------------------
# Read file
# ----------------------------
anno <- read.csv("VIZ_0024_SP_B01.hg38_multianno.csv", stringsAsFactors = FALSE)

# ----------------------------
# 1. Filter only exonic regions
# ----------------------------
# Correction: Use == "exonic" (your condition is correct)
exonic_variants <- subset(anno, Func.refGene == "exonic")

# Save to CSV
write.csv(exonic_variants, "exonic_variants.csv", row.names = FALSE)
cat("Exonic variants saved:", nrow(exonic_variants), "\n")


# ----------------------------
# 2. Filter nonsynonymous SNVs
# ----------------------------
library(dplyr)

variants <- read.csv("exonic_variants.csv", stringsAsFactors = FALSE)

# Important correction:
# "nonsynonymous SNV" is correct. No change.
nonsyn_variants <- variants %>%
  filter(ExonicFunc.refGene == "nonsynonymous SNV")

write.csv(nonsyn_variants, "nonsynonymous_variants.csv", row.names = FALSE)
cat("Nonsynonymous SNVs saved:", nrow(nonsyn_variants), "\n")


# ----------------------------
# 3. SIFT Filtering (three-step)
# ----------------------------

library(dplyr)
variants <- read.csv("nonsynonymous_variants.csv", stringsAsFactors = FALSE)

# CORRECTION: convert "." to NA before numeric conversion
variants$SIFT_score <- suppressWarnings(as.numeric(ifelse(variants$SIFT_score == ".", NA, variants$SIFT_score)))

# Remove variants with no SIFT score
variants_with_sift <- variants %>% filter(!is.na(SIFT_score))

# Apply SIFT threshold
sift_filtered <- variants_with_sift %>% filter(SIFT_score <= 0.05)

write.csv(sift_filtered, "sift_filtered_variants_nonsyno.csv", row.names = FALSE)
cat("Variants with valid SIFT score ≤ 0.05 saved:", nrow(sift_filtered), "\n")


# ----------------------------
# 4. GO / KEGG / Reactome Enrichment Analysis
# ----------------------------

library(clusterProfiler)
library(org.Hs.eg.db)
library(ReactomePA)
library(openxlsx)
library(GO.db)

# Load final filtered genes
gene_data <- read.csv("sift_filtered_variants_nonsyno.csv", stringsAsFactors = FALSE)

# CORRECTION:
# Your column name MUST exist. You used Gene.refGene everywhere correctly.
gene_symbols <- unique(gene_data$Gene.refGene)

# Convert gene symbols to ENTREZ IDs
gene_entrez <- bitr(
  gene_symbols,
  fromType = "SYMBOL",
  toType = "ENTREZID",
  OrgDb = org.Hs.eg.db
)

entrez_ids <- unique(gene_entrez$ENTREZID)

# Output folder
output_dir <- "enrichment_results_exonicvariantsnew"
if (!dir.exists(output_dir)) dir.create(output_dir)

# ----------------------------
# Helper function for saving results
# ----------------------------
save_results <- function(enrich_result, analysis_type) {
  df <- as.data.frame(enrich_result)
  
  if (!is.null(df) && nrow(df) > 0) {
    write.xlsx(df, file.path(output_dir, paste0(analysis_type, "_results.xlsx")))
    
    png(file.path(output_dir, paste0(analysis_type, "_dotplot.png")), 
        width = 1200, height = 800, res = 150)
    
    print(dotplot(enrich_result, showCategory = 20))
    dev.off()
    
    cat(analysis_type, "pathways found:", nrow(df), "\n")
  } else {
    cat("No significant", analysis_type, "pathways found.\n")
  }
}

# ----------------------------
# 5. GO enrichment
# ----------------------------
go_result <- enrichGO(
  gene          = entrez_ids,
  OrgDb         = org.Hs.eg.db,
  keyType       = "ENTREZID",
  ont           = "ALL",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.05,
  readable      = TRUE
)
save_results(go_result, "GO")

# ----------------------------
# 6. KEGG enrichment
# ----------------------------
kegg_result <- enrichKEGG(
  gene          = entrez_ids,
  organism      = "hsa",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05
)

# make readable
kegg_result <- setReadable(kegg_result, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")
save_results(kegg_result, "KEGG")

# ----------------------------
# 7. Reactome enrichment
# ----------------------------
reactome_result <- enrichPathway(
  gene          = entrez_ids,
  organism      = "human",
  pvalueCutoff  = 0.05,
  pAdjustMethod = "BH",
  readable      = TRUE
)
save_results(reactome_result, "Reactome")

cat("Analysis complete. Results saved in:", output_dir, "\n")


# ----------------------------
# 8. Extract SNP, Gene, SIFT Score
# ----------------------------

library(data.table)

dt <- fread("sift_filtered_variants_nonsyno.csv",
            select = c("avsnp151", "Gene.refGene", "SIFT_score"))

setnames(dt, 
         old = c("avsnp151", "Gene.refGene", "SIFT_score"),
         new = c("SNP", "Gene", "Sift_Score"))

# CORRECTION:
# Treat "." as NA (for numeric consistency)
dt$Sift_Score <- suppressWarnings(as.numeric(ifelse(dt$Sift_Score == ".", NA, dt$Sift_Score)))

# Filter correctly
filtered <- dt[
  SNP != "" & SNP != "." & !is.na(Sift_Score)
]

fwrite(filtered, "snp_gene.csv")
head(filtered)
