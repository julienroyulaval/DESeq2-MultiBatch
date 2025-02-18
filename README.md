# batch_correction_DESeq2
Batch Correction and Adaptive Contrast Calling in Complex RNA-seq Experimental Designs Using DESeq2

Formatting to be learned

Github content

Summary of data + design experimental as a figure

Setup with DESeq2

PCA plot before correction, explanations

Table showing DEG results with and without complex design (genotypes included for treatment) including concept shown in figure :

![image](https://github.com/user-attachments/assets/2bc47a79-397d-4e4e-8802-3de003b124a6)

Comparison with concatenation of variables/factors (a_ff_0 or a_xx_0_ctrl)

Tool for contrast following (LOCtest) if needed :

![image](https://github.com/user-attachments/assets/7b6751d0-1eeb-4dbd-844f-d2c3fd74cca2)

Contrast calling theory

including imbalance

base level genotype default result

averaging treatment effect on genotypes

Batch correction on normalised counts

"For example, if the normalized count of a gene is around 1000 for XX plants in batch A and around 2000 for XX plants in batch B, the log2 fold change (LFC) of the "batch" factor is 1. Instead of shifting both batches toward the arithmetic mean (1500), which would require different scaling factors (1.5 for batch A and 1/1.333 for batch B), a symmetric transformation is applied in log space.. Specifically, all normalized counts for XX samples from batch A are multiplied by √(batch B/batch A) = √(2) ≈ 1.414, while those from batch B are divided by the same factor (or equivalently multiplied by 0.707). This proportional adjustment in log space brings both batches closer together without distorting the relative differences between conditions."

	Quantifying basic batch effect and interaction batch effect
	Checking for consistency of batch effect across days
	Comparing basic correction with limma
	Manual data correction

PCA plot comparison

Showing bar graph with gene counts made easier

Conclusion


# Load count data
counts <- read.csv("counts.csv", sep = ",", row.names = 1, check.names = FALSE)

# Load metadata
colData <- read.csv("metadata.csv")

# Convert the condition columns to factors
colData$Treatment <- as.factor(colData$Treatment)
colData$Cultivar <- as.factor(colData$Cultivar)
colData$Day <- as.factor(colData$Day)
colData$Sex <- as.factor(colData$Sex)
colData$Batch <- as.factor(colData$Batch)
colData$Sample_ID <- as.factor(colData$Sample_ID)

# DESIGN
# Construct the design
# The interaction between treatment and cultivar was removed for simplicity in batch correction
design <- ~ Treatment + Cultivar + Sex + Day + Sex:Day + Sex:Batch + Batch

# Convert counts to matrix
count_matrix <- as.matrix(counts)

# Construct DESeqDataSet object
dds <- DESeqDataSetFromMatrix(countData = count_matrix, colData = colData, design = design)

# Set reference levels
dds$Treatment <- relevel(dds$Treatment, ref = "ctrl")
dds$Cultivar <- relevel(dds$Cultivar, ref = "DK")
dds$Sex <- relevel(dds$Sex, ref = "XX")
dds$Batch <- relevel(dds$Batch, ref = "a")
dds$Day <- relevel(dds$Day, ref = "0")

# Run DESeq
dds <- DESeq(dds)

# Filter low-count genes with HTSFilter
dds_filtered <- HTSFilter(dds, s.min = 15, s.max = 50, s.len=25, plot = TRUE)$filteredData

# Calculate batch effect and interaction term batch effect
res_batch_XX <- results(dds_filtered, name = "Batch_b_vs_a", independentFiltering = FALSE)
res_batch_XY <- results(dds_filtered, contrast=list(c("Batch_b_vs_a","SexXY.Batchb")), independentFiltering = FALSE)

# Convert to dataframes
df_batch_XX <- as.data.frame(res_batch_XX)
df_batch_XY <- as.data.frame(res_batch_XY)

# Compute scaling factors
df_batch_XX$Scaling_Factor_Batch_A <- sqrt(2^(df_batch_XX$log2FoldChange))
df_batch_XX$Scaling_Factor_Batch_B <- 1 / df_batch_XX$Scaling_Factor_Batch_A
df_batch_XY$Scaling_Factor_Batch_A <- sqrt(2^(df_batch_XY$log2FoldChange))
df_batch_XY$Scaling_Factor_Batch_B <- 1 / df_batch_XY$Scaling_Factor_Batch_A

# Extract normalized counts
normalized_counts <- counts(dds_filtered, normalized = TRUE)
normalized_counts_df <- as.data.frame(normalized_counts)
scaled_counts <- normalized_counts_df

# Apply batch correction
for (Sample_ID in colnames(normalized_counts_df)) {
  Sex <- colData[Sample_ID, "Sex"]
  Batch <- colData[Sample_ID, "Batch"]
  
  if (Sex == "XX" && Batch == "a") {
    scaling_factors <- df_batch_XX$Scaling_Factor_Batch_A
  } else if (Sex == "XX" && Batch == "b") {
    scaling_factors <- df_batch_XX$Scaling_Factor_Batch_B
  } else if (Sex == "XY" && Batch == "a") {
    scaling_factors <- df_batch_XY$Scaling_Factor_Batch_A
  } else if (Sex == "XY" && Batch == "b") {
    scaling_factors <- df_batch_XY$Scaling_Factor_Batch_B
  }
  
  scaled_counts[, Sample_ID] <- normalized_counts_df[, Sample_ID] * scaling_factors
}

# Re-run DESeq2 with corrected counts
scaled_counts_rounded <- round(scaled_counts)
dds_scaled <- DESeqDataSetFromMatrix(countData = scaled_counts_rounded, colData = colData, design = design)
dds_scaled <- DESeq(dds_scaled)

# Check batch effect after correction
res_batch_XX_corrected <- results(dds_scaled, name = "Batch_b_vs_a", independentFiltering = FALSE)
res_batch_XY_corrected <- results(dds_scaled, contrast=list(c("Batch_b_vs_a","SexXY.Batchb")), independentFiltering = FALSE)

# Convert to dataframes
df_batch_XX_corrected <- as.data.frame(res_batch_XX_corrected)
df_batch_XY_corrected <- as.data.frame(res_batch_XY_corrected)

# PCA plot function
create_pca <- function(dds, plot_title) {
  vsd <- vst(dds, blind = TRUE)
  colData_vsd <- as.data.frame(colData(vsd))
  colData_vsd$Group <- paste(colData_vsd$Treatment_Name, colData_vsd$Sex, sep = "_")
  colData(vsd) <- DataFrame(colData_vsd)
  pcaData <- plotPCA(vsd, intgroup = c("Group", "Cultivar"), returnData = TRUE)
  
  percentVar <- round(100 * attr(pcaData, "percentVar"))
  custom_colors <- c("FF_XX" = "#33cc33", "IFF_XY" = "#ff6600", "MF_XY" = "#993399", "IMF_XX" = "#0066ff")
  
  ggplot(pcaData, aes(x = PC1, y = PC2, color = Group, group = Group)) +
    geom_point(aes(shape = Cultivar), size = 2) +
    stat_ellipse(type = "t", level = 0.95, linetype = 2) +
    scale_color_manual(values = custom_colors) +
    xlab(paste0("PC1: ", percentVar[1], "% variance")) +
    ylab(paste0("PC2: ", percentVar[2], "% variance")) +
    coord_fixed() +
    ggtitle(plot_title) +
    theme(plot.title = element_text(hjust = 0.3))
}

# Subset for PCA analysis
samples_to_keep_day0 <- colData(dds_scaled)$Day == "0"
dds_scaled_0 <- dds_scaled[, samples_to_keep_day0]
dds_filtered_0 <- dds_filtered[, samples_to_keep_day0]

# Generate PCA plots
pca_plot_scaled <- create_pca(dds_scaled_0, "PCA of VSD for all genes, day 0, batch corrected")
pca_plot_filtered <- create_pca(dds_filtered_0, "PCA of VSD for all genes, day 0, before correction")

# Print PCA plots
print(pca_plot_scaled)
print(pca_plot_filtered)

# Final VST Data for WGCNA
vsd <- vst(dds_scaled, blind = FALSE)
