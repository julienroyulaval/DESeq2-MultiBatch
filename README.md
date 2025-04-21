# DESeq2-MultiBatch
Batch Correction for Multi-Factorial RNA-seq Experiments

# Context

Here, you will find the main steps to reproduce the correction method implemented in "DESeq2-MultiBatch: Batch Correction for Multi-Factorial RNA-seq Experiments".

If used in research, please cite:

Roy, J., Monthony, A. S., Torkamaneh, D. (2025). DESeq2-MultiBatch: Batch Correction for Multi-Factorial RNA-seq Experiments.

For detailed steps and comparisons with other batch correction tools, please refer to the HTML file provided in this repository.
It also explains methods of contrast calling to extract meaningful results in DESeq2 using various designs applied to the experiment outlined below, with examples and figures.
Readers will also find various ressources to better understand design matrices and apply DESeq2 to their data at the end of the HTML.

For more information about DESeq2 and the standard RNA-seq analysis steps, refer to:

- [Analyzing RNA-seq data with DESeq2](https://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html) (Love & al., 2025)

- [RNA-seq workflow: gene-level exploratory analysis and differential expression](https://bioconductor.org/packages/release/workflows/vignettes/rnaseqGene/inst/doc/rnaseqGene.html) (Love & al., 2019)

## Experimental design

The simulated dataset used in this study consists of 1,000 genes across 48 samples. The dataset is designed to be computationally efficient, with confirmed scalability to larger real-life datasets.

The experimental design involves two batches, modeling gene expression in dioecious plants (e.g., XX/XY system) undergoing sex-specific treatments at two distinct time points: Day 0 (vegetative) and Day 1 (early flowering). Due to sex-specific treatments to alter flower sex, male and female plants were separated across batches:

Batch A: Female controls (n=6) and treated males (n=6)

Batch B: Male controls (n=6) and treated females (n=6)

The treatments were applied following Day 0 measurements but before Day 1 measurements, allowing Day 0 data to serve as a baseline for controlling batch effects. Samples were collected across two genotypes, with three biological replicates per combination of sex, treatment, genotype, and time point.

# Single and double batch correction

## Load data

```
# load metadata
coldata <- read.csv(file = "coldata.csv")
coldata[,2:6] <- lapply(coldata[,2:6], as.factor)

# load counts
counts <- read.csv(file = "counts.csv", row.names = 1)
```
Samples are rows in the ColData, with factors as columns. In the count file, rows are genes with samples as columns.

## Run DESeq2

### Single batch design

```
library(DESeq2)

design_single <- ~ Day + Batch + Cultivar + Sex + Sex:Day + Treatment

dds_single <- DESeqDataSetFromMatrix(countData = counts, colData = coldata, design = design_single)

dds_single <- DESeq(dds_single)
```

### Double batch design (interaction with sex)

```
design_double <- ~ Day + Batch + Cultivar + Sex + Sex:Batch + Sex:Day + Treatment

dds_double <- DESeqDataSetFromMatrix(countData = counts, colData = coldata, design = design_double)

dds_double <- DESeq(dds_double)
```

## Call for contrasts and extract scaling factors

```
resultsNames(dds_single)
resultsNames(dds_double)

# single
batch <- results(dds_single, name = "Batch_b_vs_a", independentFiltering = FALSE)
batch <- as.data.frame(batch)

batch$Scaling_Factor_Batch_A <- sqrt(2^(batch$log2FoldChange))
batch$Scaling_Factor_Batch_B <- 1 / batch$Scaling_Factor_Batch_A

# double
batch_xx <- results(dds_double, name = "Batch_b_vs_a", independentFiltering = FALSE)
batch_xx <- as.data.frame(batch_xx)
batch_xy <- results(dds_double, contrast=list(c("Batch_b_vs_a","Batchb.Sexxy")), independentFiltering = FALSE)
batch_xy <- as.data.frame(batch_xy)

batch_xx$Scaling_Factor_Batch_A <- sqrt(2^(batch_xx$log2FoldChange))
batch_xx$Scaling_Factor_Batch_B <- 1 / batch_xx$Scaling_Factor_Batch_A

batch_xy$Scaling_Factor_Batch_A <- sqrt(2^(batch_xy$log2FoldChange))
batch_xy$Scaling_Factor_Batch_B <- 1 / batch_xy$Scaling_Factor_Batch_A
```

## Apply corrections

```
# single
normalized_counts_df <- as.data.frame(counts(dds_single, normalized = TRUE))
scaled_counts_single <- normalized_counts_df

for (Sample_ID in colnames(normalized_counts_df)) {
  
  Batch <- coldata[Sample_ID, "Batch"]
  
  if (Batch == "a") {
    scaling_factors <- batch$Scaling_Factor_Batch_A
  } else if (Batch == "b") {
    scaling_factors <- batch$Scaling_Factor_Batch_B
  }
  
  scaled_counts_single[, Sample_ID] <- normalized_counts_df[, Sample_ID] * scaling_factors
}

# double
normalized_counts_df <- as.data.frame(counts(dds_double, normalized = TRUE))
scaled_counts_double <- normalized_counts_df
for (Sample_ID in colnames(normalized_counts_df)) {

  Sex <- coldata[Sample_ID, "Sex"]
  Batch <- coldata[Sample_ID, "Batch"]

  if (Sex == "xx" && Batch == "a") {
    scaling_factors <- batch_xx$Scaling_Factor_Batch_A
  } else if (Sex == "xx" && Batch == "b") {
    scaling_factors <- batch_xx$Scaling_Factor_Batch_B
  } else if (Sex == "xy" && Batch == "a") {
    scaling_factors <- batch_xy$Scaling_Factor_Batch_A
  } else if (Sex == "xy" && Batch == "b") {
    scaling_factors <- batch_xy$Scaling_Factor_Batch_B
  }

  scaled_counts_double[, Sample_ID] <- normalized_counts_df[, Sample_ID] * scaling_factors
}

```
Scaled counts can be rounded if needed, but they are not meant to be reused for differential expression analysis with DESeq2 or other tools.

While the contrasts obtained by using DESeq2 on the scaled counts with a design that excludes the *batch* factor will provide the same log2 fold changes as the original DESeq2 analysis, the *p* values and the adjusted *p* values will be skewed, leading to an increase in false positives.

For more detail, please refer to the HTML file.
