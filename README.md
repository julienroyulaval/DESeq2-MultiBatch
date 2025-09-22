# DESeq2-MultiBatch
Batch Correction for Multi-Factorial RNA-seq Experiments

# Context

Here, you will find the main steps to reproduce the correction method implemented in "DESeq2-MultiBatch: Batch Correction for Multi-Factorial RNA-seq Experiments".

If used in research, please cite:

Roy, J., Monthony, A. S., Torkamaneh, D. (2025). DESeq2-MultiBatch: Batch Correction for Multi-Factorial RNA-seq Experiments. https://doi.org/10.1101/2025.04.20.649392

For detailed steps and comparisons with other batch correction tools, please refer to the [HTML file provided](https://julienroyulaval.github.io/DESeq2-MultiBatch/docs/JR_data_method_v5_part1_08-09-2025.html) in the docs folder of this repository.

It also explains methods of contrast calling to extract meaningful results in DESeq2 using various designs applied to the experiment outlined below, with examples and figures.
Readers will also find various ressources to better understand design matrices and apply DESeq2 to their data at the end of the HTML.

For more information about DESeq2 and the standard RNA-seq analysis steps, refer to:

- [Analyzing RNA-seq data with DESeq2](https://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html) (Love & al., 2025)

- [RNA-seq workflow: gene-level exploratory analysis and differential expression](https://bioconductor.org/packages/release/workflows/vignettes/rnaseqGene/inst/doc/rnaseqGene.html) (Love & al., 2019)

## Experimental design

The simulated dataset used in this study consists of 1,000 genes across 48 samples. The dataset is designed to be computationally efficient, with confirmed scalability to larger real-life datasets.

An HTML with the steps for creating the data is available [here](https://julienroyulaval.github.io/DESeq2-MultiBatch/docs/JR_data_creation_html_08-09-2025.html), which can be found in the docs folder. It can be used to change the experimental design and factors.

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

design_single <- ~ Day + Batch + Genotype + Sex + Sex:Day + Treatment

dds_single <- DESeqDataSetFromMatrix(countData = counts, colData = coldata, design = design_single)

dds_single <- DESeq(dds_single)
```

### Double batch design (interaction with sex)

```
design_double <- ~ Day + Batch + Genotype + Sex + Sex:Batch + Sex:Day + Treatment

dds_double <- DESeqDataSetFromMatrix(countData = counts, colData = coldata, design = design_double)

dds_double <- DESeq(dds_double)
```

## MultiBatch function

```
MultiBatch <- function(
    dds,
    batch_factor,
    interact_factor = NULL,
    mode = c("center", "reference"),
    use_normalized = TRUE,
    round_out = TRUE,
    independentFiltering = FALSE,
    verbose = TRUE
) {
  stopifnot(is(dds, "DESeqDataSet"))
  mode <- match.arg(mode)
  
  rnames <- DESeq2::resultsNames(dds)
  cd <- as.data.frame(SummarizedExperiment::colData(dds))
  if (!batch_factor %in% names(cd))
    stop("batch_factor '", batch_factor, "' not found in colData(dds).")
  if (!is.null(interact_factor) && !interact_factor %in% names(cd))
    stop("interact_factor '", interact_factor, "' not found in colData(dds).")
  
  nc <- DESeq2::counts(dds, normalized = use_normalized)
  nc <- as.matrix(nc)
  G <- nrow(nc)
  
  bf <- as.factor(cd[[batch_factor]])
  b_levels <- levels(bf)
  if (is.null(b_levels) || length(b_levels) < 2)
    stop("Batch factor must have at least 2 levels.")
  b_ref <- b_levels[1]
  
  if (is.null(interact_factor)) {
    if (verbose) message("No interact_factor provided; adjusting Batch main effect only.")
    int_levels <- "(no-int)"
    i_ref <- "(no-int)"
    int_vec <- rep("(no-int)", nrow(cd))
  } else {
    intf <- as.factor(cd[[interact_factor]])
    int_levels <- levels(intf)
    i_ref <- int_levels[1]
    int_vec <- as.character(intf)
  }
  
  safe_lfc <- function(v) {
    v[is.na(v)] <- 0
    as.numeric(v)
  }
  
  get_lfc_by_name <- function(coef_name) {
    if (!(coef_name %in% rnames)) {
      if (verbose) message("  (missing coef) ", coef_name)
      return(rep(0, G))
    }
    res <- DESeq2::results(dds, name = coef_name, independentFiltering = independentFiltering)
    safe_lfc(res$log2FoldChange)
  }

  batch_coef_name <- function(batch_level) {
    paste0(batch_factor, "_", batch_level, "_vs_", b_ref)
  }
  
  interaction_coef_name <- function(b_lev, i_lev) {
    if (is.null(interact_factor) || b_lev == b_ref || i_lev == i_ref) return(NA_character_)
    candidates <- c(
      paste0(batch_factor, b_lev, ".", interact_factor, i_lev),
      paste0(interact_factor, i_lev, ".", batch_factor, b_lev)
    )
    hit <- candidates[candidates %in% rnames]
    if (length(hit) == 0) NA_character_ else hit[1]
  }
  
  if (verbose) message("Collecting main Batch coefficients relative to '", b_ref, "'...")
  main_batch_LFC <- vector("list", length(b_levels))
  names(main_batch_LFC) <- b_levels
  main_batch_LFC[[b_ref]] <- rep(0, G)
  for (bl in b_levels[b_levels != b_ref]) {
    nm <- batch_coef_name(bl)
    if (verbose) message("  ", nm)
    main_batch_LFC[[bl]] <- get_lfc_by_name(nm)
  }
  
  delta_by_intlevel <- vector("list", length(int_levels))
  names(delta_by_intlevel) <- int_levels
  
  if (is.null(interact_factor)) {
    M <- do.call(cbind, main_batch_LFC)
    colnames(M) <- b_levels
    delta_by_intlevel[[1]] <- M
  } else {
    if (verbose) message("Collecting Batch: ", interact_factor, " interaction coefficients...")
    for (il in int_levels) {
      # start from main batch effects
      M <- do.call(cbind, main_batch_LFC)
      colnames(M) <- b_levels
      
      if (il != i_ref) {
        for (bl in b_levels[b_levels != b_ref]) {
          nm_int <- interaction_coef_name(bl, il)
          if (!is.na(nm_int)) {
            if (verbose) message("  + ", nm_int)
            M[, bl] <- M[, bl] + get_lfc_by_name(nm_int)
          } else {
            if (verbose) message("  (no explicit coef for ", bl, "×", il, ")")
          }
        }
      }
      delta_by_intlevel[[il]] <- M
    }
  }

  if (verbose) message("Computing per-gene scaling factors (mode = '", mode, "').")
  scale_by_intlevel <- lapply(delta_by_intlevel, function(M) {
    if (mode == "center") {
      row_means <- rowMeans(M, na.rm = TRUE)
      sweep_M <- sweep(M, 1, row_means, FUN = "-")
      2^(-sweep_M)
    } else {
      2^(-M)
    }
  })

  if (verbose) message("Applying scaling factors to columns of the count matrix...")
  out <- nc
  sample_batches <- as.character(bf)
  for (j in seq_len(ncol(nc))) {
    bl <- sample_batches[j]
    il <- int_vec[j]
    if (is.null(interact_factor)) il <- "(no-int)"
    
    S <- scale_by_intlevel[[il]][, bl]
    out[, j] <- nc[, j] * S
  }
  
  if (round_out) {
    out <- round(out)
    storage.mode(out) <- "integer"
  }
  
  if (verbose) {
    message("Done.")
    message("Reference levels inferred from colData(): ",
            batch_factor, "='", b_ref, "'",
            if (!is.null(interact_factor)) paste0(", ", interact_factor, "='", i_ref, "'") else "")
  }

  list(
    adjusted_counts = out,
    scaling_tables = scale_by_intlevel,
    delta_tables   = delta_by_intlevel,
    batch_levels   = b_levels,
    interact_levels = int_levels,
    reference = list(batch = b_ref, interact = i_ref),
    mode = mode,
    used_results_names = rnames
  )
}

```

## Apply corrections

```
# single

adj_s <- MultiBatch(
  dds_single,
  batch_factor    = "Batch",
  interact_factor = NULL,
  mode = "center",
  use_normalized = TRUE,
  round_out = TRUE,
  verbose = TRUE)

scaled_counts_single <- adj_s$adjusted_counts

# double
adj_d <- MultiBatch(
  dds_double,
  batch_factor    = "Batch",
  interact_factor = "Sex",
  mode = "center",
  use_normalized = TRUE,
  round_out = TRUE,
  verbose = TRUE)

scaled_counts_double <- adj_d$adjusted_counts
```
Scaled counts can be rounded if needed, but they are not meant to be reused for differential expression analysis with DESeq2 or other tools.

While the contrasts obtained by using DESeq2 on the scaled counts with a design that excludes the *batch* factor will provide the same log2 fold changes as the original DESeq2 analysis, the *p* values and the adjusted *p* values will be skewed, leading to an increase in false positives.

For more detail, please refer to the HTML file.
