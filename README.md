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
