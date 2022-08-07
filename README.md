# single_cell_Dynamic_eQTL
made by Juhyun Kim

This repository contains the code to analyze the Organoid single cell data to perform eQTL mapping through cell lineage(psuedotime). 

Inputs :
* genotype data (vcf)
* covariate data (plink covariate file )
* gene-snp pair ( ex. Tensor QTL results )
* Seurat Object to conduct slingshot

This method is comprised of 3 steps with first 1 using external tool named slignshot and the others using the code provided here.

1. Generating cell lineage and Lineage Quantile
2. Conducting eQTL analysis using mixed model with interaction term (genotype X lineage[pseudotime]). We use linear model and quadratic model to cover various biological scenario.
3. Multiple correction and plotting results. We conduct a lot of gene-SNP pairs, so before interpreting results, you have to correction the Pvalues from eQTL analysis(stpe 2). We use lose correction method, BH, and plot results by boxplot (lineage Quantile ~ gene expression value) when BH < 0.05. If you think the correction methods too lose, you can modify that methods. 

# Make Env
Download conda environment from this repository
```shell
conda env create --file scRNA_dynamiceQTL.juhyunk.20220807.yaml
```

# To run things step by step
## Generating cell lineage and Lineage Quantile

## Conducting eQTL analysis using mixed model with interaction term

## Multiple correction and plotting results




