# single_cell_Dynamic_eQTL
made by Juhyun Kim

This repository contains the code to analyze the Organoid single cell data to perform eQTL mapping through cell lineage(pseudo time). 

Inputs :
* genotype data (vcf)
* covariate data (plink covariate file )
* gene-snp pair ( ex. Tensor QTL results )
* Seurat Object to conduct slingshot

This method consists of 3 steps with the first 1 using an external tool named slingshot and the others using the code provided here.

1. Generating cell lineage and Lineage Quantile
2. Conducting eQTL analysis using mixed model with interaction term (genotype X lineage[pseudotime]). We use linear models and quadratic models to cover various biological scenarios.
3. Multiple correction and plotting results. We conduct a lot of gene-SNP pairs, so before interpreting results, you have to correct the P Values from eQTL analysis(step 2). We use lose correction method, BH, and plot results by boxplot (lineage Quantile ~ gene expression value) when BH < 0.05. If you think the correction methods too lose, you can modify that methods. 

# Make Env
Download conda environment from this repository
```shell
conda env create --file scRNA_dynamiceQTL.juhyunk.20220807.yaml
```

# To run things step by step
## Generating cell lineage and Lineage Quantile
Using jupyter notebook, 01.dynamic_eQTL-slingshot.ipynb
In this step you have to decide which lineage to choose and how many cells to divide.


## Conducting eQTL analysis using mixed model with interaction term
We conduct eQTL analysis using interaction term
And we fit the model with linear and quadratic form because there are many various biological rhythms of RNA expression by cell lineage.
So we fit a null model(means, full model including interaction term[genotype X pseudo time] for linear and [genotype X pseudo time + genotype X pseudo time ^2] for quadratic model) and compare with reduced model(without that interaction terms of each model). If there are big differences between null model and reduced model, we think that the interaction term is significant, meaning that the RNA expression level is controlled by both genotype and lineage.

## Multiple correction and plotting results
We perform association analysis for a lot of gene-SNP pairs, so you have to conduct multiple corrections to control false positives. In this code, we apply BH to correct p value, but if you want to use other methods, you can modify this code. 
And Finally, we plot the results by Box plot when the gene-SNP’s BH < 0.05. The reason for this step is to confirm that gene expression is really controlled by both genotype and lineage by eyes.
