# single_cell_Dynamic_eQTL
made by Juhyun Kim

This repository contains the code to analyze the Organoid single cell data to perform eQTL mapping through cell lineage(psuedotime). 

Inputs :
1. genotype data (vcf)
2. covariate data (plink covariate file )
3. gene-snp pair ( ex. Tensor QTL results )
4. Seurat Object to conduct slingshot

This method is comprised of 3 steps with first 1 using external tool named slignshot and the others using the code provided here.

1. Generating cell lineage and Lineage Quantile
2. Conducting eQTL analysis using mixed model with interaction term (genotype X lineage[pseudotime]). We use linear model and quadratic model to cover various biological scenario.

# Make Env
Download conda environment from this repository




