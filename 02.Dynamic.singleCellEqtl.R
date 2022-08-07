# Loading Library
print("Loading R Library")
suppressMessages(library(progress, quietly = T))
suppressMessages(library(tidyverse, quietly = T))
suppressMessages(library(slingshot, quietly = T))
suppressMessages(library(Seurat, quietly = T))
suppressMessages(library(tradeSeq, quietly = T))
suppressMessages(library(S4Vectors, quietly = T))
suppressMessages(library(SingleCellExperiment, quietly = T))
suppressMessages(library(SeuratDisk, quietly = T))
suppressMessages(library(Seurat, quietly = T))
suppressMessages(library(lme4, quietly = T))
suppressMessages(library(effects, quietly = T))
suppressMessages(library(dplyr, quietly = T))
suppressMessages(library(ggplot2, quietly = T))
suppressMessages(library('Rfast', quietly = T))
suppressMessages(library(RhpcBLASctl, quietly = T))
blas_set_num_threads(20)

# ----------------------------------------------------------------------------------- ARGUMENTS -----------------------------------------------------------------------------------------------------
print("Accept Argument from Command Line")
args = commandArgs(trailingOnly = TRUE)
sdata <- readRDS(args[1]) # slingshot RDS
eQTL_results <- read.csv(args[2], sep='\t', header = T) # eQTL results from JY
covar <- args[3]
wd <- args[4]
rnaPC <- args[5]

genotypePC <- "/mnt/gmi-analysis01/_90.User_Data/jaeyong15/singlecell/NRF_all/kchip.lungNRF.final/NARD2/IBD.PCA/extract.final.lungNRF.n68.NARD2.hg38.eigenvec"
mapFile <- "/mnt/gmi-analysis01/_90.User_Data/jaeyong15/singlecell/NRF_all/kchip.lungNRF.final/NARD2/final.lungNRF.n68.NARD2.hg38.map" 
vcf <- "/mnt/gmi-analysis01/_90.User_Data/juhyunk/project/sceQTL/dynamiceQTL/basalCell/final.lungNRF.n68.NARD2.hg38.vcf"


print(paste0("rds File : ", args[1]))
print(paste0("eQTL_results File : ", args[2]))
print(paste0("covar File : ", args[3]))
print(paste0("wd :", args[4]))

genotypePC <- read.csv(genotypePC, sep=' ', header=FALSE)[,2:5]
colnames(genotypePC) <- c("sampleid", paste0("GT_PC",1:3))
sex <- read.csv(covar, sep='\t', header=TRUE)
covar <- read.csv(rnaPC, sep='\t')
colnames(covar) <- c("barcode", paste0("RNA_PC",1:5))
#  ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

mainFunction <- function(snp,Gene){
	# SNP Genotype #
        sample_description <- system(paste0("grep -E 'CHROM|",snp,"' ",vcf), intern=TRUE)
        sample_description <- as.data.frame(strsplit(sample_description, "\t"))[-1:-9,]
        colnames(sample_description) <- c("sampleid","genotype")
        rownames(sample_description) <- NULL
        sample_description$sampleid <- gsub("_hLRO.*","",sample_description$sampleid)
	
	# RNA umi Extraction #
	sdata <- AddMetaData(sdata, metadata = sdata@assays$RNA@counts[Gene,], col.name = "Gene_umi")
        sdata_metadata <- sdata@meta.data
        sdata_metadata <- cbind(sdata_metadata, rownames(sdata@meta.data))
        colnames(sdata_metadata) <- c(colnames(sdata@meta.data),"barcode")

	# Table Joining # 
        tab <- full_join(sdata_metadata, sample_description,by='sampleid', copy = TRUE)
        tab <- full_join(tab, genotypePC, by='sampleid', copy = TRUE)
        tab <- full_join(tab, sex, by='sampleid', copy = TRUE)
        tab <- full_join(covar, tab, by='barcode', copy = TRUE)
        colnames(tab) <- gsub(pattern = "\\.", "_",colnames(tab))
        tab <- tab[!is.na(tab$DynamicLineage), ]
	
	# table formating #
        tab$nCount_RNA_x <- as.numeric(tab$nCount_RNA)
        tab$percent_mt_x <- as.numeric(tab$percent_mt)
        tab$GT_PC1 <- as.numeric(tab$GT_PC1)
        tab$GT_PC2 <- as.numeric(tab$GT_PC2)
        tab$GT_PC3 <- as.numeric(tab$GT_PC3)
        tab$RNA_PC1 <- as.numeric(tab$RNA_PC1)
        tab$RNA_PC2 <- as.numeric(tab$RNA_PC2)
        tab$RNA_PC3 <- as.numeric(tab$RNA_PC3)
        tab$RNA_PC4 <- as.numeric(tab$RNA_PC4)
        tab$RNA_PC5 <- as.numeric(tab$RNA_PC5)
        tab$DynamicLineage <- as.numeric(tab$DynamicLineage)
        tab$genotype_continue <- tab$genotype
        tab[tab$genotype_continue=="0/0","genotype_continue"]  <- 1
        tab[tab$genotype_continue=="0/1","genotype_continue"]  <- 2
        tab[tab$genotype_continue=="1/1","genotype_continue"]  <- 3
        tab$genotype_continue <- as.numeric(tab$genotype_continue)
        tab[,"poly_DynamicLineage"] <- tab$DynamicLineage**2
	tab[,"library"] <- tab$orig_ident
	
	print(paste0(snp," Ref Ref : ", table(tab$genotype)[1]))
	print(paste0(snp," Ref Alt : ", table(tab$genotype)[2]))
	print(paste0(snp," Alt Alt : ", table(tab$genotype)[3]))
	print(Gene)
	print(summary(tab$Gene_umi))
	
	# 02. Linear Poisson Univariate
	print("Linear Poisson univariate")
	## full model
	full_model_univariate <- glmer(formula = Gene_umi ~ genotype_continue + log(nCount_RNA_x) + percent_mt_x + sex + GT_PC1 +GT_PC2 +GT_PC3 + RNA_PC1 +RNA_PC2 + RNA_PC3 + RNA_PC4 + RNA_PC5  + 
                    (1|sampleid) + (1|library) +
                    DynamicLineage + DynamicLineage:genotype_continue,
                    data=tab,
                   family = "poisson", nAGQ=0, control=glmerControl(optimizer = "nloptwrap"))
	## reduced model
	reduced_model_univariate <- glmer(formula = Gene_umi ~ genotype_continue + log(nCount_RNA_x) + percent_mt_x + sex + GT_PC1 +GT_PC2 +GT_PC3 + RNA_PC1 +RNA_PC2 + RNA_PC3 + RNA_PC4 + RNA_PC5 + 
                       (1|sampleid) + (1|library)+
                       DynamicLineage,
                      data=tab,
                   family = "poisson", nAGQ=0, control=glmerControl(optimizer = "nloptwrap"))

	# 03. Quadratic Poisson
	print("Quadratic Poisson Multivariates")
	# full model
	full_model_quadratic <- glmer(formula = Gene_umi ~ genotype_continue + log(nCount_RNA_x) + percent_mt_x + sex + GT_PC1 +GT_PC2 +GT_PC3 + RNA_PC1 +RNA_PC2 + RNA_PC3 + RNA_PC4 + RNA_PC5  + 
                    (1|sampleid) + (1|library) +
                    DynamicLineage + DynamicLineage*genotype_continue + poly_DynamicLineage + genotype_continue*poly_DynamicLineage,
                    data=tab,
                   family = "poisson", nAGQ=0, control=glmerControl(optimizer = "nloptwrap"))
	# reduced model
	reduced_model_quadratic <- glmer(formula = Gene_umi ~ genotype_continue + log(nCount_RNA_x) + percent_mt_x + sex + GT_PC1 +GT_PC2 +GT_PC3 + RNA_PC1 +RNA_PC2 + RNA_PC3 + RNA_PC4 + RNA_PC5 + 
                       (1|sampleid) + (1|library)+
                       DynamicLineage +  poly_DynamicLineage,
                      data=tab,
                   family = "poisson", nAGQ=0, control=glmerControl(optimizer = "nloptwrap"))
	
	# Writing Results
	anova_results.tmp <- as.data.frame(anova(full_model_univariate, reduced_model_univariate))
        anova_results.tmp[,"SNP-eGene"] <- paste0(snp,"_",Gene)
        anova_result <- data.frame()
        anova_result <- rbind(anova_result, anova_results.tmp)
        write.csv(anova_result, paste0(wd,"/linear/",snp,"_",Gene,".linear.anovarResults.txt") ,quote = F)
	
	# Writing Results
	anova_results.tmp <- as.data.frame(anova(full_model_quadratic,reduced_model_quadratic))
	anova_results.tmp[,"SNP-eGene"] <- paste0(snp,"_",Gene)
	anova_result <- data.frame()
	anova_result <- rbind(anova_result, anova_results.tmp)
	write.csv(anova_result, paste0(wd,"/quadratic/",snp,"_",Gene,".quadratic.anovarResults.txt") ,quote = F)
}

# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

eQTL_results <- eQTL_results[eQTL_results$qval < 0.1,c("phenotype_id","variant_id","qval")]
eQTL_results <- eQTL_results[order(eQTL_results$qval), ]
eQTL_results <- unique(eQTL_results[,c("phenotype_id","variant_id")])

# Make Folder structure
print("Making Folder Structure under Working Directory - Linear / Quadratic / sig")
if (!dir.exists(paste0(wd,'/linear'))) {
	dir.create(paste0(wd,'/linear'))}

if (!dir.exists(paste0(wd,'/quadratic'))) {
        dir.create(paste0(wd,'/quadratic'))}

if (!dir.exists(paste0(wd,'/sig'))) {
        dir.create(paste0(wd,'/sig'))}

# Starting Assoc analysis
print("Starting Assoc Analysis!")


pb <- progress_bar$new(format = " Progress: [:bar] :percent, Estimated completion time: :eta", 
	total = dim(eQTL_results)[1], # totalnumber of ticks to complete (default 100)
	clear = FALSE, # whether to clear the progress bar on completion (default TRUE)
	width= 80) # width of the progress bar

for (i in 1:dim(eQTL_results)[1]) {
	snp=eQTL_results[i,"variant_id"]
	Gene=unlist(strsplit(eQTL_results[i,"phenotype_id"], split = "@"))[2]
	
	print(paste0(snp,"-",Gene))
	pb$tick()
	mainFunction(snp,Gene)
}

# v4 : Lineage name -> DynamicLineage / remove Gene expression filtering step / Accept arguments from command line / 
