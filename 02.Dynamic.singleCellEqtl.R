library(tidyverse)
library(slingshot)
library(Seurat)
library(tradeSeq)
library(S4Vectors)
library(SingleCellExperiment)
library(SeuratDisk)
library(Seurat)
library(lme4)
library(effects)
library(dplyr)
library(ggplot2)
library('Rfast')
library(RhpcBLASctl)
blas_set_num_threads(20)

sdata <- readRDS("/mnt/gmi-l1/_90.User_Data/juhyunk/project/sceQTL/dynamiceQTL/sec-deut-cil/50.NRFonly.n11.sec.to.cil.slingshot.rds")

#Gene <- "SPATA20"
#snp = "chr17_50547162_A_C"

#FeaturePlot(sdata, features = Gene ,reduction = "phate", slot = "data")+ xlim(c(-0.03,0.06)) + ylim(c(-0.02,0.03))

mainFunction <- function(snp,Gene){

    	#Gene <- "SPATA20"
	#FeaturePlot(sdata, features = Gene ,reduction = "phate", slot = "data")+ xlim(c(-0.03,0.06)) + ylim(c(-0.02,0.03))
	#sdata <- AddMetaData(sdata, metadata = sdata@assays$RNA@data[Gene,], col.name = "Gene")
	sdata <- AddMetaData(sdata, metadata = sdata@assays$RNA@counts[Gene,], col.name = "Gene_umi")
	#sdata <- AddMetaData(sdata, metadata = sdata@assays$RNA@scale.data[Gene,], col.name = "Gene_scaleData")

        lineNum = as.integer(system(paste0("grep -n ",snp," /mnt/gmi-l1/_90.User_Data/jaeyong15/singlecell/NRF_all/kchip.lungNRF.final/NARD2/final.lungNRF.n68.NARD2.hg38.map | cut -d ':' -f 1"),
                                    intern= TRUE)) + 6
        lineNum

        sample_description <- system("grep -E 'CHROM|chr1_777778_T_A' /mnt/gmi-l1/_90.User_Data/juhyunk/project/sceQTL/dynamiceQTL/basalCell/final.lungNRF.n68.NARD2.hg38.vcf", intern=TRUE)
        sample_description <- as.data.frame(strsplit(sample_description, "\t"))[-1:-9,]
        colnames(sample_description) <- c("sampleid","genotype")
        rownames(sample_description) <- NULL
        sample_description$sampleid <- gsub("_hLRO.*","",sample_description$sampleid)
        head(sample_description)

        covar <- read.csv("/mnt/gmi-l1/_90.User_Data/jaeyong15/singlecell/NRF_all/01.automated/01.NRFonly/50.NRFonly.n11.sec.to.cil.RNA.PCA.tsv", sep='\t')
        colnames(covar) <- c("barcode", paste0("RNA_PC",1:5))
	
	sex <- read.csv("/gmi-l1/_90.User_Data/juhyunk/project/sceQTL/dynamiceQTL/sec-deut-cil/sex.covar", sep='\t', header=TRUE)

        genotypePC <- read.csv("/mnt/gmi-l1/_90.User_Data/jaeyong15/singlecell/NRF_all/kchip.lungNRF.final/NARD2/IBD.PCA/extract.final.lungNRF.n68.NARD2.hg38.eigenvec", sep=' ', header=FALSE)[,2:5]
        colnames(genotypePC) <- c("sampleid", paste0("GT_PC",1:3))

        sdata_metadata <- sdata@meta.data
        sdata_metadata <- cbind(sdata_metadata, rownames(sdata@meta.data))
        colnames(sdata_metadata) <- c(colnames(sdata@meta.data),"barcode")
        tab <- full_join(sdata_metadata, sample_description,by='sampleid', copy = TRUE)
        tab <- full_join(tab, genotypePC, by='sampleid', copy = TRUE)
	tab <- full_join(tab, sex, by='sampleid', copy = TRUE)
        tab <- full_join(covar, tab, by='barcode', copy = TRUE)
        colnames(tab) <- gsub(pattern = "\\.", "_",colnames(tab))
        tab <- tab[!is.na(tab$Lineage1), ]

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
        tab$Lineage1 <- as.numeric(tab$Lineage1)
        #tab$Lineage2_Quantile <- as.numeric(tab$Lineage2_Quantile)
        #tab["Lineage2_Quantile_character"] <- as.character(tab$Lineage2_Quantile)
        #tab["Lineage2_Quantile_poly2"] <- tab$Lineage2_Quantile**2


	
        tab$genotype_continue <- tab$genotype
        tab[tab$genotype_continue=="0/0","genotype_continue"]  <- 1
        tab[tab$genotype_continue=="0/1","genotype_continue"]  <- 2
        tab[tab$genotype_continue=="1/1","genotype_continue"]  <- 3
        tab$genotype_continue <- as.numeric(tab$genotype_continue)
        tab[,"poly_Lineage1"] <- tab$Lineage1**2
	tab <- tab[tab$Lineage1_10!= "NA",]
        head(tab)

	# 01. Nostate
	## Linear Poisson Nostate
	print("Linear Poisson Nostate")
	## full model
#	full_model_nostate <- glmer(formula = Gene_umi ~ genotype_continue + log(nCount_RNA_x) + percent_mt_x + sex + GT_PC1 +GT_PC2 +GT_PC3 + RNA_PC1 +RNA_PC2 + RNA_PC3 + RNA_PC4 + RNA_PC5  + 
#                    (1|sampleid_x) + (1|library),
#                    data=tab,
#                   family = "poisson", nAGQ=0, control=glmerControl(optimizer = "nloptwrap"))
	## reduced model
#	reduced_model_nostate <- glmer(formula = Gene_umi ~ log(nCount_RNA_x) + percent_mt_x + sex + GT_PC1 +GT_PC2 +GT_PC3 + RNA_PC1 +RNA_PC2 + RNA_PC3 + RNA_PC4 + RNA_PC5 + 
#                       (1|sampleid_x) + (1|library),
#                      data=tab,
#                   family = "poisson", nAGQ=0, control=glmerControl(optimizer = "nloptwrap"))

	# 02. Linear Poisson Univariate
	print("Linear Poisson univariate")
	## full model
	full_model_univariate <- glmer(formula = Gene_umi ~ genotype_continue + log(nCount_RNA_x) + percent_mt_x + sex + GT_PC1 +GT_PC2 +GT_PC3 + RNA_PC1 +RNA_PC2 + RNA_PC3 + RNA_PC4 + RNA_PC5  + 
                    (1|sampleid) + (1|orig_ident) +
                    Lineage1 + Lineage1:genotype_continue,
                    data=tab,
                   family = "poisson", nAGQ=0, control=glmerControl(optimizer = "nloptwrap"))
	## reduced model
	reduced_model_univariate <- glmer(formula = Gene_umi ~ genotype_continue + log(nCount_RNA_x) + percent_mt_x + sex + GT_PC1 +GT_PC2 +GT_PC3 + RNA_PC1 +RNA_PC2 + RNA_PC3 + RNA_PC4 + RNA_PC5 + 
                       (1|sampleid) + (1|orig_ident)+
                       Lineage1,
                      data=tab,
                   family = "poisson", nAGQ=0, control=glmerControl(optimizer = "nloptwrap"))

	# 03. Quadratic Poisson
	print("Quadratic Poisson Multivariates")
	# full model
	full_model_quadratic <- glmer(formula = Gene_umi ~ genotype_continue + log(nCount_RNA_x) + percent_mt_x + sex + GT_PC1 +GT_PC2 +GT_PC3 + RNA_PC1 +RNA_PC2 + RNA_PC3 + RNA_PC4 + RNA_PC5  + 
                    (1|sampleid) + (1|orig_ident) +
                    Lineage1 + Lineage1*genotype_continue + poly_Lineage1 + genotype_continue*poly_Lineage1,
                    data=tab,
                   family = "poisson", nAGQ=0, control=glmerControl(optimizer = "nloptwrap"))
	# reduced model
	reduced_model_quadratic <- glmer(formula = Gene_umi ~ genotype_continue + log(nCount_RNA_x) + percent_mt_x + sex + GT_PC1 +GT_PC2 +GT_PC3 + RNA_PC1 +RNA_PC2 + RNA_PC3 + RNA_PC4 + RNA_PC5 + 
                       (1|sampleid) + (1|orig_ident)+
                       Lineage1 +  poly_Lineage1,
                      data=tab,
                   family = "poisson", nAGQ=0, control=glmerControl(optimizer = "nloptwrap"))
	
	# Writing Results
#	anova_results.tmp <- as.data.frame(anova(full_model_nostate, reduced_model_nostate))
#	anova_results.tmp[,"SNP-eGene"] <- paste0(snp,"_",Gene)
#	anova_result <- data.frame()
#	anova_result <- rbind(anova_result, anova_results.tmp)
#	saveRDS(full_model, file=paste0("/mnt/gmi-l1/_90.User_Data/juhyunk/project/sceQTL/dynamiceQTL/basalCell/tmp/lme4_results/linear/",snp,"_",Gene,".linear.fullModel.rds"))
#	saveRDS(reduced_model, file=paste0("/mnt/gmi-l1/_90.User_Data/juhyunk/project/sceQTL/dynamiceQTL/basalCell/tmp/lme4_results/linear/",snp,"_",Gene,".linear.reducedModel.rds"))
#	write.csv(anova_result, paste0("/mnt/gmi-l1/_90.User_Data/juhyunk/project/sceQTL/dynamiceQTL/basalCell/tmp/lme4_results/noState/",snp,"_",Gene,".noState.anovarResults.txt") ,quote = F)
	
	anova_results.tmp <- as.data.frame(anova(full_model_univariate, reduced_model_univariate))
        anova_results.tmp[,"SNP-eGene"] <- paste0(snp,"_",Gene)
        anova_result <- data.frame()
        anova_result <- rbind(anova_result, anova_results.tmp)
#       saveRDS(full_model, file=paste0("/mnt/gmi-l1/_90.User_Data/juhyunk/project/sceQTL/dynamiceQTL/basalCell/tmp/lme4_results/linear/",snp,"_",Gene,".linear.fullModel.rds"))
#       saveRDS(reduced_model, file=paste0("/mnt/gmi-l1/_90.User_Data/juhyunk/project/sceQTL/dynamiceQTL/basalCell/tmp/lme4_results/linear/",snp,"_",Gene,".linear.reducedModel.rds"))
        write.csv(anova_result, paste0("/gmi-l1/_90.User_Data/juhyunk/project/sceQTL/dynamiceQTL/sec-deut-cil/linear/",snp,"_",Gene,".linear.anovarResults.txt") ,quote = F)
	
	# Writing Results
	anova_results.tmp <- as.data.frame(anova(full_model_quadratic,reduced_model_quadratic))
	anova_results.tmp[,"SNP-eGene"] <- paste0(snp,"_",Gene)
	anova_result <- data.frame()
	anova_result <- rbind(anova_result, anova_results.tmp)
	
#	saveRDS(full_model, file=paste0("/mnt/gmi-l1/_90.User_Data/juhyunk/project/sceQTL/dynamiceQTL/basalCell/tmp/lme4_results/quadratic/",snp,"_",Gene,".quadratic.fullModel.rds"))
#	saveRDS(reduced_model, file=paste0("/mnt/gmi-l1/_90.User_Data/juhyunk/project/sceQTL/dynamiceQTL/basalCell/tmp/lme4_results/quadratic/",snp,"_",Gene,".quadratic.reducedModel.rds"))
	write.csv(anova_result, paste0("/gmi-l1/_90.User_Data/juhyunk/project/sceQTL/dynamiceQTL/sec-deut-cil/quadratic/",snp,"_",Gene,".quadratic.anovarResults.txt") ,quote = F)
}



eQTL_results <- read.csv("/gmi-l1/_90.User_Data/juhyunk/project/sceQTL/dynamiceQTL/sec-deut-cil/PEER1.final_celltype.sec-deut-cil.n60.cis_qtl.txt", sep='\t', header =T)
head(eQTL_results)
eQTL_results <- eQTL_results[eQTL_results$qval < 0.1,c("phenotype_id","variant_id","qval")]

eQTL_results <- eQTL_results[order(eQTL_results$qval), ]
eQTL_results[,"Gene_id"] <- sub(".*@","",eQTL_results$phenotype_id)

head(eQTL_results)
dim(eQTL_results)[1]
for (i in 1:dim(eQTL_results)[1]) {
	print(eQTL_results[i,"phenotype_id"])
	snp=eQTL_results[i,"variant_id"]
	Gene=eQTL_results[i,"Gene_id"]

	sdata <- AddMetaData(sdata, metadata = sdata@assays$RNA@data[Gene,], col.name = "Gene_data")
	sdata <- AddMetaData(sdata, metadata = sdata@assays$RNA@counts[Gene,], col.name = "Gene_umi")
	#sdata <- AddMetaData(sdata, metadata = sdata@assays$RNA@scale.data[Gene,], col.name = "Gene_scaleData")

	exp_tab <- sdata@meta.data %>%
  		group_by(Lineage1_10) %>%
		summarise_at(vars(Gene_umi), list(name = mean))
	
	if ( sum(exp_tab[1:6,"name"]>0.3, na.rm = TRUE) > 3 ){
    		print(paste0("Gene name : ",Gene, " is expressed at least 3 Quantiles"))
		mainFunction(snp,Gene)
	} else {
		print(paste0("Gene name : ",Gene, " is not expressed at least 3 Quantiles"))	
	}
}


# version 4
