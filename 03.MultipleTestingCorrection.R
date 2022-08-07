# Loading Library
print("Loading R Library")
suppressMessages(library(slingshot))
suppressMessages(library(Seurat))
suppressMessages(library(S4Vectors))
suppressMessages(library(SingleCellExperiment))
suppressMessages(library(SeuratDisk))
suppressMessages(library(Seurat))
suppressMessages(library(lme4))
suppressMessages(library(effects))
suppressMessages(library(dplyr))
suppressMessages(library(ggplot2))
suppressMessages(library('Rfast'))
suppressMessages(library(RhpcBLASctl))
blas_set_num_threads(20)

# ----------------------------------------------------------------------------------- ARGUMENTS -----------------------------------------------------------------------------------------------------
print("Accept Argument from Command Line")
args = commandArgs(trailingOnly = TRUE)
wd <- args[1]
sdata <- readRDS(args[2])
sex <- args[3]
rnaPC <- args[4]



genotypePC <- "/mnt/gmi-analysis01/_90.User_Data/jaeyong15/singlecell/NRF_all/kchip.lungNRF.final/NARD2/IBD.PCA/extract.final.lungNRF.n68.NARD2.hg38.eigenvec"
mapFile <- "/mnt/gmi-analysis01/_90.User_Data/jaeyong15/singlecell/NRF_all/kchip.lungNRF.final/NARD2/final.lungNRF.n68.NARD2.hg38.map"
vcf <- "/mnt/gmi-analysis01/_90.User_Data/juhyunk/project/sceQTL/dynamiceQTL/basalCell/final.lungNRF.n68.NARD2.hg38.vcf"

genotypePC <- read.csv(genotypePC, sep=' ', header=FALSE)[,2:5]
colnames(genotypePC) <- c("sampleid", paste0("GT_PC",1:3))
sex <- read.csv(sex, sep='\t', header=TRUE)
rnaPC <- read.csv(rnaPC, sep='\t')
colnames(rnaPC) <- c("barcode", paste0("RNA_PC",1:5))
print("All Arguments are loaded")

#  ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

system(paste0("cat ", wd, "linear/*.linear.anovarResults.txt | fgrep full_model > ",wd,"linear/linearAllResult.csv"))
system(paste0("cat ", wd, "quadratic/*.quadratic.anovarResults.txt | fgrep full_model > ",wd,"quadratic/quadraticAllResult.csv"))

linearResult <- read.csv(paste0(wd, "linear/linearAllResult.csv"), header=FALSE)
quadraticResult <- read.csv(paste0(wd, "quadratic/quadraticAllResult.csv"), header=FALSE)

colnames(linearResult) <- c("model","npar","AIC","BIC","logLik","deviance","Chisq","Df","Pval","SNP-eGene")
colnames(quadraticResult) <- c("model","npar","AIC","BIC","logLik","deviance","Chisq","Df","Pval","SNP-eGene")

linearResult[,"BH"] <- p.adjust(linearResult$Pval, 'BH')
quadraticResult[,"BH"] <- p.adjust(quadraticResult$Pval, 'BH')

linearResult <- linearResult[linearResult$BH < 0.05,]
quadraticResult <- quadraticResult[quadraticResult$BH < 0.05,]

write.table(linearResult, paste0(wd, "linear/linearAllResult.csv"), sep="\t", quote= FALSE, row.names = FALSE)
write.table(quadraticResult, paste0(wd, "quadratic/quadraticAllResult.csv"), sep="\t", quote= FALSE, row.names = FALSE)

#  ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
summaryResults <- function(snp,Gene){
        # SNP Genotype #
	print("SNP Genotype")
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
        tab <- full_join(rnaPC, tab, by='barcode', copy = TRUE)
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
	print("Step 1")
        print(paste0(snp," Ref Ref : ", table(tab$genotype)[1]))
        print(paste0(snp," Ref Alt : ", table(tab$genotype)[2]))
        print(paste0(snp," Alt Alt : ", table(tab$genotype)[3]))
        print(Gene)
	head(tab)
        print(summary(tab$Gene_umi))
	print("dynamic?")
	
	tab$DynamicLineage_Quantile <- factor(tab$DynamicLineage_Quantile, levels = seq(1:length(unique(tab$DynamicLineage_Quantile))), order = T)
    	print("Step 2")
   	 # Feature plot # 
	FeaturePlot(sdata, features = Gene ,reduction = "phate", slot = "data", pt.size = 1.5)+ 
	xlim(c(-0.03,0.06)) + 
	ylim(c(-0.02,0.03))
    	ggsave(paste0(wd, "sig/",Gene,"_",snp,"_FeaturePlot.png"))
    
    	# make metaData # 
    	sdata <- AddMetaData(sdata, metadata = sdata@assays$RNA@data[Gene,], col.name = "Gene_data")
    	sdata <- AddMetaData(sdata, metadata = sdata@assays$RNA@counts[Gene,], col.name = "Gene_umi")
    
    	# Ridge Plot #
    	# possion distribution?
    	RidgePlot(sdata, features = c("Gene_umi","Gene_data"), ncol=2, group.by ="DynamicLineage_Quantile")
    	ggsave(paste0(wd,"sig/",Gene,"_",snp,"_RidgePlot.png"))
    
	############################################################################################################

        # 02. Linear Poisson Univariate
        print("Linear Poisson univariate")
        ## full model
        full_model_linear <- glmer(formula = Gene_umi ~ genotype_continue + log(nCount_RNA_x) + percent_mt_x + sex + GT_PC1 +GT_PC2 +GT_PC3 + RNA_PC1 +RNA_PC2 + RNA_PC3 + RNA_PC4 + RNA_PC5  +
                    (1|sampleid) + (1|library) +
                    DynamicLineage + DynamicLineage:genotype_continue,
                    data=tab,
                   family = "poisson", nAGQ=0, control=glmerControl(optimizer = "nloptwrap"))
        ## reduced model
        reduced_model_linear <- glmer(formula = Gene_umi ~ genotype_continue + log(nCount_RNA_x) + percent_mt_x + sex + GT_PC1 +GT_PC2 +GT_PC3 + RNA_PC1 +RNA_PC2 + RNA_PC3 + RNA_PC4 + RNA_PC5 +
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

	###########################################################
    	# grouped boxplot
    	ggplot(tab, aes(x=DynamicLineage_Quantile, y=log(Gene_umi+1), fill=genotype)) +
        	geom_boxplot()+
        	scale_x_discrete(limits=seq(1,length(unique(sdata@meta.data[,"DynamicLineage_Quantile"])))) + 
        	stat_summary(fun=mean, colour="red", aes(group=1),
                   	geom="line", lwd=1, lty=1) +
        	ggtitle(paste0(Gene)) + theme(plot.title = element_text(hjust = 0.5), 
                              axis.title.x = element_text(size = 16),
                              axis.text.x = element_text(size = 14),
                              axis.title.y = element_text(size = 16)) +
        	theme(axis.title=element_text(size=60)) +
        	theme_bw() +
        	labs(subtitle = paste(paste0("Linear : beta geno X Lieage = ", round(full_model_linear@beta[15],digits = 3)),
                              paste0("Linear : Pval geno X Lieage = ", round(as.data.frame(summary(full_model_linear)[10])[15,4],digits = 15)),
                              paste0("Linear : Pval of model diff = ", round(anova(full_model_linear, reduced_model_linear)["full_model_linear","Pr(>Chisq)"],digits = 15)),
                              paste0("Quadratic : beta geno X Lieage = ", round(full_model_quadratic@beta[17],digits = 3)),
                              paste0("Quadratic : Pval geno X Lieage = ", round(as.data.frame(summary(full_model_quadratic)[10])[17,4],digits = 15)),
                              paste0("Quadratic : Pval of model diff = ", round(anova(full_model_quadratic, reduced_model_quadratic)["full_model_quadratic","Pr(>Chisq)"],digits = 15)),
                              sep='\n'))
    ggsave(paste0(wd, "sig/",Gene,"_",snp,"_plot03.png"))

    ###########################################################
    # grouped boxplot
    ggplot(tab, aes(x=DynamicLineage_Quantile, y=log(Gene_umi+1), fill=genotype)) +
        geom_boxplot()+
        scale_x_discrete(limits=seq(1,length(unique(sdata@meta.data[,"DynamicLineage_Quantile"])))) + 
        stat_summary(fun=mean, colour="red", aes(group=1),
                   geom="line", lwd=1, lty=1) +
        ggtitle(paste0(Gene)) + theme(plot.title = element_text(hjust = 0.5), 
                              axis.title.x = element_text(size = 16),
                              axis.text.x = element_text(size = 14),
                              axis.title.y = element_text(size = 16)) +
        theme(axis.title=element_text(size=60)) +
        theme_bw() +
        labs(subtitle = paste(paste0("Linear : beta geno X Lieage = ", round(full_model_linear@beta[15],digits = 3)),
                              paste0("Linear : Pval geno X Lieage = ", round(as.data.frame(summary(full_model_linear)[10])[15,4],digits = 15)),
                              paste0("Linear : Pval of model diff = ", round(anova(full_model_linear, reduced_model_linear)["full_model_linear","Pr(>Chisq)"],digits = 15)),
                              paste0("Quadratic : beta geno X Lieage = ", round(full_model_quadratic@beta[17],digits = 3)),
                              paste0("Quadratic : Pval geno X Lieage = ", round(as.data.frame(summary(full_model_quadratic)[10])[17,4],digits = 15)),
                              paste0("Quadratic : Pval of model diff = ", round(anova(full_model_quadratic, reduced_model_quadratic)["full_model_quadratic","Pr(>Chisq)"],digits = 15)),
                              sep='\n')) +
        facet_grid( .~genotype )
    ggsave(paste0(wd,"sig/",Gene,"_",snp,"_plot01.png"))
    ###########################################################
    ggplot(tab, aes(x=genotype, y=log(Gene_umi+1), fill=DynamicLineage_Quantile)) + 
    geom_boxplot() +
    facet_wrap(~DynamicLineage_Quantile) +
    stat_summary(fun=mean, colour="red", aes(group=2),
               geom="line", lwd=1, lty=1)+ 
    
    ggtitle(Gene) + theme(plot.title = element_text(hjust = 0.5), 
                          axis.title.x = element_text(size = 16),
                          axis.text.x = element_text(size = 14),
                          axis.title.y = element_text(size = 16)) +
    theme(axis.title=element_text(size=60)) +
    theme_bw()+ 
    labs(subtitle = paste(paste0("Linear : beta geno X Lieage = ", round(full_model_linear@beta[15],digits = 3)),
                          paste0("Linear : Pval geno X Lieage = ", round(as.data.frame(summary(full_model_linear)[10])[15,4],digits = 15)),
                          paste0("Linear : Pval of model diff = ", round(anova(full_model_linear, reduced_model_linear)["full_model_linear","Pr(>Chisq)"],digits = 15)),
                          paste0("Quadratic : beta geno X Lieage = ", round(full_model_quadratic@beta[17],digits = 3)),
                          paste0("Quadratic : Pval geno X Lieage = ", round(as.data.frame(summary(full_model_quadratic)[10])[17,4],digits = 15)),
                          paste0("Quadratic : Pval of model diff = ", round(anova(full_model_quadratic, reduced_model_quadratic)["full_model_quadratic","Pr(>Chisq)"],digits = 15)),
                          sep='\n'))
    ggsave(paste0(wd, "sig/",Gene,"_",snp,"_plot02.png"))
    ###########################################################
    resultsDB <- data.frame(modelName = c(paste0(Gene,"_",snp)),
                        model = c("Linear","Quadratic"),
                        modelPval = c(anova(full_model_linear, reduced_model_linear)["full_model_linear","Pr(>Chisq)"], anova(full_model_quadratic, reduced_model_quadratic)["full_model_quadratic","Pr(>Chisq)"]),
                        interactionTermBeta = c(full_model_linear@beta[15], full_model_quadratic@beta[17]),
                        interactionTermPval = c(as.data.frame(summary(full_model_linear)[10])[15,4], as.data.frame(summary(full_model_quadratic)[10])[17,4]))
    resultsDB
    write.csv(resultsDB, paste0(wd, "sig/",Gene,"_",snp,".modelResult.txt") ,quote = F)

    ggplot(tab, aes(x=DynamicLineage, y=log(Gene_umi+1),color = genotype)) + 
            geom_point(size = 1, alpha = .1) +
            geom_smooth(aes(group=genotype), size = 2) + 
            theme_classic()+
            ggtitle(paste(Gene,"_",snp))
   ggsave(paste0(wd, "sig/",Gene,"_",snp,"_plot04.png"))

   ggplot(tab, aes(x=DynamicLineage_Quantile, y=log(Gene_umi+1), fill=genotype)) +
   geom_boxplot()+
   scale_x_discrete(limits=seq(1,length(unique(sdata@meta.data[,"DynamicLineage_Quantile"])))) +
#        stat_summary(fun=mean, colour="red", aes(group=1), geom="line", lwd=1, lty=1, ) +
   ggtitle(paste0(Gene)) + theme(plot.title = element_text(hjust = 0.5),
                              axis.title.x = element_text(size = 16),
                              axis.text.x = element_text(size = 14),
                              axis.title.y = element_text(size = 16)) +
   theme(axis.title=element_text(size=60)) +
   theme_bw() +
   labs(subtitle = paste(paste0(Gene," : " , snp , " Linear"),
                              paste0("Linear : beta geno X Lieage = ", round(full_model_linear@beta[15],digits = 3)),
                              paste0("Linear : Pval geno X Lieage = ", round(as.data.frame(summary(full_model_linear)[10])[15,4],digits = 15)),
                              paste0("Linear : Pval of model diff = ", round(anova(full_model_linear, reduced_model_linear)["full_model_linear","Pr(>Chisq)"],digits = 15)),
                              paste0("Quadratic : beta geno X Lieage = ", round(full_model_quadratic@beta[17],digits = 3)),
                              paste0("Quadratic : Pval geno X Lieage = ", round(as.data.frame(summary(full_model_quadratic)[10])[17,4],digits = 15)),
                              paste0("Quadratic : Pval of model diff = ", round(anova(full_model_quadratic, reduced_model_quadratic)["full_model_quadratic","Pr(>Chisq)"],digits = 15)),
                              sep='\n')) +
   facet_grid( .~genotype ) +
   stat_smooth(aes(group = genotype),method = "lm", formula = y ~ x , size = 2, colour = "black") 
   ggsave(paste0(wd, "sig/",Gene,"_",snp,"_plot05.png"))	

   ggplot(tab, aes(x=DynamicLineage_Quantile, y=log(Gene_umi+1), fill=genotype)) +
        geom_boxplot()+
        scale_x_discrete(limits=seq(1,length(unique(sdata@meta.data[,"DynamicLineage_Quantile"])))) +
#        stat_summary(fun=mean, colour="red", aes(group=1), geom="line", lwd=1, lty=1, ) +
        ggtitle(paste0(Gene)) + theme(plot.title = element_text(hjust = 0.5),
                              axis.title.x = element_text(size = 16),
                              axis.text.x = element_text(size = 14),
                              axis.title.y = element_text(size = 16)) +
        theme(axis.title=element_text(size=60)) +
        theme_bw() +
        labs(subtitle = paste(paste0(Gene," : " , snp , " Quadratic"),
                              paste0("Linear : beta geno X Lieage = ", round(full_model_linear@beta[15],digits = 3)),
                              paste0("Linear : Pval geno X Lieage = ", round(as.data.frame(summary(full_model_linear)[10])[15,4],digits = 15)),
                              paste0("Linear : Pval of model diff = ", round(anova(full_model_linear, reduced_model_linear)["full_model_linear","Pr(>Chisq)"],digits = 15)),
                              paste0("Quadratic : beta geno X Lieage = ", round(full_model_quadratic@beta[17],digits = 3)),
                              paste0("Quadratic : Pval geno X Lieage = ", round(as.data.frame(summary(full_model_quadratic)[10])[17,4],digits = 15)),
                              paste0("Quadratic : Pval of model diff = ", round(anova(full_model_quadratic, reduced_model_quadratic)["full_model_quadratic","Pr(>Chisq)"],digits = 15)),
                              sep='\n')) +
        facet_grid( .~genotype ) +
        stat_smooth(aes(group = genotype),method = "lm", formula = y ~ x + I(x^2), size = 2, colour = "black") 
   ggsave(paste0(wd, "sig/",Gene,"_",snp,"_plot06.png"))
}

# -----------------------------------------------------------------------------------------------------------------------------------------------------------------------

for (pair in rbind(linearResult,quadraticResult)[,"SNP-eGene"] %>% unique() ) {
    if ( "DynamicLineage_Quantile" %in% colnames(sdata@meta.data)) {
    	print(pair)
    	snp=paste0(unlist(strsplit(pair,split="_"))[1:4],collapse = "_")
    	Gene=unlist(strsplit(pair,split="_"))[5]
    	summaryResults(snp,Gene)
    } else {
    	print("No DynamicLineage_Quantile in MetaData")
    }
}


# version 2
