#############################################################
#
# Ref to the ARTICLE 
# 
# Code to compute calculations presented in: https://www.biorxiv.org/content/10.1101/2021.12.20.472907v3
#
# Panels 'a' Figure 5 and Supplementary Figures 13-18
#  
# max.coulter@xelect.co.uk
# d.bulgarelli@dundee.ac.uk 
#
#############################################################

############################################################################
# Clean-up the memory and start a new session (note the requirement below though)
#############################################################################

rm(list=ls())
dev.off()

#############################################################
# Libraries required
#############################################################

#this code requires outputs of the code https://github.com/BulgarelliD-Lab/Microbiota_mapping/tree/main/QRMC-3HS_Fig4_SData4 

#set working directory for revision (i.e., this will be user specific)
setwd("/cluster/db/R_shared/3H_manuscript/")

#import annotation data (available at )
annotation <- read.csv("/cluster/db/mecoulter/BaRT2v18/BaRT_2_18_annotation_genes.txt",sep = "\t", header = TRUE)
rownames(annotation) <- annotation$BaRTv2.gene
colnames(annotation)

#HEB124_52 - HEB_124_17 subset in HEB124_52-Barke comparison
col0='black'
col1='#E69F00'
col2='#56B4E9'
adj_p52<-genes_3D_stat$DE.pval[,3]
lgfc52<-genes_3D_stat$DE.lfc[,3]
chromosome = "chr3H"
title='HEB_124_52-HEB_124_17_subset_2only'

merged_52 <- as.data.frame(cbind(adj_p52,lgfc52))

#Get significant genes of interest
HEB52_HEB17_genes <- rownames(DE_genes)[which(DE_genes$contrast == "TissueHEB_124_52-TissueHEB_124_17")]
HEB52_HEB17 <- DE_genes[HEB52_HEB17_genes, ]
dim(HEB52_HEB17)
HEB52_Barke_genes <- rownames(DE_genes)[which(DE_genes$contrast == "TissueHEB_124_52-TissueBarke")]
HEB52_Barke <- DE_genes[HEB52_Barke_genes, ]
dim(HEB52_Barke)

HEB52_HEB17_HEB52_Barke <-intersect(HEB52_HEB17$target,HEB52_Barke$target)
length(HEB52_HEB17_HEB52_Barke)

#Add significance column with information whether up or down regulated

significance <- list()
all_genes <-rownames(merged_52)
for(i in 1:nrow(merged_52)) {
  if (all_genes[i] %in% HEB52_HEB17_HEB52_Barke){
    if (lgfc52[i] >= 1 ){
      significance[length(significance) + 1] <- "significant up"
    }else {
      if (lgfc52[i] <= -1){
        significance[length(significance) + 1] <- "significant down"
      } else {
        print("problem!")
      }
    }
  } else{
    significance[length(significance) + 1] <- "not significant"
  }
}





merged_52$significance <- unlist(significance)

merged_52$gene_name <- rownames(merged_52)
#gene2term <- gene2term %>% filter(gene_id %in% genes_all)
merged_522 <- merge(merged_52, annotation, by.x = "gene_name", by.y = "BaRTv2.gene", all.x = TRUE, all.y = FALSE)


for (i in 1:7){
  chromosome <- paste0("chr",i,"H")
  lgfc_merged_3H <- subset(merged_522, merged_522$Chromosome==chromosome,drop=TRUE)
  #NOw plot according to start position
  ylab='Log2FC'
  xlab=paste0('Position on ', chromosome, " (Mbp)")
  b_i<-100000000#break interval
  png(filename = paste0(figure.folder,"/",chromosome,title,"_lgfc_start.png"),width = 9000, height = 4000, units = "px", pointsize = 12,bg = "white", res = 1000)
  print({
    g <- ggplot(data = lgfc_merged_3H,aes(x=Start,y=lgfc52)) + geom_point(aes(colour=significance),size=1.5) + theme_bw(base_size = 16) + labs(x=xlab,y=ylab,title = title) + scale_y_continuous(breaks = pretty(lgfc_merged_3H$lgfc52, n = 10)) + scale_x_continuous(breaks=c(0,1*b_i,2*b_i,3*b_i,4*b_i,5*b_i,6*b_i), labels=c(0,100,200,300,400,500,600)) + scale_color_manual(values=c('not significant'=col0,'significant up'=col1,'significant down'=col2))
  })
  dev.off()
}

#individual figures merged with cognate chromosomal composition panel in Illustrator