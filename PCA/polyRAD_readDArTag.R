library(VariantAnnotation)
library(polyRAD)
library(qqman)
library(dplyr)
library(data.table)

setwd("~/Desktop/Pecan-update/PCA")


# Importing the data
botloci=fread("Pecan_microhaplo.botloci",data.table=F,header = F)
botloci <- as.character(botloci)
#print(typeof(botloci))

data<-readDArTag("DPec23-8107_MADC_rename_rm10plusHaps_filter_miss.csv",
           botloci=botloci,
           sample.name.row = 1,
           possiblePloidies = list(2),
           taxaPloidy = 2L,
           n.header.rows = 0)


# Look at the RADdata object
data

# View the imported taxa names
samples <- GetTaxa(data)

# PCA plot to check sample clustering
source("/Users/sc3339-admin/Desktop/Pecan/PCA/AddPCA function.R")
#set.seed(2024)
data_pca <- AddPCA.RADdata(data)

library(tibble)
pca_df <- as.data.frame(data_pca$PCA)
pca_df <- tibble::rownames_to_column(pca_df, "Sample_ID")

pca_df_f1 <- pca_df[grep(pattern = "2016|2017|Lakota|87MX3_2.11",pca_df$Sample_ID),]
pca_df_dive <- pca_df[-grep(pattern = "2016|2017|Lakota|87MX3_2.11",pca_df$Sample_ID),]

pca_re <- rbind.data.frame(pca_df_f1,pca_df_dive)
d <-read.csv("DPec23-8107_MADC_rename_rm10plusHaps_filter_miss.csv")
d_dive <- d[,-grep(pattern = "2016|2017|Lakota|87MX3_2.11",colnames(d))]
ID_match <- colnames(d_dive[-c(1:3)])


#all pop

pca_re$Population <- rep(c("F1 (Carya illinoinensis)","Diverse"),c(190,185))
pca_re$Population <- factor(pca_re$Population,levels =c("F1 (Carya illinoinensis)","Diverse") )

p_population <- read.csv("Diversity_passport.csv")
p_population$Sample_ID <- ID_match

pa_re_p <- merge(pca_re,p_population,by="Sample_ID",all=T)
pa_re_p[which(pa_re_p$Population=="F1 (Carya illinoinensis)"),]$Species <- "Carya illinoinensis"

rownames(pa_re_p) <- pa_re_p[,1]
pa_re_p <- pa_re_p[,-1]

#View(pa_re_p[which(pa_re_p$Population=="Diverse"),])

library(ggrepel)
library(factoextra)
pa_re_p$label <- NA
# Then 'relabel' the points of interest
pa_re_p[which(grepl(pattern = "Lakota",rownames(pa_re_p))),]$label <- "Lakota"
pa_re_p[which(grepl(pattern = "87MX3_2.11",rownames(pa_re_p))),]$label <- "87MX3-2.11"

pa_re_p[which(grepl(pattern = "Carya illinoinensis × C. cordiformis",pa_re_p$Species)),]$label <- "Carya illinoinensis × C. cordiformis"
pa_re_p[which(grepl(pattern = "Carya illinoinensis × C. aquatica",pa_re_p$Species)),]$label <- "Carya illinoinensis × C. aquatica"
pa_re_p[which(grepl(pattern = "Carya cordiformis × C. ovata",pa_re_p$Species)),]$label <- "Carya cordiformis × C. ovata"
pa_re_p[which(grepl(pattern = "Carya illinoinensis × C. laciniosa",pa_re_p$Species)),]$label <- "Carya illinoinensis × C. laciniosa"

write.csv(pa_re_p,"PCA.csv")
#c("#fa9fb5","#addd8e","#fdae6b","#9ebcda","#8856a7","#e34a33","#1c9099")

png("PCA2.png",width  = 5000,height = 3400,res=300)

ggplot(pa_re_p,aes(pa_re_p$PC1,pa_re_p$PC2,color=Population,shape=Species)) + 
  geom_point(size=2) +
  scale_shape_manual(values = c(1:6,19,7:18),name="Species in diverse population")+
  scale_color_manual(values=c("#d7191c","#2c7bb6"))+
  #stat_ellipse(level = 0.95, show.legend = F)+
  #geom_text(aes(label=Population),vjust = "outward") + 
  #geom_hline(yintercept = 0,lty=2,col="red") + 
  #geom_vline(xintercept = 0,lty=2,col="blue",lwd=1) +
  theme_bw() + 
  labs(x="PC1 (24.0%)",y="PC2 (12.4%)",title="")+
  theme(text =element_text(size = 15),
        legend.position = "right",
     #   plot.title = element_text(hjust = 0.5, size = 20, face = "bold"),
        axis.title.x = element_text(face = "bold"), axis.title.y = element_text(face = "bold"),
        )+
  geom_text_repel(aes(label = pa_re_p$label),box.padding = 0.5,col="black")

dev.off()

