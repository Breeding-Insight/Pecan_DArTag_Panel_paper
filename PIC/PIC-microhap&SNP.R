library(tidyverse)
library(dplyr)
library(stringr)
library(ggpmisc)
setwd("~/Desktop/Pecan-update/PIC")
d_diplo <- read.csv("DPec23-8107_MADC_rename_rm10plusHaps_filter_miss.csv",check.names = F)
d_diplo_div <- d_diplo[,-grep(pattern = "2017|2016|Lakota|87MX3_2.11",colnames(d_diplo))]
#View(data.frame(table(d_diplo_div$CloneID)))

##rowSUM across all samples for each mHaplo
d_sum <- data.frame(CloneID=d_diplo_div$CloneID, SUM=rowSums(d_diplo_div[,-(1:3)]))
rownames(d_sum) <- d_diplo_div$AlleleID

## total reads count for each DArT locus
d_sum_total <- aggregate(d_sum[,2],by=d_sum[1],sum)

d_sum_total_merge <- merge(d_sum,d_sum_total,by="CloneID")

names(d_sum_total_merge)[3] <- c('TotalSum')
rownames(d_sum_total_merge) <- rownames(d_sum)

### Calculate frequecy for each mHaplo
d_sum_total_merge$Frequency <- d_sum_total_merge$SUM/d_sum_total_merge$TotalSum
View(data.frame(table(d_sum_total_merge$CloneID)))

d_fre_sum <- aggregate(d_sum_total_merge[,4],by=d_sum_total_merge[1],sum)

#d_sum_total_merge_1 <- d_sum_total_merge[-which(d_sum_total_merge$Frequency==1),]
#range <- aggregate(d_sum_total_merge_1[,4],by=d_sum_total_merge_1[1],FUN=mean)

### Function for PIC calculation
calc_pic <- function(x) {
  freq_squared <- x^2
  #cat("Allele Frequencies Squared: ", freq_squared, "\n")
  
  outer_matrix <- outer(freq_squared, freq_squared)
  #cat("Outer Matrix: \n")
  #print(outer_matrix)
  
  upper_tri_sum <- sum(outer_matrix[upper.tri(outer_matrix)])
  #cat("Sum of Upper Triangular: ", upper_tri_sum, "\n")
  
  pic <- 1 - sum(freq_squared) - 2*upper_tri_sum
  return(pic)
}

PIC <- aggregate(d_sum_total_merge[,4],by=d_sum_total_merge[1],FUN=calc_pic)
names(PIC)[2] <-  "PIC"
PIC_rm0 <- PIC[PIC$PIC>0,]
summary(PIC_rm0)
length(which(PIC_rm0$PIC>0.458983))

write.csv(PIC_rm0,"PIC_mHaplo.csv")

## the average microhaplotype frequencies
fre_mean <- aggregate(d_sum_total_merge[,4],by=d_sum_total_merge[1],FUN=mean)
fre_rm0 <- fre_mean[which(fre_mean$CloneID %in% PIC_rm0$CloneID),]

###SNP PIC
library(vcfR)
d_snp <- read.vcfR("/Users/sc3339-admin/Desktop/Pecan-update/Raw files/DPec23-8107_MADC_rename_rm10plusHaps_filter_miss_snps_biallelic_sorted.vcf")
Total_size <- extract.gt(d_snp, element = "DP")
ref_size <- extract.gt(d_snp, element = "RA")

class(Total_size) <- "numeric"
class(ref_size) <- "numeric"


## Diverse
Total_size_dive <- as.data.frame(Total_size[,-grep(pattern = "2016|2017|Lakota|87MX3_2.11",colnames(Total_size))])
ref_size_dive <- as.data.frame(ref_size[,-grep(pattern = "2016|2017|Lakota|87MX3_2.11",colnames(Total_size))])
identical(rownames(Total_size_dive),rownames(ref_size_dive ))
alt_size_dive  <- Total_size_dive-ref_size_dive

ref_size_dive$AlleleID <- paste0(rownames(ref_size_dive),"_","REF")
ref_size_dive$CloneID <- rownames(ref_size_dive)
alt_size_dive$AlleleID <- paste0(rownames(alt_size_dive),"_","ALT")
alt_size_dive$CloneID <- rownames(alt_size_dive)

size_dive <- rbind.data.frame(ref_size_dive,alt_size_dive)

## total reads count for each SNP
Total_size_sum <- aggregate(size_dive[,-c(186:187)],by=size_dive[187],sum)
length(unique(Total_size_sum$CloneID))
ref_size_sum <- size_dive[grep(pattern = "REF",size_dive$AlleleID),]
length(unique(ref_size_sum$CloneID))
alt_size_sum <- size_dive[grep(pattern = "ALT",size_dive$AlleleID),]
length(unique(alt_size_sum$CloneID))

## rowSUM 
Total_size_sum_t <- data.frame("CloneID"=Total_size_sum$CloneID,SUM=rowSums(Total_size_sum[,-1]))
alt_size_sum_t <- data.frame("CloneID"=alt_size_sum$CloneID,SUM=rowSums(alt_size_sum[,-c(186:187)]))
ref_size_sum_t <- data.frame("CloneID"=ref_size_sum$CloneID,SUM=rowSums(ref_size_sum[,-c(186:187)]))

## Frequecy
d_fre_ref <- data.frame("CloneID"=Total_size_sum_t$CloneID,"Frequency"=ref_size_sum_t$SUM/Total_size_sum_t$SUM)
d_fre_alt <- data.frame("CloneID"=Total_size_sum_t$CloneID,"Frequency"=alt_size_sum_t$SUM/Total_size_sum_t$SUM)


d_fre_rbind <- rbind.data.frame(d_fre_ref,d_fre_alt)

PIC_SNP <- aggregate(d_fre_rbind[,2],by=d_fre_rbind[1],FUN=calc_pic)
names(PIC_SNP)[2] <-  "PIC"
PIC_SNPrm0 <- PIC_SNP[PIC_SNP$PIC>0,]
summary(PIC_SNPrm0)
hist(PIC_SNPrm0$PIC,xlim = c(0,0.4))

#d_snp2 <-  read_csv("PIC_diversity_rmF1parents.csv")
#setdiff(PIC_SNPrm0$CloneID,d_snp2$SNP)

#fre <- d_fre_rbind[which(d_fre_rbind$CloneID %in% setdiff(PIC_SNPrm0$CloneID,d_snp2$SNP)),]
#hist(d_fre_ref[which(d_fre_ref$CloneID %in% setdiff(PIC_SNPrm0$CloneID,d_snp2$SNP)),]$Frequency)
#hist(d_fre_ref[-which(d_fre_ref$CloneID %in% setdiff(PIC_SNPrm0$CloneID,d_snp2$SNP)),]$Frequency)

write.csv(PIC_SNPrm0,"PIC_SNPreadcount_diversity.csv")

### Histogram plot for mHaplo and all SNPs

d_rbind <- rbind.data.frame(PIC_rm0,PIC_SNPrm0)
d_rbind$Type <- rep(c("Microhaplotype-based","SNP-based"),c(nrow(PIC_rm0),nrow(PIC_SNPrm0)))
  
d_rbind$Type <- factor(d_rbind$Type,levels =c("Microhaplotype-based","SNP-based"))
png("PIC_dartag_snp3.png",height = 2000,width = 2500,res = 300)
ggplot(d_rbind,aes(PIC, fill = Type)) +
  geom_histogram(alpha=0.4, position="identity",bins = 30)+
  theme_bw()+scale_fill_manual(values = c("red",  "blue"))+
  theme(text = element_text(size = 15),legend.position = "bottom")
dev.off()
### Histogram plot for mHaplo, all SNPs, and Target SNPs

TargetSNP <- read.csv("/Users/sc3339-admin/Desktop/Pecan/PIC/targetSNP_update.csv")
names(TargetSNP)[1] <- "CloneID"

SNP_Y <- merge(PIC_SNPrm0,TargetSNP,by="CloneID")
d_rbind_Y <- rbind.data.frame(PIC_rm0,PIC_SNPrm0,SNP_Y[,1:2])
d_rbind_Y$Type <- rep(c("DArTag loci","SNPs","Target SNPs"),c(nrow(PIC_rm0),nrow(PIC_SNPrm0),nrow(SNP_Y)))

d_rbind_Y$Type <- factor(d_rbind_Y$Type,levels =c("DArTag loci","SNPs","Target SNPs"))
ggplot(d_rbind_Y,aes(PIC, fill = Type)) + geom_histogram(alpha=0.3, position="identity",bins = 20)+
  theme_bw()+scale_fill_manual(values = c("red",  "blue","green"))+
  theme(text = element_text(size = 15),legend.position = "bottom")

### Scatter plot for mHaplo and Target SNPs

chr <- str_split_fixed(PIC_rm0$CloneID,pattern = "_",n=2)[,1]
pos <- str_split_fixed(PIC_rm0$CloneID,pattern = "_",n=2)[,2]
PIC_rm0$CloneID <- paste0(chr,"_",str_pad(pos,"9","left",pad = "0"))
write.csv(PIC_rm0,"PI_mHaplo.csv",row.names = F)

PIC_overlap <- PIC_rm0[which(PIC_rm0$CloneID %in% intersect(SNP_Y$CloneID,PIC_rm0$CloneID)),]
SNP_overlap <- SNP_Y[which(SNP_Y$CloneID %in% intersect(SNP_Y$CloneID,PIC_rm0$CloneID)),]

identical(PIC_overlap$CloneID,SNP_overlap$CloneID)

merge_all <- merge(PIC_overlap,SNP_overlap, by="CloneID")
names(merge_all)[2:3] <- c("microhaplotypes","TargetSNPs")

length(which(merge_all$microhaplotypes > merge_all$TargetSNPs))
length(which(merge_all$microhaplotypes == merge_all$TargetSNPs))
length(which(merge_all$microhaplotypes < merge_all$TargetSNPs))

#d2 <- merge_all[which(merge_all$microhaplotypes < merge_all$TargetSNPs),]

#chr <- str_split_fixed(d_diplo_div$CloneID,pattern = "_",n=2)[,1]
#pos <- str_split_fixed(d_diplo_div$CloneID,pattern = "_",n=2)[,2]
#d_diplo_div$CloneID2 <- paste0(chr,"_",str_pad(pos,"9","left",pad = "0"))
#d_micro <- d_diplo_div[which(d_diplo_div$CloneID2  %in% d2$CloneID),]

#d_SNPs <- size_dive[which(size_dive$CloneID  %in% d2$CloneID),]

#write.csv(d_micro,"d_micro2.csv")


png("PIC_dartag_targetsnp4.png",height = 2000,width = 2500,res = 300)
ggplot(merge_all, aes(x = microhaplotypes, y = TargetSNPs)) +
  geom_point(col = "#88000070") +
  #ylim(c(0,1))+
  # xlim(c(0,1))+
  theme_bw() +
  xlab("PIC (Microhaplotype-based)") +
  ylab("PIC (Taget SNP-based)") +
  theme(text = element_text(size = 15)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", col = "green",lwd=1.1) +
  geom_text(aes(x = 0.5, y = 0.5, label = "y = x"), color = "black", size = 6, hjust = 0.4, vjust = 0.5)
dev.off()



ggplot(merge_all, aes(x = DArTagloci, y = TargetSNPs)) +
  geom_point(col = "#88000070") +
  #ylim(c(0,1))+
 # xlim(c(0,1))+
  theme_bw() +
  xlab("PIC (DArTag loci)") +
  ylab("PIC (Target SNPs)") +
  theme(text = element_text(size = 15)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", col = "green",lwd=1.1) +
  geom_text(aes(x = 0.5, y = 0.5, label = "y = x"), color = "black", size = 6, hjust = 0.4, vjust = 0.5)+
  geom_smooth(method = "lm", se = F, col = "grey30") +
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")),
    label.x = "left", label.y = "top",
    size = 6, 
    formula = y ~ x, parse = TRUE
  )
