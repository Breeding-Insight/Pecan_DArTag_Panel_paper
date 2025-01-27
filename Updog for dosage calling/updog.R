library(vcfR)
library(updog)
setwd("~/Desktop/Pecan-update/Updog")
setwd("/local/workdir/sc3339/Pecan")
setwd("/home/sc3339/sc3339")

df <- read.vcfR("DPec23-8107_MADC_rename_rm10plusHaps_filter_miss_snps_biallelic_sorted.vcf") #Chr12_026896592 all reads=0
d <- data.table::fread("DPec23-8107_MADC_rename_rm10plusHaps_filter_miss_snps_biallelic_sorted.vcf") 


Total_size <- extract.gt(df, element = "DP")
ref_size <- extract.gt(df, element = "RA")

class(Total_size) <- "numeric"
class(ref_size) <- "numeric"

## F1
Total_size_f1 <- Total_size[,grep(pattern = "2016|2017|Lakota|87MX3_2.11",colnames(Total_size))]
ref_size_f1 <- ref_size[,grep(pattern = "2016|2017|Lakota|87MX3_2.11",colnames(Total_size))]


mout <- multidog(refmat = ref_size_f1,
                 sizemat = Total_size_f1,
                 ploidy=2,
                 model = "f1",
                 p1_id="Lakota",
                 p2_id="87MX3_2.11",
                 nc = 50)

save(mout,file = "mout_pecan_F1parents.RData")

## Diverse
Total_size_dive <- Total_size[,-grep(pattern = "2016|2017|Lakota|87MX3_2.11",colnames(Total_size))]
ref_size_dive <- ref_size[,-grep(pattern = "2016|2017|Lakota|87MX3_2.11",colnames(Total_size))]

mout <- multidog(refmat = ref_size_dive,
                 sizemat = Total_size_dive,
                 ploidy=2,
                 model = "norm",
                 nc = 50
                 )


save(mout,file = "mout_pecan_Diversity.RData")

load("~/Desktop/Pecan-update/Updog/mout_pecan_Diversity.RData")
dim(mout$snpdf)

#geno <- format_multidog(mout, varname = "geno")
#write.csv(data.frame('Lakota'=mout$snpdf$p1geno,"87MX3_2.11"=mout$snpdf$p2geno),"parent-geno.csv")

hist(mout$snpdf$prop_mis,breaks = 50, xlab = "The estimated proportion of individuals misclassified in the SNP",ylab = "Count")
abline(v=0.05,col="red",lty=2,lwd=2)
axis(1, at=0.05,labels=T)

hist(mout$snpdf$bias,breaks = 50,xlab = "The estimated allele bias of the SNP",ylab = "Count")
abline(v=2,col="red",lty=2,lwd=2)
abline(v=0.05,col="red",lty=2,lwd=2)
axis(1, at=2,labels=T)

hist(mout$snpdf$od,breaks = 50,xlab = "The estimated overdispersion parameter of the SNP",ylab = "Count")

hist(mout$inddf$maxpostprob,breaks = 50,xlab = "The maximum posterior probability",ylab = "Count")

mout_cleaned <- filter_snp(mout, prop_mis < 0.05 & bias > 0.05 & bias < 2 & od < 0.05)
dim(mout_cleaned$snpdf)

save(mout_cleaned,file = "mout_cleaned_pecan_diversity_propmis0.05&bias0.05&2&od0.05.RData")

geno <- format_multidog(mout_cleaned, varname = "geno")
#write.csv(geno,"genotype-diversity.csv")

chrom_pos_inf<- data.frame(stringr::str_split_fixed(rownames(geno),pattern = "_",n=2))
pos_inf <- chrom_pos_inf$X2
#chrom_inf <- data.frame(stringr::str_remove_all(chrom_pos_inf$X1,"chr"))
chrom_inf <- chrom_pos_inf$X1

identical(mout_cleaned$snpdf$snp,rownames(geno))

d_geno <- cbind.data.frame("Lakota"=mout_cleaned$snpdf$p1geno,"87MX3_2.11"=mout_cleaned$snpdf$p2geno,"chrom"=chrom_inf,"sequence"=pos_inf,geno)
d_geno$ID <- rownames(d_geno)

library(data.table)
ref_alt <- df@fix[,3:5]
geno_merge <- merge(d_geno, ref_alt, by = "ID")

write.csv(geno_merge,file = "pecan_filter promiss bias od for mappoly.csv")



