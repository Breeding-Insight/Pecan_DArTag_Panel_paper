library(mappoly2)
setwd("~/Desktop/Pecan-update/Linkage map")

cranberry.f1 <- read_geno_csv(file.in = "pecan_filter promiss bias od for mappoly.csv",
                              ploidy.p1 = 2,
                              ploidy.p2 = 2,
                              name.p1 = "Lakota",
                              name.p2 = "X87MX3-2.11")

#write.csv(cranberry.f1$redundant,"redundant markers.csv")
#write.csv(cranberry.f1$mrk.names,"No. unique markers.csv")

plot(cranberry.f1)
cranberry.f2 <- filter_data(cranberry.f1, mrk.thresh = 0.02, ind.thresh = 0.02) 

#write.csv(cranberry.f2$screened.data$mrk.names,"Screened mrk.csv")

cranberry.f3  <- filter_individuals(cranberry.f2,type="Gmat") #Removing individual(s)  "X2017_01_0115" "X2017_01_0251" "X2017_01_0335" "X2017_01_0623" "X2017_01_0876" "X2017_01_1032"

plot(cranberry.f3, type = "density")
cranberry.all <- pairwise_rf(cranberry.f3, mrk.scope = "all", ncpus = 8)
plot(cranberry.all)

g <- group(x = cranberry.all, expected.groups = 16, comp.mat = T, inter = T)
g
plot(g)

##AllSNP_updog with od: delete markers: Chr15_022728777

s <- make_sequence(g, 
                   lg = list(1, c(2,7), c(8,10), c(4,8), 3, 3, c(5,13), 11, 14, c(4,6), 15, 9, 2, 2, 13, 2),#LG
                   ch = list(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16))

##AllSNP_updog with od
s <- make_sequence(g, 
                   lg = list(1, c(2,7), c(8,10), c(4,8), 3, 3, c(5,13), 11, 14, c(4,6), 15, 9, 2, 2, 13, 2),#LG
                   ch = list(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16))



s <- order_sequence(s, type = "genome")

plot_rf_matrix(s, type = "genome", fact = 5)
plot_mds_vs_genome(s)

s <- pairwise_phasing(s, 
                      type = "genome",
                      thresh.LOD.ph = 3, 
                      thresh.LOD.rf = 3, 
                      thresh.rf = 0.5, 
                      max.search.expansion.p1 = 10, 
                      max.search.expansion.p2 = 10)

print(s, type = "genome")


s <- mapping(s, type = "genome", parent = "p1p2", ncpus = 8)

s <- calc_haplotypes(s, type = "genome", ncpus = 8)
plot_haplotypes(s, lg = 1, ind = "Lakota")

png("Linkage map_pecan.png",height = 2000, width = 2400,res = 300)
plot_map_list(s, type = "genome", parent = "p1p2",col = mp_pal(16))
dev.off()
H <- plot_map_list(s, type = "genome", parent = "p1p2",col = mp_pal(16))
#write.csv(H,"Markers in map_with promiss&bias&od.csv")
png("mapVsphysical_pecan.png",height = 2400, width = 2400,res = 300)
plot_genome_vs_map(s, type="genome",same.ch.lg = T)
dev.off()


#map_summary(s,type = "mds")

map_summary(s,type = "genome")
t <- map_summary(s,type = "genome")
write.csv(t,"map_summary.csv")
