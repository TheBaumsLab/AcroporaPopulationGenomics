library('vcfR')
library('plyr')
library("StAMPP")
library("adegenet")


metadata <- read.csv("../stagdb_paper_popmap.csv")
colnames(metadata) = c("id", "Region")
vcf <- read.vcfR("../allsamples_lifted_to_jaAcrPala1.3_conffiltered_indvmissfiltered_downsampled_maffiltered_ldfiltered.vcf.gz")
x <- vcfR2genlight(vcf)
ploidy(x) <- 2
inds = data.frame(indNames(x))
colnames(inds) = "id"
metadata = join(inds, metadata)
pop(x) <- as.factor(c(metadata$Region))
fst <- stamppFst(x, nboots = 100, percent = 95, nclusters = 1)
write.table(fst$Fsts, file="../stagdb_paper_fst_matrix.txt", row.names=TRUE, col.names=TRUE)
write.table(fst$Bootstraps, file="../stagdb_paper_fst_matrix_boot.txt", row.names=TRUE, col.names=TRUE)
