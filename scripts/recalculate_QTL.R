#recalculate QTL
library(qtl)
setwd("/Users/Cody_2/git.repos/brassica_reflectance/data")

ref <- read.table("lmeCoefs2011Spectra.csv", sep = ",", header = TRUE)
head(ref)[1:2]

ref_sub <- ref[c("X", "X400nm", "X415nm","X435nm")]
head(ref_sub)
ref_sub$npqi <- (ref_sub$X415nm - ref_sub$X435nm)/ (ref_sub$X415nm + ref_sub$X435nm)
ref_sub$npqi2 <- (ref_sub$X400nm - ref_sub$X435nm)/ (ref_sub$X400nm + ref_sub$X435nm)
ref_sub
ref_sub[c(ref_sub$coefs.lineR500, ref_sub$coefs.lineIMB211)] 
ref_sub <- ref_sub[-c(49, 132),] 

ref_sub$RIL <- sub("(coefs.line)(\\d+)","RIL_\\2", ref_sub$X)
head(ref_sub)

ref_sub <- ref_sub[-1]
ref_sub <- ref_sub[c(6,1:5)]
head(ref_sub)

ref_sub_t <- t(ref_sub)
dim(ref_sub_t )
head(ref_sub_t)
rownames(ref_sub_t)[6] <- "id"

setwd("/Users/Cody_2/git.repos/brassica_reflectance/data")
write.table(ref_sub_t, "reflectance_traits.csv", row.names = TRUE, col.names = FALSE, sep = ",")

library(qtl)
brassica_ref <- read.cross("csvsr", genfile ="Brassica_F8_v2.1_gen.csv", 
	                       phefile="reflectance_traits.csv", 
	                       genotypes=c("AA","BB"), na.strings = "-")

head(brassica_ref)

class(brassica_ref)[1] <- "riself"
brassica_ref <- jittermap(brassica_ref)
str(brassica_ref)
brassica_ref$pheno

scanone.imp.1 <- scanone(brassica_ref, pheno.col = 4, method = "imp", use="all.obs", chr = "A10")
plot(scanone.imp.1)

scanone.imp.2 <- scanone(brassica_ref, pheno.col = 5, method = "imp", use="all.obs")
plot(scanone.imp.2)

# nothing obvious




