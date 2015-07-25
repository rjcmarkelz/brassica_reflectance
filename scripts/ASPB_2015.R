# ASPB Poster
# also see wyoming_UV.R
library(qtl)
library(qtlbim)
library(ggplot2)
setwd("/Users/Cody_2/git.repos/brassica_eqtl_v1.5/data")
spec_traits <- read.cross("csvsr", genfile ="spec_gen.csv", 
                         phefile="spec_phe.csv", 
                         genotypes=c("AA","BB"), na.strings = c("-","NA"))
summary(spec_traits)
spec_traits

spec_traits_qb <- qb.genoprob(spec_traits, step=2, stepwidth = "variable")
summary(spec_traits_qb)
str(spec_traits_qb)
setwd("/Users/Cody_2/git.repos/brassica_reflectance/data/")

ref_2011 <- read.table("npqi_2011.csv", sep = ",", header = TRUE)
ref_2010 <- read.table("npqi2010.csv", sep = ",", header = TRUE)
head(ref_2011)
head(ref_2010)
dim(ref_2010)
dim(ref_2011)
str(ref_2011)
str(ref_2010)

ref_1011 <- merge(ref_2011, ref_2010, all.y = TRUE)
head(ref_1011)
dim(ref_1011)

ref_1011$Line
# now format for RQTL
ref_1011$Line <- sub("(\\d+)", "RIL_\\1", ref_1011$Line)
head(ref_1011$Line)

ref_1011_t <- as.data.frame(t(ref_1011))
ref_1011_t
head(ref_1011_t)
dim(ref_1011_t)
colnames(ref_1011_t)
dim(ref_1011_t)
ref_1011_t[5,] <- ref_1011_t[1,]
rownames(ref_1011_t)[5] <- "id"
rownames(ref_1011_t)
ref_1011_t <- ref_1011_t[-1,]
head(ref_1011_t)
tail(ref_1011_t)
ref_1011_t

setwd("/Users/Cody_2/git.repos/brassica_reflectance/data/")
write.table(ref_1011_t, "reflectance_2010_2011.csv", sep = ",", col.names = FALSE)

# Now pull into RQTL and do some Bayesian magic. 
library(qtl)
library(qtlbim)
ref_traits <- read.cross("csvsr", genfile ="Brassica_F8_v2.1_gen.csv", 
                         phefile="reflectance_2010_2011.csv", genotypes=c("AA","BB"), na.strings = c("-", "NA"))
head(ref_traits)
plot(ref_traits)

ref_traits_qb <- ref_traits
ref_traits_qb <- qb.genoprob(ref_traits_qb, step=2, stepwidth = "variable")
summary(ref_traits_qb)

ref_traits_qb_2010 <- qb.mcmc(ref_traits_qb, pheno.col = 2, seed = 1616, epistasis = FALSE)
#3000 iterations
#3000 iterations
summary(ref_traits_qb_2010)
str(ref_traits_qb_2010)
plot(ref_traits_qb_2010)
plot(qb.coda(ref_traits_qb_2010))

ref_2010 <- qb.scanone(ref_traits_qb_2010, type.scan = "posterior",  epistasis = FALSE)
plot(ref_2010)

ref_2010_post <- as.data.frame(ref_2010)
ref_2010_post <- subset(ref_2010_post, chr == "A10")
ref_2010_post
sum(ref_2010_post$main)
dim(ref_2010_post)

peak2010 <- max(ref_2010_post$main)
ref_2010bayes <- ggplot(ref_2010_post)
ref_2010bayes <- ref_2010bayes +  theme_bw() + geom_line(aes(x = pos, y = main), size = 2) +
                        geom_segment(aes(x = pos, xend = pos), y = (peak2010 * -0.02), yend = (peak2010 * -0.05)) +
                        theme(text = element_text(size = 20)) + ylab("Posterior") + 
                        xlab("Genetic Distance (cM)")
                        # theme(legend.position = "none",
                        #   axis.text.x = element_text(angle = 90),
                        #   axis.line=element_line())
                        #   # panel.margin = unit(0, "cm")) +
                        # ggtitle("LOD Curves for QTLs") +
                  
                        #ylab("LOD Score") 
ref_2010bayes


# 2011
ref_traits_qb_2011 <- qb.mcmc(ref_traits_qb, pheno.col = 1, seed = 1616, epistasis = FALSE)
#3000 iterations
#3000 iterations
summary(ref_traits_qb_2011)
str(ref_traits_qb_2011)
plot(ref_traits_qb_2011)
plot(qb.coda(ref_traits_qb_2011))

ref_2011 <- qb.scanone(ref_traits_qb_2011, type.scan = "posterior",  epistasis = FALSE)
plot(ref_2011)

ref_2011_post <- as.data.frame(ref_2011)
ref_2011_post <- subset(ref_2011_post, chr == "A10")
ref_2011_post
sum(ref_2011_post$main)
dim(ref_2011_post)

peak2011 <- max(ref_2011_post$main)
ref_2011bayes <- ggplot(ref_2011_post)
ref_2011bayes <- ref_2011bayes +  theme_bw() + geom_line(aes(x = pos, y = main), size = 2) +
                        geom_segment(aes(x = pos, xend = pos), y = (peak2011 * -0.02), yend = (peak2011 * -0.05)) +
                        theme(text = element_text(size = 20)) + ylab("Posterior") + 
                        xlab("Genetic Distance (cM)")
                        # theme(legend.position = "none",
                        #   axis.text.x = element_text(angle = 90),
                        #   axis.line=element_line())
                        #   # panel.margin = unit(0, "cm")) +
                        # ggtitle("LOD Curves for QTLs") +
                  
                        #ylab("LOD Score") 
ref_2011bayes



