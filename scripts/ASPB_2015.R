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
ref_traits_qb
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
?qb.hpdone
spec_pri_HDI <- qb.hpdone(ref_traits_qb_2011, effects = "cellmean")
spec_pri_HDI
warnings()



# make a merged plot for the 2 years
ref_2010_post$year <- paste("2010")
ref_2011_post$year <- paste("2011")

ref_post_merged <- merge(ref_2011_post, ref_2010_post, all = TRUE)
head(ref_post_merged)
tail(ref_post_merged)

peak2011 <- max(ref_post_merged$main)
ref_mergedbayes <- ggplot(ref_post_merged)
ref_mergedbayes <- ref_mergedbayes +  theme_bw() + geom_line(aes(x = pos, y = main, color = year), size = 1.5) +
                        geom_segment(aes(x = pos, xend = pos), y = (peak2011 * -0.02), yend = (peak2011 * -0.05)) +
                        theme(text = element_text(size = 20)) + ylab("Posterior") + 
                        xlab("Genetic Distance (cM)") + scale_colour_manual(values = c("red","black")) +
                        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

ref_mergedbayes


setwd("/Users/Cody_2/git.repos/brassica_reflectance/output/")
ggsave(file="NPQI_posterior.pdf")

setwd("/Users/Cody_2/git.repos/brassica_eqtl_v1.5/data")
#load data
load("eqtl_PAG.RData")

scanone_imp_p450 <- scanone(brass_total, pheno.col = "Bra009312", method = "imp", use="all.obs", chr = "A10")
plot(scanone_imp_p450)
scanone_imp_p450

setwd("/Users/Cody_2/git.repos/brassica_reflectance/data/")
specgenes <- read.csv("spectral_gene_candidates.csv", header = TRUE)
head(specgenes)

specgenes$Gene <- as.character(specgenes$Gene)
# QTL mapping on the genes in the reflectance QTL interval -----------
scanone_imp_p450_all <- scanone(brass_total, pheno.col = specgenes$Gene,
                         method = "imp", use = "all.obs")
#####Error in scanone(brass_total, pheno.col = specgenes$Gene, method = "imp", : 
# there are genes that are not expressed in our dataset.

# drop these non-expressed genes -----------
drops <- as.character(c("Bra009297", "Bra009304", "Bra1002808", "Bra1002809",
                        "Bra1002810", "Bra1002811", "Bra1002812",
                         "Bra1002813", "Bra1002814", "Bra009314", "Bra009324"))
drops

specgenes2 <- subset(specgenes, !(specgenes$Gene %in% drops))
specgenes2
str(specgenes2)
specgenes2$Gene <- as.character(specgenes2$Gene)

# QTL mapping on the cleaned data -----------
scanone_imp_p450_all <- scanone(brass_total, pheno.col = specgenes2$Gene,
                           method = "imp", use="all.obs", chr = "A10")
head(scanone_imp_p450_all)
tail(scanone_imp_p450_all)

# plot a few LOD traces to get a feel for things -----------
plot(scanone_imp_p450_all$Bra009331)
plot(scanone_imp_p450_all$Bra009326)
plot(scanone_imp_p450_all$Bra009312)

###bayes####
p450_traits_qb <- brass_total
class(p450_traits_qb)[1] <- "bc"
summary(p450_traits_qb)

p450_traits_qb <- qb.genoprob(p450_traits_qb, step=2, stepwidth = "variable")
summary(p450_traits_qb)


p450_pri <- qb.mcmc(p450_traits_qb, epistasis = TRUE, pheno.col = "Bra009312", seed = 1616,  n.iter = 3000)
plot(p450_pri)
p450_pri$pairs
?qb.scanone

p450_pri_post <- qb.scanone(p450_pri, type.scan = "posterior", epistasis = FALSE, smooth = 2)
p450_pri_post <- as.data.frame(p450_pri_post)
p450_pri_post <- subset(p450_pri_post, chr == "A10")
p450_pri_post
sum(p450_pri_post$main)
dim(p450_pri_post)

peak4 <- max(p450_pri_post$main)
p450_pribayes <- ggplot(p450_pri_post)
p450_pribayes <- p450_pribayes +  theme_bw() + geom_line(aes(x = pos, y = main), size = 2) +
                        geom_segment(aes(x = pos, xend = pos), y = (peak4 * -0.02), yend = (peak4 * -0.05)) +
                        theme(text = element_text(size = 20)) + ylab("Posterior") + 
                        xlab("Genetic Distance Along Chromosome")
                        # theme(legend.position = "none",
                        #   axis.text.x = element_text(angle = 90),
                        #   axis.line=element_line())
                        #   # panel.margin = unit(0, "cm")) +
                        # ggtitle("LOD Curves for QTLs") +
                  
                        #ylab("LOD Score") 
p450_pribayes
p450_HDI <- qb.hpdone(p450_pri, effects = "estimate")
p450_HDI


# other QTL

other_pri <- qb.mcmc(p450_traits_qb, epistasis = TRUE, pheno.col = "Bra009331", seed = 1616,  n.iter = 3000)
plot(other_pri)
other_pri$pairs
?qb.scanone

other_pri_post <- qb.scanone(other_pri, type.scan = "posterior", epistasis = FALSE, smooth = 2)
other_pri_post <- as.data.frame(other_pri_post)
other_pri_post <- subset(other_pri_post, chr == "A10")
other_pri_post
sum(other_pri_post$main)
dim(other_pri_post)

peak4 <- max(other_pri_post$main)
other_pribayes <- ggplot(other_pri_post)
other_pribayes <- other_pribayes +  theme_bw() + geom_line(aes(x = pos, y = main), size = 2) +
                        geom_segment(aes(x = pos, xend = pos), y = (peak4 * -0.02), yend = (peak4 * -0.05)) +
                        theme(text = element_text(size = 20)) + ylab("Posterior") + 
                        xlab("Genetic Distance Along Chromosome")
                        # theme(legend.position = "none",
                        #   axis.text.x = element_text(angle = 90),
                        #   axis.line=element_line())
                        #   # panel.margin = unit(0, "cm")) +
                        # ggtitle("LOD Curves for QTLs") +
                  
                        #ylab("LOD Score") 
other_pribayes

str(brass_total_qb)

other_HDI <- qb.hpdone(other_pri, effects = "estimate")
str(other_HDI)
other_HDI

spec_summary  <- as.data.frame(summary(spec_pri_HDI))
p450_summary  <- as.data.frame(summary(p450_HDI))
other_summary <- as.data.frame(summary(other_HDI))

p450_summary
?find.marker

mspec <- find.marker(combined, chr = spec_summary$chr, pos = spec_summary$pos)
m450 <-  find.marker(p450_traits_qb, chr = p450_summary$chr, pos = p450_summary$pos)
m_other <- find.marker(p450_traits_qb, chr = other_summary$chr, pos = other_summary$pos)
mspec
m_other
m450
markers <- list(mspec, m_other, m450)
markers

