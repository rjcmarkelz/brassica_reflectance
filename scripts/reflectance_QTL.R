##########
# Cody Markelz and Rob Baker
# markelz@gmail.com
# reflectance QTL for brassica rapa
##########
library(qtlbim)
library(ggplot2)

head(myspecdata)
class(myspecdata)
summary(myspecdata)
plot.map(myspecdata)

?scanone
plot(out.em)
?pull.pheno
plot(myspecdata)
?write.cross
setwd("/Users/Cody_2/git.repos/brassica_eqtl_v1.5/data")
write.cross(myspecdata, format = "csvsr", filestem = "spec")

out_em_2 <- scanone(myspecdata, chr = "A10", pheno.col = 6, method = "em", chr = "A10")
plot(out_em_2)
out_em_2

out_em_2 <- scanone(myspecdata, chr = "A10", pheno.col = 6, method = "em", chr = "A10")

A10_16069633
effectplot(myspecdata, mname1 = "A10_16069633", pheno.col = 6)
# effect plots
efplot <- effectplot(myspecdata, mname1 = "A10_16069633", pheno.col = 6)
efplot[2]
efplotdf <- as.data.frame(efplot[1])
efplotdf$SE <- as.data.frame(efplot[2])
efplotdf$marker <- row.names(efplotdf)
efplotdf

limits <- aes(ymax= efplotdf$Means + 0.0011021828, ymin = efplotdf$Means - 0.0011021828)
dodge <- position_dodge(width=0.9)
plotwidth5 <- ggplot(data = efplotdf,  aes(y=Means, x=marker, fill = marker))
plotwidth5 <- plotwidth5 + geom_bar(position=dodge, stat="identity") 
plotwidth5 <- plotwidth5 + geom_errorbar(limits, position = dodge, width = 0.25) 
plotwidth5 <- plotwidth5 + xlab("Alleles") + ylab("Absorption")
plotwidth5 <- plotwidth5 + ggtitle("NPQI Absorption") + guides(color = FALSE)

plotwidth5 <- plotwidth5 + theme(axis.title=element_text(face="bold",
                               size="20"), axis.text=element_text(face="bold",
                               size="20"))  
plotwidth5


peakout_em_2 <- max(out_em_2$lod)
spectral_qtl <- ggplot(out_em_2)
spectral_qtl <- spectral_qtl +  theme_bw() + geom_line(aes(x = pos, y = lod), size = 2) + facet_grid(~ chr) +
                        geom_hline(yintercept = 3.23, color = "red", size = 1) +
                        geom_segment(aes(x = pos, xend = pos), y = (peakout_em_2 * -0.02), yend = (peakout_em_2 * -0.05)) +
                        theme(text = element_text(size = 20))
                        # theme(legend.position = "none",
                        #   axis.text.x = element_text(angle = 90),
                        #   axis.line=element_line())
                        #   # panel.margin = unit(0, "cm")) +
                        # ggtitle("LOD Curves for QTLs") +
                        #xlab("Position in cM") +
                        #ylab("LOD Score") 
spectral_qtl
?bayesint
bayesint(out_em_2, "A10", prob = 0.98, expandtomarkers = TRUE)
bayesint(out_em_2, "A10", prob = 0.9999, expandtomarkers = TRUE)
?effectplot

# > bayesint(out_em_2, "A10", prob = 0.98, expandtomarkers = TRUE)
#              chr      pos      lod
# A10_16230621 A10 82.90340 13.10773
# A10_16069633 A10 83.30983 14.08732
# A10_16032467 A10 84.54993 12.45285


spec_traits
effectplot(spec_traits, pheno.col = 6, mname1 = "A10_16069633")
?plot.pxg
mname <- "A10_16069633"
plot.pxg(spec_traits, mname, pheno.col = 6,  jitter = 1)


?pull.map
chr_A10 <- pull.map(spec_traits, chr = "A10")
plot(chr_A10)
chr_A10

read.table("spectral_gene_candidates.csv", header = TRUE,  sep = ",",  stringsAsFactors = FALSE)
setwd("/Users/Cody_2/git.repos/brassica_reflectance/")
candidates <- read.table("spectral_gene_candidates.csv", header = TRUE,  sep = ",",  stringsAsFactors = FALSE)
head(candidates)

descriptions <- read.table("Brassica_rapa_final.annotation.txt", header = FALSE,  sep = ",",  stringsAsFactors = FALSE)
head(descriptions)

spec_traits <- read.cross("csvsr", genfile ="spec_gen.csv", 
                         phefile="spec_phe.csv", 
                         genotypes=c("AA","BB"), na.strings = c("-","NA"))
summary(spec_traits)
spec_traits

names(spec_traits$pheno)



spec_traits_qb <- qb.genoprob(spec_traits, step=2, stepwidth = "variable")
summary(spec_traits_qb)
str(spec_traits_qb)

spec_pri <- qb.mcmc(spec_traits_qb, epistasis = TRUE, pheno.col = 6, seed = 1616,  n.iter = 10000)
plot(spec_pri)
summary(spec_pri)
?qb.scanone

# spec_traits
# effectplot(spec_traits, pheno.col = 6, mname1 = "A10_16069633")
# ?plot.pxg
# mname <- "A10_16069633"
# plot.pxg(spec_traits, mname, pheno.col = 6,  jitter = 1)

# ?pull.geno
# geno_10 <- as.data.frame(pull.geno(spec_traits, chr = "A10"))
# str(geno_10)
# names(geno_10)

# A10_marker <- as.data.frame(geno_10[,"A10_16069633"])
# A10_marker

# A10_pheno <- as.data.frame(pull.pheno(spec_traits, pheno.col = c(1,6)))
# A10_pheno
# A10_16069633

# A10_pheno$geno <- A10_marker[,1]
# A10_pheno

# A10_pheno2 <- A10_pheno[order(A10_pheno$npqi_UN),] 
# A10_pheno2

# setwd("~/git.repos/brassica_meta_analysis/Output/")  
# ?write.table
# write.table(A10_pheno2, "A10_geno_pheno.csv")

spec_pri_post <- qb.scanone(spec_pri, type.scan = "posterior",  epistasis = FALSE, smooth = 2, chr = "A10")
spec_pri_post <- as.data.frame(spec_pri_post)
head(spec_pri_post)
spec_pri_post

spec_pri_post <- subset(spec_pri_post, chr == "A10")
spec_pri_post
sum(spec_pri_post$main)
dim(spec_pri_post)

peak4 <- max(spec_pri_post$main)
spec_pribayes <- ggplot(spec_pri_post)
spec_pribayes <- spec_pribayes +  theme_bw() + geom_line(aes(x = pos, y = main), size = 2) +
                        geom_segment(aes(x = pos, xend = pos), y = (peak4 * -0.02), yend = (peak4 * -0.05)) +
                        theme(text = element_text(size = 20)) + ylab("Posterior") + 
                        xlab("Genetic Distance Along Chromosome")
                        # theme(legend.position = "none",
                        #   axis.text.x = element_text(angle = 90),
                        #   axis.line=element_line())
                        #   # panel.margin = unit(0, "cm")) +
                        # ggtitle("LOD Curves for QTLs") +
                  
                        #ylab("LOD Score") 
spec_pribayes


###############
###############
###############
#### genes
###############
###############
###############

setwd("/Users/Cody_2/git.repos/brassica_eqtl_v1.5/data")
#load data
load("eqtl_PAG.RData")
ls()

scanone_imp_p450 <- scanone(brass_total, pheno.col = "Bra009312", method = "imp", use="all.obs", chr = "A10")
plot(scanone_imp_p450)
scanone_imp_p450

peak <- max(scanone_imp_p450$lod)
peak 

p450_plot <- ggplot(scanone_imp_p450)
p450_plot <- p450_plot +  theme_bw() + geom_line(aes(x = pos, y = lod), size = 2) + 
                        geom_hline(yintercept = 3.23, color = "red", size = 1) +
                        geom_segment(aes(x = pos, xend = pos), y = (peak * -0.02), yend = (peak * -0.05)) +
                        theme(text = element_text(size = 20))
p450_plot
head(br_phys)

###bayes####
p450_traits_qb <- brass_total
class(p450_traits_qb)[1] <- "bc"

p450_traits_qb <- qb.genoprob(p450_traits_qb, step=2, stepwidth = "variable")
summary(p450_traits_qb)

summary()

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

str(brass_total_qb)

p450_HDI <- qb.hpdone(p450_pri, effects = "estimate")
p450_HDI



spec_pri_HDI <- qb.hpdone(spec_pri, effects = "estimate")
spec_pri_HDI

spec_summary  <- as.data.frame(summary(spec_pri_HDI))
p450_summary  <- as.data.frame(summary(p450_HDI))


spec_traits_qb
summary(p450_traits_qb)
summary(spec_traits_qb)

spec_traits_qb$pheno$id
p450_traits_qb$pheno$id2 <- as.numeric(sub("(RIL_)(\\d+)", "\\2", p450_traits_qb$pheno$id))
p450_traits_qb$pheno$id2


combined <- p450_traits_qb
p450_red <- pull.pheno(p450_traits_qb, pheno.col = c("id", "id2", "Bra009331", "Bra009312"))
p450_red

spec_traits_qb$pheno
spec_red <- pull.pheno(spec_traits_qb, pheno.col = c(1, 6))
spec_red
merged_pheno <- merge(p450_red, spec_red, by.x = "id2", by.y = "id")
merged_pheno
merged_pheno <- merged_pheno[order(merged_pheno$id),]
merged_pheno$doublcheck <- as.numeric(sub("(RIL_)(\\d+)", "\\2", p450_traits_qb$pheno$id))
merged_pheno
combined$pheno <- merged_pheno[,c(2,5,4)]
combined$pheno


combined_qb <- qb.genoprob(combined, step=2, stepwidth = "variable")
summary(combined_qb)
str(combined_qb)

combined_qb_spec <- qb.mcmc(combined_qb, epistasis = TRUE, pheno.col = 2, seed = 1616,  n.iter = 10000)
plot(combined_qb_spec)
summary(combined_qb_spec)
?qb.scanone

combined_qb_gene <- qb.mcmc(combined_qb, epistasis = TRUE, pheno.col = 3, seed = 1616,  n.iter = 10000)
plot(combined_qb_gene)
summary(combined_qb_gene)

spec_HDI <- qb.hpdone(combined_qb_spec, effects = "estimate")
str(spec_HDI)
spec_HDI

spec_post <- qb.scanone(combined_qb_spec, type.scan = "posterior", epistasis = FALSE, smooth = 2, chr = "A10")
spec_post <- as.data.frame(spec_post)
spec_post
#cA10.A10x14408017

gene_post <- qb.scanone(combined_qb_gene, type.scan = "posterior", epistasis = FALSE, smooth = 2, chr = "A10")
gene_post <- as.data.frame(gene_post)
gene_post
#cA10.A10x14220630

gene_HDI <- qb.hpdone(combined_qb_gene, effects = "estimate")
str(gene_HDI)
gene_HDI







combined
effectplot(combined, pheno.col = 2, mname1 = "A10x14408017")
mname <- "A10x14408017"
plot.pxg(combined, mname, pheno.col = 2,  jitter = 1)
plot.pxg(combined, mname, pheno.col = 3,  jitter = 1)


?pull.geno
geno_10 <- as.data.frame(pull.geno(combined, chr = "A10"))
str(geno_10)
names(geno_10)

A10_marker <- as.data.frame(geno_10[,"A10x14408017"])
A10_marker

A10_pheno <- as.data.frame(pull.pheno(combined, pheno.col = c(1,2,3)))
A10_pheno
A10_16069633

A10_pheno$A10x14408017 <- A10_marker[,1]
A10_pheno

A10_pheno2 <- A10_pheno[order(A10_pheno$npqi_UN),] 
A10_pheno2

c
?write.table
write.table(A10_pheno2, "A10_geno_pheno.csv", sep = ",", row.names = FALSE)
cor(A10_pheno2$npqi_UN, A10_pheno2$A10x14408017 )

linear <- lm(npqi_UN ~ Bra009312*A10x14408017, data = A10_pheno2)
summary(linear)
str(A10_pheno)
A10_pheno$A10x14408017 <- as.factor(A10_pheno$A10x14408017)

# 2011, 2010, gene expression all in the same RQTL object
setwd("/Users/Cody_2/git.repos/brassica_eqtl_v1.5/data")
#load data
load("eqtl_PAG.RData")



setwd("~/git.repos/brassica_reflectance/")  
npqi2010 <- read.table("npqi2010.csv", sep = ",", header = TRUE)
npqi2011 <- read.table("npqi_2011.csv", sep = ",", header = TRUE)
head(npqi2011)
npqimerge <- merge(npqi2010, npqi2011, all = TRUE )
npqimerge

###bayes####
p450_traits_qb <- brass_total
class(p450_traits_qb)[1] <- "bc"

p450_traits_qb <- qb.genoprob(p450_traits_qb, step=2, stepwidth = "variable")
summary(p450_traits_qb)

summary()

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

p450_traits_qb$pheno$id2 <- as.numeric(sub("(RIL_)(\\d+)", "\\2", p450_traits_qb$pheno$id))
p450_traits_qb$pheno$id2

p450_red <- pull.pheno(p450_traits_qb, pheno.col = c("id", "id2", "Bra009312"))
p450_red
npqimerge

merged_pheno <- merge(p450_red, npqimerge, by.x = "id2", by.y = "Line")
merged_pheno
merged_pheno <- merged_pheno[order(merged_pheno$id),]
merged_pheno$doublcheck <- as.numeric(sub("(RIL_)(\\d+)", "\\2", p450_traits_qb$pheno$id))
merged_pheno

combined <- p450_traits_qb
combined$pheno <- merged_pheno[,c(2:6)]
names(combined$pheno)

str(combined)
combined_qb <- qb.genoprob(combined, step=2, stepwidth = "variable")
summary(combined_qb)


p450_pri <- qb.mcmc(combined_qb, epistasis = TRUE, pheno.col = "Bra009312", seed = 1616,  n.iter = 10000)
plot(p450_pri)
gene_post <- qb.scanone(p450_pri, type.scan = "posterior", epistasis = FALSE, smooth = 2, chr = "A10")
gene_post <- as.data.frame(gene_post)
gene_post
#A10x14408017

gene_HDI <- qb.hpdone(combined_qb_gene, effects = "estimate")
str(gene_HDI)
gene_HDI

UNnpqi2010_qb <- qb.mcmc(combined_qb, epistasis = TRUE, pheno.col = "UNnpqi2010", seed = 1616,  n.iter = 10000)
plot(UNnpqi2010_qb)
un2010_post <- qb.scanone(UNnpqi2010_qb, type.scan = "posterior", epistasis = FALSE, smooth = 2, chr = "A10")
un2010_post <- as.data.frame(un2010_post)
un2010_post
#A10x14408017 

CRnpqi2010_qb <- qb.mcmc(combined_qb, epistasis = TRUE, pheno.col = "CRnpqi2010", seed = 1616,  n.iter = 10000)
plot(CRnpqi2010_qb)
CR2010_post <- qb.scanone(CRnpqi2010_qb, type.scan = "posterior", epistasis = FALSE, smooth = 2, chr = "A10")
CR2010_post <- as.data.frame(CR2010_post)
CR2010_post
#A10x14408017 

npqi2011qb <- qb.mcmc(combined_qb, epistasis = TRUE, pheno.col = "npqi2011", seed = 1616,  n.iter = 10000)
plot(npqi2011qb)
npqi2011_post <- qb.scanone(npqi2011qb, type.scan = "posterior", epistasis = FALSE, smooth = 2, chr = "A10")
npqi2011_post <- as.data.frame(npqi2011_post)
npqi2011_post
#A10x14408017 

geno_10 <- as.data.frame(pull.geno(combined, chr = "A10"))
str(geno_10)
names(geno_10)

A10_marker <- as.data.frame(geno_10[,"A10x14408017"])
A10_marker

A10_pheno <- as.data.frame(pull.pheno(combined, pheno.col = c(1:5)))
A10_pheno


A10_pheno$A10x14408017 <- A10_marker[,1]
A10_pheno

A10_pheno2 <- A10_pheno[order(A10_pheno$npqi2011),] 
A10_pheno2
setwd("~/git.repos/brassica_reflectance/") 

write.table(A10_pheno2, "A10_geno_pheno.csv", sep = ",", row.names = FALSE)
plot()
plot()



