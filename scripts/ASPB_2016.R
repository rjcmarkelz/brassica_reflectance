library(ggplot2)
library(Rmisc)
setwd("~/git.repos/brassica_reflectance/data/")

metabolites <- read.table("k_q_mean_se_all.csv", 
                        header=TRUE, sep = ",")
metabolites

metabolites <- read.table("2015_metabolites_all.csv", 
                        header=TRUE, sep = ",")
metabolites
head(metabolites)

limits <- aes(ymax= k_dic_mean + k_dic_se, ymin = k_dic_mean - k_dic_se)
dodge <- position_dodge(width=0.9)
plotheight5 <- ggplot(data = metabolites, 
                  aes(y=k_dic_mean, x=RIL))
plotheight5 <- plotheight5 + geom_bar(position=dodge, stat="identity") 
plotheight5 <- plotheight5 + geom_errorbar(limits, position = dodge, width = 0.25) 
plotheight5 <- plotheight5 + xlab("RIL") + ylab("Keampferol Abundance")
# plotheight5 <- plotheight5 + ggtitle("Plant Height") + facet_grid(. ~ RIL)
plotheight5 <- plotheight5 + theme(axis.title=element_text(face="bold",
                               size="14"), axis.text=element_text(face="bold",
                               size="10")) + theme_bw()
plotheight5
setwd("~/git.repos/brassica_reflectance/output/")
ggsave(plotheight5, file="Keampferol-diglycoside.png", width=15, height=8)

limits <- aes(ymax= q_dic_mean + q_dic_se, ymin = q_dic_mean - q_dic_se)
dodge <- position_dodge(width=0.9)
plotheight5 <- ggplot(data = metabolites, 
                  aes( y=q_dic_mean, x=RIL))
plotheight5 <- plotheight5 + geom_bar(position=dodge, stat="identity") 
plotheight5 <- plotheight5 + geom_errorbar(limits, position = dodge, width = 0.25) 
plotheight5 <- plotheight5 + xlab("RIL") + ylab("Quercetin Abundance")
# plotheight5 <- plotheight5 + ggtitle("Plant Height") + facet_grid(. ~ RIL)
plotheight5 <- plotheight5 + theme(axis.title=element_text(face="bold",
                               size="14"), axis.text=element_text(face="bold",
                               size="10"))  + theme_bw()
plotheight5
ggsave(plotheight5, file="quercetin-diglycoside.png", width=15, height=8)

head(metabolites)
limits <- aes(ymax= Kaempferol_mean + Kaempferol_se, ymin = Kaempferol_mean - Kaempferol_se)
dodge <- position_dodge(width=0.9)
plotheight5 <- ggplot(data = metabolites, 
                  aes( y=Kaempferol_mean, x=RIL))
plotheight5 <- plotheight5 + geom_bar(position=dodge, stat="identity") 
plotheight5 <- plotheight5 + geom_errorbar(limits, position = dodge, width = 0.25) 
plotheight5 <- plotheight5 + xlab("RIL") + ylab("Keampferol Abundance")
# plotheight5 <- plotheight5 + ggtitle("Plant Height") + facet_grid(. ~ RIL)
plotheight5 <- plotheight5 + theme(axis.title=element_text(face="bold",
                               size="14"), axis.text=element_text(face="bold",
                               size="10"))  + theme_bw()
plotheight5
ggsave(plotheight5, file="keampferol.png", width=15, height=8)

limits <- aes(ymax= Quercetin_methyl_ether_mean + Quercetin_methyl_ether_se, ymin = Quercetin_methyl_ether_mean - Quercetin_methyl_ether_se)
dodge <- position_dodge(width=0.9)
plotheight5 <- ggplot(data = metabolites, 
                  aes( y=Quercetin_methyl_ether_mean, x=RIL))
plotheight5 <- plotheight5 + geom_bar(position=dodge, stat="identity") 
plotheight5 <- plotheight5 + geom_errorbar(limits, position = dodge, width = 0.25) 
plotheight5 <- plotheight5 + xlab("RIL") + ylab("Quercetin methyl ether Abundance")
# plotheight5 <- plotheight5 + ggtitle("Plant Height") + facet_grid(. ~ RIL)
plotheight5 <- plotheight5 + theme(axis.title=element_text(face="bold",
                               size="14"), axis.text=element_text(face="bold",
                               size="10"))  + theme_bw()
plotheight5
ggsave(plotheight5, file="Quercetin methyl ether Abundance.png", width=15, height=8)

head(metabolites)
limits <- aes(ymax= Kaempferol_7_O_glucoside_mean + Kaempferol_7_O_glucoside_se, ymin = Kaempferol_7_O_glucoside_mean - Kaempferol_7_O_glucoside_se)
dodge <- position_dodge(width=0.9)
plotheight5 <- ggplot(data = metabolites, 
                  aes( y=Kaempferol_7_O_glucoside_mean, x=RIL))
plotheight5 <- plotheight5 + geom_bar(position=dodge, stat="identity") 
plotheight5 <- plotheight5 + geom_errorbar(limits, position = dodge, width = 0.25) 
plotheight5 <- plotheight5 + xlab("RIL") + ylab("Kaempferol_7_O_glucoside Abundance")
# plotheight5 <- plotheight5 + ggtitle("Plant Height") + facet_grid(. ~ RIL)
plotheight5 <- plotheight5 + theme(axis.title=element_text(face="bold",
                               size="14"), axis.text=element_text(face="bold",
                               size="10"))  + theme_bw()
plotheight5
ggsave(plotheight5, file="Kaempferol_7_O_glucoside Abundance.png", width=15, height=8)

head(metabolites)
limits <- aes(ymax= Quercetin_methyl_ether_monoglucoside_end_mean + Quercetin_methyl_ether_monoglucoside_end_se, ymin = Quercetin_methyl_ether_monoglucoside_end_mean - Quercetin_methyl_ether_monoglucoside_end_se)
dodge <- position_dodge(width=0.9)
plotheight5 <- ggplot(data = metabolites, 
                  aes( y=Quercetin_methyl_ether_monoglucoside_end_mean, x=RIL))
plotheight5 <- plotheight5 + geom_bar(position=dodge, stat="identity") 
plotheight5 <- plotheight5 + geom_errorbar(limits, position = dodge, width = 0.25) 
plotheight5 <- plotheight5 + xlab("RIL") + ylab("Quercetin_methyl_ether_monoglucoside_end Abundance")
# plotheight5 <- plotheight5 + ggtitle("Plant Height") + facet_grid(. ~ RIL)
plotheight5 <- plotheight5 + theme(axis.title=element_text(face="bold",
                               size="14"), axis.text=element_text(face="bold",
                               size="10"))  + theme_bw()
plotheight5
ggsave(plotheight5, file="Quercetin_methyl_ether_monoglucoside_end Abundance.png", width=15, height=8)

library(tidyr)
library(dplyr)
setwd("~/git.repos/brassica_reflectance/data/")
reflectance <- read.table("Brapa2015_leaf_indices_clean.csv", 
                        header=TRUE, sep = ",")
head(reflectance)


reflect_long <- gather(reflectance, index, reflect, mcari1:nkfi)
head(reflect_long)
reflect_long <- reflect_long[,-c(4:7)]

rfl_sum <- summarySE(reflect_long, measurevar = "reflect", groupvars = c("RIL", "index"))
head(rfl_sum)
rfl_npqi <- subset(rfl_sum, index == "npqi")
rfl_npqi
limits <- aes(ymax= reflect + se, ymin = reflect - se)
dodge <- position_dodge(width=0.9)
reflect_plot <- ggplot(subset(rfl_sum, index == "npqi"), aes( y=reflect, x=RIL))
reflect_plot <- reflect_plot + geom_bar(position=dodge, stat="identity") 
reflect_plot <- reflect_plot + geom_errorbar(limits, position = dodge, width = 0.25) 
reflect_plot <- reflect_plot + xlab("RIL") + ylab("NPQI")
reflect_plot <- reflect_plot + facet_grid(. ~ index)
reflect_plot <- reflect_plot + theme(axis.title=element_text(face="bold",
                               size="14"), axis.text=element_text(face="bold",
                               size="10"))  + theme_bw()
reflect_plot
setwd("~/git.repos/brassica_reflectance/output/")
ggsave(reflect_plot, file="npqi_2015.png", width=15, height=8)


head(reflectance)
out <- split(metabolites, metabolites$Molecule)
head(out)
names(out)
kaemp <- split(metabolites, metabolites$Molecule)[[2]]
kaemp_dic <- split(metabolites, metabolites$Molecule)[[2]]
querc_me <- split(metabolites, metabolites$Molecule)[[6]]
head(kaemp)
head(querc_me)
rfl_sum <- summarySE(querc_me, measurevar = "area", groupvars = c("Genotype"))
head(rfl_sum)

rfl_npqi <- subset(rfl_sum, index == "npqi")
rfl_npqi

dim(kaemp)
dim(reflectance)
kaemp_reflect <- merge(kaemp, reflectance, by.x = "Individual", by.y = "Plant")
querc_reflect <- merge(querc_me, reflectance, by.x = "Individual", by.y = "Plant")
head(keamp_reflect)
head(querc_reflect)
plot(kaemp_reflect$npqi, kaemp_reflect$area)
cor(kaemp_reflect$npqi, kaemp_reflect$area)
plot(querc_reflect$npqi, querc_reflect$area)
cor(querc_reflect$npqi, querc_reflect$area)
?qplot

qplot(data = querc_reflect, npqi, area, color = RIL)
head(querc_me)
rfl_sum <- summarySE(querc_me, measurevar = "area", groupvars = c("Genotype"))
head(rfl_sum)

querc_sum <- summarySE(querc_reflect, measurevar = c("area"), groupvars = c("RIL"))
head(querc_sum)
rfl_sum <- summarySE(querc_reflect, measurevar = c("npqi"), groupvars = c("RIL"))
head(rfl_sum)

out <- merge(querc_sum, rfl_sum, by = "RIL")
head(out)

cor(out$npqi, out$area)
qplot(data = out, npqi, area, color = RIL)
?lm()

dim(querc_reflect)
model <- lm(area ~ RIL + npqi, data = querc_reflect)
summary(model)

expression <- read.table("A10_geno_pheno.csv", 
                        header=TRUE, sep = ",")
expression
str(expression)
plot(expression$Bra009312)
qplot(data = expression, Bra009312 )
ggplot(expression, UNnpqi2010, color = "A10x14408017")
ggplot(data = expression, CRnpqi2010 )

head(expression)
plotheight5 <- ggplot(data = expression, y = npqi2011  )
plotheight5 <- plotheight5 + geom_histogram() 
plotheight5 <- plotheight5 + geom_errorbar(limits, position = dodge, width = 0.25) 
# plotheight5 <- plotheight5 + xlab("RIL") + ylab("Quercetin methyl ether Abundance")
# # plotheight5 <- plotheight5 + ggtitle("Plant Height") + facet_grid(. ~ RIL)
# plotheight5 <- plotheight5 + theme(axis.title=element_text(face="bold",
#                                size="14"), axis.text=element_text(face="bold",
#                                size="10"))  + theme_bw()
plotheight5
