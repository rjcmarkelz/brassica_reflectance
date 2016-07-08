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
querc_reflect
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



#########
head(out)
str(out)
out$area2 <- out$area/100000

npqi_met <- ggplot(data = out, aes(x = npqi, y = area2)) 
npqi_met <- npqi_met + geom_point(aes(size = 1, color = RIL)) + geom_smooth(data = out,aes(x = npqi, y = area2), method = "lm")
npqi_met <- npqi_met + xlab("NPQI")  + ylab("Quercetin Diglycoside") + theme_bw()
npqi_met <- npqi_met + theme(axis.title=element_text(face="bold", size="30"),
                        axis.text=element_text(face="bold",size="25"), legend.position="none")
npqi_met
setwd("~/git.repos/brassica_reflectance/output/")
ggsave("npqi-met.pdf")






#########


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
head(out)

ex_out <- merge(expression, out, by = "RIL" )
head(ex_out)
dim(ex_out)
#########
str(ex_out)
npqi_exp <- ggplot(data = ex_out, aes(x = npqi, y = area2)) 
npqi_exp <- npqi_exp + geom_smooth(data = out,aes(x = npqi, y = area2), method = "lm") + geom_point(aes(color = Bra009312), size = 8) + scale_colour_gradient(low = "black", high = "red", name = "Relative \nExpression")
npqi_exp <- npqi_exp + xlab("NPQI")  + ylab("Quercetin Diglycoside Abundance (10e5)") + theme_bw()
npqi_exp <- npqi_exp + theme(axis.title.x=element_text(face="bold", size="25"), axis.title.y=element_text(face="bold", size="18"),
                        axis.text=element_text(face="bold",size="20"),
                        legend.position = c(0.8,0.80), legend.title = element_text(size = 15, face = "bold"),
                        legend.title.align = 0.0)
npqi_exp
setwd("~/git.repos/brassica_reflectance/output/")
ggsave("npqi-expression.pdf")

?geom_point
p <- ggplot(mtcars, aes(wt, mpg))
p + geom_point(aes(colour = cyl)) + scale_colour_gradient(low = "blue")
p + geom_point(aes(colour = cyl)) + scale_colour_gradient(low = "blue")
p

p <- ggplot(mtcars, aes(wt, mpg))
p + geom_point()
#########

expression$A10x14408017 <- as.factor(expression$A10x14408017)
npqi2011 <- ggplot(data = expression, aes(x = npqi2011))
# npqi2011 <- npqi2011 + geom_histogram(binwidth = 0.01, aes(color = expression$A10x14408017)) 
# npqi2011 <- npqi2011 + geom_density(aes(color = expression$A10x14408017)) 
npqi2011 <- npqi2011 + geom_histogram(binwidth = 0.001)
# npqi2011 <- npqi2011 + geom_errorbar(limits, position = dodge, width = 0.25) 
npqi2011 <- npqi2011 + xlab("NPQI")  + ylab("Count") + theme_bw() + xlim(-0.09, -0.03) + ylim(0, 10)
# # npqi2011 <- npqi2011 + ggtitle("Plant Height") + facet_grid(. ~ RIL)
npqi2011 <- npqi2011 + theme(axis.title=element_text(face="bold", size="30"),axis.text=element_text(face="bold",size="25"))
npqi2011 <- npqi2011  
npqi2011

setwd("~/git.repos/brassica_reflectance/output/")
?ggsave
ggsave("npqi2011-histogram.pdf")

npqi2010 <- ggplot(data = expression, aes(x = UNnpqi2010))
# npqi2010 <- npqi2010 + geom_histogram(binwidth = 0.01, aes(color = expression$A10x14408017)) 
# npqi2010 <- npqi2010 + geom_density(aes(color = expression$A10x14408017)) 
npqi2010 <- npqi2010 + geom_histogram(binwidth = 0.001)
# npqi2010 <- npqi2010 + geom_errorbar(limits, position = dodge, width = 0.25) 
npqi2010 <- npqi2010 + xlab("NPQI") + ylab("Count")
# # npqi2010 <- npqi2010 + ggtitle("Plant Height") + facet_grid(. ~ RIL)
npqi2010 <- npqi2010 + theme(axis.title=element_text(face="bold",
                               size="20"), axis.text=element_text(face="bold",
                               size="20"))  + theme_bw() + xlim(-0.095, -0.04) + ylim(0, 10)
npqi2010


npqi2010 <- ggplot(data = expression, aes(x = UNnpqi2010 ))
# npqi2010 <- npqi2010 + geom_histogram(binwidth = 0.01, aes(color = expression$A10x14408017)) 
npqi2010 <- npqi2010 + geom_density() 

# npqi2010 <- npqi2010 + geom_errorbar(limits, position = dodge, width = 0.25) 
# npqi2010 <- npqi2010 + xlab("RIL") + ylab("Quercetin methyl ether Abundance")
# # npqi2010 <- npqi2010 + ggtitle("Plant Height") + facet_grid(. ~ RIL)
# npqi2010 <- npqi2010 + theme(axis.title=element_text(face="bold",
#                                size="14"), axis.text=element_text(face="bold",
#                                size="10"))  + theme_bw()
npqi2010

Bra009312 <- ggplot(data = expression, aes(x = Bra009312 ))
# Bra009312 <- Bra009312 + geom_histogram(binwidth = 0.01, aes(color = expression$A10x14408017)) 
Bra009312 <- Bra009312 + geom_density(aes(color = expression$A10x14408017)) 
Bra009312 <- Bra009312 + geom_density() 
# Bra009312 <- Bra009312 + geom_errorbar(limits, position = dodge, width = 0.25) 
# Bra009312 <- Bra009312 + xlab("RIL") + ylab("Quercetin methyl ether Abundance")
# # Bra009312 <- Bra009312 + ggtitle("Plant Height") + facet_grid(. ~ RIL)
# Bra009312 <- Bra009312 + theme(axis.title=element_text(face="bold",
#                                size="14"), axis.text=element_text(face="bold",
#                                size="10"))  + theme_bw()
Bra009312

area <- ggplot(data = querc_sum, aes(x = area ))
# area <- area + geom_histogram(binwidth = 0.01, aes(color = expression$A10x14408017)) 
area <- area + geom_density(aes(color = querc_sum$A10x14408017)) 
area <- area + geom_density() 
# area <- area + geom_errorbar(limits, position = dodge, width = 0.25) 
# area <- area + xlab("RIL") + ylab("Quercetin methyl ether Abundance")
# # area <- area + ggtitle("Plant Height") + facet_grid(. ~ RIL)
# area <- area + theme(axis.title=element_text(face="bold",
#                                size="14"), axis.text=element_text(face="bold",
#                                size="10"))  + theme_bw()
area

setwd("~/git.repos/brassica_meta_analysis/Cleaned_data/")

traits <- read.csv("all_traits_clean.csv", header = TRUE)
head(traits)
head(expression)
head(querc_sum)

expression$RIL <- sub("(RIL_)(\\d+)", "\\2", expression$id)

ex_tr <- merge(traits, expression, by.x = "Line", by.y = "id")
head(ex_tr)
ex_tr$RIL <- sub("(RIL_)(\\d+)", "\\2", ex_tr$Line)

ex_tr_q <- merge(ex_tr, out, by = "RIL") 
dim(ex_tr_q)
head(ex_tr_q)

seeds <- ggplot(data = ex_tr_q, aes(x = SLA  ))
# seeds <- seeds + geom_histogram(binwidth = 0.01, aes(color = expression$A10x14408017)) 
seeds <- seeds + geom_density(aes(color = ex_tr_q$A10x14408017)) 


# seeds <- seeds + geom_errorbar(limits, position = dodge, width = 0.25) 
# seeds <- seeds + xlab("RIL") + ylab("Quercetin methyl ether Abundance")
# # seeds <- seeds + ggtitle("Plant Height") + facet_grid(. ~ RIL)
# seeds <- seeds + theme(axis.title=element_text(face="bold",
#                                size="14"), axis.text=element_text(face="bold",
#                                size="10"))  + theme_bw()
seeds

head(ex_tr_q)
out3 <- lm(SLA ~ Bra009312 , data = ex_tr_q)
out4 <- lm(SLA ~ A10x14408017 + Bra009312 + npqi, data = ex_tr_q)
out5 <- lm(SLA ~ A10x14408017, data = ex_tr_q)

summary(out3)
summary(out4)
summary(out5)
anova(out3, out4)
anova(out4, out5)
plot(ex_tr_q$seed_per_plant_UN2012)
ex_tr_q$seed_per_plant_UN2012

cor(ex_tr_q$SLA, ex_tr_q$LMA, na.rm = TRUE)
?cor

# density plot expression
# density plot npqi
plot(ex_tr_q$SLA)