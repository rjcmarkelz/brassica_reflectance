
library(ggplot2)
setwd("~/git.repos/brassica_reflectance/data/")

metabolites <- read.table("k_q_mean_se_all.csv", 
                        header=TRUE, sep = ",")
metabolites




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




