# remake some plots for ASPB 2015 poster
setwd("/Users/Cody_2/git.repos/brassica_reflectance/data/")

UV2010 <- read.table("2010_UV_Wyoming.csv", header = FALSE,  sep = ",",  stringsAsFactors = FALSE)
UV2011 <- read.table("2011_UV_Wyoming.csv", header = FALSE,  sep = ",",  stringsAsFactors = FALSE)

head(UV2010)
head(UV2011)

UV2010$V2 <- as.Date(as.character(UV2010$V2), format = "%Y%m%d")
UV2010$V2

UV2011$V2 <- as.Date(as.character(UV2011$V2), format = "%Y%m%d")
UV2011$V2

plot(UV2010$V2, UV2010$V3)
plot(UV2011$V2, UV2011$V3)


#181 to 211
# subset based on July
july2010 <- UV2010[181:211,2:3]
july2010

july2011 <- UV2011[181:211,2:3]
july2011

plot(july2010$V2, july2010$V3)
plot(july2011$V2, july2011$V3)

UVtotal <- rbind(july2010, july2011)
head(UVtotal)
tail(UVtotal)

UVtotal$Date <- format(UVtotal$V2, "%d") 
head(UVtotal)
UVtotal$Year <- format(UVtotal$V2, "%Y")

UVtotal$date <- factor(UVtotal$Date, levels=unique(UVtotal$Date), ordered=TRUE)

UVtotal

library(ggplot2)
set.seed(12345)
UV <- ggplot(data = UVtotal, mapping=aes(x=date, y=V3)) 
UV <- UV + theme_bw() + facet_grid(facets = Year ~ ., margins = FALSE)
UV <- UV + xlab("July") + ylab("UV Index") + geom_hline(yintercept = 11, color = "purple")
UV <- UV + geom_hline(yintercept = 8, color = "red") 
UV <- UV + scale_y_continuous(limits=c(0, 13)) 
UV <- UV + geom_point() + geom_line(aes(group = 1)) 
UV <- UV + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
	        axis.text.x = element_text(angle = 90))
UV
setwd("/Users/Cody_2/git.repos/brassica_reflectance/output/")
ggsave(file="UV_plot_wyoming.pdf")

UV <- ggplot(data = UVtotal, mapping=aes(x=date, y=V3, group = Year)) 
UV <- UV + theme_bw() #+ facet_grid(facets = Year ~ ., margins = FALSE)
UV <- UV + geom_line(aes(color = Year), size = 2)
UV <- UV + xlab("July") + ylab("UV Index") + geom_hline(yintercept = 11, color = "purple", size = 2)
UV <- UV + geom_hline(yintercept = 8, color = "red", size = 2) 
UV <- UV + scale_y_continuous(limits=c(6, 13))  
UV <- UV + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
	        axis.text.x = element_text(angle = 90, size = 14, face = "bold"),
	        axis.title.x = element_text(size = 16, face = "bold"),
	        axis.text.y = element_text(size = 14, face = "bold"),
	        axis.title.y = element_text(size = 16, face = "bold"))
UV
setwd("/Users/Cody_2/git.repos/brassica_reflectance/output/")
ggsave(file="UV_plot_wyoming.pdf")

