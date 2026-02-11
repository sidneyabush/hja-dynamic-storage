require(vegan)
require(ggplot2)
require(ggrepel)
require(dplyr)
require(tibble)
require(reshape2)

##rda for SiSyn Data - updated with new data on 1/18/22
setwd("/Users/keirajohnson/Box Sync/Keira_Johnson/SiSyn")

new_RDA<-read.csv("RDA_data_080422.csv")

#assign to "trends tot"
trends_tot<-new_RDA

##select conc, yield, and Q - will change for percent vs absolute
trends_tot<-trends_tot[,c(1,2,3,4,6,8)]

#set up dataframe for centered and scaled data
final_center<-trends_tot

site_chars<-new_RDA[,c(1:3,9:27)]

#scale data - this will change if including yield (2:6 vs 2:8)
final_center[c(4:6)]<-data.frame(sapply(final_center[c(4:6)], scale))

#keep only columns to run in RDA
final_trends<-final_center

#unique column for later PCA formatting
site_chars$unique<-paste0(site_chars$LTER, "_", site_chars$site)

#remove rows with missing data
site_chars<-site_chars[complete.cases(site_chars),]

#center and scale
final_center<-site_chars

#center and scale - this will not change depending on yield
final_center[c(4:22)]<-data.frame(sapply(final_center[c(4:22)], scale))

final_chars<-final_center

#run rda - this will change with yield data
my.rda<-rda(final_trends[,c(4:6)], final_chars[,c(4:22)]) #trends = storage metrics, chars = watershed charactericts

#summarize RDA
rda_sum<-summary(my.rda)

#extract data from rda summary - from example code from kjo
st=as.data.frame(rda_sum$sites[,1:2])
sp=as.data.frame(rda_sum$species[,1:2])*2 #alter this scaling factor to align the arrows properly
yz=as.data.frame(rda_sum$biplot[,1:2])

#set up df for loadings plot
rda_loadings<-rda_sum$biplot

rda_loadings<-rda_loadings[rownames(rda_loadings) %in% rownames(yz),]

#melt
rda_loadings_melt<-melt(rda_loadings)
names(rda_loadings_melt)<-c("Variable", "RDA_axis", "Eigenvalue")
#only keep RDA 1 and 2
rda_loadings_melt<-subset(rda_loadings_melt, rda_loadings_melt$RDA_axis=="RDA1"|rda_loadings_melt$RDA_axis=="RDA2")

ggplot() +
  geom_point(data = st,aes(RDA1,RDA2),size=4)+ #can add color variable here using call col="XX" inside aes()) function
  theme_bw()+
  theme(legend.position = "right")+
  geom_segment(data = sp,aes(x = 0, y = 0, xend = RDA1, yend = RDA2), 
               arrow = arrow(angle=22.5,length = unit(0.35,"cm"),
                             type = "closed"),linetype=1, size=0.6,colour = "grey")+
  geom_text_repel(data = sp,aes(RDA1,RDA2,label=row.names(sp)),colour="grey")+
  geom_segment(data = yz,aes(x = 0, y = 0, xend = RDA1, yend = RDA2), 
               arrow = arrow(angle=22.5,length = unit(0.35,"cm"),
                             type = "closed"), size=0.6)+
  geom_text_repel(data = yz,aes(RDA1,RDA2,label=rownames(yz)))
