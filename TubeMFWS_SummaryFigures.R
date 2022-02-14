# Rapid Floristic Assessment of the MFWS 2021 Project

# Luke O'Loughlin
# date started: 2021-11-22

#Script 1: Summary raw data figures

library(dplyr)
library(ggplot2)
library(readxl)
library(ggpubr)

# 1. READ DATA ####

inputDir <- "J:/Policy/Nature Conservation Policy/CPR/Govt/Initiatives and Projects/Mulligans Flat Woodlands Sanctuary/Monitoring Plan continued - Luke OLoughlin/Bulb_rapid_monitoring/"
dat <- read_excel(paste0(inputDir, "TuberousPlants_MFWS_Data.xlsx"))
str(dat)
View(dat)

# 2. INDICATOR RICHNESS ####

#Spatial XY figure for "Indicator Richness"


dat1 <- dat %>% mutate(INDI.RICH = as.factor(INDI.RICH))
BT <- dat1 %>% filter(Exp.Treat %in% c("Bettong fence LG", "Bettong fence HG", "Bettong fence FT"))
plot1<- dat %>% 
  ggplot(aes(x=x))+
  geom_point(aes(y=y, colour = Area, size = INDI.RICH))+
  theme_bw()+
  ylab("Latitude")+
  xlab("Longitude")+
  labs(size = "Indicator richness")
  scale_size_manual(values = c(1,3,6,9,12,15))
plot1

ggsave("INDIRICH.tiff",
       plot=plot1,
       width = 7.4, height = 5.7, units = "in",
       dpi = 300)  


# 3. SPECIES RICHNESS ####

dat2 <- dat %>% mutate(SPP.RICH = as.factor(SPP.RICH))
plot2<- dat2 %>% 
  ggplot(aes(x=x))+
  geom_point(aes(y=y, colour = Area, size = SPP.RICH))+
  theme_bw()+
  ylab("Latitude")+
  xlab("Longitude")+
  labs(size = "Species richness")+
  scale_size_manual(values = c(1,2,4,6,8,10,12,14,16))
plot2

ggsave("SPPRICH.tiff",
       plot=plot2,
       width = 7.4, height = 5.7, units = "in",
       dpi = 300)  



# 4. INDICATOR ABUNDANCE ####



dat3 <- dat %>% mutate(MAX.ABUND = as.factor(MAX.ABUND))
plot3<- dat3 %>% 
  ggplot(aes(x=x))+
  geom_point(aes(y=y, colour = Area, size = MAX.ABUND))+
  theme_bw()+
  ylab("Latitude")+
  xlab("Longitude")+
  labs(size = "Indicator abundance")+
  scale_size_manual(values = c(1,3,6,9,12))
plot3

ggsave("INDIABUND.tiff",
       plot=plot3,
       width = 7.4, height = 5.7, units = "in",
       dpi = 300) 

dat4 <- dat %>% mutate(Lilies = as.factor(Lilies))
plot4<- dat4 %>% 
  ggplot(aes(x=x))+
  geom_point(aes(y=y, colour = Area, size = Lilies))+
  theme_bw()+
  theme(axis.title = element_blank())+
  theme(axis.text = element_blank())+
  theme(legend.position = "none")+
  ggtitle("Lily abundance")+
  scale_size_manual(values = c(1,3,6,9,12))
plot4

dat5 <- dat %>% mutate(Orchids = as.factor(Orchids))
plot5<- dat5 %>% 
  ggplot(aes(x=x))+
  geom_point(aes(y=y, colour = Area, size = Orchids))+
  theme_bw()+
  theme(axis.title = element_blank())+
  theme(axis.text = element_blank())+
  theme(legend.position = "none")+
  ggtitle("Orchid abundance")+
  scale_size_manual(values = c(1,3,6,9,12))
plot5

dat6 <- dat %>% mutate(Daisies.T = as.factor(Daisies.T))
plot6<- dat6 %>% 
  ggplot(aes(x=x))+
  geom_point(aes(y=y, colour = Area, size = Daisies.T))+
  theme_bw()+
  theme(axis.title = element_blank())+
  theme(axis.text = element_blank())+
  theme(legend.position = "none")+
  ggtitle("Daisies (tuberous) abundance")+
  scale_size_manual(values = c(1,3,6,9,12))
plot6

dat7 <- dat %>% mutate(Daisies.C = as.factor(Daisies.C))
plot7<- dat7 %>% 
  ggplot(aes(x=x))+
  geom_point(aes(y=y, colour = Area, size = Daisies.C))+
  theme_bw()+
  theme(axis.title = element_blank())+
  theme(axis.text = element_blank())+
  theme(legend.position = "none")+
  ggtitle("Daisies (control) abundance")+
  scale_size_manual(values = c(1,3,6,9,12))
plot7

INDI.FIG <- ggarrange(plot4, plot5, plot6, plot7,
                     ncol = 2,  nrow = 2)
INDI.FIG

ggsave("IndiALL.tiff",
       plot=INDI.FIG,
       width = 6, height = 6, units = "in",
       dpi = 300)  

#### With Bettong Fences Labelled

datBT <- dat %>% mutate(AreaX = ifelse(Exp.Treat  %in% c("Bettong fence LG", "Bettong fence HG", "Bettong fence FT"),
                                       "MULL_BettongFence", Area))

datBT4 <- datBT %>% mutate(Lilies = as.factor(Lilies))
plotBT4<- datBT4 %>% 
  ggplot(aes(x=x))+
  geom_point(aes(y=y, colour = AreaX, size = Lilies))+
  theme_bw()+
  theme(axis.title = element_blank())+
  theme(axis.text = element_blank())+
  scale_colour_discrete(labels = c("Mulligans \n(inside)", "Goorooyarroo", "Mulligans \n(outside)", "Offsets", "Mulligans \n(bettong fence)"))+
  labs(colour = "Area", size = "Count")+
  ggtitle("Lilies")+
  theme(plot.title = element_text(size = 10, face = "bold"))+
  scale_size_manual(values = c(1,3,6,9,12),
                    labels = c("0", "<10", "<100", "<1000", "1000+"))
plotBT4

datBT5 <- datBT %>% mutate(Orchids = as.factor(Orchids))
plotBT5<- datBT5 %>% 
  ggplot(aes(x=x))+
  geom_point(aes(y=y, colour = AreaX, size = Orchids))+
  theme_bw()+
  theme(axis.title = element_blank())+
  theme(axis.text = element_blank())+
  scale_colour_discrete(labels = c("Mulligans \n(inside)", "Goorooyarroo", "Mulligans \n(outside)", "Offsets", "Mulligans \n(bettong fence)"))+
  labs(colour = "Area", size = "Count")+
  ggtitle("Orchids")+
  theme(plot.title = element_text(size = 10, face = "bold"))+
  scale_size_manual(values = c(1,3,6,9,12),
                    labels = c("0", "<10", "<100", "<1000", "1000+"))
plotBT5

datBT6 <- datBT %>% mutate(Daisies.T = as.factor(Daisies.T))
plotBT6<- datBT6 %>% 
  ggplot(aes(x=x))+
  geom_point(aes(y=y, colour = AreaX, size = Daisies.T))+
  theme_bw()+
  theme(axis.title = element_blank())+
  theme(axis.text = element_blank())+
  scale_colour_discrete(labels = c("Mulligans \n(inside)", "Goorooyarroo", "Mulligans \n(outside)", "Offsets", "Mulligans \n(bettong fence)"))+
  labs(colour = "Area", size = "Count")+
  ggtitle("Daisies (tuberous)")+
  theme(plot.title = element_text(size = 10, face = "bold"))+
  scale_size_manual(values = c(1,3,6,9,12),
                    labels = c("0", "<10", "<100", "<1000", "1000+"))
plotBT6

datBT7 <- datBT %>% mutate(Daisies.C = as.factor(Daisies.C))
plotBT7<- datBT7 %>% 
  ggplot(aes(x=x))+
  geom_point(aes(y=y, colour = AreaX, size = Daisies.C))+
  theme_bw()+
  theme(axis.title = element_blank())+
  theme(axis.text = element_blank())+
  scale_colour_discrete(labels = c("Mulligans \n(inside)", "Goorooyarroo", "Mulligans \n(outside)", "Offsets", "Mulligans \n(bettong fence)"))+
  labs(colour = "Area", size = "Count")+
  ggtitle("Daisies (control)")+
  theme(plot.title = element_text(size = 10, face = "bold"))+
  scale_size_manual(values = c(1,3,6,9,12),
                    labels = c("0", "<10", "<100", "<1000", "1000+"))
plotBT7

INDI.FIG.BT <- ggarrange(plotBT4, plotBT5, plotBT6, plotBT7,
                      ncol = 2,  nrow = 2,
                      common.legend = TRUE, legend = "left")
INDI.FIG.BT

ggsave("IndiALLBT.tiff",
       plot=INDI.FIG.BT,
       width = 9, height = 6.5, units = "in",
       dpi = 300)  


# 5. SPECIES PRESENCE ####

plot8<- dat %>% 
  ggplot(aes(x=x))+
  geom_point(aes(y=y, colour = Area, size = Arth.spp))+
  theme_bw()+
  theme(axis.title = element_blank())+
  theme(axis.text = element_blank())+
  ggtitle("Arthropodium spp.")+
  theme(legend.position = "none")
plot8

plot9<- dat %>% 
  ggplot(aes(x=x))+
  geom_point(aes(y=y, colour = Area, size = Bulb.spp))+
  theme_bw()+
  theme(axis.title = element_blank())+
  theme(axis.text = element_blank())+
  ggtitle("Bulbine spp.")+
  theme(legend.position = "none")
plot9

plot10<- dat %>% 
  ggplot(aes(x=x))+
  geom_point(aes(y=y, colour = Area, size = Wurm.dio))+
  theme_bw()+
  theme(axis.title = element_blank())+
  theme(axis.text = element_blank())+
  ggtitle("Wurmbia dioicia")+
  theme(legend.position = "none")
plot10

plot11<- dat %>% 
  ggplot(aes(x=x))+
  geom_point(aes(y=y, colour = Area, size = Microt.spp))+
  theme_bw()+
  theme(axis.title = element_blank())+
  theme(axis.text = element_blank())+
  ggtitle("Microtis spp.")+
  theme(legend.position = "none")
plot11

plot12<- dat %>% 
  ggplot(aes(x=x))+
  geom_point(aes(y=y, colour = Area, size = Thely.spp))+
  theme_bw()+
  theme(axis.title = element_blank())+
  theme(axis.text = element_blank())+
  ggtitle("Thelymitra spp.")+
  theme(legend.position = "none")
plot12

plot13<- dat %>% 
  ggplot(aes(x=x))+
  geom_point(aes(y=y, colour = Area, size = Diur.spp))+
  theme_bw()+
  theme(axis.title = element_blank())+
  theme(axis.text = element_blank())+
  ggtitle("Diuris spp.")+
  theme(legend.position = "none")
plot13

plot14<- dat %>% 
  ggplot(aes(x=x))+
  geom_point(aes(y=y, colour = Area, size = Micros.spp))+
  theme_bw()+
  theme(axis.title = element_blank())+
  theme(axis.text = element_blank())+
  ggtitle("Microseris spp.")+
  theme(legend.position = "none")
plot14

plot15<- dat %>% 
  ggplot(aes(x=x))+
  geom_point(aes(y=y, colour = Area, size = Cras.spp))+
  theme_bw()+
  theme(axis.title = element_blank())+
  theme(axis.text = element_blank())+
  ggtitle("Craspedia spp.")+
  theme(legend.position = "none")
plot15

plot16<- dat %>% 
  ggplot(aes(x=x))+
  geom_point(aes(y=y, colour = Area, size = Lept.squ))+
  theme_bw()+
  theme(axis.title = element_blank())+
  theme(axis.text = element_blank())+
  ggtitle("Leptorhynchos squamatus")+
  theme(legend.position = "none")
plot16

plot17<- dat %>% 
  ggplot(aes(x=x))+
  geom_point(aes(y=y, colour = Area, size = Chry.api))+
  theme_bw()+
  theme(axis.title = element_blank())+
  theme(axis.text = element_blank())+
  ggtitle("Chrysocephalum apiculatum")+
  theme(legend.position = "none")
plot17

#Join all spp. figs together
SPP.FIG <- ggarrange(plot8, plot9, plot10,
                     plot11, plot12, plot13,
                     plot14, plot15, plot16,
                     plot17,
                     ncol = 3,  nrow = 4)
SPP.FIG

ggsave("SppALL.tiff",
       plot=SPP.FIG,
       width = 5.4, height = 8, units = "in",
       dpi = 300)  
