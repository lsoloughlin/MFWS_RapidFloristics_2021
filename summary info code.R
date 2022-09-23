a<-dat88 %>% select(Area, Exp.Treat, Average.Grass.Height)%>%
  group_by(Area, Exp.Treat) %>% 
  summarise(mean = mean(Average.Grass.Height),
            n = n(),
            sd = sd(Average.Grass.Height))%>%
  mutate(se = sd/n,
         ci_l=mean-(se*1.96),
         ci_u=mean+(se*1.96))

#Make a species occurrence matrix
b<-dat %>%
  filter(Area != "MULL_OUT")%>%
  group_by(Area)%>%
  summarise(Arthropodium = sum(Arth.spp)/length(Area),
            Bulbine = sum(Bulb.spp)/length(Area),
            Wurmbea = sum(Wurm.dio)/length(Area),
            Burchardia = sum(Burch.umb)/length(Area),
            Microtis = sum(Microt.spp)/length(Area),
            Thelymitra = sum(Thely.spp)/length(Area),
            Diuris = sum(Diur.spp)/length(Area),
            Hymenochilus = sum(Hymen.spp)/length(Area),
            Microseris = sum(Micros.spp)/length(Area),
            Craspedia = sum(Cras.spp)/length(Area),
            Leptorynchos = sum(Lept.squ)/length(Area),
            Chrysocephalum = sum(Chry.api)/length(Area),
            Stackhousia = sum(Stac.mon)/length(Area))%>%
  pivot_longer(!Area, names_to = "Species", values_to = "Occurrence")
b$Species <- as.factor(b$Species) 
b <- as.data.frame(b)
b <- b%>%mutate(Species = fct_relevel(Species, "Stackhousia", "Chrysocephalum", "Leptorynchos",
                            "Craspedia", "Microseris", "Hymenochilus", "Diuris",
                            "Thelymitra", "Microtis", "Burchardia", "Bulbine",
                            "Wurmbea", "Arthopodium"))
Occ.Plot <- b%>%
           ggplot(aes(Area))+
           geom_bin2d(aes(y = Species, fill = Occurrence))+
           theme_bw()+
           scale_fill_continuous(high = 'dark blue', low = 'white')+
           theme(axis.text.x = element_text(angle = 45,hjust=0.95, colour='black', size = 14))+
           theme(axis.text.y = element_text(colour='black'))+
           scale_x_discrete(labels = c( "Mulligans Flat", "Goorooyarroo", "Offsets"))
          
Occ.Plot


ggsave("OccSppMatrix1.tiff",
       plot=Occ.Plot,
       width = 5.5, height = 5.5, units = "in",
       dpi = 300)


#Make a species occurrence matrix without exotic dom sites
c<-dat %>%
  filter(Area != "MULL_OUT")%>%
  filter(PCT.STATE2 != "Exotic Woodland",
         PCT.STATE2 != "Forest")%>%
  group_by(Area)%>%
  summarise(Arthropodium = sum(Arth.spp)/length(Area),
            Bulbine = sum(Bulb.spp)/length(Area),
            Wurmbea = sum(Wurm.dio)/length(Area),
            Burchardia = sum(Burch.umb)/length(Area),
            Microtis = sum(Microt.spp)/length(Area),
            Thelymitra = sum(Thely.spp)/length(Area),
            Diuris = sum(Diur.spp)/length(Area),
            Hymenochilus = sum(Hymen.spp)/length(Area),
            Microseris = sum(Micros.spp)/length(Area),
            Craspedia = sum(Cras.spp)/length(Area),
            Leptorynchos = sum(Lept.squ)/length(Area),
            Chrysocephalum = sum(Chry.api)/length(Area),
            Stackhousia = sum(Stac.mon)/length(Area))%>%
  pivot_longer(!Area, names_to = "Species", values_to = "Occurrence")
c$Species <- as.factor(c$Species) 
c <- as.data.frame(c)
c <- c%>%mutate(Species = fct_relevel(Species, "Stackhousia", "Chrysocephalum", "Leptorynchos",
                                      "Craspedia", "Microseris", "Hymenochilus", "Diuris",
                                      "Thelymitra", "Microtis", "Burchardia", "Bulbine",
                                      "Wurmbea", "Arthopodium"))
Occ.Plot2 <- c%>%
  ggplot(aes(Area))+
  geom_bin2d(aes(y = Species, fill = Occurrence))+
  theme_bw()+
  scale_fill_continuous(high = 'darkblue', low = 'white')+
  theme(axis.text.x = element_text(angle = 45,hjust=0.95, colour='black', size = 14))+
  theme(axis.text.y = element_text(colour='black'))+
  scale_x_discrete(labels = c( "Mulligans Flat", "Goorooyarroo", "Offsets"))

Occ.Plot2


ggsave("OccSppMatrix2.tiff",
       plot=Occ.Plot2,
       width = 5.5, height = 5.5, units = "in",
       dpi = 300)

#Make a species occurrence matrix without exotic dom sites or Bettong FT in MULL
d<-dat %>%
  filter(Area != "MULL_OUT")%>%
  filter(PCT.STATE2 != "Exotic Woodland",
         PCT.STATE2 != "Forest",
         Exp.Treat != "Bettong fence FT")%>%
  group_by(Area)%>%
  summarise(Arthropodium = sum(Arth.spp)/length(Area),
            Bulbine = sum(Bulb.spp)/length(Area),
            Wurmbea = sum(Wurm.dio)/length(Area),
            Burchardia = sum(Burch.umb)/length(Area),
            Microtis = sum(Microt.spp)/length(Area),
            Thelymitra = sum(Thely.spp)/length(Area),
            Diuris = sum(Diur.spp)/length(Area),
            Hymenochilus = sum(Hymen.spp)/length(Area),
            Microseris = sum(Micros.spp)/length(Area),
            Craspedia = sum(Cras.spp)/length(Area),
            Leptorynchos = sum(Lept.squ)/length(Area),
            Chrysocephalum = sum(Chry.api)/length(Area),
            Stackhousia = sum(Stac.mon)/length(Area))%>%
  pivot_longer(!Area, names_to = "Species", values_to = "Occurrence")
d$Species <- as.factor(d$Species) 
d <- as.data.frame(d)
d <- d%>%mutate(Species = fct_relevel(Species, "Stackhousia", "Chrysocephalum", "Leptorynchos",
                                      "Craspedia", "Microseris", "Hymenochilus", "Diuris",
                                      "Thelymitra", "Microtis", "Burchardia", "Bulbine",
                                      "Wurmbea", "Arthopodium"))
Occ.Plot3 <- d%>%
  ggplot(aes(Area))+
  geom_bin2d(aes(y = Species, fill = Occurrence))+
  theme_bw()+
  scale_fill_continuous(high = 'darkblue', low = 'white')+
  theme(axis.text.x = element_text(angle = 45,hjust=0.95, colour='black', size = 14))+
  theme(axis.text.y = element_text(colour='black'))+
  scale_x_discrete(labels = c( "Mulligans Flat", "Goorooyarroo", "Offsets"))

Occ.Plot3


ggsave("OccSppMatrix3.tiff",
       plot=Occ.Plot3,
       width = 5.5, height = 5.5, units = "in",
       dpi = 300)
c<-dat %>% 
  count(DomGrass.sum)
        