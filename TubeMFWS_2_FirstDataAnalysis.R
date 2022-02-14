# Rapid Floristic Assessment of the MFWS 2021 Project

# Luke O'Loughlin
# dates written: 2022-02-01 --> 2022-02-14


#Script 2: First data analysis

# 1. PREPARATION_________________________________####

## total clean up 
rm(list = ls()) 
invisible(gc())

## Libraries 
# for this whole script
library(dplyr)
library(ggplot2)
library(tidyr)
library(glmmTMB)
library(DHARMa)
library(ggpubr)
library(performance)
library(MuMIn)
library(readr)
library(broom)
library(stringr)
library(DataExplorer)
library(corrplot)
library(forcats)

#Read in (1) monitoring data that includes Veg Comm and PCT info, and 
# (2) the 2021 veg structure monitoring data from MFWS
dat1 <- read_csv("MFWS-RapidFloristicData-2021.csv", 
                                         col_types = cols(Date = col_date(format = "%Y-%m-%d"), 
                                         Start.time = col_time(format = "%H:%M"), 
                                         End.time = col_time(format = "%H:%M")))
str(dat1)

dat2 <- read_csv("MFWS-VegStructData-2021.csv")
colnames(dat2) <- make.names(colnames(dat2), unique=TRUE) #fix col names so they dont have spaces
str(dat2)

#---#
#Edit the veg struct data frame
#---#

#Remove CR veg team observations
dat2 <- dat2 %>%
  filter(Observer != "Felicity Grant",
         Observer != "Emma Cook")

#calculate proportional cover (0-1) of variables from count/total steps
dat2 <- dat2 %>%
  mutate(NatGrass = (C4.Native+C4.Native..Themeda.+
                       C3.Native..Other.+C3.Native..Rytidosperma.sp..+
                       C3.Native..Stipa.grass.) / Total,
         NatShrub = Native.Ground..Shrubs. / Total,
         NatOther = Native.Ground..Other. / Total,
         NatTotal = NatGrass + NatShrub + NatOther,
         Rock = Rock.Cover / Total,
         BareGround = Bare.Ground.Cover / Total,
         Litter = (Leaf.Litter.Cover + Dead.Plant) / Total,
         Thatch = Thatch.Cover / Total,
         ExBroadleaf = (Exotic.Broadleaf..Annual.+
                          Exotic.Broadleaf..Prennial.) / Total,
         ExGrass = (Wild.Oats+Other.Annual.Grasses+Phalaris+
                      Other.Exotic.Perennial.Grass+
                      Chilean.Needle.Grass+Love.Grass+
                      Serrated.Tussock)/Total,
         ExTotal = ExBroadleaf + ExGrass + (Exotic.Cover..Clover./Total),
         ExPTotal = ExTotal - ((Wild.Oats + Other.Annual.Grasses+
                       Exotic.Broadleaf..Annual.)/Total), #Exotic perennial species only
         Nativeness = NatTotal / (NatTotal + ExPTotal) #Nativeness of perennial veg
  )

#select only cols of interest to join
dat3 <- dat2 %>%
  select(Plot.ID,#location identifier as in other dataframe
         NatGrass,NatShrub,NatOther,NatTotal,#Native veg cover
         Rock,BareGround,Litter,Thatch,#biophysical cover
         ExBroadleaf,ExGrass,ExTotal,ExPTotal,#Exotic veg cover
         Nativeness,Dom.Grass,Average.Grass.Height#Key structure indicators
         )%>%
  rename(Plot = Plot.ID)

#Join dat3 (structure of interest) to dat1 (monitoring data)
dat <- left_join(dat1, dat3, by = "Plot") #dat = master data for analysis



# 2. EXPLORE EXPLANATORY DATA____________________####

#---#
# Replications of Stratification Factors
#---#

# PCT Condition Zones
dat <- dat %>% mutate(PCT.ZONE = as.factor(PCT_ZONE))
a <- dat %>% group_by(Area) %>% count(PCT.ZONE)
dat <- dat %>% mutate(PCT.STATE = if_else(PCT.ZONE == "16.1" | #Group Zones into 4 broader states
                                            PCT.ZONE == "16.2" |
                                            PCT.ZONE == "16.3" |
                                            PCT.ZONE == "16.4", "Woodland",
                                          if_else(PCT.ZONE == "16.5" |
                                                    PCT.ZONE == "16.6" |
                                                    PCT.ZONE == "16.7" |
                                                    PCT.ZONE == "16.8", "Derived Grassland",
                                                  if_else(PCT.ZONE == "25.1" |
                                                            PCT.ZONE == "25.2"|
                                                            PCT.ZONE == "25.3"|
                                                            PCT.ZONE == "25.8"|
                                                            PCT.ZONE == "99", "Forest", 
                                                            "Exotic Woodland")))) #Exotic Woodland and Exotic Derived Grassland
                                                                                  #lumped together as there is only a little of the latter

a <- dat %>% group_by(Area) %>% count(PCT.STATE)
#                 MF    MF_OUT    GOO   OFF
#Woodland         34      4       18    11      #Mature Tree, native dominate ground
#Derived Grass    10      2       3     8       #No trees, native dominate ground
#Exotic Wood      2       -       -     11      #Maybe trees, exotic dominate ground
#Forest           -       -       9     -       #Mature trees not YB and/or BRRG

a <- dat %>% select(Nativeness, PCT.STATE) #look for missmatch b/w PCT mapped and observed
#10 "woodlands" and 2 "derived grasslands" are less than 50% nativeness
dat$PCT.STATE <- as.character(dat$PCT.STATE)
dat <- dat %>%
  mutate(nat2 = replace_na(Nativeness, 1),
         PCT.STATE2 = ifelse(nat2 < 0.5 & PCT.STATE == "Woodland", "Exotic Woodland",
                             ifelse(nat2 < 0.5 & PCT.STATE == "Derived Grassland", "Exotic Woodland", 
                                    PCT.STATE)))
a <- dat %>% select(nat2, PCT.STATE, PCT.STATE2)
a <- dat %>% group_by(Area) %>% count(PCT.STATE2)
#                 MF    MF_OUT    GOO   OFF     #MULL didnt change, 9 OFF and 4 GOO
#Woodland         34      4       14    4       #Mature Tree, native dominate ground
#Derived Grass    10      2       3     6       #No trees, native dominate ground
#Exotic Wood      2       -       4     20      #Maybe trees, exotic dominate ground
#Forest           -       -       9     -       #Mature trees not YB and/or BRRG

# Vegetation Communities
a <- dat %>% group_by(Area) %>% count(VEG_COMM)
#                 MF    MF_OUT    GOO   OFF
#u19 BGGGW        28      4       14    15      #Yellow box / Blakleyi
#u178 YBAB        3       -       -     -       #Yellow box / Applebox
#p14 DryScl For   -       -       3     -       #Dry schlorphyll grassy forest
#NG               9       2       6     14      #"Native grassland"
#DNS              1       -       -     -       #"Derived native scrubland"
#EPN              5       -       7     -       #"Environmental Planting"

# Grazing treatment
a <- dat %>% group_by(Area) %>% count(Exp.Treat)
#                 MF    MF_OUT    GOO   OFF
#H Grazing        12      -       20    30      #"High"
#L Grazing        10      -       10    -       #"Low"
#BettongF HG      12      -       -     -       # Fence in a High G context
#BettongF LG      8       -       -     -       # Fence in a Low G context
#BettongF FT      4       -       -     -       # Bettong exclusion

# C3 or C4 grass mapped
a <- dat %>% group_by(Area) %>% count(C3_C4)
#                 MF    MF_OUT    GOO   OFF
#Native C3         8      -       15    17      #"High"
#Native C4        38      6       14    2       #"Low"
#Exotic C3         -      -       1     11       # Fence in a High G context

# dominant grass observed 
a <- dat %>% group_by(Area) %>% count(Dom.Grass)
a<-dat%>%distinct(Dom.Grass)
dat <- dat %>% mutate(DomGrass.sum = if_else(Dom.Grass == "Themeda australis", "Themeda", #Group exotics and natives except T,R,A
                                             if_else(Dom.Grass == "Rytidosperma spp.", "Rytid",
                                                     if_else(str_starts(Dom.Grass, "Aust"), "Stipa",
                                                             if_else(str_starts(Dom.Grass, "Both") |
                                                                       str_starts(Dom.Grass, "Mic")|
                                                                       str_starts(Dom.Grass, "Ari")|
                                                                       str_starts(Dom.Grass, "Ann")|
                                                                       str_starts(Dom.Grass, "Era"), "Other native", "Exotic"
                                                             )))))
                                                                    
                                                             if_else(str_starts(Dom.Grass, "Aust"), "Stipa", "Other native"
                                                                     )))))

a <- dat %>% group_by(Area) %>% count(DomGrass.sum)
#                 MF    MF_OUT    GOO   OFF
#Themeda          8       3        9     1       #Themeda australis
#Rytid            5       1        1     1       #Rytidosperma spp.
#Stipa            -       -        7     8       #multiple Austrostipa species
#Other native     5       1        6     4       #5 other natives
#Exotic           4       1        7     16      #8 exotics




#---#
# Scale numerical variables
#---#

dat <- dat %>%
  mutate(zCanopy.Cov  = scale(Canopy.Cov),
         zUnderS.Cov  = scale(UnderS.Cov),
         zNatGrass    = scale(NatGrass),
         zNatShrub    = scale(NatShrub),
         zNatOther    = scale(NatOther),
         zNatTotal    = scale(NatTotal),
         zRock        = scale(Rock),
         zBareGround  = scale(BareGround),
         zLitter      = scale(Litter),
         zThatch      = scale(Thatch),
         zExBroadleaf = scale(ExBroadleaf),
         zExGrass     = scale(ExGrass),
         zExTotal     = scale(ExTotal),
         zExPTotal    = scale(ExPTotal), #Exotic perennial species only
         zNativeness  = scale(Nativeness),
         zAvGrassHt   = scale(Average.Grass.Height)#Nativeness of perennial veg
  )

#---#
# Check correlation among full list of potential explanatory variables
#---#

dat.corr <- dat %>% select(zCanopy.Cov, zUnderS.Cov, zNatGrass, zNatShrub,
                           zNatOther, zNatTotal, zRock, zBareGround, zLitter,
                           zThatch, zExBroadleaf, zExGrass, zExTotal, zExPTotal,
                           zNativeness, zAvGrassHt)

cor_matrix<-cor(dat.corr, method = "spearman", use="complete.obs") #use spearmans to assess both linear and monotonic relationships
corrplot(cor_matrix,
         method = "number",
         type = "upper",
)
# using r < 0.7 as cut off
#
# NatGrass (+0.81 w NatTotal) & (-0.73 w ExTotal)
# NatOther (+0.71 w NatTotal)
# NatTotal (-0.70 w ExGrass) & (-0.90 e ExTotal) & (+0.76 w Nativeness)
# ExTotal (+0.81 w ExPTotal) & (-0.86 w Nativeness)
#
# retain: Nativeness, NatGrass, NatOther, ExGrass
# drop:   NatTotal, ExTotal, ExPTotal
#


dat.corr2 <- dat %>% select(zCanopy.Cov, zUnderS.Cov, zNatGrass, zNatShrub,
                           zNatOther, zRock, zBareGround, zLitter,
                           zThatch, zExBroadleaf, zExGrass,
                           zNativeness, zAvGrassHt)

cor_matrix2<-cor(dat.corr2, method = "spearman", use="complete.obs") #use spearmans to assess both linear and monotonic relationships
corrplot(cor_matrix2,
         method = "number",
         type = "upper",
)

# using r < 0.6 as cut off
#
# Canopy.Cov (+0.63 w Litter)
# NatGrass (-0.64 w ExGrass)
# NatOther (0.60 w Nativeness)
# Nativeness (-0.65 w AvGrHieght)

# 3. FULL MODEL SELECTION________________________####

#reorder factors in the dataset
dat<-dat%>%
  mutate(Area = fct_relevel(Area, "MULL"),
         PCT.STATE = fct_relevel(PCT.STATE, "Woodland"),
         UMC_ID = fct_relevel(UMC_ID, "u19"),
         Exp.Treat = fct_relevel(Exp.Treat, "High grazing"),
         C3_C4 = fct_relevel(C3_C4, "Native_C4"),
         DomGrass.sum = fct_relevel(DomGrass.sum, "Themeda"),
         PCT.STATE2 = fct_relevel(PCT.STATE, "Woodland"))


# 3.1  INDICATOR RICHNESS________________________####

# “Indicator richness” is the count of indicator groups present at a site. 
# There were 5 groups (Lilies, Orchids, Daisies [tuberous], Daisies [control], and Stachhousia).

ggplot(dat, aes(INDI.RICH))+  #check distribution of response variable
  geom_histogram()            #All should poisson (b/c it's count data)
                              #unless there's overdispersion

####
# Modelling
####

# Full model

ir.f <- glmmTMB(INDI.RICH ~ Area + 
                PCT.STATE2 +
                UMC_ID +
                Exp.Treat +
                C3_C4 +
                DomGrass.sum +
                zCanopy.Cov + zUnderS.Cov + zNatGrass + zNatShrub+
                zNatOther+ zRock+ zBareGround+ zLitter+
                zThatch+ zExBroadleaf+ zExGrass+
                zNativeness+ zAvGrassHt +
                  (1|Polygon),
                data = dat,
                family = poisson)


summary(ir.f)


#check model
simulationOutput <- simulateResiduals(fittedModel = ir.f, n = 250)
plotSimulatedResiduals(simulationOutput = simulationOutput) #Overdispersed. 
testUniformity(simulationOutput = simulationOutput)
testDispersion(ir.f)
testZeroInflation(simulationOutput)
testSpatialAutocorrelation(simulationOutput = simulationOutput, x = dat$x, y= dat$y)


#Model Selection  on Full model
options(na.action = "na.fail")
#ir.d <- dredge(ir.f) #Dont run again unless global model changes
ir.d
ir.sub<-subset(ir.d, delta<=2,  recalc.weights = FALSE)
ir.sub #What are the supported (delta AICc <2) models?
options(na.action = "na.omit") #change NA default back


#There are four supported models, for the first cut lets just look at the "best"
ir1<-glmmTMB(INDI.RICH ~ Area + C3_C4 + PCT.STATE2 + 
               (1|Polygon), data = dat, family = poisson)
summary(ir1)
r2_nakagawa(ir1)

#check model
simulationOutput <- simulateResiduals(fittedModel = ir1, n = 250)
plotSimulatedResiduals(simulationOutput = simulationOutput) #NO issues. 

####
#Predicitions 
####

#PREDICTIONS
#Make Dataframes
ir.predA <- data.frame(Area = unique(dat$Area),         #For predicting based on Area
                       Model = "Indicator Richness",
                       Predicting = "Area",
                       Polygon = NA,
                       C3_C4 = "Native_C4",
                       PCT.STATE2 = "Woodland")
ir.predB1 <- data.frame(Area = "MULL",                  #For predicting based on Mapped C3/C4
                       Model = "Indicator Richness",
                       Predicting = "C3_C4",
                       Polygon = NA,
                       C3_C4 = unique(dat$C3_C4),
                       PCT.STATE2 = "Woodland")
ir.predB2 <- data.frame(Area = "GOO",                   #For predicting based on Mapped C3/C4
                        Model = "Indicator Richness",
                        Predicting = "C3_C4",
                        Polygon = NA,
                        C3_C4 = unique(dat$C3_C4),
                        PCT.STATE2 = "Woodland")
ir.predC1 <- data.frame(Area = "MULL",                   #For predicting based on PCT State
                       Model = "Indicator Richness",
                       Predicting = "PCT.STATE",
                       Polygon = NA,
                       C3_C4 = "Native_C4",
                       PCT.STATE2 = unique(dat$PCT.STATE2))
ir.predC2 <- data.frame(Area = "GOO",                   #For predicting based on PCT State
                       Model = "Indicator Richness",
                       Predicting = "PCT.STATE",
                       Polygon = NA,
                       C3_C4 = "Native_C4",
                       PCT.STATE2 = unique(dat$PCT.STATE2))


# Generate predictions from best model
{
P1 <- predict(ir1, newdata = ir.predA, type = "response", se.fit = TRUE, re.form = NA)  #Predictions for Area
ir.predA$pred <- P1$fit 
ir.predA$SE <- P1$se.fit 
ir.predA$up <- (ir.predA$pred+(1.96 * ir.predA$SE)) 
ir.predA$low <- (ir.predA$pred-(1.96 * ir.predA$SE))
ir.predA$resp <- (ir.predA$pred)

P1 <- predict(ir1, newdata = ir.predB1, type = "response", se.fit = TRUE, re.form = NA)  #Predictions for C3_C4 "MULL"
ir.predB1$pred <- P1$fit 
ir.predB1$SE <- P1$se.fit 
ir.predB1$up <- (ir.predB1$pred+(1.96 * ir.predB1$SE)) 
ir.predB1$low <- (ir.predB1$pred-(1.96 * ir.predB1$SE))
ir.predB1$resp <- (ir.predB1$pred)

P1 <- predict(ir1, newdata = ir.predB2, type = "response", se.fit = TRUE, re.form = NA)  #Predictions for C3_C4 "GOO"
ir.predB2$pred <- P1$fit 
ir.predB2$SE <- P1$se.fit 
ir.predB2$up <- (ir.predB2$pred+(1.96 * ir.predB2$SE)) 
ir.predB2$low <- (ir.predB2$pred-(1.96 * ir.predB2$SE))
ir.predB2$resp <- (ir.predB2$pred)

P1 <- predict(ir1, newdata = ir.predC1, type = "response", se.fit = TRUE, re.form = NA)  #Predictions for PCT STATE "MULL"
ir.predC1$pred <- P1$fit 
ir.predC1$SE <- P1$se.fit 
ir.predC1$up <- (ir.predC1$pred+(1.96 * ir.predC1$SE)) 
ir.predC1$low <- (ir.predC1$pred-(1.96 * ir.predC1$SE))
ir.predC1$resp <- (ir.predC1$pred)

P1 <- predict(ir1, newdata = ir.predC2, type = "response", se.fit = TRUE, re.form = NA)  #Predictions for PCT STATE "GOO"
ir.predC2$pred <- P1$fit 
ir.predC2$SE <- P1$se.fit 
ir.predC2$up <- (ir.predC2$pred+(1.96 * ir.predC2$SE)) 
ir.predC2$low <- (ir.predC2$pred-(1.96 * ir.predC2$SE))
ir.predC2$resp <- (ir.predC2$pred)

}
# Attach all prediction df's together and export

ir.pred.df <- rbind(ir.predA,
                    ir.predB1, ir.predB2,
                    ir.predC1, ir.predC2)  

write.csv(ir.pred.df, "Indic.Richness_predictions.csv")


####
# Making Figures of Predictions
####

#labelbg = paste0("m~italic(R)^2 == 0.04:~italic(P) < 0.01")

#PLOT Indicator Richness for Area

Fig.IR.Area <- ir.pred.df %>% 
  filter(Predicting == "Area") %>%
  mutate(Area = fct_relevel(Area, "OFF", "MULL_OUT", "GOO", "MULL"))%>%
  
  ggplot(aes(y=Area))+
  scale_y_discrete(labels = c("Offsets", "Mulligans (outside)", "Goorooyarroo", "Mulligans (inside)"))+
  
  geom_pointrange(aes(x=resp, xmin=low, xmax=up)) +
  
  theme_bw()+xlab("Indicator richness")+ylab("Area of the MFWS")
Fig.IR.Area


Fig.IR.C3C4 <- ir.pred.df %>% 
  filter(Predicting == "C3_C4") %>%
  mutate(C3_C4 = fct_relevel(C3_C4, "Exotic_C3", "Native_C3", "Native_C4"))%>%
  mutate(low = if_else(low < 0, 0, low))%>%
  
  ggplot(aes(y=C3_C4))+
  scale_y_discrete(labels = c("Exotic C3", "Native C3", "Native C4"))+
  
  geom_pointrange(aes(x=resp, xmin=low, xmax=up, colour = Area),
                      position=position_dodge(width=0.3)) +
  
  theme_bw()+xlab("Indicator richness")+ylab("Mapped grass type")
Fig.IR.C3C4


Fig.IR.PCT <- ir.pred.df %>% 
  filter(Predicting == "PCT.STATE")%>%
  filter(!row_number()%in%c(2))%>%
  mutate(PCT.STATE2 = fct_relevel(PCT.STATE2, "Forest", "Exotic Woodland", "Derived Grassland", "Woodland"))%>%
  mutate(low = if_else(low < 0, 0, low))%>%
  
  ggplot(aes(y=PCT.STATE2))+
  scale_y_discrete(labels = c("Forest", "Exotic \nWoodland", "Derived \nGrassland", "Woodland"))+
  
  geom_pointrange(aes(x=resp, xmin=low, xmax=up, colour = Area),
                  position=position_dodge(width=0.3)) +
  
  theme_bw()+xlab("Indicator richness")+ylab("PCT Condition State")
Fig.IR.PCT


# 3.2  SPECIES RICHNESS__________________________####

# “Species richness” is the count of targeted species present at a site. 
# Four different orchid species were observed meaning the total species pool considered was 13.

ggplot(dat, aes(SPP.RICH))+  #check distribution of response variable
  geom_histogram()            #All should poisson (b/c it's count data)
#unless there's overdispersion


# Full model

sr.f <- glmmTMB(SPP.RICH ~ Area + 
                  PCT.STATE2 +
                  UMC_ID +
                  Exp.Treat +
                  C3_C4 +
                  DomGrass.sum +
                  zCanopy.Cov + zUnderS.Cov + zNatGrass + zNatShrub+
                  zNatOther+ zRock+ zBareGround+ zLitter+
                  zThatch+ zExBroadleaf+ zExGrass+
                  zNativeness+ zAvGrassHt +
                  (1|Polygon),
                data = dat,
                family = poisson)


summary(sr.f)


#check model
simulationOutput <- simulateResiduals(fittedModel = sr.f, n = 250)
plotSimulatedResiduals(simulationOutput = simulationOutput) #Overdispersed. 
testUniformity(simulationOutput = simulationOutput)
testDispersion(sr.f)
testZeroInflation(simulationOutput)
testSpatialAutocorrelation(simulationOutput = simulationOutput, x = dat$x, y= dat$y)


#Model Selection on Full model
options(na.action = "na.fail")
#sr.d <- dredge(sr.f) #Dont run again unless global model changes
sr.d
sr.sub<-subset(sr.d, delta<=2,  recalc.weights = FALSE)
sr.sub #What are the supported (delta AICc <2) models?
options(na.action = "na.omit")

#Just one best model
sr1<-glmmTMB(SPP.RICH ~ Area + PCT.STATE2 + 
               (1|Polygon), data = dat, family = poisson)
summary(sr1)
r2_nakagawa(sr1)

#check model
simulationOutput <- simulateResiduals(fittedModel = sr1, n = 250)
plotSimulatedResiduals(simulationOutput = simulationOutput) #NO issues. 

####
#Predicitions 
####

#PREDICTIONS
#Make Dataframes
sr.predA <- data.frame(Area = unique(dat$Area),         #For predicting based on Area
                       Model = "Species Richness",
                       Predicting = "Area",
                       Polygon = NA,
                       PCT.STATE2 = "Woodland")
sr.predC1 <- data.frame(Area = "MULL",                   #For predicting based on PCT State (MULL)
                        Model = "Species Richness",
                        Predicting = "PCT.STATE",
                        Polygon = NA,
                        PCT.STATE2 = unique(dat$PCT.STATE2))
sr.predC2 <- data.frame(Area = "GOO",                   #For predicting based on PCT State (GOO)
                        Model = "Species Richness",
                        Predicting = "PCT.STATE",
                        Polygon = NA,
                        PCT.STATE2 = unique(dat$PCT.STATE2))


# Generate predictions from best model
{
  P1 <- predict(sr1, newdata = sr.predA, type = "response", se.fit = TRUE, re.form = NA)  #Predictions for Area
  sr.predA$pred <- P1$fit 
  sr.predA$SE <- P1$se.fit 
  sr.predA$up <- (sr.predA$pred+(1.96 * sr.predA$SE)) 
  sr.predA$low <- (sr.predA$pred-(1.96 * sr.predA$SE))
  sr.predA$resp <- (sr.predA$pred)
  
  P1 <- predict(sr1, newdata = sr.predC1, type = "response", se.fit = TRUE, re.form = NA)  #Predictions for PCT STATE "MULL"
  sr.predC1$pred <- P1$fit 
  sr.predC1$SE <- P1$se.fit 
  sr.predC1$up <- (sr.predC1$pred+(1.96 * sr.predC1$SE)) 
  sr.predC1$low <- (sr.predC1$pred-(1.96 * sr.predC1$SE))
  sr.predC1$resp <- (sr.predC1$pred)
  
  P1 <- predict(sr1, newdata = sr.predC2, type = "response", se.fit = TRUE, re.form = NA)  #Predictions for PCT STATE "GOO"
  sr.predC2$pred <- P1$fit 
  sr.predC2$SE <- P1$se.fit 
  sr.predC2$up <- (sr.predC2$pred+(1.96 * sr.predC2$SE)) 
  sr.predC2$low <- (sr.predC2$pred-(1.96 * sr.predC2$SE))
  sr.predC2$resp <- (sr.predC2$pred)
  
}
# Attach all prediction df's together and export

sr.pred.df <- rbind(sr.predA,
                    sr.predC1, sr.predC2)  

write.csv(sr.pred.df, "Spp.Richness_predictions.csv")


####
# Making Figures of Predictions
####

#labelbg = paste0("m~italic(R)^2 == 0.04:~italic(P) < 0.01")

#PLOT Indicator Richness for Area

Fig.sr.Area <- sr.pred.df %>% 
  filter(Predicting == "Area") %>%
  mutate(Area = fct_relevel(Area, "OFF", "MULL_OUT", "GOO", "MULL"))%>%
  
  ggplot(aes(y=Area))+
  scale_y_discrete(labels = c("Offsets", "Mulligans (outside)", "Goorooyarroo", "Mulligans (inside)"))+
  
  geom_pointrange(aes(x=resp, xmin=low, xmax=up)) +
  
  theme_bw()+xlab("Species richness")+ylab("Area of the MFWS")
Fig.sr.Area

Fig.sr.PCT <- sr.pred.df %>% 
  filter(Predicting == "PCT.STATE")%>%
  filter(!row_number()%in%c(2))%>%
  mutate(PCT.STATE2 = fct_relevel(PCT.STATE2, "Forest", "Exotic Woodland", "Derived Grassland", "Woodland"))%>%
  mutate(low = if_else(low < 0, 0, low))%>%
  
  ggplot(aes(y=PCT.STATE2))+
  scale_y_discrete(labels = c("Forest", "Exotic \nWoodland", "Derived \nGrassland", "Woodland"))+
  
  geom_pointrange(aes(x=resp, xmin=low, xmax=up, colour = Area),
                  position=position_dodge(width=0.3)) +
  
  theme_bw()+xlab("Species richness")+ylab("PCT Condition State")
Fig.sr.PCT

# 3.3  MAX INDICATOR ABUNDANCE___________________####

#“Indicator abundance” is the “number of individuals” count for the most abundant indicator at a site. 
# The count scale is in orders of magnitude (1 = at least 1 individual but <10, 2 = at least 10 individuals but <100, 
# 3 = at least 100 individuals but <1000, 4 = 1000 or more individuals

ggplot(dat, aes(MAX.ABUND))+  #check distribution of response variable
  geom_histogram()            #All should poisson (b/c it's count data)
                              #unless there's overdispersion


# Full model

ma.f <- glmmTMB(MAX.ABUND ~ Area + 
                  PCT.STATE2 +
                  UMC_ID +
                  Exp.Treat +
                  C3_C4 +
                  DomGrass.sum +
                  zCanopy.Cov + zUnderS.Cov + zNatGrass + zNatShrub+
                  zNatOther+ zRock+ zBareGround+ zLitter+
                  zThatch+ zExBroadleaf+ zExGrass+
                  zNativeness+ zAvGrassHt +
                  (1|Polygon),
                data = dat,
                family = poisson)


summary(ma.f)


#check model
simulationOutput <- simulateResiduals(fittedModel = ma.f, n = 250)
plotSimulatedResiduals(simulationOutput = simulationOutput) #Overdispersed. 
testUniformity(simulationOutput = simulationOutput)
testDispersion(ma.f)
testZeroInflation(simulationOutput)
testSpatialAutocorrelation(simulationOutput = simulationOutput, x = dat$x, y= dat$y)


#Model Selection on Full model
options(na.action = "na.fail")
#ma.d <- dredge(ma.f) #Dont run again unless global model changes
ma.d
ma.sub<-subset(ma.d, delta<=2,  recalc.weights = FALSE)
ma.sub #What are the supported (delta AICc <2) models?
options(na.action = "na.omit")

#Just one best model
ma1<-glmmTMB(MAX.ABUND ~ Area + PCT.STATE2 + 
               (1|Polygon), data = dat, family = poisson)
summary(ma1)
r2_nakagawa(ma1)

#check model
simulationOutput <- simulateResiduals(fittedModel = ma1, n = 250)
plotSimulatedResiduals(simulationOutput = simulationOutput) #very minor issues. 

####
#Predicitions 
####

#PREDICTIONS
#Make Dataframes
ma.predA <- data.frame(Area = unique(dat$Area),         #For predicting based on Area
                       Model = "Max Abundance",
                       Predicting = "Area",
                       Polygon = NA,
                       PCT.STATE2 = "Woodland")
ma.predC1 <- data.frame(Area = "MULL",                   #For predicting based on PCT State (MULL)
                        Model = "Max Abundance",
                        Predicting = "PCT.STATE",
                        Polygon = NA,
                        PCT.STATE2 = unique(dat$PCT.STATE2))
ma.predC2 <- data.frame(Area = "GOO",                   #For predicting based on PCT State (GOO)
                        Model = "Max Abundance",
                        Predicting = "PCT.STATE",
                        Polygon = NA,
                        PCT.STATE2 = unique(dat$PCT.STATE2))


# Generate predictions from best model
{
  P1 <- predict(ma1, newdata = ma.predA, type = "response", se.fit = TRUE, re.form = NA)  #Predictions for Area
  ma.predA$pred <- P1$fit 
  ma.predA$SE <- P1$se.fit 
  ma.predA$up <- (ma.predA$pred+(1.96 * ma.predA$SE)) 
  ma.predA$low <- (ma.predA$pred-(1.96 * ma.predA$SE))
  ma.predA$resp <- (ma.predA$pred)
  
  P1 <- predict(ma1, newdata = ma.predC1, type = "response", se.fit = TRUE, re.form = NA)  #Predictions for PCT STATE "MULL"
  ma.predC1$pred <- P1$fit 
  ma.predC1$SE <- P1$se.fit 
  ma.predC1$up <- (ma.predC1$pred+(1.96 * ma.predC1$SE)) 
  ma.predC1$low <- (ma.predC1$pred-(1.96 * ma.predC1$SE))
  ma.predC1$resp <- (ma.predC1$pred)
  
  P1 <- predict(ma1, newdata = ma.predC2, type = "response", se.fit = TRUE, re.form = NA)  #Predictions for PCT STATE "GOO"
  ma.predC2$pred <- P1$fit 
  ma.predC2$SE <- P1$se.fit 
  ma.predC2$up <- (ma.predC2$pred+(1.96 * ma.predC2$SE)) 
  ma.predC2$low <- (ma.predC2$pred-(1.96 * ma.predC2$SE))
  ma.predC2$resp <- (ma.predC2$pred)
  
}
# Attach all prediction df's together and export

ma.pred.df <- rbind(ma.predA,
                    ma.predC1, ma.predC2)  

write.csv(ma.pred.df, "Max.Abundance_predictions.csv")


####
# Making Figures of Predictions
####

#labelbg = paste0("m~italic(R)^2 == 0.04:~italic(P) < 0.01")

#PLOT Max ABundance for Area

Fig.ma.Area <- ma.pred.df %>% 
  filter(Predicting == "Area") %>%
  mutate(Area = fct_relevel(Area, "OFF", "MULL_OUT", "GOO", "MULL"))%>%
  
  ggplot(aes(y=Area))+
  scale_y_discrete(labels = c("Offsets", "Mulligans (outside)", "Goorooyarroo", "Mulligans (inside)"))+
  
  geom_pointrange(aes(x=resp, xmin=low, xmax=up)) +
  scale_x_continuous(limits = c(0,NA),
                     labels = c("0", "<10", "<100", "<1000"))+
  
  theme_bw()+xlab("Count of most abundant indicator group")+ylab("Area of the MFWS")
Fig.ma.Area

Fig.ma.PCT <- ma.pred.df %>% 
  filter(Predicting == "PCT.STATE")%>%
  filter(!row_number()%in%c(2))%>%
  mutate(PCT.STATE2 = fct_relevel(PCT.STATE2, "Forest", "Exotic Woodland", "Derived Grassland", "Woodland"))%>%
  mutate(low = if_else(low < 0, 0, low))%>%
  
  ggplot(aes(y=PCT.STATE2))+
  scale_y_discrete(labels = c("Forest", "Exotic \nWoodland", "Derived \nGrassland", "Woodland"))+
  
  geom_pointrange(aes(x=resp, xmin=low, xmax=up, colour = Area),
                  position=position_dodge(width=0.3)) +
  scale_x_continuous(limits = c(0,NA),
                     labels = c("0", "<10", "<100", "<1000"))+
  
  theme_bw()+xlab("Count of most abundant indicator group")+ylab("PCT Condition State")
Fig.ma.PCT


# 3.4  LILIES ABUNDANCE__________________________####

# Abundance the Lilies indicator group. Count is the scale previously as described 

ggplot(dat, aes(Lilies))+  #check distribution of response variable
  geom_histogram()            #All should poisson (b/c it's count data)
#unless there's overdispersion


# Full model

ll.f <- glmmTMB(Lilies ~ Area + 
                  PCT.STATE2 +
                  UMC_ID +
                  Exp.Treat +
                  C3_C4 +
                  DomGrass.sum +
                  zCanopy.Cov + zUnderS.Cov + zNatGrass + zNatShrub+
                  zNatOther+ zRock+ zBareGround+ zLitter+
                  zThatch+ zExBroadleaf+ zExGrass+
                  zNativeness+ zAvGrassHt +
                  (1|Polygon),
                data = dat,
                family = poisson)


summary(ll.f)


#check model
simulationOutput <- simulateResiduals(fittedModel = ll.f, n = 250)
plotSimulatedResiduals(simulationOutput = simulationOutput) #Overdispersed. 
testUniformity(simulationOutput = simulationOutput)
testDispersion(ll.f)
testZeroInflation(simulationOutput)
testSpatialAutocorrelation(simulationOutput = simulationOutput, x = dat$x, y= dat$y)


#Model Selection on Full model
options(na.action = "na.fail")
#ll.d <- dredge(ll.f) #Dont run again unless global model changes
ll.d
ll.sub<-subset(ll.d, delta<=2,  recalc.weights = FALSE)
ll.sub #What are the supported (delta AICc <2) models?
options(na.action = "na.omit")

#Just one best model
ll1<-glmmTMB(Lilies ~ Area + C3_C4 + PCT.STATE2 + 
               (1|Polygon), data = dat, family = poisson)
summary(ll1)
r2_nakagawa(ll1)

#check model
simulationOutput <- simulateResiduals(fittedModel = ll1, n = 250)
plotSimulatedResiduals(simulationOutput = simulationOutput) #NO issues. 

####
#Predicitions 
####

#PREDICTIONS
#Make Dataframes
ll.predA <- data.frame(Area = unique(dat$Area),         #For predicting based on Area
                       Model = "Lilies count",
                       Predicting = "Area",
                       Polygon = NA,
                       C3_C4 = "Native_C4",
                       PCT.STATE2 = "Woodland")
ll.predB1 <- data.frame(Area = "MULL",                  #For predicting based on Mapped C3/C4
                        Model = "Lilies count",
                        Predicting = "C3_C4",
                        Polygon = NA,
                        C3_C4 = unique(dat$C3_C4),
                        PCT.STATE2 = "Woodland")
ll.predB2 <- data.frame(Area = "GOO",                   #For predicting based on Mapped C3/C4
                        Model = "Lilies count",
                        Predicting = "C3_C4",
                        Polygon = NA,
                        C3_C4 = unique(dat$C3_C4),
                        PCT.STATE2 = "Woodland")
ll.predC1 <- data.frame(Area = "MULL",                   #For predicting based on PCT State (MULL)
                        Model = "Lilies count",
                        Predicting = "PCT.STATE",
                        Polygon = NA,
                        C3_C4 = "Native_C4",
                        PCT.STATE2 = unique(dat$PCT.STATE2))
ll.predC2 <- data.frame(Area = "GOO",                   #For predicting based on PCT State (GOO)
                        Model = "Lilies count",
                        Predicting = "PCT.STATE",
                        Polygon = NA,
                        C3_C4 = "Native_C4",
                        PCT.STATE2 = unique(dat$PCT.STATE2))


# Generate predictions from best model
{
  P1 <- predict(ll1, newdata = ll.predA, type = "response", se.fit = TRUE, re.form = NA)  #Predictions for Area
  ll.predA$pred <- P1$fit 
  ll.predA$SE <- P1$se.fit 
  ll.predA$up <- (ll.predA$pred+(1.96 * ll.predA$SE)) 
  ll.predA$low <- (ll.predA$pred-(1.96 * ll.predA$SE))
  ll.predA$resp <- (ll.predA$pred)
  
  P1 <- predict(ll1, newdata = ll.predB1, type = "response", se.fit = TRUE, re.form = NA)  #Predictions for PCT STATE "MULL"
  ll.predB1$pred <- P1$fit 
  ll.predB1$SE <- P1$se.fit 
  ll.predB1$up <- (ll.predB1$pred+(1.96 * ll.predB1$SE)) 
  ll.predB1$low <- (ll.predB1$pred-(1.96 * ll.predB1$SE))
  ll.predB1$resp <- (ll.predB1$pred)
  
  P1 <- predict(ll1, newdata = ll.predB2, type = "response", se.fit = TRUE, re.form = NA)  #Predictions for PCT STATE "GOO"
  ll.predB2$pred <- P1$fit 
  ll.predB2$SE <- P1$se.fit 
  ll.predB2$up <- (ll.predB2$pred+(1.96 * ll.predB2$SE)) 
  ll.predB2$low <- (ll.predB2$pred-(1.96 * ll.predB2$SE))
  ll.predB2$resp <- (ll.predB2$pred)
  
  P1 <- predict(ll1, newdata = ll.predC1, type = "response", se.fit = TRUE, re.form = NA)  #Predictions for PCT STATE "MULL"
  ll.predC1$pred <- P1$fit 
  ll.predC1$SE <- P1$se.fit 
  ll.predC1$up <- (ll.predC1$pred+(1.96 * ll.predC1$SE)) 
  ll.predC1$low <- (ll.predC1$pred-(1.96 * ll.predC1$SE))
  ll.predC1$resp <- (ll.predC1$pred)
  
  P1 <- predict(ll1, newdata = ll.predC2, type = "response", se.fit = TRUE, re.form = NA)  #Predictions for PCT STATE "GOO"
  ll.predC2$pred <- P1$fit 
  ll.predC2$SE <- P1$se.fit 
  ll.predC2$up <- (ll.predC2$pred+(1.96 * ll.predC2$SE)) 
  ll.predC2$low <- (ll.predC2$pred-(1.96 * ll.predC2$SE))
  ll.predC2$resp <- (ll.predC2$pred)
  
}
# Attach all prediction df's together and export

ll.pred.df <- rbind(ll.predA,
                    ll.predB1, ll.predB2,
                    ll.predC1, ll.predC2)  

write.csv(ll.pred.df, "Lilies_predictions.csv")


####
# Making Figures of Predictions
####

#labelbg = paste0("m~italic(R)^2 == 0.04:~italic(P) < 0.01")

#PLOT Max ABundance for Area

Fig.ll.Area <- ll.pred.df %>% 
  filter(Predicting == "Area") %>%
  mutate(Area = fct_relevel(Area, "OFF", "MULL_OUT", "GOO", "MULL"))%>%
  
  ggplot(aes(y=Area))+
  scale_y_discrete(labels = c("Offsets", "Mulligans (outside)", "Goorooyarroo", "Mulligans (inside)"))+
  
  geom_pointrange(aes(x=resp, xmin=low, xmax=up)) +
  scale_x_continuous(limits = c(0,NA),
                     labels = c("0", "<10", "<100", "<1000"))+
  
  theme_bw()+xlab("Count of individuals")+ylab("Area of the MFWS")+
  ggtitle("Lilies Indicator Group")
Fig.ll.Area

Fig.ll.C3C4 <- ll.pred.df %>% 
  filter(Predicting == "C3_C4") %>%
  mutate(C3_C4 = fct_relevel(C3_C4, "Exotic_C3", "Native_C3", "Native_C4"))%>%
  mutate(low = if_else(low < 0, 0, low))%>%
  
  ggplot(aes(y=C3_C4))+
  scale_y_discrete(labels = c("Exotic C3", "Native C3", "Native C4"))+
  
  geom_pointrange(aes(x=resp, xmin=low, xmax=up, colour = Area),
                  position=position_dodge(width=0.3)) +
  scale_x_continuous(limits = c(0,NA),
                     labels = c("0", "<10", "<100", "<1000"))+
  
  theme_bw()+xlab("Count of individuals")+ylab("Mapped grass type")+
  ggtitle("Lilies Indicator Group")
Fig.ll.C3C4


Fig.ll.PCT <- ll.pred.df %>% 
  filter(Predicting == "PCT.STATE")%>%
  filter(!row_number()%in%c(2))%>%
  mutate(PCT.STATE2 = fct_relevel(PCT.STATE2, "Forest", "Exotic Woodland", "Derived Grassland", "Woodland"))%>%
  mutate(low = if_else(low < 0, 0, low))%>%
  
  ggplot(aes(y=PCT.STATE2))+
  scale_y_discrete(labels = c("Forest", "Exotic \nWoodland", "Derived \nGrassland", "Woodland"))+
  
  geom_pointrange(aes(x=resp, xmin=low, xmax=up, colour = Area),
                  position=position_dodge(width=0.3)) +
  scale_x_continuous(limits = c(0,NA),
                     labels = c("0", "<10", "<100", "<1000"))+
  
  theme_bw()+xlab("Count of individuals")+ylab("PCT Condition State")+
  ggtitle("Lilies Indicator Group")
Fig.ll.PCT



# 3.5  ORCHIDS ABUNDANCE_________________________####

# Abundance the Orchids indicator group. Count is the scale previously as described 

ggplot(dat, aes(Orchids))+  #check distribution of response variable
  geom_histogram()            #All should poisson (b/c it's count data)
#unless there's overdispersion


# Full model

or.f <- glmmTMB(Orchids ~ Area + 
                  PCT.STATE +
                  UMC_ID +
                  Exp.Treat +
                  C3_C4 +
                  DomGrass.sum +
                  zCanopy.Cov + zUnderS.Cov + zNatGrass + zNatShrub+
                  zNatOther+ zRock+ zBareGround+ zLitter+
                  zThatch+ zExBroadleaf+ zExGrass+
                  zNativeness+ zAvGrassHt +
                  (1|Polygon),
                data = dat,
                family = poisson)


summary(or.f)


#check model
simulationOutput <- simulateResiduals(fittedModel = or.f, n = 250)
plotSimulatedResiduals(simulationOutput = simulationOutput) #Overdispersed. 
testUniformity(simulationOutput = simulationOutput)
testDispersion(or.f)
testZeroInflation(simulationOutput)
testSpatialAutocorrelation(simulationOutput = simulationOutput, x = dat$x, y= dat$y)


#Model Selection on Full model
options(na.action = "na.fail")
#or.d <- dredge(or.f) #only run again if global model changes
or.d
or.sub<-subset(or.d, delta<=2,  recalc.weights = FALSE)
or.sub #What are the supported (delta AICc <2) models?
options(na.action = "na.omit")

#Just one best model (NO AREA)
or1<-glmmTMB(Orchids ~ C3_C4 + PCT.STATE2 + UMC_ID + 
               (1|Polygon) + (1|Obs), data = dat, family = poisson)
summary(or1)
r2_nakagawa(or1)

#check model
simulationOutput <- simulateResiduals(fittedModel = or1, n = 250)
plotSimulatedResiduals(simulationOutput = simulationOutput) #NO issues. 

####
#Predicitions 
####

#PREDICTIONS
#Make Dataframes
or.predA <- data.frame(UMC_ID = unique(dat$UMC_ID),         #For predicting based on Area
                       Model = "Orchids count",
                       Predicting = "UMC_ID",
                       Polygon = NA,
                       Obs = NA,
                       C3_C4 = "Native_C4",
                       PCT.STATE2 = "Woodland")
or.predB1 <- data.frame(UMC_ID = "u19",                  #For predicting based on Mapped C3/C4
                        Model = "Orchids count",
                        Predicting = "C3_C4",
                        Polygon = NA,
                        Obs = NA,
                        C3_C4 = unique(dat$C3_C4),
                        PCT.STATE2 = "Woodland")
or.predC1 <- data.frame(UMC_ID = "u19",                   #For predicting based on PCT State
                        Model = "Orchids count",
                        Predicting = "PCT.STATE",
                        Polygon = NA,
                        Obs = NA,
                        C3_C4 = "Native_C4",
                        PCT.STATE2 = unique(dat$PCT.STATE2))


# Generate predictions from best model
{
  P1 <- predict(or1, newdata = or.predA, type = "response", se.fit = TRUE, re.form = NA)  #Predictions for Area
  or.predA$pred <- P1$fit 
  or.predA$SE <- P1$se.fit 
  or.predA$up <- (or.predA$pred+(1.96 * or.predA$SE)) 
  or.predA$low <- (or.predA$pred-(1.96 * or.predA$SE))
  or.predA$resp <- (or.predA$pred)
  
  P1 <- predict(or1, newdata = or.predB1, type = "response", se.fit = TRUE, re.form = NA)  #Predictions for PCT STATE "MULL"
  or.predB1$pred <- P1$fit 
  or.predB1$SE <- P1$se.fit 
  or.predB1$up <- (or.predB1$pred+(1.96 * or.predB1$SE)) 
  or.predB1$low <- (or.predB1$pred-(1.96 * or.predB1$SE))
  or.predB1$resp <- (or.predB1$pred)
  
  
  P1 <- predict(or1, newdata = or.predC1, type = "response", se.fit = TRUE, re.form = NA)  #Predictions for PCT STATE "MULL"
  or.predC1$pred <- P1$fit 
  or.predC1$SE <- P1$se.fit 
  or.predC1$up <- (or.predC1$pred+(1.96 * or.predC1$SE)) 
  or.predC1$low <- (or.predC1$pred-(1.96 * or.predC1$SE))
  or.predC1$resp <- (or.predC1$pred)
  

  
}
# Attach all prediction df's together and export

or.pred.df <- rbind(or.predA,
                    or.predB1,
                    or.predC1)  

write.csv(or.pred.df, "Orchids_predictions.csv")


####
# Making Figures of Predictions
####

#labelbg = paste0("m~italic(R)^2 == 0.04:~italic(P) < 0.01")

#PLOT Orchids for PCT.STATE2
Fig.or.C3C4 <- or.pred.df %>% 
  filter(Predicting == "C3_C4") %>%
  mutate(C3_C4 = fct_relevel(C3_C4, "Exotic_C3", "Native_C3", "Native_C4"))%>%
  mutate(low = if_else(low < 0, 0, low))%>%
  
  ggplot(aes(y=C3_C4))+
  scale_y_discrete(labels = c("Exotic C3", "Native C3", "Native C4"))+
  
  geom_pointrange(aes(x=resp, xmin=low, xmax=up)) +
  scale_x_continuous(limits = c(0,NA),
                     labels = c("0", "<10", "<100", "<1000"))+
  
  theme_bw()+xlab("Count of individuals")+ylab("Mapped grass type")+
  ggtitle("Orchids Indicator Group")
Fig.or.C3C4


Fig.or.PCT <- or.pred.df %>% 
  filter(Predicting == "PCT.STATE")%>%
  mutate(PCT.STATE2 = fct_relevel(PCT.STATE2, "Forest", "Exotic Woodland", "Derived Grassland", "Woodland"))%>%
  mutate(low = if_else(low < 0, 0, low))%>%
  
  ggplot(aes(y=PCT.STATE2))+
  scale_y_discrete(labels = c("Forest", "Exotic \nWoodland", "Derived \nGrassland", "Woodland"))+
  
  geom_pointrange(aes(x=resp, xmin=low, xmax=up)) +
  scale_x_continuous(limits = c(0,NA),
                     labels = c("0", "<10", "<100", "<1000"))+
  
  theme_bw()+xlab("Count of individuals")+ylab("PCT Condition State")+
  ggtitle("Orchids Indicator Group")
Fig.or.PCT


#UMC_ID doesn't work

#   3.5.1  ORCHID SPP RICHNESS___________________####

# Create a "orchid species richness" response variable

dat <- dat %>%
  mutate(ORCH.SPPRICH = Microt.spp + Thely.spp + Diur.spp + Hymen.spp) #sum of occurrences of 4 genera


ggplot(dat, aes(ORCH.SPPRICH))+  #check distribution of response variable
  geom_histogram()            #All should poisson (b/c it's count data)
#unless there's overdispersion


# Full model

orSR.f <- glmmTMB(ORCH.SPPRICH ~ Area + 
                  PCT.STATE +
                  UMC_ID +
                  Exp.Treat +
                  C3_C4 +
                  DomGrass.sum +
                  zCanopy.Cov + zUnderS.Cov + zNatGrass + zNatShrub+
                  zNatOther+ zRock+ zBareGround+ zLitter+
                  zThatch+ zExBroadleaf+ zExGrass+
                  zNativeness+ zAvGrassHt +
                  (1|Polygon),
                data = dat,
                family = poisson)


summary(orSR.f)


#check model
simulationOutput <- simulateResiduals(fittedModel = orSR.f, n = 250)
plotSimulatedResiduals(simulationOutput = simulationOutput) #Overdispersed. 
testUniformity(simulationOutput = simulationOutput)
testDispersion(orSR.f)
testZeroInflation(simulationOutput)
testSpatialAutocorrelation(simulationOutput = simulationOutput, x = dat$x, y= dat$y)


#Model Selection on Full model
options(na.action = "na.fail")
orSR.d <- dredge(orSR.f) #only run again if global model changes
orSR.d
orSR.sub<-subset(orSR.d, delta<=2,  recalc.weights = FALSE)
orSR.sub #What are the supported (delta AICc <2) models?
options(na.action = "na.omit")

# 5 supported models, lets just look at the best model
orSR1<-glmmTMB(ORCH.SPPRICH ~ Area + PCT.STATE2 +
               (1|Polygon), data = dat, family = poisson)
summary(orSR1)
r2_nakagawa(orSR1)

#check model
simulationOutput <- simulateResiduals(fittedModel = orSR1, n = 250)
plotSimulatedResiduals(simulationOutput = simulationOutput) #NO issues. 


####
#Predicitions 
####

#PREDICTIONS
#Make Dataframes
orSR.predA <- data.frame(Area = unique(dat$Area),         #For predicting based on Area
                       Model = "Orch Spp Rich",
                       Predicting = "Area",
                       Polygon = NA,
                       PCT.STATE2 = "Woodland")



# Generate predictions from best model
{
  P1 <- predict(orSR1, newdata = orSR.predA, type = "response", se.fit = TRUE, re.form = NA)  #Predictions for Area
  orSR.predA$pred <- P1$fit 
  orSR.predA$SE <- P1$se.fit 
  orSR.predA$up <- (orSR.predA$pred+(1.96 * orSR.predA$SE)) 
  orSR.predA$low <- (orSR.predA$pred-(1.96 * orSR.predA$SE))
  orSR.predA$resp <- (orSR.predA$pred)
  
  
}
# Attach all prediction df's together and export

orSR.pred.df <- rbind(orSR.predA)


write.csv(orSR.pred.df, "OrchSppRich_predictions.csv")

####
# Making Figures of Predictions
####

#PLOT ANTTHRO for Area

Fig.orSR.Area <- orSR.pred.df %>% 
  filter(Predicting == "Area") %>%
  mutate(Area = fct_relevel(Area, "OFF", "MULL_OUT", "GOO", "MULL"))%>%
  mutate(low = if_else(low < 0, 0, low))%>%
  
  ggplot(aes(y=Area))+
  scale_y_discrete(labels = c("Offsets", "Mulligans (outside)", "Goorooyarroo", "Mulligans (inside)"))+
  
  geom_pointrange(aes(x=resp, xmin=low, xmax=up)) +
  
  theme_bw()+xlab("Orchid Species Richness")+ylab("Area of the MFWS")
Fig.orSR.Area


# 3.6  DAISIES (TUBEROUS) ABUNDANCE______________####

# Abundance the Daisies (Tuberous) indicator group. Count is the scale previously as described 

ggplot(dat, aes(Daisies.T))+  #check distribution of response variable
  geom_histogram()            #All should poisson (b/c it's count data)
#unless there's overdispersion


# Full model

dT.f <- glmmTMB(Daisies.T ~ Area + 
                  PCT.STATE +
                  UMC_ID +
                  Exp.Treat +
                  C3_C4 +
                  DomGrass.sum +
                  zCanopy.Cov + zUnderS.Cov + zNatGrass + zNatShrub+
                  zNatOther+ zRock+ zBareGround+ zLitter+
                  zThatch+ zExBroadleaf+ zExGrass+
                  zNativeness+ zAvGrassHt +
                  (1|Polygon),
                data = dat,
                family = poisson)


summary(dT.f)


#check model
simulationOutput <- simulateResiduals(fittedModel = dT.f, n = 250)
plotSimulatedResiduals(simulationOutput = simulationOutput) #Overdispersed. 
testUniformity(simulationOutput = simulationOutput)
testDispersion(dT.f)
testZeroInflation(simulationOutput)
testSpatialAutocorrelation(simulationOutput = simulationOutput, x = dat$x, y= dat$y)


#Model Selection on Full model
options(na.action = "na.fail")
dT.d <- dredge(dT.f)
dT.d
dT.sub<-subset(dT.d, delta<=2,  recalc.weights = FALSE)
dT.sub #What are the supported (delta AICc <2) models?
options(na.action = "na.omit")

#Just one supported model that is no better than the null model
dT1<-glmmTMB(Daisies.T ~ C3_C4 + 
               (1|Polygon), data = dat, family = poisson)
summary(dT1)    #No sig effect of mapped grass type
r2_nakagawa(dT1)

#check model
simulationOutput <- simulateResiduals(fittedModel = dT1, n = 250)
plotSimulatedResiduals(simulationOutput = simulationOutput) #NO issues. 

####
#Predicitions 
####

#PREDICTIONS
#Make Dataframes
dT.predB1 <- data.frame(                  #For predicting based on Mapped C3/C4
                        Model = "DaisiesT count",
                        Predicting = "C3_C4",
                        Polygon = NA,
                        C3_C4 = unique(dat$C3_C4)
                        )



# Generate predictions from best model
{
  
  P1 <- predict(dT1, newdata = dT.predB1, type = "response", se.fit = TRUE, re.form = NA)  #Predictions for PCT STATE "MULL"
  dT.predB1$pred <- P1$fit 
  dT.predB1$SE <- P1$se.fit 
  dT.predB1$up <- (dT.predB1$pred+(1.96 * dT.predB1$SE)) 
  dT.predB1$low <- (dT.predB1$pred-(1.96 * dT.predB1$SE))
  dT.predB1$resp <- (dT.predB1$pred)
  
  
}
# Attach all prediction df's together and export

dT.pred.df <- rbind(
                    dT.predB1
                    )  

write.csv(dT.pred.df, "Daisies.T_predictions.csv")


####
# Making Figures of Predictions
####

#labelbg = paste0("m~italic(R)^2 == 0.04:~italic(P) < 0.01")

#PLOT Daisies T for Area

Fig.dT.C3C4 <- dT.pred.df %>% 
  filter(Predicting == "C3_C4") %>%
  mutate(Area = fct_relevel(C3_C4, "Exotic_C3", "Native_C3", "Native_C4"))%>%
  mutate(low = if_else(low < 0, 0, low))%>%
  
  ggplot(aes(y=C3_C4))+
  scale_y_discrete(labels = c("Exotic C3", "Native C3", "Native C4"))+
  
  geom_pointrange(aes(x=resp, xmin=low, xmax=up)) +
  scale_x_continuous(breaks = c(0,1),
                     labels = c("0", "<10"))+
  coord_cartesian(xlim = c(0, 1))+
  
  theme_bw()+xlab("Count of individuals")+ylab("Mapped grass type")+
  ggtitle("Daisies (tuberous) Indicator Group")

Fig.dT.C3C4

# 3.7  DAISIES (CONTROL) ABUNDANCE_______________####

# Abundance the Daisies (Tuberous) indicator group. Count is the scale previously as described 

ggplot(dat, aes(Daisies.C))+  #check distribution of response variable
  geom_histogram()            #All should poisson (b/c it's count data)
#unless there's overdispersion


# Full model

dC.f <- glmmTMB(Daisies.C ~ Area + 
                  PCT.STATE +
                  UMC_ID +
                  Exp.Treat +
                  C3_C4 +
                  DomGrass.sum +
                  zCanopy.Cov + zUnderS.Cov + zNatGrass + zNatShrub+
                  zNatOther+ zRock+ zBareGround+ zLitter+
                  zThatch+ zExBroadleaf+ zExGrass+
                  zNativeness+ zAvGrassHt +
                  (1|Polygon),
                data = dat,
                family = poisson)


summary(dC.f)


#check model
simulationOutput <- simulateResiduals(fittedModel = dC.f, n = 250)
plotSimulatedResiduals(simulationOutput = simulationOutput) #Overdispersed. 
testUniformity(simulationOutput = simulationOutput)
testDispersion(dC.f)
testZeroInflation(simulationOutput)
testSpatialAutocorrelation(simulationOutput = simulationOutput, x = dat$x, y= dat$y)


#Model Selection on Full model
options(na.action = "na.fail")
#dC.d <- dredge(dC.f)
dC.d
dC.sub<-subset(dC.d, delta<=2,  recalc.weights = FALSE)
dC.sub #What are the supported (delta AICc <2) models?
options(na.action = "na.omit")

#3 supported models, lets just look at the best model
dC1<-glmmTMB(Daisies.C ~ Area +  
               (1|Polygon), data = dat, family = poisson)
summary(dC1) #No sig effects
r2_nakagawa(dC1)

#check model
simulationOutput <- simulateResiduals(fittedModel = dC1, n = 250)
plotSimulatedResiduals(simulationOutput = simulationOutput) #NO issues. 

####
#Predicitions 
####

#PREDICTIONS
#Make Dataframes
dC.predA <- data.frame(Area = unique(dat$Area),         #For predicting based on Area
                       Model = "DaisiesC count",
                       Predicting = "Area",
                       Polygon = NA,
                       Obs = NA)


# Generate predictions from best model
{
  P1 <- predict(dC1, newdata = dC.predA, type = "response", se.fit = TRUE, re.form = NA)  #Predictions for Area
  dC.predA$pred <- P1$fit 
  dC.predA$SE <- P1$se.fit 
  dC.predA$up <- (dC.predA$pred+(1.96 * dC.predA$SE)) 
  dC.predA$low <- (dC.predA$pred-(1.96 * dC.predA$SE))
  dC.predA$resp <- (dC.predA$pred)
  
}
# Attach all prediction df's together and export

dC.pred.df <- rbind(dC.predA)
                    

write.csv(dC.pred.df, "Daisies.C_predictions.csv")


####
# Making Figures of Predictions
####

#labelbg = paste0("m~italic(R)^2 == 0.04:~italic(P) < 0.01")

#PLOT Daisies T for Area

Fig.dC.Area <- dC.pred.df %>% 
  filter(Predicting == "Area") %>%
  mutate(Area = fct_relevel(Area, "OFF", "MULL_OUT", "GOO", "MULL"))%>%
  mutate(low = if_else(low < 0, 0, low))%>%
  
  ggplot(aes(y=Area))+
  scale_y_discrete(labels = c("Offsets", "Mulligans (outside)", "Goorooyarroo", "Mulligans (inside)"))+
  
  geom_pointrange(aes(x=resp, xmin=low, xmax=up)) +
  scale_x_continuous(limits = c(0,NA),
                     labels = c("0", "<10", "<100", "<1000"))+
  
  theme_bw()+xlab("Count of individuals")+ylab("Area of the MFWS")+
  ggtitle("Daisies (control) Indicator Group")
Fig.dC.Area


# 4. SOME DESCRIPTIVE RESULTS____________________####

#Proportion of plots that didn't have anything?
dat%>%count(SPP.RICH == 0)
(76/112)*100#%had something
(36/112)*100#had nothing

#Proportion of plots that didn't have anything by Area?
datBT%>%group_by(AreaX)%>%count(SPP.RICH == 0)
19/22 #MULL in === 86%
22/30 #GOO         73%
5/6   #MULL out    83%
9/30  #off         30%
21/24 #Mull BT     88%


#Proportion of plots that had Lilies?
dat%>%group_by(Area)%>%count(Lilies == 0)
24/46 #MULL in === 52%
18/30 #GOO         60%
4/6   #MULL out    66%
6/30  #off         20%
dat%>%count(Lilies == 0)
52/112 #overall    46%

#Proportion of plots that had Orchids?
dat%>%group_by(Area)%>%count(Orchids == 0)
29/46 #MULL in === 63%
13/30 #GOO         43%
3/6   #MULL out    50%
4/30  #off         13%
dat%>%count(Orchids == 0)
49/112 #overall    46%

#Proportion of plots that had Daisies T?
dat%>%group_by(Area)%>%count(Daisies.T == 0)
14/46 #MULL in === 30%
6/30 #GOO          20%
1/6   #MULL out    16%
4/30  #off         13%
dat%>%count(Daisies.T == 0)
25/112 #overall    22%

#Proportion of plots that had Daisies T?
dat%>%group_by(Area)%>%count(Daisies.C == 0)
11/46 #MULL in === 23%
18/30 #GOO         60%
2/6   #MULL out    33%
1/30  #off         3%
dat%>%count(Daisies.C == 0)
32/112 #overall    22%

#Proportion of plots that had Daisies T?
dat%>%group_by(Area)%>%count(Stackhousia == 0)
0/46 #MULL in ===  0%
5/30 #GOO          16%
0/6   #MULL out    0%
1/30  #off         3%


# 5. INDIVIDUAL SPECIES MODELS __________________####

# What % of observations for each species

Occ.Spp <- dat %>%
  summarise(Arth.spp = sum(Arth.spp)/112,       #30%
            Bulb.spp = sum(Bulb.spp)/112,       #13%
            Wurm.dio = sum(Wurm.dio)/112,       #27%
            Burch.umb = sum(Burch.umb)/112,     # 5%
            Microt.spp = sum(Microt.spp)/112,   #42%
            Thely.spp = sum(Thely.spp)/112,     #13%
            Diur.spp = sum(Diur.spp)/112,       # 4%
            Hymen.spp = sum(Hymen.spp)/112,     # 2%
            Micros.spp = sum(Micros.spp)/112,   # 5%
            Cras.spp = sum(Cras.spp)/112,       #16%
            Lept.squ = sum(Lept.squ)/112,       #20%
            Chry.api = sum(Chry.api)/112,       #15%
            Stac.mon = sum(Stac.mon)/112)       # 1%

# Binomial GLMM Full Models for 8 of 13 species present in >10% of plots
AS.f <- glmmTMB(Arth.spp ~ Area + 
                  PCT.STATE +
                  UMC_ID +
                  Exp.Treat +
                  C3_C4 +
                  DomGrass.sum +
                  zCanopy.Cov + zUnderS.Cov + zNatGrass + zNatShrub+
                  zNatOther+ zRock+ zBareGround+ zLitter+
                  zThatch+ zExBroadleaf+ zExGrass+
                  zNativeness+ zAvGrassHt +
                  (1|Polygon),
                data = dat,
                family = binomial)
BS.f <- glmmTMB(Bulb.spp ~ Area + 
                  PCT.STATE +
                  UMC_ID +
                  Exp.Treat +
                  C3_C4 +
                  DomGrass.sum +
                  zCanopy.Cov + zUnderS.Cov + zNatGrass + zNatShrub+
                  zNatOther+ zRock+ zBareGround+ zLitter+
                  zThatch+ zExBroadleaf+ zExGrass+
                  zNativeness+ zAvGrassHt +
                  (1|Polygon),
                data = dat,
                family = binomial)
WD.f <- glmmTMB(Wurm.dio ~ Area + 
                  PCT.STATE +
                  UMC_ID +
                  Exp.Treat +
                  C3_C4 +
                  DomGrass.sum +
                  zCanopy.Cov + zUnderS.Cov + zNatGrass + zNatShrub+
                  zNatOther+ zRock+ zBareGround+ zLitter+
                  zThatch+ zExBroadleaf+ zExGrass+
                  zNativeness+ zAvGrassHt +
                  (1|Polygon),
                data = dat,
                family = binomial)
MS.f <- glmmTMB(Microt.spp ~ Area + 
                  PCT.STATE +
                  UMC_ID +
                  Exp.Treat +
                  C3_C4 +
                  DomGrass.sum +
                  zCanopy.Cov + zUnderS.Cov + zNatGrass + zNatShrub+
                  zNatOther+ zRock+ zBareGround+ zLitter+
                  zThatch+ zExBroadleaf+ zExGrass+
                  zNativeness+ zAvGrassHt +
                  (1|Polygon),
                data = dat,
                family = binomial)
TS.f <- glmmTMB(Thely.spp ~ Area + 
                  PCT.STATE +
                  UMC_ID +
                  Exp.Treat +
                  C3_C4 +
                  DomGrass.sum +
                  zCanopy.Cov + zUnderS.Cov + zNatGrass + zNatShrub+
                  zNatOther+ zRock+ zBareGround+ zLitter+
                  zThatch+ zExBroadleaf+ zExGrass+
                  zNativeness+ zAvGrassHt +
                  (1|Polygon),
                data = dat,
                family = binomial)
CS.f <- glmmTMB(Cras.spp ~ Area + 
                  PCT.STATE +
                  UMC_ID +
                  Exp.Treat +
                  C3_C4 +
                  DomGrass.sum +
                  zCanopy.Cov + zUnderS.Cov + zNatGrass + zNatShrub+
                  zNatOther+ zRock+ zBareGround+ zLitter+
                  zThatch+ zExBroadleaf+ zExGrass+
                  zNativeness+ zAvGrassHt +
                  (1|Polygon),
                data = dat,
                family = binomial)
LS.f <- glmmTMB(Lept.squ ~ Area + 
                  PCT.STATE +
                  UMC_ID +
                  Exp.Treat +
                  C3_C4 +
                  DomGrass.sum +
                  zCanopy.Cov + zUnderS.Cov + zNatGrass + zNatShrub+
                  zNatOther+ zRock+ zBareGround+ zLitter+
                  zThatch+ zExBroadleaf+ zExGrass+
                  zNativeness+ zAvGrassHt +
                  (1|Polygon),
                data = dat,
                family = binomial)
CA.f <- glmmTMB(Chry.api ~ Area + 
                  PCT.STATE +
                  UMC_ID +
                  Exp.Treat +
                  C3_C4 +
                  DomGrass.sum +
                  zCanopy.Cov + zUnderS.Cov + zNatGrass + zNatShrub+
                  zNatOther+ zRock+ zBareGround+ zLitter+
                  zThatch+ zExBroadleaf+ zExGrass+
                  zNativeness+ zAvGrassHt +
                  (1|Polygon),
                data = dat,
                family = binomial)

# model selection on 8 species models
options(na.action = "na.fail")
AS.d <- dredge(AS.f)
BS.d <- dredge(BS.f)
WD.d <- dredge(WD.f)
MS.d <- dredge(MS.f)
TS.d <- dredge(TS.f)
CS.d <- dredge(CS.f)
LS.d <- dredge(LS.f)
CA.d <- dredge(CA.f)

# subset best models on 8 species

AS.sub<-subset(AS.d, delta<=2,  recalc.weights = FALSE)
BS.sub<-subset(BS.d, delta<=2,  recalc.weights = FALSE)
WD.sub<-subset(WD.d, delta<=2,  recalc.weights = FALSE)
MS.sub<-subset(MS.d, delta<=2,  recalc.weights = FALSE)
TS.sub<-subset(TS.d, delta<=2,  recalc.weights = FALSE)
CS.sub<-subset(CS.d, delta<=2,  recalc.weights = FALSE)
LS.sub<-subset(LS.d, delta<=2,  recalc.weights = FALSE)
CA.sub<-subset(CA.d, delta<=2,  recalc.weights = FALSE)
options(na.action = "na.omit")

# what are the best models on 8 species?

AS.sub  # 1 model. ~ Area + PCT
BS.sub  # best = null, ~ Area also supported
WD.sub  # 5 supported, best = ~ PCT
MS.sub  # best C3_C4 + UMC, also supported Area + PCT
TS.sub  # best = null, ~ UMC also supported
CS.sub  # best null, 3 also supported, ~ C3   UMC   C3+UMC
LS.sub  # 1 model. ~ Area + C3_C4
CA.sub  # 1 model. ~ Area



# 5.1. ARTHROPODIUM SPP. ________________________####


#BEST MODEL
AS1<-glmmTMB(Arth.spp ~ Area + PCT.STATE2 + 
               (1|Polygon), data = dat, family = binomial)
summary(AS1) 
r2_nakagawa(AS1)

#check model
simulationOutput <- simulateResiduals(fittedModel = AS1, n = 250)
plotSimulatedResiduals(simulationOutput = simulationOutput) #NO issues. 

####
#Predicitions 
####

#PREDICTIONS
#Make Dataframes
AS.predA <- data.frame(Area = unique(dat$Area),         #For predicting based on Area
                       Model = "Arthro",
                       Predicting = "Area",
                       Polygon = NA,
                       PCT.STATE2 = "Woodland")
AS.predC1 <- data.frame(Area = "MULL",                   #For predicting based on PCT State (MULL)
                        Model = "Arthro",
                        Predicting = "PCT.STATE",
                        Polygon = NA,
                        PCT.STATE2 = unique(dat$PCT.STATE2))
AS.predC2 <- data.frame(Area = "GOO",                   #For predicting based on PCT State (GOO)
                        Model = "Arthro",
                        Predicting = "PCT.STATE",
                        Polygon = NA,
                        PCT.STATE2 = unique(dat$PCT.STATE2))


# Generate predictions from best model
{
  P1 <- predict(AS1, newdata = AS.predA, type = "response", se.fit = TRUE, re.form = NA)  #Predictions for Area
  AS.predA$pred <- P1$fit 
  AS.predA$SE <- P1$se.fit 
  AS.predA$up <- (AS.predA$pred+(1.96 * AS.predA$SE)) 
  AS.predA$low <- (AS.predA$pred-(1.96 * AS.predA$SE))
  AS.predA$resp <- (AS.predA$pred)
  
  P1 <- predict(AS1, newdata = AS.predC1, type = "response", se.fit = TRUE, re.form = NA)  #Predictions for PCT STATE "MULL"
  AS.predC1$pred <- P1$fit 
  AS.predC1$SE <- P1$se.fit 
  AS.predC1$up <- (AS.predC1$pred+(1.96 * AS.predC1$SE)) 
  AS.predC1$low <- (AS.predC1$pred-(1.96 * AS.predC1$SE))
  AS.predC1$resp <- (AS.predC1$pred)
  
  P1 <- predict(AS1, newdata = AS.predC2, type = "response", se.fit = TRUE, re.form = NA)  #Predictions for PCT STATE "GOO"
  AS.predC2$pred <- P1$fit 
  AS.predC2$SE <- P1$se.fit 
  AS.predC2$up <- (AS.predC2$pred+(1.96 * AS.predC2$SE)) 
  AS.predC2$low <- (AS.predC2$pred-(1.96 * AS.predC2$SE))
  AS.predC2$resp <- (AS.predC2$pred)
  
}
# Attach all prediction df's together and export

AS.pred.df <- rbind(AS.predA,
                    AS.predC1,AS.predC2)


write.csv(AS.pred.df, "Anthro_predictions.csv")


####
# Making Figures of Predictions
####

#PLOT ANTTHRO for Area

Fig.AS.Area <- AS.pred.df %>% 
  filter(Predicting == "Area") %>%
  mutate(Area = fct_relevel(Area, "OFF", "MULL_OUT", "GOO", "MULL"))%>%
  mutate(low = if_else(low < 0, 0, low),
         up = if_else(up >1, 1, up))%>%
  
  ggplot(aes(y=Area))+
  scale_y_discrete(labels = c("Offsets", "Mulligans (outside)", "Goorooyarroo", "Mulligans (inside)"))+
  
  geom_pointrange(aes(x=resp, xmin=low, xmax=up)) +
  
  theme_bw()+xlab("Occurence")+ylab("Area of the MFWS")+
  ggtitle("Arthropodium spp.")
Fig.AS.Area

Fig.AS.PCT <- AS.pred.df %>% 
  filter(Predicting == "PCT.STATE")%>%
  filter(!row_number()%in%c(2))%>%
  mutate(PCT.STATE2 = fct_relevel(PCT.STATE2, "Forest", "Exotic Woodland", "Derived Grassland", "Woodland"))%>%
  mutate(low = if_else(low < 0, 0, low),
         up = if_else(up >1, 1, up))%>%
  
  ggplot(aes(y=PCT.STATE2))+
  scale_y_discrete(labels = c("Forest", "Exotic \nWoodland", "Derived \nGrassland", "Woodland"))+
  
  geom_pointrange(aes(x=resp, xmin=low, xmax=up, colour = Area),
                  position=position_dodge(width=0.3)) +
  
  theme_bw()+xlab("Occurence")+ylab("PCT Condition State")+
  ggtitle("Arthropodium spp.")
Fig.AS.PCT

# 5.2. BULBINE SPP. _____________________________####

#BEST MODEL
BS1<-glmmTMB(Bulb.spp ~ Area + 
               (1|Polygon), data = dat, family = binomial)
summary(BS1) 
r2_nakagawa(BS1)

#check model
simulationOutput <- simulateResiduals(fittedModel = BS1, n = 250)
plotSimulatedResiduals(simulationOutput = simulationOutput) #NO issues. 

####
#Predicitions 
####

#PREDICTIONS
#Make Dataframes
BS.predA <- data.frame(Area = unique(dat$Area),         #For predicting based on Area
                       Model = "Bulb",
                       Predicting = "Area",
                       Polygon = NA)


# Generate predictions from best model
{
  P1 <- predict(BS1, newdata = BS.predA, type = "response", se.fit = TRUE, re.form = NA)  #Predictions for Area
  BS.predA$pred <- P1$fit 
  BS.predA$SE <- P1$se.fit 
  BS.predA$up <- (BS.predA$pred+(1.96 * BS.predA$SE)) 
  BS.predA$low <- (BS.predA$pred-(1.96 * BS.predA$SE))
  BS.predA$resp <- (BS.predA$pred)

  
}
# Attach all prediction df's together and export

BS.pred.df <- rbind(BS.predA)


write.csv(BS.pred.df, "Bulbine_predictions.csv")


####
# Making Figures of Predictions
####

#PLOT ANTTHRO for Area

Fig.BS.Area <- BS.pred.df %>% 
  filter(Predicting == "Area") %>%
  mutate(Area = fct_relevel(Area, "OFF", "MULL_OUT", "GOO", "MULL"))%>%
  mutate(low = if_else(low < 0, 0, low),
         up = if_else(up >1, 1, up))%>%
  
  ggplot(aes(y=Area))+
  scale_y_discrete(labels = c("Offsets", "Mulligans (outside)", "Goorooyarroo", "Mulligans (inside)"))+
  
  geom_pointrange(aes(x=resp, xmin=low, xmax=up)) +
  #scale_x_continuous(limits = c(0,NA),
  # labels = c("0", "<10", "<100", "<1000"))+
  
  theme_bw()+xlab("Occurence")+ylab("Area of the MFWS")+
  ggtitle("Bulbine spp.")
Fig.BS.Area

# 5.3. WURMBEA DIOICA ___________________________####


#BEST MODEL
WD1<-glmmTMB(Wurm.dio ~ PCT.STATE2 + 
               (1|Polygon), data = dat, family = binomial)
summary(WD1) 
r2_nakagawa(WD1)

#check model
simulationOutput <- simulateResiduals(fittedModel = WD1, n = 250)
plotSimulatedResiduals(simulationOutput = simulationOutput) #NO issues. 

####
#Predicitions 
####

#PREDICTIONS
#Make Dataframes
WD.predA <- data.frame(PCT.STATE2 = unique(dat$PCT.STATE2),         #For predicting based on Area
                       Model = "Bulb",
                       Predicting = "PCT.STATE",
                       Polygon = NA)


# Generate predictions from best model
{
  P1 <- predict(WD1, newdata = WD.predA, type = "response", se.fit = TRUE, re.form = NA)  #Predictions for Area
  WD.predA$pred <- P1$fit 
  WD.predA$SE <- P1$se.fit 
  WD.predA$up <- (WD.predA$pred+(1.96 * WD.predA$SE)) 
  WD.predA$low <- (WD.predA$pred-(1.96 * WD.predA$SE))
  WD.predA$resp <- (WD.predA$pred)
  
  
}
# Attach all prediction df's together and export

WD.pred.df <- rbind(WD.predA)


write.csv(WD.pred.df, "Wurmbia_predictions.csv")


####
# Making Figures of Predictions
####

#PLOT ANTTHRO for Area

Fig.WD.PCT <- WD.pred.df %>% 
  filter(Predicting == "PCT.STATE")%>%
  mutate(PCT.STATE2 = fct_relevel(PCT.STATE2, "Forest", "Exotic Woodland", "Derived Grassland", "Woodland"))%>%
  mutate(low = if_else(low < 0, 0, low),
         up = if_else(up >1, 1, up))%>%
  
  ggplot(aes(y=PCT.STATE2))+
  scale_y_discrete(labels = c("Forest", "Exotic \nWoodland", "Derived \nGrassland", "Woodland"))+
  
  geom_pointrange(aes(x=resp, xmin=low, xmax=up))+
  
  theme_bw()+xlab("Occurence")+ylab("PCT Condition State")+
  ggtitle("Wurmbia dioica")
Fig.WD.PCT



# 5.4. MICROTIS SPP. ____________________________####


#BEST MODEL & SUPPORTED MODEL
MS1.1<-glmmTMB(Microt.spp ~ C3_C4 + UMC_ID + 
               (1|Polygon), data = dat, family = binomial)
MS1.2<-glmmTMB(Microt.spp ~ Area + PCT.STATE2 + 
                 (1|Polygon), data = dat, family = binomial)
summary(MS1.1) 
r2_nakagawa(MS1.1)
summary(MS1.2) 
r2_nakagawa(MS1.2)

#check model
simulationOutput <- simulateResiduals(fittedModel = MS1.1, n = 250)
plotSimulatedResiduals(simulationOutput = simulationOutput) #NO issues. 

simulationOutput <- simulateResiduals(fittedModel = MS1.2, n = 250)
plotSimulatedResiduals(simulationOutput = simulationOutput) #NO issues. 

####
#Predicitions 
####

#PREDICTIONS
#Make Dataframes
MS.predA <- data.frame(Area = unique(dat$Area),         #For predicting based on Area
                       Model = "Microt",
                       Predicting = "Area",
                       Polygon = NA,
                       C3_C4 = NA,
                       UMC_ID = NA,
                       PCT.STATE2 = "Woodland")
MS.predB1 <- data.frame(Area = NA,                   #For predicting based on PCT State (MULL)
                        Model = "Microt",
                        Predicting = "C3_C4",
                        Polygon = NA,
                        C3_C4 = unique(dat$C3_C4),
                        UMC_ID = "u19",
                        PCT.STATE2 = NA)
MS.predC1 <- data.frame(Area = "MULL",                   #For predicting based on PCT State (MULL)
                        Model = "Arthro",
                        Predicting = "PCT.STATE",
                        Polygon = NA,
                        C3_C4 = NA,
                        UMC_ID = NA,
                        PCT.STATE2 = unique(dat$PCT.STATE2))
MS.predC2 <- data.frame(Area = "GOO",                   #For predicting based on PCT State (GOO)
                        Model = "Arthro",
                        Predicting = "PCT.STATE",
                        Polygon = NA,
                        C3_C4 = NA,
                        UMC_ID = NA,
                        PCT.STATE2 = unique(dat$PCT.STATE2))


# Generate predictions from best model
{
  P1 <- predict(MS1.2, newdata = MS.predA, type = "response", se.fit = TRUE, re.form = NA)  #Predictions for Area
  MS.predA$pred <- P1$fit 
  MS.predA$SE <- P1$se.fit 
  MS.predA$up <- (MS.predA$pred+(1.96 * MS.predA$SE)) 
  MS.predA$low <- (MS.predA$pred-(1.96 * MS.predA$SE))
  MS.predA$resp <- (MS.predA$pred)
  
  P1 <- predict(MS1.1, newdata = MS.predB1, type = "response", se.fit = TRUE, re.form = NA)  #Predictions for PCT STATE "MULL"
  MS.predB1$pred <- P1$fit 
  MS.predB1$SE <- P1$se.fit 
  MS.predB1$up <- (MS.predB1$pred+(1.96 * MS.predB1$SE)) 
  MS.predB1$low <- (MS.predB1$pred-(1.96 * MS.predB1$SE))
  MS.predB1$resp <- (MS.predB1$pred)
  

}
# Attach all prediction df's together and export

MS.pred.df <- rbind(MS.predA,
                    MS.predB1)


write.csv(MS.pred.df, "Microtis_predictions.csv")

####
# Making Figures of Predictions
####

#PLOT ANTTHRO for Area

Fig.MS.Area <- MS.pred.df %>% 
  filter(Predicting == "Area") %>%
  mutate(Area = fct_relevel(Area, "OFF", "MULL_OUT", "GOO", "MULL"))%>%
  mutate(low = if_else(low < 0, 0, low),
         up = if_else(up >1, 1, up))%>%
  
  ggplot(aes(y=Area))+
  scale_y_discrete(labels = c("Offsets", "Mulligans (outside)", "Goorooyarroo", "Mulligans (inside)"))+
  
  geom_pointrange(aes(x=resp, xmin=low, xmax=up)) +
  
  theme_bw()+xlab("Occurence")+ylab("Area of the MFWS")+
  ggtitle("Microtis spp.")
Fig.MS.Area


Fig.MS.C3C4 <- MS.pred.df %>% 
  filter(Predicting == "C3_C4") %>%
  mutate(Area = fct_relevel(C3_C4, "Exotic_C3", "Native_C3", "Native_C4"))%>%
  mutate(low = if_else(low < 0, 0, low))%>%
  
  ggplot(aes(y=C3_C4))+
  scale_y_discrete(labels = c("Exotic C3", "Native C3", "Native C4"))+
  
  geom_pointrange(aes(x=resp, xmin=low, xmax=up)) +
  
  theme_bw()+xlab("Occurence")+ylab("Mapped grass type")+
  ggtitle("Microtis spp.")

Fig.MS.C3C4

# 5.6. CRASPEDIA SPP. ___________________________####


#BEST MODEL & SUPPORTED MODEL
CS1<-glmmTMB(Cras.spp ~ C3_C4 +
                 (1|Polygon), data = dat, family = binomial)

summary(CS1) 
r2_nakagawa(CS1)


#check model
simulationOutput <- simulateResiduals(fittedModel = CS1, n = 250)
plotSimulatedResiduals(simulationOutput = simulationOutput) #NO issues. 


####
#Predicitions 
####

#PREDICTIONS
#Make Dataframes

CS.predB1 <- data.frame(Area = NA,                   #For predicting based on PCT State (MULL)
                        Model = "Crasp",
                        Predicting = "C3_C4",
                        Polygon = NA,
                        C3_C4 = unique(dat$C3_C4),
                        PCT.STATE2 = NA)


# Generate predictions from best model
{
  
  
  P1 <- predict(CS1, newdata = CS.predB1, type = "response", se.fit = TRUE, re.form = NA)  #Predictions for PCT STATE "MULL"
  CS.predB1$pred <- P1$fit 
  CS.predB1$SE <- P1$se.fit 
  CS.predB1$up <- (CS.predB1$pred+(1.96 * CS.predB1$SE)) 
  CS.predB1$low <- (CS.predB1$pred-(1.96 * CS.predB1$SE))
  CS.predB1$resp <- (CS.predB1$pred)
  
  
}
# Attach all prediction df's together and export

CS.pred.df <- rbind(
                    CS.predB1)


write.csv(CS.pred.df, "Craspedia_predictions.csv")

####
# Making Figures of Predictions
####

#PLOT 


Fig.CS.C3C4 <- CS.pred.df %>% 
  filter(Predicting == "C3_C4") %>%
  mutate(Area = fct_relevel(C3_C4, "Exotic_C3", "Native_C3", "Native_C4"))%>%
  mutate(low = if_else(low < 0, 0, low))%>%
  
  ggplot(aes(y=C3_C4))+
  scale_y_discrete(labels = c("Exotic C3", "Native C3", "Native C4"))+
  
  geom_pointrange(aes(x=resp, xmin=low, xmax=up)) +
  
  theme_bw()+xlab("Occurence")+ylab("Mapped grass type")+
  ggtitle("Craspedia spp.")

Fig.CS.C3C4


# 5.7. LEPTORHYNCHOS SQUAMATUS __________________####


#BEST MODEL 
LS1<-glmmTMB(Lept.squ ~ Area + C3_C4 + 
                 (1|Polygon), data = dat, family = binomial)

summary(LS1) 
r2_nakagawa(LS1)


#check model
simulationOutput <- simulateResiduals(fittedModel = LS1, n = 250)
plotSimulatedResiduals(simulationOutput = simulationOutput) #NO issues. 



####
#Predicitions 
####

#PREDICTIONS
#Make Dataframes
LS.predA <- data.frame(Area = unique(dat$Area),         #For predicting based on Area
                       Model = "Lepto",
                       Predicting = "Area",
                       Polygon = NA,
                       C3_C4 = "Native_C4")
LS.predB1 <- data.frame(Area = "MULL",                   #For predicting based on grass type (MULL)
                        Model = "Lepto",
                        Predicting = "C3_C4",
                        Polygon = NA,
                        C3_C4 = unique(dat$C3_C4))
LS.predB2 <- data.frame(Area = "GOO",                   #For predicting based on grass type (GOO)
                        Model = "Lepto",
                        Predicting = "C3_C4",
                        Polygon = NA,
                        C3_C4 = unique(dat$C3_C4))


# Generate predictions from best model
{
  P1 <- predict(LS1, newdata = LS.predA, type = "response", se.fit = TRUE, re.form = NA)  #Predictions for Area
  LS.predA$pred <- P1$fit 
  LS.predA$SE <- P1$se.fit 
  LS.predA$up <- (LS.predA$pred+(1.96 * LS.predA$SE)) 
  LS.predA$low <- (LS.predA$pred-(1.96 * LS.predA$SE))
  LS.predA$resp <- (LS.predA$pred)
  
  P1 <- predict(LS1, newdata = LS.predB1, type = "response", se.fit = TRUE, re.form = NA)  #Predictions for PCT STATE "MULL"
  LS.predB1$pred <- P1$fit 
  LS.predB1$SE <- P1$se.fit 
  LS.predB1$up <- (LS.predB1$pred+(1.96 * LS.predB1$SE)) 
  LS.predB1$low <- (LS.predB1$pred-(1.96 * LS.predB1$SE))
  LS.predB1$resp <- (LS.predB1$pred)
  
  P1 <- predict(LS1, newdata = LS.predB2, type = "response", se.fit = TRUE, re.form = NA)  #Predictions for PCT STATE "MULL"
  LS.predB2$pred <- P1$fit 
  LS.predB2$SE <- P1$se.fit 
  LS.predB2$up <- (LS.predB2$pred+(1.96 * LS.predB2$SE)) 
  LS.predB2$low <- (LS.predB2$pred-(1.96 * LS.predB2$SE))
  LS.predB2$resp <- (LS.predB2$pred)
  
  
}
# Attach all prediction df's together and export

LS.pred.df <- rbind(LS.predA,
                    LS.predB1, LS.predB2)


write.csv(LS.pred.df, "Lepto_predictions.csv")

####
# Making Figures of Predictions
####

#PLOT ANTTHRO for Area

Fig.LS.Area <- LS.pred.df %>% 
  filter(Predicting == "Area") %>%
  mutate(Area = fct_relevel(Area, "OFF", "MULL_OUT", "GOO", "MULL"))%>%
  mutate(low = if_else(low < 0, 0, low),
         up = if_else(up >1, 1, up))%>%
  
  ggplot(aes(y=Area))+
  scale_y_discrete(labels = c("Offsets", "Mulligans (outside)", "Goorooyarroo", "Mulligans (inside)"))+
  
  geom_pointrange(aes(x=resp, xmin=low, xmax=up)) +
  
  theme_bw()+xlab("Occurence")+ylab("Area of the MFWS")+
  ggtitle("Leptorhynchos squamatus")
Fig.LS.Area


Fig.LS.C3C4 <- LS.pred.df %>% 
  filter(Predicting == "C3_C4") %>%
  mutate(C3_C4 = fct_relevel(C3_C4, "Exotic_C3", "Native_C3", "Native_C4"))%>%
  mutate(low = if_else(low < 0, 0, low))%>%
  
  ggplot(aes(y=C3_C4))+
  scale_y_discrete(labels = c("Exotic C3", "Native C3", "Native C4"))+
  
  geom_pointrange(aes(x=resp, xmin=low, xmax=up, colour = Area),
                  position=position_dodge(width=0.3)) +
  
  theme_bw()+xlab("Occurence")+ylab("Mapped grass type")+
  ggtitle("Leptorhynchos squamatus")

Fig.LS.C3C4


# 5.8. CHRYSOCEPHALUM APICULATUM ________________####


#BEST MODEL 
CA1<-glmmTMB(Chry.api ~ Area +
               (1|Polygon), data = dat, family = binomial)

summary(CA1) 
r2_nakagawa(CA1)


#check model
simulationOutput <- simulateResiduals(fittedModel = CA1, n = 250)
plotSimulatedResiduals(simulationOutput = simulationOutput) #NO issues. 



####
#Predicitions 
####

#PREDICTIONS
#Make Dataframes
CA.predA <- data.frame(Area = unique(dat$Area),         #For predicting based on Area
                       Model = "Chrys",
                       Predicting = "Area",
                       Polygon = NA)



# Generate predictions from best model
{
  P1 <- predict(CA1, newdata = CA.predA, type = "response", se.fit = TRUE, re.form = NA)  #Predictions for Area
  CA.predA$pred <- P1$fit 
  CA.predA$SE <- P1$se.fit 
  CA.predA$up <- (CA.predA$pred+(1.96 * CA.predA$SE)) 
  CA.predA$low <- (CA.predA$pred-(1.96 * CA.predA$SE))
  CA.predA$resp <- (CA.predA$pred)
  
  
}
# Attach all prediction df's together and export

CA.pred.df <- rbind(CA.predA)


write.csv(CA.pred.df, "Chrys_predictions.csv")

####
# Making Figures of Predictions
####

#PLOT ANTTHRO for Area

Fig.CA.Area <- CA.pred.df %>% 
  filter(Predicting == "Area") %>%
  mutate(Area = fct_relevel(Area, "OFF", "MULL_OUT", "GOO", "MULL"))%>%
  mutate(low = if_else(low < 0, 0, low),
         up = if_else(up >1, 1, up))%>%
  
  ggplot(aes(y=Area))+
  scale_y_discrete(labels = c("Offsets", "Mulligans (outside)", "Goorooyarroo", "Mulligans (inside)"))+
  
  geom_pointrange(aes(x=resp, xmin=low, xmax=up)) +
  
  theme_bw()+xlab("Occurrence")+ylab("Area of the MFWS")+
  ggtitle("Chrysocephalum apiculatum")
Fig.CA.Area


# 6. JUST MULLIGANS DATA

# Are there experimental treatment effects (LG, HG, BF) just in the Mulligans data?

# Load in MULLonly data that has more "treatment" variables to explore

mdat <- read_csv("MULLonly-RapidFloristicData-2021.csv", 
                 col_types = cols(Date = col_date(format = "%d/%m/%Y"), 
                                  Start.time = col_time(format = "%H:%M"), 
                                  End.time = col_time(format = "%H:%M")))
str(mdat)


#6.1  MULLonly INDICATOR RICHNESS________________####

# “Indicator richness” is the count of indicator groups present at a site. 
# There were 5 groups (Lilies, Orchids, Daisies [tuberous], Daisies [control], and Stachhousia).

ggplot(mdat, aes(INDI.RICH))+  #check distribution of response variable
  geom_histogram()            #All should poisson (b/c it's count data)
#unless there's overdispersion

####
# Modelling
####

# Manual comparison of a few model structures

M.ir.1 <- glmmTMB(INDI.RICH ~ Exp.Treat + (1|Polygon), data = mdat, family = poisson)
M.ir.2 <- glmmTMB(INDI.RICH ~ Grazing + (1|Polygon), data = mdat, family = poisson)
M.ir.3 <- glmmTMB(INDI.RICH ~ Bet.Fence + (1|Polygon), data = mdat, family = poisson)
M.ir.4 <- glmmTMB(INDI.RICH ~ Grazing * Bet.Fence + (1|Polygon), data = mdat, family = poisson)
M.ir.5 <- glmmTMB(INDI.RICH ~ Flop.Fence + (1|Polygon), data = mdat, family = poisson)
M.ir.6 <- glmmTMB(INDI.RICH ~ Grazing * Flop.Fence + (1|Polygon), data = mdat, family = poisson)
M.ir.0 <- glmmTMB(INDI.RICH ~ 1 + (1|Polygon), data = mdat, family = poisson)

M.ir.AICc <- AICc(M.ir.1,M.ir.2,M.ir.3,M.ir.4,M.ir.5,M.ir.6,M.ir.0)
# Null is best, models 3, 3, and 4 are supported.
# Look at 4 as it's got both terms

summary(M.ir.4)
r2_nakagawa(M.ir.4)

#check model
simulationOutput <- simulateResiduals(fittedModel = M.ir.4, n = 250)
plotSimulatedResiduals(simulationOutput = simulationOutput) #NO issues. 

####
#Predicitions 
####

#PREDICTIONS
#Make Dataframes
m.ir.predA <- data.frame(Grazing = unique(mdat$Grazing), #For predicting based on Grazing
                       Bet.Fence = "Bettong fence",         
                       Model = "Indicator Richness",
                       Predicting = "Both",
                       Polygon = NA)
m.ir.predB <- data.frame(Grazing = unique(mdat$Grazing), #For predicting based on Grazing
                         Bet.Fence = "No fence",         
                         Model = "Indicator Richness",
                         Predicting = "Both",
                         Polygon = NA)
m.ir.pred <- rbind(m.ir.predA, m.ir.predB) #Combine the two frames

# Make predictions

P1 <- predict(M.ir.4, newdata = m.ir.pred, type = "response", se.fit = TRUE, re.form = NA)  #Predictions for Area
m.ir.pred$pred <- P1$fit 
m.ir.pred$SE <- P1$se.fit 
m.ir.pred$up <- (m.ir.pred$pred+(1.96 * m.ir.pred$SE)) 
m.ir.pred$low <- (m.ir.pred$pred-(1.96 * m.ir.pred$SE))
m.ir.pred$resp <- (m.ir.pred$pred)

# Plot predictions

Fig.m.ir <- m.ir.pred %>% 
  
  ggplot(aes(y=Grazing))+
  
  geom_pointrange(aes(x=resp, xmin=low, xmax=up, colour = Bet.Fence),
                  position=position_dodge(width=0.3)) +
  
  theme_bw()+xlab("Indicator Group Richness")+ylab("Grazing Treatment")+
  ggtitle("Mulligans only")
Fig.m.ir

# 
#6.2  MULLonly SPECIES RICHNESS__________________####


ggplot(mdat, aes(SPP.RICH))+  #check distribution of response variable
  geom_histogram()            #All should poisson (b/c it's count data)
#unless there's overdispersion

####
# Modelling
####

# Manual comparison of a few model structures

M.sr.1 <- glmmTMB(SPP.RICH ~ Exp.Treat + (1|Polygon), data = mdat, family = poisson)
M.sr.2 <- glmmTMB(SPP.RICH ~ Grazing + (1|Polygon), data = mdat, family = poisson)
M.sr.3 <- glmmTMB(SPP.RICH ~ Bet.Fence + (1|Polygon), data = mdat, family = poisson)
M.sr.4 <- glmmTMB(SPP.RICH ~ Grazing * Bet.Fence + (1|Polygon), data = mdat, family = poisson)
M.sr.5 <- glmmTMB(SPP.RICH ~ Flop.Fence + (1|Polygon), data = mdat, family = poisson)
M.sr.6 <- glmmTMB(SPP.RICH ~ Grazing * Flop.Fence + (1|Polygon), data = mdat, family = poisson)
M.sr.0 <- glmmTMB(SPP.RICH ~ 1 + (1|Polygon), data = mdat, family = poisson)

M.sr.AICc <- AICc(M.sr.1,M.sr.2,M.sr.3,M.sr.4,M.sr.5,M.sr.6,M.sr.0)
# Null is best, models 2, 3, and 4 are supported.
# Look at 4 as it's got both terms

summary(M.sr.4)
r2_nakagawa(M.sr.4)

#check model
simulationOutput <- simulateResiduals(fittedModel = M.sr.4, n = 250)
plotSimulatedResiduals(simulationOutput = simulationOutput) #NO issues. 

####
#Predicitions 
####

#PREDICTIONS
#Make Dataframes
m.sr.predA <- data.frame(Grazing = unique(mdat$Grazing), #For predicting based on Grazing
                         Bet.Fence = "Bettong fence",         
                         Model = "Indicator Richness",
                         Predicting = "Both",
                         Polygon = NA)
m.sr.predB <- data.frame(Grazing = unique(mdat$Grazing), #For predicting based on Grazing
                         Bet.Fence = "No fence",         
                         Model = "Indicator Richness",
                         Predicting = "Both",
                         Polygon = NA)
m.sr.pred <- rbind(m.sr.predA, m.sr.predB) #Combine the two frames

# Make predictions

P1 <- predict(M.sr.4, newdata = m.sr.pred, type = "response", se.fit = TRUE, re.form = NA)  #Predictions for Area
m.sr.pred$pred <- P1$fit 
m.sr.pred$SE <- P1$se.fit 
m.sr.pred$up <- (m.sr.pred$pred+(1.96 * m.sr.pred$SE)) 
m.sr.pred$low <- (m.sr.pred$pred-(1.96 * m.sr.pred$SE))
m.sr.pred$resp <- (m.sr.pred$pred)

# Plot predictions

Fig.m.sr <- m.sr.pred %>% 
  
  ggplot(aes(y=Grazing))+
  
  geom_pointrange(aes(x=resp, xmin=low, xmax=up, colour = Bet.Fence),
                  position=position_dodge(width=0.3)) +
  
  theme_bw()+xlab("Species Richness")+ylab("Grazing Treatment")+
  ggtitle("Mulligans only")
Fig.m.sr

#6.3  MULLonly MAX INDICATOR ABUNDANCE___________####


ggplot(mdat, aes(MAX.ABUND))+  #check distribution of response variable
  geom_histogram()            #All should poisson (b/c it's count data)
#unless there's overdispersion

####
# Modelling
####

# Manual comparison of a few model structures

M.ma.1 <- glmmTMB(MAX.ABUND ~ Exp.Treat + (1|Polygon), data = mdat, family = poisson)
M.ma.2 <- glmmTMB(MAX.ABUND ~ Grazing + (1|Polygon), data = mdat, family = poisson)
M.ma.3 <- glmmTMB(MAX.ABUND ~ Bet.Fence + (1|Polygon), data = mdat, family = poisson)
M.ma.4 <- glmmTMB(MAX.ABUND ~ Grazing * Bet.Fence + (1|Polygon), data = mdat, family = poisson)
M.ma.5 <- glmmTMB(MAX.ABUND ~ Flop.Fence + (1|Polygon), data = mdat, family = poisson)
M.ma.6 <- glmmTMB(MAX.ABUND ~ Grazing * Flop.Fence + (1|Polygon), data = mdat, family = poisson)
M.ma.0 <- glmmTMB(MAX.ABUND ~ 1 + (1|Polygon), data = mdat, family = poisson)

M.ma.AICc <- AICc(M.ma.1,M.ma.2,M.ma.3,M.ma.4,M.ma.5,M.ma.6,M.ma.0)
# Null is best, models 2, 3, and 4 are supported.
# Look at 4 as it's got both terms

summary(M.ma.2)
r2_nakagawa(M.ma.2)

#check model
simulationOutput <- simulateResiduals(fittedModel = M.ma.4, n = 250)
plotSimulatedResiduals(simulationOutput = simulationOutput) #Disperso. 

#No sig effects 

#6.4  MULLonly LILIES ABUNDANCE__________________####

ggplot(mdat, aes(Lilies))+  #check distribution of response variable
  geom_histogram()            #All should poisson (b/c it's count data)
#unless there's overdispersion

####
# Modelling
####

# Manual comparison of a few model structures

M.ll.1 <- glmmTMB(Lilies ~ Exp.Treat + (1|Polygon), data = mdat, family = poisson)
M.ll.2 <- glmmTMB(Lilies ~ Grazing + (1|Polygon), data = mdat, family = poisson)
M.ll.3 <- glmmTMB(Lilies ~ Bet.Fence + (1|Polygon), data = mdat, family = poisson)
M.ll.4 <- glmmTMB(Lilies ~ Grazing * Bet.Fence + (1|Polygon), data = mdat, family = poisson)
M.ll.5 <- glmmTMB(Lilies ~ Flop.Fence + (1|Polygon), data = mdat, family = poisson)
M.ll.6 <- glmmTMB(Lilies ~ Grazing * Flop.Fence + (1|Polygon), data = mdat, family = poisson)
M.ll.0 <- glmmTMB(Lilies ~ 1 + (1|Polygon), data = mdat, family = poisson)

M.ll.AICc <- AICc(M.ll.1,M.ll.2,M.ll.3,M.ll.4,M.ll.5,M.ll.6,M.ll.0)
# Null is best, models 3 also supported.
# 

summary(M.ll.3)
r2_nakagawa(M.ll.3)

#check model
simulationOutput <- simulateResiduals(fittedModel = M.ll.3, n = 250)
plotSimulatedResiduals(simulationOutput = simulationOutput) #NO issues. 
#No sig effects 

#6.5  MULLonly ORCHIDS ABUNDANCE_________________####

ggplot(mdat, aes(Orchids))+  #check distribution of response variable
  geom_histogram()            #All should poisson (b/c it's count data)
#unless there's overdispersion

####
# Modelling
####

mdat<-mdat%>%                  # actually sig effects here so lets relevel to "floppy" as the base so it's being compared to all
  mutate(Flop.Fence = fct_relevel(Flop.Fence, "Floppy bettong fence"))

# Manual comparison of a few model structures

M.or.1 <- glmmTMB(Orchids ~ Exp.Treat + (1|Polygon), data = mdat, family = poisson)
M.or.2 <- glmmTMB(Orchids ~ Grazing + (1|Polygon), data = mdat, family = poisson)
M.or.3 <- glmmTMB(Orchids ~ Bet.Fence + (1|Polygon), data = mdat, family = poisson)
M.or.4 <- glmmTMB(Orchids ~ Grazing * Bet.Fence + (1|Polygon), data = mdat, family = poisson)
M.or.5 <- glmmTMB(Orchids ~ Flop.Fence + (1|Polygon), data = mdat, family = poisson)
M.or.6 <- glmmTMB(Orchids ~ Grazing * Flop.Fence + (1|Polygon), data = mdat, family = poisson)
M.or.0 <- glmmTMB(Orchids ~ 1 + (1|Polygon), data = mdat, family = poisson)

M.or.AICc <- AICc(M.or.1,M.or.2,M.or.3,M.or.4,M.or.5,M.or.6,M.or.0)
# Model  1 is best, Model 5 supported
# 

summary(M.or.1)
r2_nakagawa(M.or.1)

summary(M.or.5)
r2_nakagawa(M.or.5)

#check model
simulationOutput <- simulateResiduals(fittedModel = M.or.1, n = 250)
plotSimulatedResiduals(simulationOutput = simulationOutput) #NO issues. 
simulationOutput <- simulateResiduals(fittedModel = M.or.5, n = 250)
plotSimulatedResiduals(simulationOutput = simulationOutput) #NO issues. 


####
#Predicitions 
####

#PREDICTIONS
#Make Dataframes
m.or.predA <- data.frame(Exp.Treat = unique(mdat$Exp.Treat), #For predicting based on Grazing
                         Flop.Fence = NA,         
                         Model = "Mull Orchid",
                         Predicting = "Exp.Treat",
                         Polygon = NA)
m.or.predB <- data.frame(Exp.Treat = NA, #For predicting based on Grazing
                         Flop.Fence = unique(mdat$Flop.Fence),         
                         Model = "Mull Orchid",
                         Predicting = "Flop.Fence",
                         Polygon = NA)


# Make predictions

P1 <- predict(M.or.1, newdata = m.or.predA, type = "response", se.fit = TRUE, re.form = NA)  #Predictions for Area
m.or.predA$pred <- P1$fit 
m.or.predA$SE <- P1$se.fit 
m.or.predA$up <- (m.or.predA$pred+(1.96 * m.or.predA$SE)) 
m.or.predA$low <- (m.or.predA$pred-(1.96 * m.or.predA$SE))
m.or.predA$resp <- (m.or.predA$pred)

P1 <- predict(M.or.5, newdata = m.or.predB, type = "response", se.fit = TRUE, re.form = NA)  #Predictions for Area
m.or.predB$pred <- P1$fit 
m.or.predB$SE <- P1$se.fit 
m.or.predB$up <- (m.or.predB$pred+(1.96 * m.or.predB$SE)) 
m.or.predB$low <- (m.or.predB$pred-(1.96 * m.or.predB$SE))
m.or.predB$resp <- (m.or.predB$pred)


m.or.pred <- rbind(m.or.predA, m.or.predB) #Combine the two frames
# Plot predictions

Fig.m.or1 <- m.or.pred %>%
  filter(Predicting == "Exp.Treat")%>%
  
  ggplot(aes(y=Exp.Treat))+
  
  geom_pointrange(aes(x=resp, xmin=low, xmax=up))+
  scale_x_continuous(limits = c(0,NA),
                       labels = c("0", "<10", "<100", "<1000", "1000+"))+
  
  theme_bw()+xlab("Count of individuals")+ylab("Experimental Treatment")+
  ggtitle("Mulligans only - Orchids Indicator Group")
Fig.m.or1

Fig.m.or2 <- m.or.pred %>%
  filter(Predicting == "Flop.Fence")%>%
  
  ggplot(aes(y=Flop.Fence))+
  
  geom_pointrange(aes(x=resp, xmin=low, xmax=up))+
  scale_x_continuous(limits = c(0,NA),
                     labels = c("0", "<10", "<100", "<1000", "1000+"))+
  
  theme_bw()+xlab("Count of individuals")+ylab("Fencing")+
  ggtitle("Mulligans only - Orchids Indicator Group")
Fig.m.or2

#6.5  MULLonly DAISIES (tuberous) ABUNDANCE______####

ggplot(mdat, aes(Daisies.T))+  #check distribution of response variable
  geom_histogram()            #All should poisson (b/c it's count data)
#unless there's overdispersion

####
# Modelling
####


# Manual comparison of a few model structures

M.dT.1 <- glmmTMB(Daisies.T ~ Exp.Treat + (1|Polygon), data = mdat, family = poisson)
M.dT.2 <- glmmTMB(Daisies.T ~ Grazing + (1|Polygon), data = mdat, family = poisson)
M.dT.3 <- glmmTMB(Daisies.T ~ Bet.Fence + (1|Polygon), data = mdat, family = poisson)
M.dT.4 <- glmmTMB(Daisies.T ~ Grazing * Bet.Fence + (1|Polygon), data = mdat, family = poisson)
M.dT.5 <- glmmTMB(Daisies.T ~ Flop.Fence + (1|Polygon), data = mdat, family = poisson)
M.dT.6 <- glmmTMB(Daisies.T ~ Grazing * Flop.Fence + (1|Polygon), data = mdat, family = poisson)
M.dT.0 <- glmmTMB(Daisies.T ~ 1 + (1|Polygon), data = mdat, family = poisson)

M.dT.AICc <- AICc(M.dT.1,M.dT.2,M.dT.3,M.dT.4,M.dT.5,M.dT.6,M.dT.0)
# Only the null model is supported
# 
#6.5  MULLonly DAISIES (control) ABUNDANCE_______####

ggplot(mdat, aes(Daisies.C))+  #check distribution of response variable
  geom_histogram()            #All should poisson (b/c it's count data)
#unless there's overdispersion

####
# Modelling
####


# Manual comparison of a few model structures

M.dC.1 <- glmmTMB(Daisies.C ~ Exp.Treat + (1|Polygon), data = mdat, family = poisson)
M.dC.2 <- glmmTMB(Daisies.C ~ Grazing + (1|Polygon), data = mdat, family = poisson)
M.dC.3 <- glmmTMB(Daisies.C ~ Bet.Fence + (1|Polygon), data = mdat, family = poisson)
M.dC.4 <- glmmTMB(Daisies.C ~ Grazing * Bet.Fence + (1|Polygon), data = mdat, family = poisson)
M.dC.5 <- glmmTMB(Daisies.C ~ Flop.Fence + (1|Polygon), data = mdat, family = poisson)
M.dC.6 <- glmmTMB(Daisies.C ~ Grazing * Flop.Fence + (1|Polygon), data = mdat, family = poisson)
M.dC.0 <- glmmTMB(Daisies.C ~ 1 + (1|Polygon), data = mdat, family = poisson)

M.dC.AICc <- AICc(M.dC.1,M.dC.2,M.dC.3,M.dC.4,M.dC.5,M.dC.6,M.dC.0)
#  null model is best, but 2,3 and 4 are also supported
# Look at 4 as it's got both terms

summary(M.dC.4)
r2_nakagawa(M.dC.4)
# no sig effects 



#______________________________________#
#                                      #  
#      ###   ##   #   ###     #        #
#      #     # #  #   #  #    #        #
#      ###   #  # #   #   #   #        #
#      #     #   ##   #  #             #
#      ###   #    #   ###     #        #
#                                      #
#                                      #
#___________2022-02-14_________________#