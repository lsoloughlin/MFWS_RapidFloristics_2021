# Rapid Floristic Assessment of the MFWS 2021 Project

# FINAL ANALYSIS

# Luke O'Loughlin
# dates written: 2022-05-19 --> 

# 1. PREPARATION_________________________________####

## total clean up 
rm(list = ls()) 
invisible(gc())

## Load Libraries 
# for this whole script. Maybe.
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
library(vegan)
library(REdaS) #this is for checking the adequacy of your data for PCA
library(corrplot)
library(see)


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
dat <- dat %>%
  mutate(Dom.Grass = ifelse(is.na(Dom.Grass), GROUNDSP1, Dom.Grass)) #Fill in Dom Grass for bettong fences

####
# Summarise PCT ZONES into STATES
###

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

####
# Summarise Dominant Grass into GROUPS
###
a <- dat %>% group_by(Area) %>% count(Dom.Grass)
a<-dat%>%distinct(Dom.Grass)
dat <- dat %>% mutate(DomGrass.sum = if_else(Dom.Grass == "Themeda australis" | Dom.Grass == "Themeda triandra", "Themeda", #Group exotics and natives except T,R,A
                                             if_else(Dom.Grass == "Rytidosperma spp.", "Rytid",
                                                     if_else(str_starts(Dom.Grass, "Aust"), "Stipa",
                                                             if_else(str_starts(Dom.Grass, "Both") |
                                                                       str_starts(Dom.Grass, "Mic")|
                                                                       str_starts(Dom.Grass, "Ari")|
                                                                       str_starts(Dom.Grass, "Ann")|
                                                                       str_starts(Dom.Grass, "Era"), "Other native", "Exotic"
                                                             )))))

a <- dat %>% select(Area, Plot, Dom.Grass, DomGrass.sum)
a <- dat %>% group_by(Area) %>% count(DomGrass.sum)
#                 MF    MF_OUT    GOO   OFF
#Themeda          28      3        9     1       #Themeda australis
#Rytid            5       1        1     1       #Rytidosperma spp.
#Stipa            -       -        7     8       #multiple Austrostipa species
#Other native     5       1        6     4       #5 other natives
#Exotic           8       1        7     16      #8 exotics

#Calculate response variable "SUM.ABUND" an "obs.rand"

dat <- dat %>%
  mutate(SUM.ABUND = Lilies + Orchids + Daisies.T + Daisies.C + Stackhousia,
         obs.rand = 1:112)


# 2. PRINCIPLE COMPONENTS________________________####

# make a dataframe that is just the numerical data that PCA works on

dat.2 <- dat %>% # "dat" is the already loaded full dataframe
  select(Plot, #keep identifyer for linking later
         NatGrass, NatShrub, NatOther,          #columns selected from name.
         Rock, BareGround, Litter, Thatch,      #could also select by position %>% select(35:46)
         ExBroadleaf, ExGrass, Nativeness,      #if easier
         Average.Grass.Height)

# PCA doesn't do NA, and dat.2 should have 24 obs that are no good. Omit NAs

dat.3 <- na.omit(dat.2)
nrow(dat.2) #n = 112
nrow(dat.3) #n = 88

#Just have a squizz at the data, looking at correlationss
dat.cor <- dat.3 %>% select(!Plot) #remove the non-numerical col for this
cor_matrix<-cor(dat.cor, method = "spearman", use="complete.obs") #use spearmans to assess both linear and monotonic relationships
corrplot.mixed(cor(cor_matrix),
               upper = "number",
               lower = "ellipse")

# check the adequacy of the data for component analysis

KMOS(dat.cor) #all good if KMO-Criterion is > 0.5 ........ # JUST!!!
bart_spher(dat.cor) #all goof if P < 0.05 i guess .........# Good 


# do the PCA

dat.pca <- rda(dat.cor, scale = T) #on a dataframe that is only the numeric variables


# Check the explained variance 

summary(dat.pca)$cont   # Eigenvalue should be at least 1 .... we have 4 principle components that meeet that
# Cumulative proportion of variaiton explained should be 80%  ..... not the case, only 73%

# Check what the commponents are

biplot(dat.pca, choices = c(1, 2), scaling = "symmetric", type = c("text", "points"))
# PCA 1 = increasing Av.Grass,Height in the negative direction
#         increasing nativeness / Native Other in the positive direction

# PCA 2 = increasing Nat Grass / thatch in (-) 
#         increasing Bare Ground / Litter / Shrub in (+)

biplot(dat.pca, choices = c(1, 3), scaling = "symmetric", type = c("text", "points"))
# PCA 3 = increasing ExBraodleaf in the negative direction
#         increasing Thatch in the positive direction

biplot(dat.pca, choices = c(1, 4), scaling = "symmetric", type = c("text", "points"))
# PCA 4 = increasing Rock in the negative direction
#         increasing Ex.Broadleaf in the positive direction


# Extract components
scrs <- scores(dat.pca,display=c("sites","species"),choices=c(1,2,3,4));
dat.3 <- cbind(dat.3,scrs$sites)
head(dat.3)

# Join these Principle Components to the response variables dataframe
pc.scores <- dat.3 %>% select(Plot, PC1, PC2, PC3, PC4) 
dat <- left_join(dat, pc.scores)


# 3. QUESTION 1__________________________________####
#                Does the abundance and occurrence of plant indicators
#                differ among the different areas of the Sanctuary?

# Full GLMM models
# To address Question 1 all observations were included (n = 112) 
# with each response fitted with a GLMM that included the additive effects of 
# sanctuary area, condition state, experimental grazing treatment and dominant grass species. 

# Abundance scores and richness response variables were modelled with Poisson error distributions and log-link functions, 
# while occurrence variables were modelled with binomial error distributions and logit link functions. 
# The random effect of “polygon” was included in all models to account for any spatial autocorrelation associated with the paired nature of plots.

# Set reference level for each of the FACTORS
dat<-dat%>%
  mutate(Area = fct_relevel(Area, "MULL"),
         PCT.STATE = fct_relevel(PCT.STATE, "Woodland"),
         UMC_ID = fct_relevel(UMC_ID, "u19"),
         Exp.Treat = fct_relevel(Exp.Treat, "High grazing"),
         C3_C4 = fct_relevel(C3_C4, "Native_C4"),
         DomGrass.sum = fct_relevel(DomGrass.sum, "Themeda"),
         PCT.STATE2 = fct_relevel(PCT.STATE, "Woodland"))

# Full models for all 15 response variables
ir.q1F <- glmmTMB(INDI.RICH ~ Area + PCT.STATE2 + Exp.Treat + DomGrass.sum +
                  (1|Polygon), data = dat, family = poisson) #Indicator Richness
sr.q1F <- glmmTMB(SPP.RICH ~ Area + PCT.STATE2 + Exp.Treat + DomGrass.sum +
                    (1|Polygon), data = dat, family = poisson) #Species Richness
ma.q1F <- glmmTMB(MAX.ABUND ~ Area + PCT.STATE2 + Exp.Treat + DomGrass.sum +
                    (1|Polygon), data = dat, family = poisson) #Highest Indicator Group Abund
sa.q1F <- glmmTMB(SUM.ABUND ~ Area + PCT.STATE2 + Exp.Treat + DomGrass.sum +
                    (1|Polygon), data = dat, family = poisson) #Sum of Indicator Group Abund
ll.q1F <- glmmTMB(Lilies ~ Area + PCT.STATE2 + Exp.Treat + DomGrass.sum +
                    (1|Polygon), data = dat, family = poisson) #LILIES count
or.q1F <- glmmTMB(Orchids ~ Area + PCT.STATE2 + Exp.Treat + DomGrass.sum +
                    (1|Polygon), data = dat, family = poisson) #ORCHIDS count
dT.q1F <- glmmTMB(Daisies.T ~ Area + PCT.STATE2 + Exp.Treat + DomGrass.sum +
                    (1|Polygon), data = dat, family = poisson) #DAISIES (Tuberous) count
dC.q1F <- glmmTMB(Daisies.C ~ Area + PCT.STATE2 + Exp.Treat + DomGrass.sum +
                    (1|Polygon), data = dat, family = poisson) #Daisies (Control) count
as.q1F <- glmmTMB(Arth.spp ~ Area + PCT.STATE2 + Exp.Treat + DomGrass.sum +
                    (1|Polygon), data = dat, family = binomial) #Arthro OCC
bs.q1F <- glmmTMB(Bulb.spp ~ Area + PCT.STATE2 + Exp.Treat + DomGrass.sum +
                    (1|Polygon), data = dat, family = binomial) #Bulbine OCC
wd.q1F <- glmmTMB(Wurm.dio ~ Area + PCT.STATE2 + Exp.Treat + DomGrass.sum +
                    (1|Polygon), data = dat, family = binomial) #Wurmbia OCC
ms.q1F <- glmmTMB(Microt.spp ~ Area + PCT.STATE2 + Exp.Treat + DomGrass.sum +
                    (1|Polygon), data = dat, family = binomial) #Microtis OCC
cs.q1F <- glmmTMB(Cras.spp ~ Area + PCT.STATE2 + Exp.Treat + DomGrass.sum +
                    (1|Polygon), data = dat, family = binomial) #Craspedia OCC
ls.q1F <- glmmTMB(Lept.squ ~ Area + PCT.STATE2 + Exp.Treat + DomGrass.sum +
                    (1|Polygon), data = dat, family = binomial) #Lepto OCC
ca.q1F <- glmmTMB(Chry.api ~ Area + PCT.STATE2 + Exp.Treat + DomGrass.sum +
                    (1|Polygon), data = dat, family = binomial) #Chrys OCC

#Model Selection  on all 15 Full model
options(na.action = "na.fail")

ir.q1D <- dredge(ir.q1F)
sr.q1D <- dredge(sr.q1F)
ma.q1D <- dredge(ma.q1F)
sa.q1D <- dredge(sa.q1F)
ll.q1D <- dredge(ll.q1F)
or.q1D <- dredge(or.q1F)
dT.q1D <- dredge(dT.q1F)
dC.q1D <- dredge(dC.q1F)
as.q1D <- dredge(as.q1F)
bs.q1D <- dredge(bs.q1F)
wd.q1D <- dredge(wd.q1F)
ms.q1D <- dredge(ms.q1F)
cs.q1D <- dredge(cs.q1F)
ls.q1D <- dredge(ls.q1F)
ca.q1D <- dredge(ca.q1F)

# subset best models on 15 response variables

ir.q1S <- subset(ir.q1D, delta<=2,  recalc.weights = FALSE)
sr.q1S <- subset(sr.q1D, delta<=2,  recalc.weights = FALSE)
ma.q1S <- subset(ma.q1D, delta<=2,  recalc.weights = FALSE)
sa.q1S <- subset(sa.q1D, delta<=2,  recalc.weights = FALSE)
ll.q1S <- subset(ll.q1D, delta<=2,  recalc.weights = FALSE)
or.q1S <- subset(or.q1D, delta<=2,  recalc.weights = FALSE)
dT.q1S <- subset(dT.q1D, delta<=2,  recalc.weights = FALSE)
dC.q1S <- subset(dC.q1D, delta<=2,  recalc.weights = FALSE)
as.q1S <- subset(as.q1D, delta<=2,  recalc.weights = FALSE)
bs.q1S <- subset(bs.q1D, delta<=2,  recalc.weights = FALSE)
wd.q1S <- subset(wd.q1D, delta<=2,  recalc.weights = FALSE)
ms.q1S <- subset(ms.q1D, delta<=2,  recalc.weights = FALSE)
cs.q1S <- subset(cs.q1D, delta<=2,  recalc.weights = FALSE)
ls.q1S <- subset(ls.q1D, delta<=2,  recalc.weights = FALSE)
ca.q1S <- subset(ca.q1D, delta<=2,  recalc.weights = FALSE)

# what are the best / supported models (delta AICc < 2) on 15 response variables?

ir.q1S # Area + DomGrass + PCT
sr.q1S # Area + DomGrass + PCT
ma.q1S # Area + PCT
sa.q1S # Area + DomGrass + PCT
ll.q1S # Area + DomGrass + PCT
or.q1S # Area + DomGrass + PCT
dT.q1S # DomGrass
dC.q1S # Area
as.q1S # Area + PCT
bs.q1S # Area (null also supported) 
wd.q1S # Area + DomGrass
ms.q1S # DomGrass + PCT
cs.q1S # DomGrass + PCT
ls.q1S # Area + DomGrass
ca.q1S # DomGrass  

options(na.action = "na.omit") #change NA default back

# Best models for all 15 response variables
ir.q1B <- glmmTMB(INDI.RICH ~ Area + PCT.STATE2 +  DomGrass.sum +
                    (1|Polygon) + (1|obs.rand), data = dat, family = poisson) #Indicator Richness
sr.q1B <- glmmTMB(SPP.RICH ~ Area + PCT.STATE2 +  DomGrass.sum +
                    (1|Polygon) + (1|obs.rand), data = dat, family = poisson) #Species Richness
ma.q1B <- glmmTMB(MAX.ABUND ~ Area + PCT.STATE2 +  
                    (1|Polygon) + (1|obs.rand), data = dat, family = poisson) #Highest Indicator Group Abund
sa.q1B <- glmmTMB(SUM.ABUND ~ Area + PCT.STATE2 +  DomGrass.sum +
                    (1|Polygon) + (1|obs.rand), data = dat, family = poisson) #Sum of Indicator Group Abund
ll.q1B <- glmmTMB(Lilies ~ Area + PCT.STATE2 +  DomGrass.sum +
                    (1|Polygon), data = dat, family = poisson) #LILIES count
or.q1B <- glmmTMB(Orchids ~ Area + PCT.STATE2 +  DomGrass.sum +
                    (1|Polygon), data = dat, family = poisson) #ORCHIDS count
dT.q1B <- glmmTMB(Daisies.T ~    DomGrass.sum +
                    (1|Polygon), data = dat, family = poisson) #DAISIES (Tuberous) count
dC.q1B <- glmmTMB(Daisies.C ~ Area +  
                    (1|Polygon), data = dat, family = poisson) #Daisies (Control) count
as.q1B <- glmmTMB(Arth.spp ~ Area + PCT.STATE2 +  
                    (1|Polygon), data = dat, family = binomial) #Arthro OCC
bs.q1B <- glmmTMB(Bulb.spp ~ Area +   
                    (1|Polygon), data = dat, family = binomial) #Bulbine OCC
wd.q1B <- glmmTMB(Wurm.dio ~ Area + DomGrass.sum +
                    (1|Polygon), data = dat, family = binomial) #Wurmbia OCC
ms.q1B <- glmmTMB(Microt.spp ~  PCT.STATE2 +  DomGrass.sum +
                    (1|Polygon), data = dat, family = binomial) #Microtis OCC
cs.q1B <- glmmTMB(Cras.spp ~  PCT.STATE2 +  DomGrass.sum +
                    (1|Polygon), data = dat, family = binomial) #Craspedia OCC
ls.q1B <- glmmTMB(Lept.squ ~ Area +  DomGrass.sum +
                    (1|Polygon), data = dat, family = binomial) #Lepto OCC
ca.q1B <- glmmTMB(Chry.api ~    DomGrass.sum +
                    (1|Polygon), data = dat, family = binomial) #Chrys OCC

# Summary of effects

summary(ir.q1B) # Goo +, PCT and Dom exotic = negs
summary(sr.q1B) # Goo +, PCT and Dom exotic = negs
summary(ma.q1B) # OFF -, PCT exotic = negs
summary(sa.q1B) # GOO +, PCT and Dom exotic = negs
summary(ll.q1B) # GOO +, PCT and Dom exotic = negs
summary(or.q1B) # Just neg effects of dom grass
summary(dT.q1B) # Just neg effects of dom grass
summary(dC.q1B) # Goo +. OFF -
summary(as.q1B) # Goo +, PCT negs
summary(bs.q1B) # none
summary(wd.q1B) # GOO +, dom grass negatives
summary(ms.q1B) # DomGrass + PCT both with some negs
summary(cs.q1B) # derived grassland +
summary(ls.q1B) # model error
summary(ca.q1B) # none

# R2 of best models

r2_nakagawa(ir.q1B) # 66
r2_nakagawa(sr.q1B) # 71
r2_nakagawa(ma.q1B) # 62
r2_nakagawa(sa.q1B) # 75
r2_nakagawa(ll.q1B) # NA, 98
r2_nakagawa(or.q1B) # 98, 97
r2_nakagawa(dT.q1B) # 96, 95
r2_nakagawa(dC.q1B) # 62, 32
r2_nakagawa(as.q1B) # 94, 93
r2_nakagawa(bs.q1B) # 86, 81
r2_nakagawa(wd.q1B) # 70, 43
r2_nakagawa(ms.q1B) # 99, 45
r2_nakagawa(cs.q1B) # 99, 53
r2_nakagawa(ls.q1B) # nnnnnn
r2_nakagawa(ca.q1B) # 99, 05

#confindence intervals
# R2 of best models

confint(ir.q1B) # 66
confint(sr.q1B) # 71
confint(ma.q1B) # 62
confint(sa.q1B) # 75
confint(ll.q1B) # NA, 98
confint(or.q1B) # 98, 97
confint(dT.q1B) # 96, 95
confint(dC.q1B) # 62, 32
confint(as.q1B) # 94, 93
confint(bs.q1B) # 86, 81
confint(wd.q1B) # 70, 43
confint(ms.q1B) # 99, 45
confint(cs.q1B) # 99, 53
confint(ls.q1B) # nnnnnn
confint(ca.q1B) # 99, 05

# AICc of best models

AICc(ir.q1B) # 66
AICc(sr.q1B) # 71
AICc(ma.q1B) # 62
AICc(sa.q1B) # 75
AICc(ll.q1B) # NA, 98
AICc(or.q1B) # 98, 97
AICc(dT.q1B) # 96, 95
AICc(dC.q1B) # 62, 32
AICc(as.q1B) # 94, 93
AICc(bs.q1B) # 86, 81
AICc(wd.q1B) # 70, 43
AICc(ms.q1B) # 99, 45
AICc(cs.q1B) # 99, 53
AICc(ls.q1B) # nnnnnn
AICc(ca.q1B) # 99, 05


# check and correct models as needed

check_model(ir.q1B)
simulationOutput <- simulateResiduals(fittedModel = ir.q1B, n = 250)
plotSimulatedResiduals(simulationOutput = simulationOutput) 

check_model(sr.q1B) 
simulationOutput <- simulateResiduals(fittedModel = sr.q1B, n = 250)
plotSimulatedResiduals(simulationOutput = simulationOutput) 

check_model(ma.q1B)
simulationOutput <- simulateResiduals(fittedModel = ma.q1B, n = 250)
plotSimulatedResiduals(simulationOutput = simulationOutput) 


check_model(sa.q1B)
simulationOutput <- simulateResiduals(fittedModel = sa.q1B, n = 250)
plotSimulatedResiduals(simulationOutput = simulationOutput) 

check_model(ll.q1B)
simulationOutput <- simulateResiduals(fittedModel = ll.q1B, n = 250)
plotSimulatedResiduals(simulationOutput = simulationOutput) # NA, 98

check_model(or.q1B)
simulationOutput <- simulateResiduals(fittedModel = or.q1B, n = 250)
plotSimulatedResiduals(simulationOutput = simulationOutput) # NA, 98

check_model(dT.q1B)
simulationOutput <- simulateResiduals(fittedModel = dT.q1B, n = 250)
plotSimulatedResiduals(simulationOutput = simulationOutput) # 

check_model(dC.q1B)
simulationOutput <- simulateResiduals(fittedModel = dC.q1B, n = 250)
plotSimulatedResiduals(simulationOutput = simulationOutput)

check_model(as.q1B)
simulationOutput <- simulateResiduals(fittedModel = as.q1B, n = 250)
plotSimulatedResiduals(simulationOutput = simulationOutput)

check_model(bs.q1B)
simulationOutput <- simulateResiduals(fittedModel = bs.q1B, n = 250)
plotSimulatedResiduals(simulationOutput = simulationOutput)

check_model(wd.q1B)
simulationOutput <- simulateResiduals(fittedModel = wd.q1B, n = 250)
plotSimulatedResiduals(simulationOutput = simulationOutput)

check_model(ms.q1B) 
simulationOutput <- simulateResiduals(fittedModel = ms.q1B, n = 250)
plotSimulatedResiduals(simulationOutput = simulationOutput)

check_model(cs.q1B)
simulationOutput <- simulateResiduals(fittedModel = cs.q1B, n = 250)
plotSimulatedResiduals(simulationOutput = simulationOutput)

check_model(ls.q1B)
simulationOutput <- simulateResiduals(fittedModel = ls.q1B, n = 250)
plotSimulatedResiduals(simulationOutput = simulationOutput)

check_model(ca.q1B)
simulationOutput <- simulateResiduals(fittedModel = ca.q1B, n = 250)
plotSimulatedResiduals(simulationOutput = simulationOutput)





#   3.1 PREDICTIONS AND FIGURE (Q1)______________####

####
# PANAL 1: INDICATOR RICHNESS
#
#Make Dataframes

ir.pred1 <- data.frame(Area = "MULL",                  
                        Model = "Indi Rich",
                        Polygon = NA,
                        obs.rand = NA,
                        PCT.STATE2 = unique(dat$PCT.STATE2),  
                        DomGrass.sum = "Themeda")
ir.pred2 <- data.frame(Area = "GOO",                   
                        Model = "Indi Rich",
                        Polygon = NA,
                        obs.rand = NA,
                        PCT.STATE2 = unique(dat$PCT.STATE2),  
                        DomGrass.sum = "Themeda")
ir.pred3 <- data.frame(Area = "OFF",                   
                        Model = "Indi Rich",
                        Polygon = NA,
                        obs.rand = NA,
                        PCT.STATE2 = unique(dat$PCT.STATE2),  
                        DomGrass.sum = "Themeda")
ir.pred <- rbind(ir.pred1, ir.pred2, ir.pred3)
ir.pred <- ir.pred%>%filter(PCT.STATE2 != "Forest")

# Make predictions

P1 <- predict(ir.q1B, newdata = ir.pred, type = "response", se.fit = TRUE, re.form = NA)  #Predictions for Area
ir.pred$resp <- P1$fit 
ir.pred$SE <- P1$se.fit 
ir.pred$up <- (ir.pred$resp+(1.96 * ir.pred$SE)) 
ir.pred$low <- (ir.pred$resp-(1.96 * ir.pred$SE))

# Make figure

ir.fig <- ir.pred %>% 
  mutate(PCT.STATE2 = fct_relevel(PCT.STATE2, "Exotic Woodland", "Derived Grassland", "Woodland"),
         Area = fct_relevel(Area, "OFF", "GOO", "MULL"))%>%
  mutate(low = if_else(low < 0, 0, low),
         up = if_else(up >5, 5, up))%>%
  
  ggplot(aes(y=Area))+
  scale_y_discrete(labels = c( "Offsets", "Goorooyarroo", "Mulligans Flat"))+
  
  geom_pointrange(aes(x=resp, xmin=low, xmax=up, colour = PCT.STATE2),
                  position=position_dodge(width=0.3)) +
  
  theme_bw()+xlab("Indicator group richness")+
  theme(axis.text.x = element_text(colour = "black"),
        axis.text.y = element_text(colour = "black"),
        legend.position = "none",
        axis.title.y = element_blank(),
        legend.title = element_blank())

ir.fig

####
# PANAL 2: SPECIES RICHNESS
#
#Make Dataframes

sr.pred1 <- data.frame(Area = "MULL",                  
                       Model = "Spp Rich",
                       Polygon = NA,
                       obs.rand = NA,
                       PCT.STATE2 = unique(dat$PCT.STATE2),  
                       DomGrass.sum = "Themeda")
sr.pred2 <- data.frame(Area = "GOO",                   
                       Model = "Spp Rich",
                       Polygon = NA,
                       obs.rand = NA,
                       PCT.STATE2 = unique(dat$PCT.STATE2),  
                       DomGrass.sum = "Themeda")
sr.pred3 <- data.frame(Area = "OFF",                   
                       Model = "Spp Rich",
                       Polygon = NA,
                       obs.rand = NA,
                       PCT.STATE2 = unique(dat$PCT.STATE2),  
                       DomGrass.sum = "Themeda")
sr.pred <- rbind(sr.pred1, sr.pred2, sr.pred3)
sr.pred <- sr.pred%>%filter(PCT.STATE2 != "Forest")

# Make predictions

P1 <- predict(sr.q1B, newdata = sr.pred, type = "response", se.fit = TRUE, re.form = NA)  #Predictions for Area
sr.pred$resp <- P1$fit 
sr.pred$SE <- P1$se.fit 
sr.pred$up <- (sr.pred$resp+(1.96 * sr.pred$SE)) 
sr.pred$low <- (sr.pred$resp-(1.96 * sr.pred$SE))

# Make figure

sr.fig <- sr.pred %>% 
  mutate(PCT.STATE2 = fct_relevel(PCT.STATE2, "Exotic Woodland", "Derived Grassland", "Woodland"),
         Area = fct_relevel(Area, "OFF", "GOO", "MULL"))%>%
  mutate(low = if_else(low < 0, 0, low),
         up = if_else(up >13, 13, up))%>%
  
  ggplot(aes(y=Area))+
  scale_y_discrete(labels = c( "Offsets", "Goorooyarroo", "Mulligans Flat"))+
  
  geom_pointrange(aes(x=resp, xmin=low, xmax=up, colour = PCT.STATE2),
                  position=position_dodge(width=0.4)) +
  
  theme_bw()+xlab("Target species richness")+
  theme(axis.text.x = element_text(colour = "black"),
        axis.text.y = element_text(colour = "black"),
        legend.position = "none",
        axis.title.y = element_blank(),
        legend.title = element_blank())
sr.fig

####
# PANAL 3: MAX ABUND
#
#Make Dataframes

ma.pred1 <- data.frame(Area = "MULL",                  
                       Model = "Max Abund",
                       Polygon = NA,
                       obs.rand = NA,
                       PCT.STATE2 = unique(dat$PCT.STATE2),  
                       DomGrass.sum = NA)
ma.pred2 <- data.frame(Area = "GOO",                   
                       Model = "Max Abund",
                       Polygon = NA,
                       obs.rand = NA,
                       PCT.STATE2 = unique(dat$PCT.STATE2),  
                       DomGrass.sum = NA)
ma.pred3 <- data.frame(Area = "OFF",                   
                       Model = "Max Abund",
                       Polygon = NA,
                       obs.rand = NA,
                       PCT.STATE2 = unique(dat$PCT.STATE2),  
                       DomGrass.sum = NA)
ma.pred <- rbind(ma.pred1, ma.pred2, ma.pred3)
ma.pred <- ma.pred%>%filter(PCT.STATE2 != "Forest")

# Make predictions

P1 <- predict(ma.q1B, newdata = ma.pred, type = "response", se.fit = TRUE, re.form = NA)  #Predictions for Area
ma.pred$resp <- P1$fit 
ma.pred$SE <- P1$se.fit 
ma.pred$up <- (ma.pred$resp+(1.96 * ma.pred$SE)) 
ma.pred$low <- (ma.pred$resp-(1.96 * ma.pred$SE))

# Make figure

ma.fig <- ma.pred %>% 
  mutate(PCT.STATE2 = fct_relevel(PCT.STATE2, "Exotic Woodland", "Derived Grassland", "Woodland"),
         Area = fct_relevel(Area, "OFF", "GOO", "MULL"))%>%
  mutate(low = if_else(low < 0, 0, low),
         up = if_else(up >4, 4, up))%>%
  
  ggplot(aes(y=Area))+
  scale_y_discrete(labels = c( "Offsets", "Goorooyarroo", "Mulligans Flat"))+
  
  geom_pointrange(aes(x=resp, xmin=low, xmax=up, colour = PCT.STATE2),
                  position=position_dodge(width=0.4)) +
  scale_x_continuous(limits = c(0,NA),
                     labels = c("0", "1 \n1-9", "2 \n10-99", "3 \n100-999"))+
  
  theme_bw()+xlab("Highest indicator abundance score")+
  theme(axis.text.x = element_text(colour = "black"),
        axis.text.y = element_text(colour = "black"),
        legend.position = "none",
        axis.title.y = element_blank(),
        legend.title = element_blank())
ma.fig

####
# PANAL 4: SUM ABUND
#
#Make Dataframes

sa.pred1 <- data.frame(Area = "MULL",                  
                       Model = "Sum Abund",
                       Polygon = NA,
                       obs.rand = NA,
                       PCT.STATE2 = unique(dat$PCT.STATE2),  
                       DomGrass.sum = "Themeda")
sa.pred2 <- data.frame(Area = "GOO",                   
                       Model = "Sum Abund",
                       Polygon = NA,
                       obs.rand = NA,
                       PCT.STATE2 = unique(dat$PCT.STATE2),  
                       DomGrass.sum = "Themeda")
sa.pred3 <- data.frame(Area = "OFF",                   
                       Model = "Sum Abund",
                       Polygon = NA,
                       obs.rand = NA,
                       PCT.STATE2 = unique(dat$PCT.STATE2),  
                       DomGrass.sum = "Themeda")
sa.pred <- rbind(sa.pred1, sa.pred2, sa.pred3)
sa.pred <- sa.pred%>%filter(PCT.STATE2 != "Forest")

# Make predictions

P1 <- predict(sa.q1B, newdata = sa.pred, type = "response", se.fit = TRUE, re.form = NA)  #Predictions for Area
sa.pred$resp <- P1$fit 
sa.pred$SE <- P1$se.fit 
sa.pred$up <- (sa.pred$resp+(1.96 * sa.pred$SE)) 
sa.pred$low <- (sa.pred$resp-(1.96 * sa.pred$SE))

# Make figure

sa.fig <- sa.pred %>% 
  mutate(PCT.STATE2 = fct_relevel(PCT.STATE2, "Exotic Woodland", "Derived Grassland", "Woodland"),
         Area = fct_relevel(Area, "OFF", "GOO", "MULL"))%>%
  mutate(low = if_else(low < 0, 0, low),
         up = if_else(up >13, 13, up))%>%
  
  ggplot(aes(y=Area))+
  scale_y_discrete(labels = c( "Offsets", "Goorooyarroo", "Mulligans Flat"))+
  
  geom_pointrange(aes(x=resp, xmin=low, xmax=up, colour = PCT.STATE2),
                  position=position_dodge(width=0.4)) +
  
  theme_bw()+xlab("Sum of indicator abundance scores")+
  theme(axis.text.x = element_text(colour = "black"),
        axis.text.y = element_text(colour = "black"),
        legend.position = "none",
        axis.title.y = element_blank(),
        legend.title = element_blank())
sa.fig

####
# PANAL 5: LILIES ABUND
#
#Make Dataframes

ll.pred1 <- data.frame(Area = "MULL",                  
                       Model = "Lilies",
                       Polygon = NA,
                       obs.rand = NA,
                       PCT.STATE2 = unique(dat$PCT.STATE2),  
                       DomGrass.sum = "Themeda")
ll.pred2 <- data.frame(Area = "GOO",                   
                       Model = "Lilies",
                       Polygon = NA,
                       obs.rand = NA,
                       PCT.STATE2 = unique(dat$PCT.STATE2),  
                       DomGrass.sum = "Themeda")
ll.pred3 <- data.frame(Area = "OFF",                   
                       Model = "Lilies",
                       Polygon = NA,
                       obs.rand = NA,
                       PCT.STATE2 = unique(dat$PCT.STATE2),  
                       DomGrass.sum = "Themeda")
ll.pred <- rbind(ll.pred1, ll.pred2, ll.pred3)
ll.pred <- ll.pred%>%filter(PCT.STATE2 != "Forest")

# Make predictions

P1 <- predict(ll.q1B, newdata = ll.pred, type = "response", se.fit = TRUE, re.form = NA)  #Predictions for Area
ll.pred$resp <- P1$fit 
ll.pred$SE <- P1$se.fit 
ll.pred$up <- (ll.pred$resp+(1.96 * ll.pred$SE)) 
ll.pred$low <- (ll.pred$resp-(1.96 * ll.pred$SE))

# Make figure

ll.fig <- ll.pred %>% 
  mutate(PCT.STATE2 = fct_relevel(PCT.STATE2, "Exotic Woodland", "Derived Grassland", "Woodland"),
         Area = fct_relevel(Area, "OFF", "GOO", "MULL"))%>%
  mutate(low = if_else(low < 0, 0, low),
         up = if_else(up >4, 4, up))%>%
  
  ggplot(aes(y=Area))+
  scale_y_discrete(labels = c( "Offsets", "Goorooyarroo", "Mulligans Flat"))+
  
  geom_pointrange(aes(x=resp, xmin=low, xmax=up, colour = PCT.STATE2),
                  position=position_dodge(width=0.4)) +
  scale_x_continuous(limits = c(0,NA),
                     labels = c("0", "1 \n1-9", "2 \n10-99", "3 \n100-999",  "4 \n1000+"))+
  
  theme_bw()+xlab("Lilies indicator abundance score")+
  theme(axis.text.x = element_text(colour = "black"),
        axis.text.y = element_text(colour = "black"),
        legend.position = "none",
        axis.title.y = element_blank(),
        legend.title = element_blank())
ll.fig

####
# PANAL 6: DAISIES CONTROL ABUND
#
#Make Dataframes

dC.pred1 <- data.frame(Area = "MULL",                  
                       Model = "Daisies Cont",
                       Polygon = NA,
                       obs.rand = NA,
                       PCT.STATE2 = NA,  
                       DomGrass.sum = NA)
dC.pred2 <- data.frame(Area = "GOO",                   
                       Model = "Daisies Cont",
                       Polygon = NA,
                       obs.rand = NA,
                       PCT.STATE2 = NA,  
                       DomGrass.sum = NA)
dC.pred3 <- data.frame(Area = "OFF",                   
                       Model = "Daisies Cont",
                       Polygon = NA,
                       obs.rand = NA,
                       PCT.STATE2 = NA,  
                       DomGrass.sum = NA)
dC.pred <- rbind(dC.pred1, dC.pred2, dC.pred3)


# Make predictions

P1 <- predict(dC.q1B, newdata = dC.pred, type = "response", se.fit = TRUE, re.form = NA)  #Predictions for Area
dC.pred$resp <- P1$fit 
dC.pred$SE <- P1$se.fit 
dC.pred$up <- (dC.pred$resp+(1.96 * dC.pred$SE)) 
dC.pred$low <- (dC.pred$resp-(1.96 * dC.pred$SE))

# Make figure

dC.fig <- dC.pred %>% 
  mutate(
         Area = fct_relevel(Area, "OFF", "GOO", "MULL"))%>%
  mutate(low = if_else(low < 0, 0, low),
         up = if_else(up >4, 4, up))%>%
  
  ggplot(aes(y=Area))+
  scale_y_discrete(labels = c( "Offsets", "Goorooyarroo", "Mulligans Flat"))+
  
  geom_pointrange(aes(x=resp, xmin=low, xmax=up)) +
  scale_x_continuous(limits = c(0,NA),
                     labels = c("0", "1 \n1-9", "2 \n10-99", "3 \n100-999"))+
  
  theme_bw()+xlab("Daisies (control) indicator abundance score")+
  theme(axis.text.x = element_text(colour = "black"),
        axis.text.y = element_text(colour = "black"),
        legend.position = "none",
        axis.title.y = element_blank(),
        legend.title = element_blank())
dC.fig

####
# JOIN PANNELS INTO 1 FIGURE
####

Q1.FIG <- ggarrange(ir.fig, sr.fig, ma.fig, sa.fig, ll.fig, dC.fig,
                    ncol = 2, nrow = 3,
                    labels = "AUTO",
                    common.legend = TRUE,
                    vjust = 0.5)
Q1.FIG

ggsave("Q1-Figure.tiff",
       plot=Q1.FIG,
       width = 7.7, height = 5.5, units = "in",
       dpi = 300)  

# 4. QUESTION 2__________________________________####
#                does ground cover and vegetation structural attributes strongly 
#               influence plant indicator groups or species? 

# Full GLMM models
# To address Question 2 only observations with veg covariate data were included (n = 88) 
# with each response fitted with a GLMM that included the additive effects of 
# sanctuary area, and the four principle components 

# Abundance scores and richness response variables were modelled with Poisson error distributions and log-link functions, 
# while occurrence variables were modelled with binomial error distributions and logit link functions. 
# The random effect of “polygon” was included in all models to account for any spatial autocorrelation associated with the paired nature of plots.

#subset the data
dat <- dat%>%mutate(Exp.Treat = replace_na(Exp.Treat, "High grazing"))
dat88 <- dat%>%
  filter(Exp.Treat != "Bettong fence FT",
         Exp.Treat != "Bettong fence LG",
         Exp.Treat != "Bettong fence HG")

#scale PC variables
dat88$zPC1 <- scale(dat88$PC1)
dat88$zPC2 <- scale(dat88$PC2)
dat88$zPC3 <- scale(dat88$PC3)
dat88$zPC4 <- scale(dat88$PC4)

# Full models for all 15 response variables
ir.q2F <- glmmTMB(INDI.RICH ~ Area + zPC1 + zPC2 + zPC3 + zPC4 +
                    (1|Polygon), data = dat88, family = poisson) #Indicator Richness
sr.q2F <- glmmTMB(SPP.RICH ~ Area + zPC1 + zPC2 + zPC3 + zPC4 +
                    (1|Polygon), data = dat88, family = poisson) #Species Richness
ma.q2F <- glmmTMB(MAX.ABUND ~ Area + zPC1 + zPC2 + zPC3 + zPC4 +
                    (1|Polygon), data = dat88, family = poisson) #Highest Indicator Group Abund
sa.q2F <- glmmTMB(SUM.ABUND ~ Area + zPC1 + zPC2 + zPC3 + zPC4 +
                    (1|Polygon), data = dat88, family = poisson) #Sum of Indicator Group Abund
ll.q2F <- glmmTMB(Lilies ~ Area + zPC1 + zPC2 + zPC3 + zPC4 +
                    (1|Polygon), data = dat88, family = poisson) #LILIES count
or.q2F <- glmmTMB(Orchids ~ Area + zPC1 + zPC2 + zPC3 + zPC4 +
                    (1|Polygon), data = dat88, family = poisson) #ORCHIDS count
dT.q2F <- glmmTMB(Daisies.T ~ Area + zPC1 + zPC2 + zPC3 + zPC4 +
                    (1|Polygon), data = dat88, family = poisson) #DAISIES (Tuberous) count
dC.q2F <- glmmTMB(Daisies.C ~ Area + zPC1 + zPC2 + zPC3 + zPC4 +
                    (1|Polygon), data = dat88, family = poisson) #Daisies (Control) count
as.q2F <- glmmTMB(Arth.spp ~ Area + zPC1 + zPC2 + zPC3 + zPC4 +
                    (1|Polygon), data = dat88, family = binomial) #Arthro OCC
bs.q2F <- glmmTMB(Bulb.spp ~ Area + zPC1 + zPC2 + zPC3 + zPC4 +
                    (1|Polygon), data = dat88, family = binomial) #Bulbine OCC
wd.q2F <- glmmTMB(Wurm.dio ~ Area + zPC1 + zPC2 + zPC3 + zPC4 +
                    (1|Polygon), data = dat88, family = binomial) #Wurmbia OCC
ms.q2F <- glmmTMB(Microt.spp ~ Area + zPC1 + zPC2 + zPC3 + zPC4 +
                    (1|Polygon), data = dat88, family = binomial) #Microtis OCC
cs.q2F <- glmmTMB(Cras.spp ~ Area + zPC1 + zPC2 + zPC3 + zPC4 +
                    (1|Polygon), data = dat88, family = binomial) #Craspedia OCC
ls.q2F <- glmmTMB(Lept.squ ~ Area + zPC1 + zPC2 + zPC3 + zPC4 +
                    (1|Polygon), data = dat88, family = binomial) #Lepto OCC
ca.q2F <- glmmTMB(Chry.api ~ Area + zPC1 + zPC2 + zPC3 + zPC4 +
                    (1|Polygon), data = dat88, family = binomial) #Chrys OCC

#Model Selection  on all 15 Full model
options(na.action = "na.fail")

ir.q2D <- dredge(ir.q2F)
sr.q2D <- dredge(sr.q2F)
ma.q2D <- dredge(ma.q2F)
sa.q2D <- dredge(sa.q2F)
ll.q2D <- dredge(ll.q2F)
or.q2D <- dredge(or.q2F)
dT.q2D <- dredge(dT.q2F)
dC.q2D <- dredge(dC.q2F)
as.q2D <- dredge(as.q2F)
bs.q2D <- dredge(bs.q2F)
wd.q2D <- dredge(wd.q2F)
ms.q2D <- dredge(ms.q2F)
cs.q2D <- dredge(cs.q2F)
ls.q2D <- dredge(ls.q2F)
ca.q2D <- dredge(ca.q2F)

# subset best models on 15 response variables

ir.q2S <- subset(ir.q2D, delta<=2,  recalc.weights = FALSE)
sr.q2S <- subset(sr.q2D, delta<=2,  recalc.weights = FALSE)
ma.q2S <- subset(ma.q2D, delta<=2,  recalc.weights = FALSE)
sa.q2S <- subset(sa.q2D, delta<=2,  recalc.weights = FALSE)
ll.q2S <- subset(ll.q2D, delta<=2,  recalc.weights = FALSE)
or.q2S <- subset(or.q2D, delta<=2,  recalc.weights = FALSE)
dT.q2S <- subset(dT.q2D, delta<=2,  recalc.weights = FALSE)
dC.q2S <- subset(dC.q2D, delta<=2,  recalc.weights = FALSE)
as.q2S <- subset(as.q2D, delta<=2,  recalc.weights = FALSE)
bs.q2S <- subset(bs.q2D, delta<=2,  recalc.weights = FALSE)
wd.q2S <- subset(wd.q2D, delta<=2,  recalc.weights = FALSE)
ms.q2S <- subset(ms.q2D, delta<=2,  recalc.weights = FALSE)
cs.q2S <- subset(cs.q2D, delta<=2,  recalc.weights = FALSE)
ls.q2S <- subset(ls.q2D, delta<=2,  recalc.weights = FALSE)
ca.q2S <- subset(ca.q2D, delta<=2,  recalc.weights = FALSE)

# what are the best / supported models (delta AICc < 2) on 15 response variables?
# * = OTHER SUPPORTED MODELS TOO

ir.q2S # Area + PC1 + PC3 *
sr.q2S # Area + PC1 + PC3 *
ma.q2S # PC1 + PC2
sa.q2S # Area + PC1 + PC2 + PC3 *
ll.q2S # Area + PC1 *
or.q2S # PC1 + PC2 *
dT.q2S # PC1 *
dC.q2S # Area + PC1 *
as.q2S # Area + PC1 *
bs.q2S # Area + PC1 * 
wd.q2S # Area + PC1 + PC2 + PC3
ms.q2S # PC1 + PC2 *
cs.q2S # PC2 (null also supported)
ls.q2S # Area + PC1 *
ca.q2S # Area + PC1 *

options(na.action = "na.omit") #change NA default back

# Best models for all 15 response variables
ir.q2B <- glmmTMB(INDI.RICH ~ Area + zPC1 + zPC3 +
                    (1|Polygon), data = dat88, family = poisson) #Indicator Richness
sr.q2B <- glmmTMB(SPP.RICH ~ Area + zPC1 + zPC3 +
                    (1|Polygon), data = dat88, family = poisson) #Species Richness
ma.q2B <- glmmTMB(MAX.ABUND ~ zPC1 + zPC2 +  
                    (1|Polygon) + (1|obs.rand), data = dat88, family = poisson) #Highest Indicator Group Abund
sa.q2B <- glmmTMB(SUM.ABUND ~ Area + zPC1 + zPC2 + zPC3 +
                    (1|Polygon), data = dat88, family = poisson) #Sum of Indicator Group Abund
ll.q2B <- glmmTMB(Lilies ~ Area + zPC1 +
                    (1|Polygon), data = dat88, family = poisson) #LILIES count
or.q2B <- glmmTMB(Orchids ~ zPC1 + zPC2 +
                    (1|Polygon), data = dat88, family = poisson) #ORCHIDS count
dT.q2B <- glmmTMB(Daisies.T ~ zPC1 +
                    (1|Polygon), data = dat88, family = poisson) #DAISIES (Tuberous) count
dC.q2B <- glmmTMB(Daisies.C ~ Area + zPC1 +
                    (1|Polygon), data = dat88, family = poisson) #Daisies (Control) count
as.q2B <- glmmTMB(Arth.spp ~ Area + zPC1 +  
                    (1|Polygon), data = dat88, family = binomial) #Arthro OCC
bs.q2B <- glmmTMB(Bulb.spp ~ Area + zPC1  +
                    (1|Polygon), data = dat88, family = binomial) #Bulbine OCC
wd.q2B <- glmmTMB(Wurm.dio ~ Area + zPC1 + zPC2 + zPC3 +
                    (1|Polygon) + (1|obs.rand), data = dat88, family = binomial) #Wurmbia OCC
ms.q2B <- glmmTMB(Microt.spp ~  zPC1 + zPC2 +
                    (1|Polygon), data = dat88, family = binomial) #Microtis OCC
cs.q2B <- glmmTMB(Cras.spp ~  zPC2 +
                    (1|Polygon), data = dat88, family = binomial) #Craspedia OCC
ls.q2B <- glmmTMB(Lept.squ ~ Area + zPC1 +
                    (1|Polygon), data = dat88, family = binomial) #Lepto OCC
ca.q2B <- glmmTMB(Chry.api ~  Area + zPC1 +
                    (1|Polygon), data = dat88, family = binomial) #Chrys OCC

# Summary of effects

summary(ir.q2B) # Goo +, PC1 +
summary(sr.q2B) # Goo +, PC1 +, PC3 +
summary(ma.q2B) # PC1 +, PC2 -
summary(sa.q2B) # GOO +, PC1 +, PC2 -, PC3+
summary(ll.q2B) # GOO +, PC1 +
summary(or.q2B) # PC1 +, PC2 -
summary(dT.q2B) # PC1 +
summary(dC.q2B) # Goo +. PC1 +
summary(as.q2B) # Goo +, PC1 +
summary(bs.q2B) # none
summary(wd.q2B) # GOO +, PC1 +, PC2 -, PC3+
summary(ms.q2B) # PC1 +, PC2 -
summary(cs.q2B) # 
summary(ls.q2B) # Goo +. PC1 +
summary(ca.q2B) # Goo+

# R2 of best models

r2_nakagawa(ir.q2B) # 59
r2_nakagawa(sr.q2B) # 63
r2_nakagawa(ma.q2B) # 62
r2_nakagawa(sa.q2B) # 75
r2_nakagawa(ll.q2B) # NA, 98
r2_nakagawa(or.q2B) # 98, 97
r2_nakagawa(dT.q2B) # 96, 95
r2_nakagawa(dC.q2B) # 62, 32
r2_nakagawa(as.q2B) # 94, 93
r2_nakagawa(bs.q2B) # 86, 81
r2_nakagawa(wd.q2B) # 70, 43
r2_nakagawa(ms.q2B) # 99, 45
r2_nakagawa(cs.q2B) # 99, 53
r2_nakagawa(ls.q2B) # nnnnnn
r2_nakagawa(ca.q2B) # 99, 05


# AICc of best models

AICc(ir.q2B) # 59
AICc(sr.q2B) # 63
AICc(ma.q2B) # 62
AICc(sa.q2B) # 75
AICc(ll.q2B) # NA, 98
AICc(or.q2B) # 98, 97
AICc(dT.q2B) # 96, 95
AICc(dC.q2B) # 62, 32
AICc(as.q2B) # 94, 93
AICc(bs.q2B) # 86, 81
AICc(wd.q2B) # 70, 43
AICc(ms.q2B) # 99, 45
AICc(cs.q2B) # 99, 53
AICc(ls.q2B) # nnnnnn
AICc(ca.q2B) # 99, 05

# check and correct models as needed

check_model(ir.q2B)
simulationOutput <- simulateResiduals(fittedModel = ir.q2B, n = 250)
plotSimulatedResiduals(simulationOutput = simulationOutput) 

check_model(sr.q2B) 
simulationOutput <- simulateResiduals(fittedModel = sr.q2B, n = 250)
plotSimulatedResiduals(simulationOutput = simulationOutput) 

check_model(ma.q2B)
simulationOutput <- simulateResiduals(fittedModel = ma.q2B, n = 250)
plotSimulatedResiduals(simulationOutput = simulationOutput) 


check_model(sa.q2B)
simulationOutput <- simulateResiduals(fittedModel = sa.q2B, n = 250)
plotSimulatedResiduals(simulationOutput = simulationOutput) 

check_model(ll.q2B)
simulationOutput <- simulateResiduals(fittedModel = ll.q2B, n = 250)
plotSimulatedResiduals(simulationOutput = simulationOutput) # NA, 98

check_model(or.q2B)
simulationOutput <- simulateResiduals(fittedModel = or.q2B, n = 250)
plotSimulatedResiduals(simulationOutput = simulationOutput) # NA, 98

check_model(dT.q2B)
simulationOutput <- simulateResiduals(fittedModel = dT.q2B, n = 250)
plotSimulatedResiduals(simulationOutput = simulationOutput) # 

check_model(dC.q2B)
simulationOutput <- simulateResiduals(fittedModel = dC.q2B, n = 250)
plotSimulatedResiduals(simulationOutput = simulationOutput)

check_model(as.q2B)
simulationOutput <- simulateResiduals(fittedModel = as.q2B, n = 250)
plotSimulatedResiduals(simulationOutput = simulationOutput)

check_model(bs.q2B)
simulationOutput <- simulateResiduals(fittedModel = bs.q2B, n = 250)
plotSimulatedResiduals(simulationOutput = simulationOutput)

check_model(wd.q2B)
simulationOutput <- simulateResiduals(fittedModel = wd.q2B, n = 250)
plotSimulatedResiduals(simulationOutput = simulationOutput)

check_model(ms.q2B) 
simulationOutput <- simulateResiduals(fittedModel = ms.q2B, n = 250)
plotSimulatedResiduals(simulationOutput = simulationOutput)

check_model(cs.q2B)
simulationOutput <- simulateResiduals(fittedModel = cs.q2B, n = 250)
plotSimulatedResiduals(simulationOutput = simulationOutput)

check_model(ls.q2B)
simulationOutput <- simulateResiduals(fittedModel = ls.q2B, n = 250)
plotSimulatedResiduals(simulationOutput = simulationOutput)

check_model(ca.q2B)
simulationOutput <- simulateResiduals(fittedModel = ca.q2B, n = 250)
plotSimulatedResiduals(simulationOutput = simulationOutput)



#   4.1 PREDICTIONS AND FIGURE (Q2)______________####

# Determine the scale of PC1 that each area can predict on

x<-dat88%>%group_by(Area)%>%
  summarise(minPC1 = min(PC1),
            maxPC1 = max(PC1))

#Unscaled best model for predictions
ir.q2Bx <- glmmTMB(INDI.RICH ~ Area + PC1 + PC3 +
                    (1|Polygon), data = dat88, family = poisson) #Indicator Richness
sr.q2Bx <- glmmTMB(SPP.RICH ~ Area + PC1 + PC3 +
                    (1|Polygon), data = dat88, family = poisson) #Species Richness
ma.q2Bx <- glmmTMB(MAX.ABUND ~ PC1 + PC2 +  
                    (1|Polygon) + (1|obs.rand), data = dat88, family = poisson) #Highest Indicator Group Abund
sa.q2Bx <- glmmTMB(SUM.ABUND ~ Area + PC1 + PC2 + PC3 +
                    (1|Polygon), data = dat88, family = poisson) #Sum of Indicator Group Abund
ll.q2Bx <- glmmTMB(Lilies ~ Area + PC1 +
                    (1|Polygon), data = dat88, family = poisson) #LILIES count
or.q2Bx <- glmmTMB(Orchids ~ PC1 + PC2 +
                    (1|Polygon), data = dat88, family = poisson) #ORCHIDS count
dT.q2Bx <- glmmTMB(Daisies.T ~ PC1 +
                    (1|Polygon), data = dat88, family = poisson) #DAISIES (Tuberous) count
dC.q2Bx <- glmmTMB(Daisies.C ~ Area + PC1 +
                    (1|Polygon), data = dat88, family = poisson) #Daisies (Control) count
####
# PANAL 1: INDICATOR RICHNESS
#
#Make Dataframes
ir2.pred1 <- data.frame(Area = "MULL",                  
                       Model = "Indicator richness",
                       Polygon = NA,
                       obs.rand = NA,
                       PC1 = seq(-1.28, 0.86, by = 0.02),  
                       PC2 = NA,
                       PC3 = mean(dat88$PC3),
                       order = "A")
ir2.pred2 <- data.frame(Area = "GOO",                   
                        Model = "Indicator richness",
                        Polygon = NA,
                        obs.rand = NA,
                        PC1 = seq(-0.80, 0.89, by = 0.02),  
                        PC2 = NA,
                        PC3 = mean(dat88$PC3),
                        order = "A")  
ir2.pred3 <- data.frame(Area = "OFF",                   
                         Model = "Indicator richness",
                         Polygon = NA,
                         obs.rand = NA,
                         PC1 = seq(-1.19, 1.05, by = 0.02),  
                         PC2 = NA,
                         PC3 = mean(dat88$PC3),
                        order = "A") 
ir2.pred <- rbind(ir2.pred1, ir2.pred2, ir2.pred3)

# Make predictions

P1 <- predict(ir.q2Bx, newdata = ir2.pred, type = "response", se.fit = TRUE, re.form = NA)  #Predictions for Area
ir2.pred$resp <- P1$fit 
ir2.pred$SE <- P1$se.fit 
ir2.pred$up <- (ir2.pred$resp+(1.96 * ir2.pred$SE)) 
ir2.pred$low <- (ir2.pred$resp-(1.96 * ir2.pred$SE))

####
# PANAL 2: SPECIES RICHNESS
#
#Make Dataframes
sr2.pred1 <- data.frame(Area = "MULL",                  
                        Model = "Target species richness",
                        Polygon = NA,
                        obs.rand = NA,
                        PC1 = seq(-1.28, 0.86, by = 0.02),  
                        PC2 = NA,
                        PC3 = mean(dat88$PC3),
                        order = "B")
sr2.pred2 <- data.frame(Area = "GOO",                   
                        Model = "Target species richness",
                        Polygon = NA,
                        obs.rand = NA,
                        PC1 = seq(-0.80, 0.89, by = 0.02),  
                        PC2 = NA,
                        PC3 = mean(dat88$PC3),
                        order = "B")  
sr2.pred3 <- data.frame(Area = "OFF",                   
                        Model = "Target species richness",
                        Polygon = NA,
                        obs.rand = NA,
                        PC1 = seq(-1.19, 1.05, by = 0.02),  
                        PC2 = NA,
                        PC3 = mean(dat88$PC3),
                        order = "B") 
sr2.pred <- rbind(sr2.pred1, sr2.pred2, sr2.pred3)

# Make predictions

P1 <- predict(sr.q2Bx, newdata = sr2.pred, type = "response", se.fit = TRUE, re.form = NA)  #Predictions for Area
sr2.pred$resp <- P1$fit 
sr2.pred$SE <- P1$se.fit 
sr2.pred$up <- (sr2.pred$resp+(1.96 * sr2.pred$SE)) 
sr2.pred$low <- (sr2.pred$resp-(1.96 * sr2.pred$SE))

####
# PANAL 3: MAX ABUND
#
#Make Dataframes
ma2.pred <- data.frame(Area = NA,                  
                        Model = "Highest indocator abundance score",
                        Polygon = NA,
                        obs.rand = NA,
                        PC1 = seq(-1.28, 1.05, by = 0.02),  
                        PC2 = mean(dat88$PC2),
                        PC3 = NA,
                       order = "C")

# Make predictions

P1 <- predict(ma.q2Bx, newdata = ma2.pred, type = "response", se.fit = TRUE, re.form = NA)  #Predictions for Area
ma2.pred$resp <- P1$fit 
ma2.pred$SE <- P1$se.fit 
ma2.pred$up <- (ma2.pred$resp+(1.96 * ma2.pred$SE)) 
ma2.pred$low <- (ma2.pred$resp-(1.96 * ma2.pred$SE))

####
# PANAL 4: SUM ABUND
#
#Make Dataframes
sa2.pred1 <- data.frame(Area = "MULL",                  
                        Model = "Sum of indocator abundance scores",
                        Polygon = NA,
                        obs.rand = NA,
                        PC1 = seq(-1.28, 0.86, by = 0.02),  
                        PC2 = mean(dat88$PC2),
                        PC3 = mean(dat88$PC3),
                        order = "D")
sa2.pred2 <- data.frame(Area = "GOO",                   
                        Model = "Sum of indocator abundance scores",
                        Polygon = NA,
                        obs.rand = NA,
                        PC1 = seq(-0.80, 0.89, by = 0.02),  
                        PC2 = mean(dat88$PC2),
                        PC3 = mean(dat88$PC3),
                        order = "D")  
sa2.pred3 <- data.frame(Area = "OFF",                   
                        Model = "Sum of indocator abundance scores",
                        Polygon = NA,
                        obs.rand = NA,
                        PC1 = seq(-1.19, 1.05, by = 0.02),  
                        PC2 = mean(dat88$PC2),
                        PC3 = mean(dat88$PC3),
                        order = "D") 
sa2.pred <- rbind(sa2.pred1, sa2.pred2, sa2.pred3)

# Make predictions

P1 <- predict(sa.q2Bx, newdata = sa2.pred, type = "response", se.fit = TRUE, re.form = NA)  #Predictions for Area
sa2.pred$resp <- P1$fit 
sa2.pred$SE <- P1$se.fit 
sa2.pred$up <- (sa2.pred$resp+(1.96 * sa2.pred$SE)) 
sa2.pred$low <- (sa2.pred$resp-(1.96 * sa2.pred$SE))

####
# PANAL 5: Lilies
#
#Make Dataframes
ll2.pred1 <- data.frame(Area = "MULL",                  
                        Model = "Lilies",
                        Polygon = NA,
                        obs.rand = NA,
                        PC1 = seq(-1.28, 0.86, by = 0.02),  
                        PC2 = NA,
                        PC3 = NA,
                        order = "E")
ll2.pred2 <- data.frame(Area = "GOO",                   
                        Model = "Lilies",
                        Polygon = NA,
                        obs.rand = NA,
                        PC1 = seq(-0.80, 0.89, by = 0.02),  
                        PC2 = NA,
                        PC3 = NA,
                        order = "E") 
ll2.pred3 <- data.frame(Area = "OFF",                   
                        Model = "Lilies",
                        Polygon = NA,
                        obs.rand = NA,
                        PC1 = seq(-1.19, 1.05, by = 0.02),  
                        PC2 = NA,
                        PC3 = NA,
                        order = "E") 
ll2.pred <- rbind(ll2.pred1, ll2.pred2, ll2.pred3)

# Make predictions

P1 <- predict(ll.q2Bx, newdata = ll2.pred, type = "response", se.fit = TRUE, re.form = NA)  #Predictions for Area
ll2.pred$resp <- P1$fit 
ll2.pred$SE <- P1$se.fit 
ll2.pred$up <- (ll2.pred$resp+(1.96 * ll2.pred$SE)) 
ll2.pred$low <- (ll2.pred$resp-(1.96 * ll2.pred$SE))

####
# PANAL 6: Orchids
#
#Make Dataframes
or2.pred <- data.frame(Area = NA,                  
                       Model = "Orchids",
                       Polygon = NA,
                       obs.rand = NA,
                       PC1 = seq(-1.28, 1.05, by = 0.02),  
                       PC2 = mean(dat88$PC2),
                       PC3 = NA,
                       order = "F")

# Make predictions

P1 <- predict(or.q2Bx, newdata = or2.pred, type = "response", se.fit = TRUE, re.form = NA)  #Predictions for Area
or2.pred$resp <- P1$fit 
or2.pred$SE <- P1$se.fit 
or2.pred$up <- (or2.pred$resp+(1.96 * or2.pred$SE)) 
or2.pred$low <- (or2.pred$resp-(1.96 * or2.pred$SE))

####
# PANAL 7: Daisies (tuberous)
#
#Make Dataframes
dT2.pred <- data.frame(Area = NA,                  
                       Model = "Daisies T",
                       Polygon = NA,
                       obs.rand = NA,
                       PC1 = seq(-1.28, 1.05, by = 0.02),  
                       PC2 = NA,
                       PC3 = NA,
                       order = "G")

# Make predictions

P1 <- predict(dT.q2Bx, newdata = dT2.pred, type = "response", se.fit = TRUE, re.form = NA)  #Predictions for Area
dT2.pred$resp <- P1$fit 
dT2.pred$SE <- P1$se.fit 
dT2.pred$up <- (dT2.pred$resp+(1.96 * dT2.pred$SE)) 
dT2.pred$low <- (dT2.pred$resp-(1.96 * dT2.pred$SE))

####
# PANAL 8: Daisies (control)
#
#Make Dataframes
dC2.pred1 <- data.frame(Area = "MULL",                  
                        Model = "Daisies C",
                        Polygon = NA,
                        obs.rand = NA,
                        PC1 = seq(-1.28, 0.86, by = 0.02),  
                        PC2 = NA,
                        PC3 = NA,
                        order = "H")
dC2.pred2 <- data.frame(Area = "GOO",                   
                        Model = "Daisies C",
                        Polygon = NA,
                        obs.rand = NA,
                        PC1 = seq(-0.80, 0.89, by = 0.02),  
                        PC2 = NA,
                        PC3 = NA,
                        order = "H") 
dC2.pred3 <- data.frame(Area = "OFF",                   
                        Model = "Daisies C",
                        Polygon = NA,
                        obs.rand = NA,
                        PC1 = seq(-1.19, 1.05, by = 0.02),  
                        PC2 = NA,
                        PC3 = NA,
                        order = "H") 
dC2.pred <- rbind(dC2.pred1, dC2.pred2, dC2.pred3)

# Make predictions

P1 <- predict(dC.q2Bx, newdata = dC2.pred, type = "response", se.fit = TRUE, re.form = NA)  #Predictions for Area
dC2.pred$resp <- P1$fit 
dC2.pred$SE <- P1$se.fit 
dC2.pred$up <- (dC2.pred$resp+(1.96 * dC2.pred$SE)) 
dC2.pred$low <- (dC2.pred$resp-(1.96 * dC2.pred$SE))

# Make figure

#combine predicition data
Q2.pred <- rbind(ir2.pred, sr2.pred, ma2.pred, sa2.pred,
                 ll2.pred, or2.pred, dT2.pred, dC2.pred)
#new facet labels
new_labels <- c("A" = "Indicator richness",
                "B" = "Target species richness", 
                "C" = "Highest indicator abundance score",
                "D" = "Sum of indicator abundance scores",
                "E" = "Lilies indicator abundance score",
                "F" = "Orchids indicator abundance score",
                "G" = "Daisies (tuberous) indicator abundance score",
                "H" = "Daisies (control) indicator abundance score")

Q2.FIG <- Q2.pred %>% 
  mutate(Area = ifelse(Area=="OFF", "Offsets",
                        ifelse(Area=="GOO", "Goorooyarroo",
                               ifelse(Area=="MULL", "Mulligans Flat", ""))))%>%
  mutate(Area = fct_relevel(Area, "Offsets", "Goorooyarroo", "Mulligans Flat"))%>%
  #mutate(low = if_else(low < 0, 0, low),
   #      up = if_else(up >5, 5, up))%>%
  
  ggplot(aes(x=PC1))+
  geom_ribbon(aes(ymin=low, ymax=up, fill = Area), alpha = 0.6)+
  geom_line(aes(y=resp, colour = Area), size =1)+
  facet_wrap(~order, labeller = labeller(order = new_labels), scales = "free_y", ncol = 2)+
  #scale_y_discrete(labels = c( "Offsets", "Goorooyarroo", "Mulligans Flat"))+
  
  theme_bw()+xlab("PC1")+ylab("Predicted reponse")+labs(fill = "MFWS Area", colour = "MFWS Area")+
  theme(axis.text.x = element_text(colour = "black"),
        axis.text.y = element_text(colour = "black"),
        axis.title = element_text(face = "bold"),
        legend.title = element_text(face = "bold"))

Q2.FIG

ggsave("Q2-Figure.tiff",
       plot=Q2.FIG,
       width = 7.5, height = 5.0, units = "in",
       dpi = 300)  

# 5. QUESTION 3__________________________________####
#                               Within the Mulligans Flat area where Ngaluda are present, 
#                               do experimental grazing treatments influence plant indicator groups or species?

# Full GLMM models
# To address Question 3 only observations from within MULLIGANS were included (n = 46) 
# with each response fitted with a GLMM that included just the single variable
# "Experimental Treatment" (factor with 5 levels)

# Abundance scores and richness response variables were modelled with Poisson error distributions and log-link functions, 
# while occurrence variables were modelled with binomial error distributions and logit link functions. 
# The random effect of “polygon” was included in all models to account for any spatial autocorrelation associated with the paired nature of plots.

#subset the data

datMULL <- dat%>%
  filter(Area == "MULL")

# Full models for all 15 response variables
ir.q3B <- glmmTMB(INDI.RICH ~ Exp.Treat +
                    (1|Polygon), data = datMULL, family = poisson) #Indicator Richness
ma.q3B <- glmmTMB(SPP.RICH ~ Exp.Treat +
                    (1|Polygon), data = datMULL, family = poisson) #Species Richness
ma.q3B <- glmmTMB(MAX.ABUND ~ Exp.Treat +
                    (1|Polygon), data = datMULL, family = poisson) #Highest Indicator Group Abund
sa.q3B <- glmmTMB(SUM.ABUND ~ Exp.Treat +
                    (1|Polygon), data = datMULL, family = poisson) #Sum of Indicator Group Abund
ll.q3B <- glmmTMB(Lilies ~ Exp.Treat +
                    (1|Polygon), data = datMULL, family = poisson) #LILIES count
or.q3B <- glmmTMB(Orchids ~ Exp.Treat +
                    (1|Polygon), data = datMULL, family = poisson) #ORCHIDS count
dT.q3B <- glmmTMB(Daisies.T ~ Exp.Treat +
                    (1|Polygon), data = datMULL, family = poisson) #DAISIES (Tuberous) count
dC.q3B <- glmmTMB(Daisies.C ~ Exp.Treat +
                    (1|Polygon), data = datMULL, family = poisson) #Daisies (Control) count
as.q3B <- glmmTMB(Arth.spp ~ Exp.Treat +
                    (1|Polygon), data = datMULL, family = binomial) #Arthro OCC
bs.q3B <- glmmTMB(Bulb.spp ~ Exp.Treat +
                    (1|Polygon), data = datMULL, family = binomial) #Bulbine OCC
wd.q3B <- glmmTMB(Wurm.dio ~ Exp.Treat +
                    (1|Polygon), data = datMULL, family = binomial) #Wurmbia OCC
ms.q3B <- glmmTMB(Microt.spp ~ Exp.Treat +
                    (1|Polygon), data = datMULL, family = binomial) #Microtis OCC
cs.q3B <- glmmTMB(Cras.spp ~ Exp.Treat +
                    (1|Polygon), data = datMULL, family = binomial) #Craspedia OCC
ls.q3B <- glmmTMB(Lept.squ ~ Exp.Treat +
                    (1|Polygon), data = datMULL, family = binomial) #Lepto OCC
ca.q3B <- glmmTMB(Chry.api ~ Exp.Treat +
                    (1|Polygon), data = datMULL, family = binomial) #Chrys OCC


# Summary of effects

summary(ir.q3B) #lower in low grazing
summary(sr.q3B) #lower in low grazing
summary(ma.q3B) # 
summary(sa.q3B) #lower in low grazing
summary(ll.q3B) # 
summary(or.q3B) # lower in bettong fences and low grazing
summary(dT.q3B) # 
summary(dC.q3B) # 
summary(as.q3B) # 
summary(bs.q3B) #
summary(wd.q3B) # 
summary(ms.q3B) # lower in low grazing / bettong fence LG
summary(cs.q3B) # 
summary(ls.q3B) # 
summary(ca.q3B) # 

# R2 of best models

r2_nakagawa(ir.q3B) # 59
r2_nakagawa(sr.q3B) # 63
r2_nakagawa(ma.q3B) # 62
r2_nakagawa(sa.q3B) # 75
r2_nakagawa(ll.q3B) # NA, 98
r2_nakagawa(or.q3B) # 98, 97
r2_nakagawa(dT.q3B) # 96, 95
r2_nakagawa(dC.q3B) # 62, 32
r2_nakagawa(as.q3B) # 94, 93
r2_nakagawa(bs.q3B) # 86, 81
r2_nakagawa(wd.q3B) # 70, 43
r2_nakagawa(ms.q3B) # 99, 45
r2_nakagawa(cs.q3B) # 99, 53
r2_nakagawa(ls.q3B) # nnnnnn
r2_nakagawa(ca.q3B) # 99, 05


# check and correct models as needed

check_model(ir.q3B)
simulationOutput <- simulateResiduals(fittedModel = ir.q3B, n = 250)
plotSimulatedResiduals(simulationOutput = simulationOutput) 

check_model(sr.q3B) 
simulationOutput <- simulateResiduals(fittedModel = sr.q3B, n = 250)
plotSimulatedResiduals(simulationOutput = simulationOutput) 

check_model(ma.q3B)
simulationOutput <- simulateResiduals(fittedModel = ma.q3B, n = 250)
plotSimulatedResiduals(simulationOutput = simulationOutput) 


check_model(sa.q3B)
simulationOutput <- simulateResiduals(fittedModel = sa.q3B, n = 250)
plotSimulatedResiduals(simulationOutput = simulationOutput) 

check_model(ll.q3B)
simulationOutput <- simulateResiduals(fittedModel = ll.q3B, n = 250)
plotSimulatedResiduals(simulationOutput = simulationOutput) # NA, 98

check_model(or.q3B)
simulationOutput <- simulateResiduals(fittedModel = or.q3B, n = 250)
plotSimulatedResiduals(simulationOutput = simulationOutput) # NA, 98

check_model(dT.q3B)
simulationOutput <- simulateResiduals(fittedModel = dT.q3B, n = 250)
plotSimulatedResiduals(simulationOutput = simulationOutput) # 

check_model(dC.q3B)
simulationOutput <- simulateResiduals(fittedModel = dC.q3B, n = 250)
plotSimulatedResiduals(simulationOutput = simulationOutput)

check_model(as.q3B)
simulationOutput <- simulateResiduals(fittedModel = as.q3B, n = 250)
plotSimulatedResiduals(simulationOutput = simulationOutput)

check_model(bs.q3B)
simulationOutput <- simulateResiduals(fittedModel = bs.q3B, n = 250)
plotSimulatedResiduals(simulationOutput = simulationOutput)

check_model(wd.q3B)
simulationOutput <- simulateResiduals(fittedModel = wd.q3B, n = 250)
plotSimulatedResiduals(simulationOutput = simulationOutput)

check_model(ms.q3B) 
simulationOutput <- simulateResiduals(fittedModel = ms.q3B, n = 250)
plotSimulatedResiduals(simulationOutput = simulationOutput)

check_model(cs.q3B)
simulationOutput <- simulateResiduals(fittedModel = cs.q3B, n = 250)
plotSimulatedResiduals(simulationOutput = simulationOutput)

check_model(ls.q3B)
simulationOutput <- simulateResiduals(fittedModel = ls.q3B, n = 250)
plotSimulatedResiduals(simulationOutput = simulationOutput)

check_model(ca.q3B)
simulationOutput <- simulateResiduals(fittedModel = ca.q3B, n = 250)
plotSimulatedResiduals(simulationOutput = simulationOutput)

#   5.1 PREDICTIONS AND FIGURE (Q3)______________####

####
# PANAL 1: INDICATOR RICHNESS
#
#Make Dataframes

ir3.pred <- data.frame(Exp.Treat = unique(datMULL$Exp.Treat),                  
                       Model = "Indi Rich",
                       Polygon = NA,
                       obs.rand = NA)

# Make predictions

P1 <- predict(ir.q3B, newdata = ir3.pred, type = "response", se.fit = TRUE, re.form = NA)  #Predictions for Area
ir3.pred$resp <- P1$fit 
ir3.pred$SE <- P1$se.fit 
ir3.pred$up <- (ir3.pred$resp+(1.96 * ir3.pred$SE)) 
ir3.pred$low <- (ir3.pred$resp-(1.96 * ir3.pred$SE))

# Make figure


ir3.fig <- ir3.pred %>% 
  mutate(Exp.Treat = fct_relevel(Exp.Treat, "Low grazing", "Bettong fence LG", "Bettong fence HG", "Bettong fence FT", "High grazing"))%>%
  mutate(low = if_else(low < 0, 0, low),
         up = if_else(up >5, 5, up))%>%
  
  ggplot(aes(y=Exp.Treat))+
  geom_vline(xintercept = 2, linetype = "dashed", size = 1, color = "black")+
  scale_y_discrete(labels = c( "High Grazing", "Bettong Fence \n(Floppy Top)", "Bettong Fence \n(High Grazing)", "Bettong Fence \n(Low Grazing)", "Low Grazing"))+
  geom_jitter(data = datMULL, aes(x=INDI.RICH), colour = "#619CFF", size = 2, height = 0.1, width = 0.1)+
  geom_pointrange(aes(x=resp, xmin=low, xmax=up), size = 1) +
  
  
  theme_bw()+xlab("Indicator group richness")+

  theme(axis.text.x = element_text(colour = "black"),
        axis.text.y = element_text(colour = "black"),
        legend.position = "none",
        axis.title.y = element_blank(),
        legend.title = element_blank())

ir3.fig

####
# PANAL 2: SPECIES RICHNESS
#
sr3.pred <- data.frame(Exp.Treat = unique(datMULL$Exp.Treat),                  
                       Model = "Species Rich",
                       Polygon = NA,
                       obs.rand = NA)

# Make predictions

P1 <- predict(sr.q3B, newdata = sr3.pred, type = "response", se.fit = TRUE, re.form = NA)  #Predictions for Area
sr3.pred$resp <- P1$fit 
sr3.pred$SE <- P1$se.fit 
sr3.pred$up <- (sr3.pred$resp+(1.96 * sr3.pred$SE)) 
sr3.pred$low <- (sr3.pred$resp-(1.96 * sr3.pred$SE))

# Make figure


sr3.fig <- sr3.pred %>% 
  mutate(Exp.Treat = fct_relevel(Exp.Treat, "Low grazing", "Bettong fence LG", "Bettong fence HG", "Bettong fence FT", "High grazing"))%>%
  mutate(low = if_else(low < 0, 0, low),
         up = if_else(up >13, 13, up))%>%
  
  ggplot(aes(y=Exp.Treat))+
  geom_vline(xintercept = 2.416667, linetype = "dashed", size = 1, color = "black")+
  scale_y_discrete(labels = c( "High Grazing", "Bettong Fence \n(Floppy Top)", "Bettong Fence \n(High Grazing)", "Bettong Fence \n(Low Grazing)", "Low Grazing"))+
  geom_jitter(data = datMULL, aes(x=SPP.RICH), colour = "#619CFF", size = 2, height = 0.1, width = 0.1)+
  geom_pointrange(aes(x=resp, xmin=low, xmax=up), size = 1) +
  
  
  theme_bw()+xlab("Target species richness")+
  
  theme(axis.text.x = element_text(colour = "black"),
        axis.text.y = element_text(colour = "black"),
        legend.position = "none",
        axis.title.y = element_blank(),
        legend.title = element_blank())

sr3.fig

####
# PANAL 3: MAX ABUND
#
ma3.pred <- data.frame(Exp.Treat = unique(datMULL$Exp.Treat),                  
                       Model = "max abund",
                       Polygon = NA,
                       obs.rand = NA)

# Make predictions

P1 <- predict(ma.q3B, newdata = ma3.pred, type = "response", se.fit = TRUE, re.form = NA)  #Predictions for Area
ma3.pred$resp <- P1$fit 
ma3.pred$SE <- P1$se.fit 
ma3.pred$up <- (ma3.pred$resp+(1.96 * ma3.pred$SE)) 
ma3.pred$low <- (ma3.pred$resp-(1.96 * ma3.pred$SE))

# Make figure


ma3.fig <- ma3.pred %>% 
  mutate(Exp.Treat = fct_relevel(Exp.Treat, "Low grazing", "Bettong fence LG", "Bettong fence HG", "Bettong fence FT", "High grazing"))%>%
  mutate(low = if_else(low < 0, 0, low),
         up = if_else(up >4, 4, up))%>%
  
  ggplot(aes(y=Exp.Treat))+
  geom_vline(xintercept = 2, linetype = "dashed", size = 1, color = "black")+
  scale_y_discrete(labels = c( "High Grazing", "Bettong Fence \n(Floppy Top)", "Bettong Fence \n(High Grazing)", "Bettong Fence \n(Low Grazing)", "Low Grazing"))+
  geom_jitter(data = datMULL, aes(x=MAX.ABUND), colour = "#619CFF", size = 2, height = 0.1, width = 0.1)+
  geom_pointrange(aes(x=resp, xmin=low, xmax=up), size = 1) +
  scale_x_continuous(limits = c(0,NA),
                     labels = c("0", "1 \n1-9", "2 \n10-99", "3 \n100-999", "4 \n1000+"))+
  
  theme_bw()+xlab("Highest indicator abundance score")+
  
  theme(axis.text.x = element_text(colour = "black"),
        axis.text.y = element_text(colour = "black"),
        legend.position = "none",
        axis.title.y = element_blank(),
        legend.title = element_blank())

ma3.fig

####
# PANAL 4: ORCHIDS
#
#Make Dataframes

or3.pred <- data.frame(Exp.Treat = unique(datMULL$Exp.Treat),                  
                       Model = "orchid",
                       Polygon = NA,
                       obs.rand = NA)

# Make predictions

P1 <- predict(or.q3B, newdata = or3.pred, type = "response", se.fit = TRUE, re.form = NA)  #Predictions for Area
or3.pred$resp <- P1$fit 
or3.pred$SE <- P1$se.fit 
or3.pred$up <- (or3.pred$resp+(1.96 * or3.pred$SE)) 
or3.pred$low <- (or3.pred$resp-(1.96 * or3.pred$SE))

# Make figure


or3.fig <- or3.pred %>% 
  mutate(Exp.Treat = fct_relevel(Exp.Treat, "Low grazing", "Bettong fence LG", "Bettong fence HG", "Bettong fence FT", "High grazing"))%>%
  mutate(low = if_else(low < 0, 0, low),
         up = if_else(up >4, 4, up))%>%
  
  ggplot(aes(y=Exp.Treat))+
  geom_vline(xintercept =  1.7499999, linetype = "dashed", size = 1, color = "black")+
  scale_y_discrete(labels = c( "High Grazing", "Bettong Fence \n(Floppy Top)", "Bettong Fence \n(High Grazing)", "Bettong Fence \n(Low Grazing)", "Low Grazing"))+
  geom_jitter(data = datMULL, aes(x=Orchids), colour = "#619CFF", size = 2, height = 0.1, width = 0.1)+
  geom_pointrange(aes(x=resp, xmin=low, xmax=up), size = 1) +
  scale_x_continuous(limits = c(0,NA),
                     labels = c("0", "1 \n1-9", "2 \n10-99", "3 \n100-999", "4 \n1000+"))+
  
  theme_bw()+xlab("Orchids indicator abundance score")+
  
  theme(axis.text.x = element_text(colour = "black"),
        axis.text.y = element_text(colour = "black"),
        legend.position = "none",
        axis.title.y = element_blank(),
        legend.title = element_blank())

or3.fig


####
# JOIN PANNELS INTO 1 FIGURE
####

Q3.FIG <- ggarrange(ir3.fig, sr3.fig, ma3.fig, or3.fig,
                    ncol = 2, nrow = 2,
                    labels = "AUTO")
Q3.FIG

ggsave("Q3-Figure.tiff",
       plot=Q3.FIG,
       width = 7.3, height = 5.1, units = "in",
       dpi = 300)  

