#### BEE FUNCTIONAL TRAIT ANALYSIS #####
#### Monika Egerer and Hamutahl
#### October 2020

### game plan: 
### 1) RLQ analysis
### 2) Fourth Corner Analysis

setwd("~/Documents/ResearchProjects/BeeTraitProject/DataForRevisions")

# Getting biodiversity data into R
library(picante)
library(ade4)
library(vegan)


# load data needed:
beediv <- read.csv("BeeDiversity2015.csv", header = TRUE)
beetraitsNoRare <- read.csv("BeeSpeciesTraitsFoundNoRareSpp.csv", header=TRUE)
beecomm <- read.csv("BeeCommunity2015.csv", header = TRUE)
metadata <- read.csv("Garden2015.csv", header = TRUE) # metadata (garden site data)
metadataT <- metadata[,1:10]
metadataNT <- metadata[,11:19]


#############################
## RLQ Analysis for Bees
#############################


#### READ IN THE DATA ####
# Need species table (rows=sites; columns=species), environmental table (rows=sites; columns=env variables[continous and/or categorical])
#and trait table  (rows=species; columns=traits [can be continous and/or categorical])
traits <- read.table('BeeSpeciesTraitsFoundNoRareSpp.csv',
                     header=TRUE,  sep=",", na.strings="NA", dec=".", strip.white=TRUE)
env <- read.table('BeeGardenVariables2015.csv',
                  header=TRUE,  sep=",", na.strings="NA", dec=".", strip.white=TRUE)
species <- read.table('BeeCommunity2015.csv',
                      header=TRUE,  sep=",", na.strings="NA", dec=".", strip.white=TRUE)
envf<- env[,2:8]
# pick: LNGardenSize+Bare1m+LNFlowers+HerbPlantSpp+Urban2kmSQRT
library(Hmisc)
rcorr(as.matrix(envf), type="pearson")
## Must pick either bare or mulch, we picked bare soil
envNT <- envf[,c(1:3,5:7)]


### DOUBLE CHECKED AND CORRECTED ON 07 OCT 2020
# Need to have matching data sets for traits and species.
# the community composition data has more species because it includes "genus spp." which we excluded from the trait dataset
row.names(species)<- species$Site
species$Site <-NULL
speciesNoRare <- species # remove singletons and species we don't have data for (13? species)
speciesNoRare$Andrenasp <-NULL
speciesNoRare$Anthophorasp <-NULL
speciesNoRare$Colletesp <-NULL
speciesNoRare$Eucerasp <-NULL
speciesNoRare$Halictusrubicundus <-NULL
speciesNoRare$Halictussp <-NULL
speciesNoRare$HylaeusspB <-NULL
speciesNoRare$HylaeusspA <-NULL
speciesNoRare$Hylaeussp <-NULL
speciesNoRare$Holcopasitessp <- NULL
speciesNoRare$Lasioglossumsp <- NULL
speciesNoRare$LasioglossumDialictusincompletum <-NULL
speciesNoRare$Megachileapicalis <-NULL
speciesNoRare$Megachilelatimanus <-NULL
speciesNoRare$Megachilerelativa <-NULL
speciesNoRare$MegachilespA <-NULL
speciesNoRare$Megachilesp <- NULL
speciesNoRare$MelissodesspA <-NULL
speciesNoRare$unknownunknown <-NULL
speciesNoRare$Nomadasp <- NULL

### DOUBLE CHECKED AND CORRECTED ON 07 OCT 2020


### FOURTH CORNER ANALYSES: 
##final traits are interteg (mm), sociality, nest location, nest behavior, pollenstructure, lecty


# remove species name from trait data
row.names(traits) <- traits$beespecies
traits$beespecies <- NULL

#i didnt remove parasitsim, but if i need to, use this code
#traits$parasitism <- NULL

# check if we need to transform
spec.dca <- decorana(speciesNoRare)
summary(spec.dca, display='none')
#we dont need to transform; axes are <3


# environmental
envNT$LNGardenSize <- with(envNT, (LNGardenSize+3)) ## put on positive scale
env.dca <- decorana(envNT)
summary(env.dca, display='none') 
#we dont need to transform; axes are <3

#### Fourth corner analysis:

four.comb.beefdfbody <- fourthcorner(envNT, speciesNoRare,
                                 traits, modeltype = 6, p.adjust.method.G = "fdr",
                                 p.adjust.method.D = "fdr", nrepet = 999) ## with adjustment

four.comb.beeTbody <- fourthcorner(envNT, speciesNoRare,
                               traits, modeltype = 6, p.adjust.method.G = "none",
                               p.adjust.method.D = "none", nrepet = 999) ## without adjustment

#stat="D2": the association is measured between the quantitative variable and each category sepa-
#rately. A correlation coefcient is used to indicate the strength of the association between the given
#category and the small or large values of the quantitative variable.
#model 6 with FDR corrects for Type I error
summary(four.comb.beefdfbody)
summary(four.comb.beeTbody)
png('4corner-bee-G-BodySize.png', width=5, height=5, units="in", res=300)
plot(four.comb.beeTbody, alpha = 0.05, stat = "G")
dev.off()
png('4corner-bee-D2-BodySize.png', width=5, height=10, units="in", res=300)
plot(four.comb.beeTbody, alpha = 0.05, stat = "D2")
dev.off()


L.bee <-dudi.coa(speciesNoRare, scannf = FALSE) # correspondence analysis for species
R.bee <-dudi.pca(envNT, row.w = L.bee$lw,
                 scannf = FALSE) # pca for environmental variables (they are all quantitative)
Q.bee <-dudi.hillsmith(traits, row.w = L.bee$cw,
                       scannf = FALSE)  # Hills smith for traits since there are quantitative and qualitative
rlq.bee <- rlq(R.bee, L.bee, Q.bee,
               scannf = FALSE)  # rlq analysis

png('RLQ-bee-Full.png', width=8, height=8, units="in", res=300)
plot(rlq.bee)
dev.off()

#extract scores
rlq.bee$c1 
rlq.bee$lQ
rlq.bee$l1
rlq.bee$lR

summary(Q.bee)
summary(R.bee)
summary(L.bee)
summary(rlq.bee)

## Percentage of co-Inertia for each axis
100*rlq.bee$eig/sum(rlq.bee$eig)

#The diferent figures can be obtained separately by plotting the different elements contained in the
#rlq.bee object:

png('RLQ-bee-nmds.png', width=9, height=6, units="in", res=300)
par(mfrow = c(1, 3))
s.arrow(rlq.bee$l1)
s.arrow(rlq.bee$c1)
s.label(rlq.bee$lQ, boxes = FALSE)
dev.off()

#A biplot representing traits and environmental variables 

par(mfrow = c(1, 1))
s.arrow(rlq.bee$c1, xlim=c(-2,2), boxes = FALSE, grid=FALSE, clab=0.8)
s.label(rlq.bee$li, add.plot=T, clab=1)

#Species scores on the frst two axes of RLQ analysis:
s.label(rlq.bee$lQ, clabel = 0)
par(mar = c(0.1, 0.1, 0.1, 0.1))
pointLabel(rlq.bee$lQ,row.names(rlq.bee$lQ), cex=0.7)

#combining both approaches

#First, a multivariate test can be applied to evaluate the global significance of the traits-environment
#relationships:
testrlq.bee <- randtest(rlq.bee, modeltype = 6, nrepet = 999)
testrlq.bee

#The total inertia of RLQ analysis is equal to the SRLQ multivariate statistic dened in Dray and
#Legendre (2008). This statistic is returned by the fourthcorner2 function:
nrepet=999
Srlq <- fourthcorner2(envNT, speciesNoRare, traits,
                      modeltype = 6, p.adjust.method.G = "none", nrepet = nrepet)
Srlq$trRLQ

summary(Srlq)

#Both approaches can be combined if RLQ scores are used to represent traits and environmental variables
#on a biplot. Then, signicant associations revealed by the fourthcorner approach can be represented using
#segments (blue lines for negative associations, red lines for positive associations, see the argument col).
#Only traits and environmental variables that have at least one significant association are represented. Here,
#we apply this method using adjusted pvalues for multiple comparisons and a signicant level  = 0:05.

#Another approach is provided by the fourthcorner.rlq function and consists in testing directly the
#links between RLQ axes and traits (typetest="Q.axes") or environmental variables (typetest="R.axes").

testQaxes.comb.bee <- fourthcorner.rlq(rlq.bee, modeltype = 6,
                                       typetest = "Q.axes", nrepet = 999, p.adjust.method.G = "none",
                                       p.adjust.method.D = "none")
testRaxes.comb.bee <- fourthcorner.rlq(rlq.bee, modeltype = 6,
                                       typetest = "R.axes", nrepet = 999, p.adjust.method.G = "none",
                                       p.adjust.method.D = "none")
print(testQaxes.comb.bee, stat = "D")
print(testRaxes.comb.bee, stat = "D")

#Results can be represented using a table with colors indicating signicance :

png('RLQ-biplot-bee-bodysize1.png', width=5, height=7, units="in", res=300)
par(mfrow = c(1, 1))
plot(testQaxes.comb.bee, alpha = 0.05, type = "table",
     stat = "D2")
dev.off()
png('RLQ-biplot-bee-bodysize2.png', width=5, height=3, units="in", res=300)
par(mfrow = c(1, 1))
plot(testRaxes.comb.bee, alpha = 0.05, type = "table",
     stat = "D2")
dev.off()

# with ITD
png('RLQ-biplot-bee-ITD1.png', width=5, height=7, units="in", res=300)
par(mfrow = c(1, 1))
plot(testQaxes.comb.bee, alpha = 0.05, type = "table",
     stat = "D2")
dev.off()
png('RLQ-biplot-bee-ITD2.png', width=5, height=3, units="in", res=300)
par(mfrow = c(1, 1))
plot(testRaxes.comb.bee, alpha = 0.05, type = "table",
     stat = "D2")
dev.off()


#Significance with axes can also be reported on the factorial map of RLQ analysis. Here, significant
#associations with the frst axis are represented in blue, with the second axis in orange, with both axes in
#green (variables with no signicant association are in black):

png('RLQ-factorialmap-bee.png', width=7, height=4, units="in", res=300)
par(mfrow = c(1, 2))
plot(testQaxes.comb.bee, alpha = 0.05, type = "biplot",
     stat = "D2", col = c("black", "blue", "orange", "green"))
plot(testRaxes.comb.bee, alpha = 0.05, type = "biplot",
     stat = "D2", col = c("black", "blue", "orange", "green"))
dev.off()

########################################################################

# with ITD

four.comb.beefdfbodyITD <- fourthcorner(envNT, speciesNoRareITD,
                                     traitsITD, modeltype = 6, p.adjust.method.G = "fdr",
                                     p.adjust.method.D = "fdr", nrepet = 999) ## with adjustment

four.comb.beeTbodyITD <- fourthcorner(envNT, speciesNoRareITD,
                                      traitsITD, modeltype = 6, p.adjust.method.G = "none",
                                   p.adjust.method.D = "none", nrepet = 999) ## without adjustment

#stat="D2": the association is measured between the quantitative variable and each category sepa-
#rately. A correlation coefcient is used to indicate the strength of the association between the given
#category and the small or large values of the quantitative variable.
#model 6 with FDR corrects for Type I error
summary(four.comb.beefdfbodyITD)
summary(four.comb.beeTbodyITD)
png('4corner-bee-G-ITD.png', width=5, height=5, units="in", res=300)
plot(four.comb.beeT, alpha = 0.05, stat = "G")
dev.off()
png('4corner-bee-D2-ITD.png', width=5, height=10, units="in", res=300)
plot(four.comb.beeT, alpha = 0.05, stat = "D2")
dev.off()


L.bee <-dudi.coa(speciesNoRareITD, scannf = FALSE) # correspondence analysis for species
R.bee <-dudi.pca(envNT, row.w = L.bee$lw,
                 scannf = FALSE) # pca for environmental variables (they are all quantitative)
Q.bee <-dudi.hillsmith(traitsITD, row.w = L.bee$cw,
                       scannf = FALSE)  # Hills smith for traits since there are quantitative and qualitative
rlq.bee <- rlq(R.bee, L.bee, Q.bee,
               scannf = FALSE)  # rlq analysis

png('RLQ-bee-Full.png', width=8, height=8, units="in", res=300)
plot(rlq.bee)
dev.off()


################################################v########################
########################################################################

# BEE RICHNESS AND ABUNDANCE MODELS

# GLMulti for each of the biodiv. measures
library(glmulti)
library(visreg)
library(RVAideMemoire)
library(MuMIn)
library(mvabund)

#check data for overdisperson if using poission
global.model<-glm(totalabundallperiod ~ GardenSizeLN+TreesShrubs+Bare1m+FlowersLN+HerbPlantSpp+Urban2kmSQRT, family=quasipoisson, data=beediv)
summary(global.model)
#is the residual deviance is greater than the degrees of freedom? if yes, there is overdispersion, 
#with some variance not accounted for by the model or error structure

#global model selection
global.model<-glm(avgabundper ~ GardenSizeLN+TreesShrubs+Bare1m+FlowersLN+HerbPlantSpp+Urban2kmSQRT, family=gaussian, data=beediv)
TotBeeAbun<- glmulti(global.model, level = 1, crit = aicc)
summary(TotBeeAbun)
weightable(TotBeeAbun)

global.model<-glm(totalabundallperiod ~ GardenSizeLN+TreesShrubs+Bare1m+FlowersLN+HerbPlantSpp+Urban2kmSQRT, family=poisson, data=beediv)
TotBeeAbun<- glmulti(global.model, level = 1, crit = aicc)
summary(TotBeeAbun)
weightable(TotBeeAbun)



###################################################
####### Figures #########################
###################################################

## plots for abundance
ba1<-ggplot(data=dfglm, aes(x=LNGardenSize, y=totalabundallperiod)) +
  geom_point(size=2) +
  geom_smooth(method ="glm", se=TRUE, fullrange=TRUE,
              method.args = list(family = "poisson"), col="gray")+
  geom_point() +
  theme_classic() +
  xlab("Garden size (LN)") +
  ylab("Bee abundance")
ba2<-ggplot(data=dfglm, aes(x=Urban2kmSQRT, y=totalabundallperiod)) +
  geom_point(size=2) +
  geom_smooth(method ="glm", se=TRUE, fullrange=TRUE,
              method.args = list(family = "poisson"), col="gray")+
  geom_point() +
  theme_classic() + 
  xlab("% Urban (2 km)") +
  ylab("")
ba3<-ggplot(data=dfglm, aes(x=LNFlowers, y=totalabundallperiod)) +
  geom_point(size=2) +
  geom_smooth(method ="glm", se=TRUE, fullrange=TRUE,
              method.args = list(family = "poisson"), col="gray")+
  geom_point() +
  theme_classic()+ 
  xlab("No. flowers (LN)") +
  ylab("")

library(gridExtra)
panel.act<-grid.arrange(ba1,ba2,ba3,ncol=3,nrow=1)
ggsave("BeeAbundance2015-untransformed.png", panel.act, dpi=300, height=2.5, width=6)
ggsave("BeeAbundance2015-transformed.png", panel.act, dpi=300, height=2.5, width=6)

## plots for species richness
br1<-ggplot(data=dfglm, aes(x=GardenSize, y=totalrichnessallperiod)) +
  geom_point(size=2) +
  geom_smooth(method ="glm", se=TRUE, fullrange=TRUE,
              method.args = list(family = "poisson"), col="gray")+
  geom_point() +
  theme_classic() +
  ylim(5,30) +
  xlab("Garden size") +
  ylab("Bee spp. richness")
br2<-ggplot(data=dfglm, aes(x=Urban2km, y=totalrichnessallperiod)) +
  geom_point(size=2) +
  geom_smooth(method ="glm", se=TRUE, fullrange=TRUE,
              method.args = list(family = "poisson"), col="gray")+
  geom_point() +
  ylim(5,30) +
  theme_classic() + 
  xlab("% Urban (2 km)") +
  ylab("")

library(gridExtra)
panel.act<-grid.arrange(br1,br2,ncol=2,nrow=1)
ggsave("BeeSRichness2015-untransformed.png", panel.act, dpi=300, height=2.5, width=6)
ggsave("BeeSRichness2015-transformed.png", panel.act, dpi=300, height=2.5, width=6)



### notes






