#### BEE Models  #####
### general diversity analysis with local and landscape factors
s

### 1) general diversity analysis with local and landscape factors

# Getting biodiversity data into R
rm(list=ls())
library(picante)
library(ade4)
library(vegan)
library(visreg)
library(RVAideMemoire)
library(MuMIn)
library(mvabund)
library(lme4)
library(car)


setwd("~/Dropbox/traits")

by.site <- read.csv("bySite.csv")

#Average Bee Richness, Average Bee Abundance

###-TEST###
ys <- c("Richness",
        "TotalAbundance")

xvars <-      c("GardenSizeLN",
                "scale(TreesShrubs)",
                "scale(Bare1m)",
                "FlowersLN",
                "scale(HerbPlantSpp)",
                "Urban2kmSQRT",
                "(1|Site)")


formulas <-lapply(ys, function(y) {
  as.formula(paste(y, "~",
                   paste(paste(xvars,
                               collapse="+"))))
})

names(formulas) <- ys

#bee abundance model
#use a negative binomial error model to account for overdispersion. (a quasipoisson will achieve the same result)
bee.abund.mod <- glmer.nb(formulas[["TotalAbundance"]],
                          na.action = "na.fail",
                          glmerControl(optimizer="bobyqa"),
                          data=by.site)

vif(bee.abund.mod)

#exclude mulch and bare soil from same model
#ms.bee.abund <- dredge(bee.abund.mod,
 #                      subset =
   #                      !("scale(Bare1m)" && "scale(Mulch1m)")) 
                      

##OR TRY THIS, without exluding
ms.bee.abund <- dredge(bee.abund.mod) 

## model average within 2 AICc of the min
ma.bee.abund <- model.avg(ms.bee.abund, subset= delta < 2,
                          revised.var = TRUE)



## full model of wild bee richness
#we fit a gaussian error model 
bee.rich.mod <- lmer(formulas[["Richness"]],
                     na.action = "na.fail",
                     data=by.site)

#exclude mulch and bare soil from same model
#ms.bee.rich <- dredge(bee.rich.mod,
#                       subset =
 #                        !("scale(Bare1m)" && "scale(Mulch1m)")) 

#try this
ms.bee.rich <- dredge(bee.rich.mod) 


ma.bee.rich <- model.avg(ms.bee.rich, subset= delta < 2,
                         revised.var = TRUE)

mods <- list(ma.bee.abund,
             ma.bee.rich)

inv.logit <- function(x){
  exp(x)/(1+exp(x))
}

sumMSdredge <- function(res){
  res <- summary(res)
  mmi <- as.data.frame(res$coefmat.subset)
  
  intercept <- mmi[ "(Intercept)", "Estimate"]
  next.levs <- c(0, mmi[, "Estimate"][-1])
  SEs <- mmi[, "Std. Error"]
  
  mmi$P1 <- inv.logit(intercept + next.levs)
  mmi$P1.ci.ub <-     mmi$P1 +  qnorm(0.975)*SEs
  mmi$P1.ci.lb <-     mmi$P1 -  qnorm(0.975)*SEs
  
  mmi$OR  <- exp(mmi[, "Estimate"])
  mmi$OR.ci.ub <-   exp(mmi[, "Estimate"] +  qnorm(0.975)*SEs)
  mmi$OR.ci.lb <-   exp(mmi[, "Estimate"] -  qnorm(0.975)*SEs)
  
  mmi$ORdelta  <- (mmi$OR-1)*100
  mmi$ORdelta.ci.ub <-    (mmi$OR.ci.ub -1)*100
  mmi$ORdelta.ci.lb <-   (mmi$OR.ci.lb -1)*100
  mmi <- round(mmi, 3)
  return(mmi)
}

coeffs <- lapply(mods, sumMSdredge)

args <- commandArgs(trailingOnly=TRUE)

if(length(args) != 0){
  focal.bee <- args[1]
} else{
  focal.bee <- "all"
}

save(ma.bee.abund,
     ms.bee.abund,
     ma.bee.rich,
     ms.bee.rich,
     file=sprintf("%s_beeMods.RData",
                  gsub(" ", "", focal.bee)))


mapply(function(x, y){
  write.csv(x,
            file=sprintf("beeMods_%s.csv",
                         y))
  write.table(x,
              file=sprintf("beeMods_%s.txt",
                           y), sep="&")
},
x=coeffs,
y=ys[1:2]
)


###PLOTTING RESULTS
rm(list=ls())
source("src/initialize.R")
source("src/misc.R")
source("src/predictIntervals.R")
source("src/plotPanels.R")
source("src/CIplotting.R")
source("src/diagnostics.R")
library(viridis)
library(boot)
library(car)
library(ggplot2)
library(gridExtra)

load(sprintf('all_BeeMods.RData'))

#ABUNDANCE RESULTS AND PLOTS
summary(ma.bee.abund)

top.mod.abund <- get.models(ms.bee.abund, 1)[[1]]
summary(top.mod.abund)
r.squaredGLMM(top.mod.abund)
vif(top.mod.abund)


plotDiagAbund <- function(){
  plotDiagnostics(top.mod.abund, by.site)
}

pdf.f(plotDiagAbund,
      file=file.path('figures/beeAbund.pdf'),
      height=7, width=3)


cols.var <- add.alpha(viridis(3), 0.5)[c(1,3,2)]
names(cols.var) <- levels(by.site$sampleround)

#plot abund by survey round by site type
by.site <- read.csv("bySite.csv",header=TRUE)
ba1<-ggplot(data=by.site, aes(x=GardenSizeLN, y=TotalAbundance)) +
  geom_point(size=2) +
  geom_smooth(method ="glm", se=TRUE, fullrange=TRUE,
              method.args = list(family = "poisson"), col="gray")+
  geom_point() +
  theme_classic() +
  xlab("Garden Size") +
  ylab("Bee Abundance")

panel.act<-grid.arrange(ba1,ncol=1,nrow=1)
ggsave("BeeAbundance2015Perid.png", panel.act, dpi=300, height=2, width=5)

#plot abund, collapsing total abundance from each sampling round for each site (adding them together)
site.avg <- read.csv("BeeDiversity2015.csv")
ba1<-ggplot(data=site.avg, aes(x=GardenSizeLN, y=totalabundallperiod)) +
  geom_point(size=2) +
  geom_smooth(method ="glm", se=TRUE, fullrange=TRUE,
              method.args = list(family = "poisson"), col="gray")+
  geom_point() +
  theme_classic() +
  xlab("Garden Size") +
  ylab("Bee Abundance")


panel.act<-grid.arrange(ba1,ncol=1,nrow=1)
ggsave("BeeAbundance2015Total.png", panel.act, dpi=300, height=2, width=5)

#plot abund, averaging total abundance from each sampling round for a site
ba1<-ggplot(data=site.avg, aes(x=GardenSizeLN, y=avgabundper)) +
  geom_point(size=2) +
  geom_smooth(method ="glm", se=TRUE, fullrange=TRUE,
              method.args = list(family = "gaussian"), col="gray")+
  geom_point() +
  theme_classic() +
  xlab("Garden Size") +
  ylab("Bee Abundance")

panel.act<-grid.arrange(ba1,ncol=1,nrow=1)
ggsave("BeeAbundance2015Avg.png", panel.act, dpi=300, height=2, width=5)


#RICHNESS RESULTS AND PLOT
summary(ma.bee.rich)

top.mod.abund <- get.models(ms.bee.rich, 1)[[1]]
summary(top.mod.rich)
r.squaredGLMM(top.mod.rich)
vif(top.mod.rich)

plotDiagRich <- function(){
  plotDiagnostics(top.mod.rich, by.site)
}

pdf.f(plotDiagRich,
      file=file.path('figures/beeRich.pdf'),
      height=7, width=3)

#avg bee richness at a site
ba1<-ggplot(data=site.avg, aes(x=GardenSizeLN, y=avgrichper)) +
  geom_point(size=2) +
  geom_smooth(method ="glm", se=TRUE, fullrange=TRUE,
              method.args = list(family = "gaussian"), col="gray")+
  geom_point() +
  theme_classic() +
  xlab("Garden Size") +
  ylab("Bee Richness")
ba2<-ggplot(data=site.avg, aes(x=FlowersLN, y=avgrichper)) +
  geom_point(size=2) +
  geom_smooth(method ="glm", se=TRUE, fullrange=TRUE,
              method.args = list(family = "gaussian"), col="gray")+
  geom_point() +
  theme_classic() + 
  xlab("No. Flowers") +
  ylab("")
ba3<-ggplot(data=site.avg, aes(x=Bare1m, y=avgrichper)) +
  geom_point(size=2) +
  geom_smooth(method ="glm", se=TRUE, fullrange=TRUE,
              method.args = list(family = "gaussian"), col="gray")+
  geom_point() +
  theme_classic()+ 
  xlab("Bare Soil") +
  ylab("")
ba4<-ggplot(data=site.avg, aes(x=Urban2kmSQRT, y=avgrichper)) +
  geom_point(size=2) +
  geom_smooth(method ="glm", se=TRUE, fullrange=TRUE,
              method.args = list(family = "gaussian"), col="gray")+
  geom_point() +
  theme_classic()+ 
  xlab("Urban (2km)") +
  ylab("")

panel.act<-grid.arrange(ba1,ba2,ba3,ba4,ncol=4,nrow=1)
ggsave("BeeRichnessAvg.png", panel.act, dpi=300, height=2, width=10)




#replot, with total richness per site

ba1<-ggplot(data=site.avg, aes(x=GardenSizeLN, y=totalrichnessallperiod)) +
  geom_point(size=2) +
  geom_smooth(method ="glm", se=TRUE, fullrange=TRUE,
              method.args = list(family = "poisson"), col="gray")+
  geom_point() +
  theme_classic() +
  xlab("Garden Size") +
  ylab("Bee Richness")
ba2<-ggplot(data=site.avg, aes(x=FlowersLN, y=totalrichnessallperiod)) +
  geom_point(size=2) +
  geom_smooth(method ="glm", se=TRUE, fullrange=TRUE,
              method.args = list(family = "poisson"), col="gray")+
  geom_point() +
  theme_classic() + 
  xlab("No. Flowers") +
  ylab("")
ba3<-ggplot(data=site.avg, aes(x=Bare1m, y=totalrichnessallperiod)) +
  geom_point(size=2) +
  geom_smooth(method ="glm", se=TRUE, fullrange=TRUE,
              method.args = list(family = "poisson"), col="gray")+
  geom_point() +
  theme_classic()+ 
  xlab("Bare Soil") +
  ylab("")
ba4<-ggplot(data=site.avg, aes(x=Urban2kmSQRT, y=totalrichnessallperiod)) +
  geom_point(size=2) +
  geom_smooth(method ="glm", se=TRUE, fullrange=TRUE,
              method.args = list(family = "poisson"), col="gray")+
  geom_point() +
  theme_classic()+ 
  xlab("Urban (2km)") +
  ylab("")

panel.act<-grid.arrange(ba1,ba2,ba3,ba4,ncol=4,nrow=1)
ggsave("BeeRichnessTotal.png", panel.act, dpi=300, height=2, width=10)

#replot, with richness at site at each survey period
by.site <- read.csv("bySite.csv",header=TRUE)

ba1<-ggplot(data=by.site, aes(x=GardenSizeLN, y=Richness)) +
  geom_point(size=2) +
  geom_smooth(method ="glm", se=TRUE, fullrange=TRUE,
              method.args = list(family = "gaussian"), col="gray")+
  geom_point() +
  theme_classic() +
  xlab("Garden Size") +
  ylab("Bee Richness")
ba2<-ggplot(data=by.site, aes(x=FlowersLN, y=Richness)) +
  geom_point(size=2) +
  geom_smooth(method ="glm", se=TRUE, fullrange=TRUE,
              method.args = list(family = "gaussian"), col="gray")+
  geom_point() +
  theme_classic() + 
  xlab("No. Flowers") +
  ylab("")
ba3<-ggplot(data=by.site, aes(x=Bare1m, y=Richness)) +
  geom_point(size=2) +
  geom_smooth(method ="glm", se=TRUE, fullrange=TRUE,
              method.args = list(family = "gaussian"), col="gray")+
  geom_point() +
  theme_classic()+ 
  xlab("Bare Soil") +
  ylab("")
ba4<-ggplot(data=by.site, aes(x=Urban2kmSQRT, y=Richness)) +
  geom_point(size=2) +
  geom_smooth(method ="glm", se=TRUE, fullrange=TRUE,
              method.args = list(family = "gaussian"), col="gray")+
  geom_point() +
  theme_classic()+ 
  xlab("Urban (2km)") +
  ylab("")

panel.act<-grid.arrange(ba1,ba2,ba3,ba4,ncol=4,nrow=1)
ggsave("BeeRichnessPeriod.png", panel.act, dpi=300, height=2, width=10)

