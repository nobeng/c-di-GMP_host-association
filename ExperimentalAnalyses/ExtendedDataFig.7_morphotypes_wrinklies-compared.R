################################
# Wrinkly phenotypes compared  #
#                              #
#    Open, plot and analyze    #
#                              #
################################

#### = 1. Set dependencies = ####
# Set working directory to source location
currentDirectory <- dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(currentDirectory)
library(ggplot2); library(multcomp); library(gridExtra); library(DescTools)

#### = 2. Read in data = ####
biofilms <- read.csv("morphotypes_biofilm.csv", header=T)
dispersal <- read.csv("morphotypes_dispersal.csv", header=T)

biofilms <- biofilms[biofilms$morphotype%in%c("MT48", "MT12", "MT14", "MT22"), ]
biofilms <- biofilms[biofilms$pop%in%c(3,5), ]

dispersal <- dispersal[dispersal$morphotype%in%c("MT48", "MT12", "MT14", "MT22"), ]
dispersal <- dispersal[dispersal$pop%in%c(3,5), ]

#### = 3. Plot data = ####

### Biofilm formation in vitro ###
plotBiofilms <- ggplot(biofilms, aes(x=morphotype, y=od))+
  geom_boxplot()+
  geom_abline(intercept = median(biofilms[biofilms$morphotype == "MT48", 'od']), slope=0, lty= 2)+
  geom_jitter(width=0.05, height=0)+
  labs(x="Evolved", y="Biofilm formation (OD600)")+
  theme_bw()
  
### Swarming (24h) ###

ggplot(dispersal, aes(x=morphotype, y=swarm24h, col = factor(pop)))+
  geom_boxplot()+
  geom_abline(intercept = median(dispersal[dispersal$pop == 3 & dispersal$morphotype == "MT48", 'swarm24h']), slope=0, lty= 2)+
  geom_jitter(width=0.05, height=0)+
  labs(x="Evolved", y="Swarming at 24h (cm)")+
  theme_bw()

swarm3 <- ggplot(dispersal[dispersal$pop == 3, ], aes(x=morphotype, y=swarm24h))+
  geom_boxplot()+
  geom_abline(intercept = median(dispersal[dispersal$pop == 3 & dispersal$morphotype == "MT48", 'swarm24h']), slope=0, lty= 2)+
  geom_jitter(width=0.05, height=0)+
  labs(x="Evolved", y="Swarming at 24h (cm)")+
  theme_bw()

swarm5 <- ggplot(dispersal[dispersal$pop == 5, ], aes(x=morphotype, y=swarm24h))+
  geom_boxplot()+
  geom_abline(intercept = median(dispersal[dispersal$pop == 5 & dispersal$morphotype == "MT48", 'swarm24h']), slope=0, lty= 2)+
  geom_jitter(width=0.05, height=0)+
  labs(x="Evolved", y="Swarming at 24h (cm)")+
  theme_bw()

plotSwarming <- grid.arrange(swarm3,swarm5, nrow=1)

### Colony expansion (3d) ###
plotColExpand <- ggplot(dispersal, aes(x=morphotype, y=expand3d))+
  geom_boxplot()+
  geom_abline(intercept = median(dispersal[dispersal$morphotype == "MT48", 'expand3d']), slope=0, lty= 2)+
  geom_jitter(width=0.05, height=0)+
  labs(x="Evolved", y="Colony expansion at 3d (cm)")+
  theme_bw()

#### Plot all same sizes 
svg("wrinklies-compared_biofilm-swarming-colexpansion.svg", 
    width = 5, height = 7)

grid.arrange(plotBiofilms, plotSwarming, plotColExpand)
dev.off()

#### = 4. Analyse data = ####

### Biofilm formation in vitro ###
biofilms$morphotype <- factor(biofilms$morphotype, levels=c("MT48", "MT12", "MT14", "MT22"))
m <- lm(od ~ morphotype+pop, data=biofilms)
m <- lm(od ~ morphotype, data=biofilms)
anova(m)

DunnettTest(x=biofilms$od, g=biofilms$morphotype) # clear differences btw. wrinkly isolates
TukeyHSD(aov(m), which = 'morphotype')

cld(glht(m, mcp(morphotype="Tukey")))

### Swarming (24h) ###
dispersal$morphotype <- factor(dispersal$morphotype, levels=c("MT48", "MT12", "MT14", "MT22"))

# Differences between populations assayed in different runs? --> yes, thus, plot and analyze separately
m <- lm(swarm24h ~ morphotype+pop, data=dispersal)
anova(m)

# Swarm pop3 
m <- lm(swarm24h ~ morphotype, data=dispersal[dispersal$pop ==3, ])
anova(m)
DunnettTest(x=dispersal[dispersal$pop ==3, ]$swarm24h,
            g=dispersal[dispersal$pop ==3, ]$morphotype) # Differences only visible in comparison to ancestral reference
TukeyHSD(aov(m), which = 'morphotype')
cld(glht(m, mcp(morphotype="Tukey")))

#Swarm pop5
t.test(dispersal[dispersal$pop ==5 & dispersal$morphotype == "MT22", 'swarm24h'], 
       dispersal[dispersal$pop ==5 & dispersal$morphotype == "MT48", 'swarm24h'])

### Colony expansion (3d) ###
dispersal$morphotype <- factor(dispersal$morphotype, levels=c("MT48", "MT12", "MT14", "MT22"))
m <- lm(expand3d ~ morphotype+pop, data=dispersal) # no differences detectable
anova(m)

m <- lm(expand3d ~ morphotype, data=dispersal) # no differences detectable
anova(m)

DunnettTest(x=dispersal$expand3d, g=dispersal$morphotype) # Differences only visible in comparison to ancestral reference
TukeyHSD(aov(m), which = 'morphotype')

cld(glht(m, mcp(morphotype="Tukey")))