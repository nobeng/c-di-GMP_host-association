###################################
#                                 #
# Quantify morphotype proportions #
#         per population          #
#                                 #
#   = bacteria per worm pop.  =   #
#                                 #
###################################

#### 1. Set dependencies ####
# Set working directory to source location
currentDirectory <- dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(currentDirectory)
source("my_theme.R")
library(plyr); library(ggplot2); library(reshape); library(stringr); library(lme4); library(multcomp);
library(dplyr)

colorsMYb11 <- c("#9E163C","#CF1332","#2e4057","#4A678C","#66947C")

#### 2. Import data ####
#data <- read.csv("cfuPerWormInPop_absData.csv", header = T)
d <- read.csv("cfuPerWormInPop_by-morphotype_MYb11.csv", header = T)
d$pop <- paste(toupper(substring(d$treatE,1,1)),d$pop, sep="")

#### 3. Plot data ####
# Focus on evolved
dProps <- d[d$treatE != "anc", ]

# Calculate proportion of morphotypes in population
dTotals <- dProps %>%
  group_by(assay, treatE, pop, sample, rep) %>%
  summarize(cfuTotal = sum(cfu_WormByMT))

# Attach population totals + calculate proportion
dProps <- merge(dProps, dTotals, by = c("assay", "treatE", "pop", "sample", "rep"))
dProps$prop <- dProps$cfu_WormByMT/dProps$cfuTotal

# Keep zeros where no colonies where detected per type
dProps[dProps$cfu_WormByMT == 0, 'prop'] <- 0

dProps <- subset(dProps, select=-strain)


# Attach zero entries to morphotypes that were not detected
absentTypes <- dProps[dProps$treatE == "mono" & dProps$morphotype == "smooth",]
absentTypes$morphotype <- rep("wrinkly", length(absentTypes$assay)) # no wrinklies emerged in free-living treat.
absentTypes$prop <- 0

# Based on missing entries for absent morphotypes in table dProps:
absentTypes <- rbind(absentTypes, 
                     data.frame(assay=rep("load",24), 
                                treatE=c(rep("bi",17),rep("mono", 7)), 
                                pop=c("B1",rep("B2",3), "B3", rep("B4",4), rep("B5",4),
                                      rep("B6",4),rep("M2",3), "M4", "M5", "M5", "M6"),
                                sample=c(23,rep(17,3), 2, rep(15,4), rep(2,4),rep(12,4),
                                         rep(6,3),25,28,28,26),
                                rep=c(1,1,2,3,1,1,2,3,3,1,1,2,3,1,2,3,3,1,2,3,2,2,3,3),
                                morphotype=c("smooth", rep("fuzzy", 7), rep("smooth",2),
                                             rep("fuzzy",3),rep("smooth", 3),rep("fuzzy",8)), 
                                cfu_WormByMT = rep(0,24),
                                cfuTotal=rep(0,24),
                                prop= rep(0,24)))
dProps <- rbind(dProps, absentTypes)

# Average of replicates
dPropsReps <- dProps
dProps <- dProps %>%
  group_by(assay, treatE,morphotype, pop, sample) %>%
  summarize(prop = mean(prop))

# Reorder morphotypes
dProps$morphotype <- factor(dProps$morphotype, levels = c("wrinkly","smooth", "fuzzy"))

p<-ggplot(dProps, aes(x=morphotype, y=prop))+
  geom_boxplot(aes(fill= morphotype), outlier.shape=NA)+
  geom_jitter(size = 1, alpha = 0.7, width = 0.1)+
  geom_point(aes(x=1, y=1.1), alpha = 0)+ # buffer for significance indication
  scale_y_continuous(breaks=c(0,0.5,1))+
  facet_wrap(~treatE)+
  labs(x = "Morphology", y = "Proportion in host")+
  scale_fill_manual(values = c("#ee9b33ff", "#7f7f7fff", "#006111ff"))+
  my_theme(d)+
  theme(axis.text.x = element_text(hjust=1, angle =45)); p

# plot without wrinklies in mono
p<-ggplot(dProps[!(dProps$morphotype == "wrinkly" & dProps$treatE == "mono"), ], aes(x=morphotype, y=prop))+
  geom_boxplot(aes(fill= morphotype), outlier.shape=NA)+
  geom_jitter(size = 1, alpha = 0.7, width = 0.1)+
  geom_point(aes(x=1, y=1.1), alpha = 0)+ # buffer for significance indication
  scale_y_continuous(breaks=c(0,0.5,1))+
  facet_wrap(~treatE)+
  labs(x = "Morphology", y = "Proportion in host")+
  scale_fill_manual(values = c("#ee9b33ff", "#7f7f7fff", "#006111ff"))+
  my_theme(d)+
  theme(axis.text.x = element_text(hjust=1, angle =45)); p

# Save plot to file
dirAbsPlots <- "C:/Users/Nancy/Documents/phd/data/experimental evolution/C. post-ee experiments part II/E. data/analysis/morphotypes in population"
fileExtensions <- c(".svg", ".png")
for(i in 1:length(fileExtensions)){
  ggsave(paste(dirAbsPlots, "/propsMT~Pop_MYb11_boxplot", fileExtensions[i], sep=""),
         dpi = 300, height = 3.5, width = 3.5)  }

# Statistical test
# For mono/control populations: exclude absent wrinklies
dPropsReps$morphotype <- as.factor(dPropsReps$morphotype)

m <- glm(prop ~ morphotype, family = quasibinomial, data=dPropsReps[dPropsReps$treatE == "bi",])

# run post-hoc (Tukey)
tuk <- glht(m, linfct = mcp(morphotype="Tukey"))
summary(tuk, adjusted("fdr"))

m <- glm(prop ~ morphotype, family = quasibinomial, data=dPropsReps[dPropsReps$treatE == "mono" & dPropsReps$morphotype != "wrinkly",])
summary(m)

# run post-hoc (Tukey)
tuk <- glht(m, linfct = mcp(morphotype="Tukey"))
summary(tuk, adjusted("fdr"))