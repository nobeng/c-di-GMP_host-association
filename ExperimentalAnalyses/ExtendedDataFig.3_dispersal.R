################################
#     Morphotype dispersal     #
################################

#### = 1.Set dependencies + open data ####
# Set working directory to source location
currentDirectory <- dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(currentDirectory)

# Specify libraries & source files
library(plyr); library(ggplot2); library(gridExtra); library(lawstat); library(multcomp); library(multcompView)
library(lme4); library(lmerTest); library(emmeans)
library(finalfit)
source("my_theme.R")

# General variables
colsTreat <- c("#8d163cff", "#2e4057ff")
colsTyps <- c()

#### = 2. Import data = ####
data <- read.csv("morphotypes_dispersal_acrossMTs.csv")

#### = 3. Function definitions = ####
plotExpansionAcrossMorphologies <- function(data, agarConc, lifeCycle){
  # Subset data from respective agar concentration & life cycle
  ifelse(lifeCycle == "mono",
         { lifeC <- "bi"
           colsTreat <- colsTreat[2]  },
         { lifeC <- "mono"
           colsTreat <- colsTreat[1] }  )
  d <- data[data$agar == agarConc & data$treatE != lifeC, ]
  ifelse(agarConc == "5 g/l",
         {agar <- "5gAgar"
         d <- d[d$tp == "24h", ]},
         {agar <- "34gAgar"
         d <- d[d$tp == "3d", ]})

  # Reverse MT order to ensure consistent shape for ancestral MT48
  d$morphotype <- factor(d$morphotype, levels = rev(levels(d$morphotype)))

  # Prepare plot
  plots <- list(length = length(unique(unique(data$tp))))
  counter <- 1

  for (p in 1:length(pops)){

      dPlot <- d[d$pop == pops[p],]
      spaceOnTop <- max(dPlot$Length)*1.25

      # Collect ancestral mean and SE for plotting as dashed line + exclude from data to be plotted
      meanANC <- mean(dPlot[dPlot$treatE == "anc", 'Length'])
      seANC <- sd(dPlot[dPlot$treatE == "anc", 'Length'])/sqrt(length(dPlot[dPlot$treatE == "anc", 'Length']))

      dPlot <- dPlot[dPlot$treatE != "anc", ]

      ifelse(agarConc == "5 g/l",
             assaySpecificMax <- 6,
             assaySpecificMax <- 2.6)

      # Set population label within plot and set of shapes for the population
      ifelse(lifeCycle == "bi",
             labelPop <- paste("Pop. B", dPlot$pop[1], sep =""),
             labelPop <- paste("Pop. M", dPlot$pop[1], sep =""))

      plots[[counter]] <- ggplot(dPlot, aes(x = morph.simple, y = Length, col = treatE)) +
                          geom_rect(xmin = 0 , xmax = 6, ymin = meanANC - seANC, ymax = meanANC + seANC, col = NA, fill = "lightgrey", alpha = 0.5)+
                          geom_hline(yintercept = meanANC, lty = 2, size = 1, alpha = 0.5)+
                          #geom_hline(yintercept = spaceOnTop, alpha = 0)+
                          geom_boxplot()+
        stat_summary(data = dPlot,
                     fun = mean,
                     fun.min = function(x) mean(x) - sd(x),
                     fun.max = function(x) mean(x) + sd(x),
                     geom = "pointrange",alpha = 0.5, size = 1,
                     aes(group = morphotype),
                     position = position_jitter(height=0, width = 0.1))+
                          labs(x = "Morphology", y = "Colony diam. (cm)")+
                          scale_color_manual(values = "black")+
        scale_y_continuous(limits=c(0, assaySpecificMax))+

        my_theme(dPlot)+
                          theme(
                            legend.position = "none",
                            axis.text.x = element_text(hjust = 1, angle = 45))

      # Advance counter to cycle through time points (tps)
      counter <- counter + 1
    }

  # Save plot to file
  ggsave(paste("dispersal_diameter~morphologySIMPLE_",toupper(lifeCycle),"_",agar ,"_2022.png",sep =""),
         do.call("grid.arrange", c(plots, nrow=6)),
         dpi = 300, height = 20, width = 3.5)

  ggsave(paste("dispersal_diameter~morphologySIMPLE_",toupper(lifeCycle),"_",agar ,"_2022.svg",sep =""),
         do.call("grid.arrange", c(plots, nrow=6)),
         dpi = 300, height = 20, width = 3.5)

  return(plots)
}

#### = 4. Run plotting functions = ####
for (i in 1:length(agars)){
  for(t in 1:length(treats)){
  plotExpansionAcrossMorphologies(data, agars[i], treats[t]) }}

allPlots1 <- c(plotExpansionAcrossMorphologies(data, agars[1], treats[2]), plotExpansionAcrossMorphologies(data, agars[1], treats[1]))
allPlots2 <- c(plotExpansionAcrossMorphologies(data, agars[2], treats[2]), plotExpansionAcrossMorphologies(data, agars[2], treats[1]))

# Save plot to file
fileExtensions <- c(".png", ".pdf")
for(i in 1:length(fileExtensions)){
  allPlots <- allPlots1
  ggsave(paste("swarm~morphologySIMPLE_all_2022_simple", fileExtensions[i], sep =""),
         do.call("grid.arrange", 
                 c(allPlots[7], allPlots[8],allPlots[9], allPlots[10],allPlots[11], allPlots[12], 
                   allPlots[1], allPlots[2], allPlots[3], allPlots[4], allPlots[5], allPlots[6],
                   ncol=6)), 
         dpi = 300, height = 5, width = 9)

  allPlots <- allPlots2
  ggsave(paste("colExpand~morphologySIMPLE_all_2022_simple", fileExtensions[i], sep =""),
         do.call("grid.arrange", 
                 c(allPlots[7], allPlots[8],allPlots[9], allPlots[10],allPlots[11], allPlots[12], 
                   allPlots[1], allPlots[2], allPlots[3], allPlots[4], allPlots[5], allPlots[6],
                   ncol=6)), 
         dpi = 300, height = 5, width = 9)}

# Clean up
rm(plotExpansionAcrossMorphologies, my_theme, colsTreat, colsTyps, i, mts, tps, t)


#### = 5. Statistical testing: max. expansion ~ morphology = ####

# Check for differences across populations
runLMM_Morphology <- function(data, agarConc, lifeCycle){
  # Subset life cycle and agar concentration
  ifelse(lifeCycle == "mono", 
         lifeC <- "bi",
         lifeC <- "mono")
  
  d <- data[data$agar == agarConc & data$treatE != lifeC, ]  
  
  # Focus on 24h of swarming (5g/l agar) and colony expansion at 3d (34g/l agar)
  ifelse(agarConc == "5 g/l", 
         
         {agar <- "5gAgar"
         tp <- "24h"
         d <- d[d$tp == "24h", ]},
         
         {agar <- "34gAgar"
         tp <- "3d"
         d <- d[d$tp == "3d", ]})
  
  # Export summary
  fileName <- paste("statsOut/effectMorphology_acrossPops/dispersal_LMM.effectMorphWithinTreatE_", lifeCycle, "_", agar, "_", tp, ".txt", sep = "")
  
  # Fit model and check for normality of residuals
  m <- lmer(Length ~ morph.simple + (1|pop), data = d)
  
  marginal <- emmeans(m, ~ morph.simple)
  
  # Export summary
  #fileName <- paste("stats/biofilm_Tukey_effectAcrossPops",toupper(lifeCycle),".txt", sep ="")
  
  # Save results to file
  sink(fileName)
  print(anova(m))
  print("")
  print(as.data.frame(pairs(marginal)))
  sink()
  # Return model for Tukey post-hoc contrasts
  return(m)
  }

# Post hoc contrasts of morphologies across populations
runTukeyPostHoc_Morphology_AcrossPops <- function(m, agarConc, lifeCycle){
  marginal <- emmeans(m, ~ morphology)
  
  # Export summary
  ifelse(agarConc == "5 g/l",
         {agar <- "5gAgar"
         tp <- "24h"},
         {agar <- "34gAgar"
         tp <- "3d"})
  
  fileName <- paste("statsOut/effectMorphology_acrossPops/dispersal_Tukey.effectMorphWithinTreatE_", agar, "_", lifeCycle,"_", tp, ".txt", sep = "")
  
  # Save to file
  write.csv(pairs(marginal), fileName, row.names = F)
  sink(fileName)
  print(as.data.frame(pairs(marginal)))
  sink()
  
  return("Tukey saved to file.")
}

# Check for differences within populations
runPairwiseTtest <- function(data, agarConc, lifeCycle, pop){
  # Subset life cycle and agar concentration
  ifelse(lifeCycle == "mono", 
         lifeC <- "bi",
         lifeC <- "mono")
  
  d <- data[data$agar == agarConc & data$treatE != lifeC & data$pop == pop, ]  
  
  # Focus on 24h of swarming (5g/l agar) and colony expansion at 3d (34g/l agar)
  ifelse(agarConc == "5 g/l", 
         
         {agar <- "5gAgar"
         tp <- "24h"
         d <- d[d$tp == "24h", ]},
         
         {agar <- "34gAgar"
         tp <- "3d"
         d <- d[d$tp == "3d", ]})
  
  # Export summary
  fileName <- paste("statsOut/dispersal_pairwiseTtestWithinPops_agar", agar,"_", toupper(lifeCycle),"_pop", pop,".txt", sep ="")
  
  # Run pairwise t-test
  pw <- pairwise.t.test(d$Length, d$morphology, p.adjust.method = "fdr", alternative = "two.sided")
  
  # Transform matrix of p-values
  tri.to.squ<-function(x) # by Fabio Marroni (https://fabiomarroni.wordpress.com/2017/03/25/perform-pairwise-wilcoxon-test-classify-groups-by-significance-and-plot-results/)
  {
    rn<-row.names(x)
    cn<-colnames(x)
    an<-unique(c(cn,rn))
    myval<-x[!is.na(x)]
    mymat<-matrix(1,nrow=length(an),ncol=length(an),dimnames=list(an,an))
    for(ext in 1:length(cn))
    {
      for(int in 1:length(rn))
      {
        if(is.na(x[row.names(x)==rn[int],colnames(x)==cn[ext]])) next
        mymat[row.names(mymat)==rn[int],colnames(mymat)==cn[ext]]<-x[row.names(x)==rn[int],colnames(x)==cn[ext]]
        mymat[row.names(mymat)==cn[ext],colnames(mymat)==rn[int]]<-x[row.names(x)==rn[int],colnames(x)==cn[ext]]
      }
      
    }
    return(mymat)
  }
  
  mymat <- tri.to.squ(pw$p.value) 
  
  # Get letters for grouping & save to file
  sink(fileName)
  print(pw)
  print(multcompLetters(mymat,compare="<=",threshold=0.05,Letters=letters))
  sink()
  
  return("Pairwise Wilcox saved to file.")
}
runLM_Tukey_WithinPops <- function(data, agarConc, lifeCycle, pop){
  # Subset life cycle and agar concentration
  ifelse(lifeCycle == "mono", 
         lifeC <- "bi",
         lifeC <- "mono")
  
  d <- data[data$agar == agarConc & data$treatE != lifeC & data$pop == pop, ]  
  d$morph.simple <- as.factor(d$morph.simple)
  
  # Focus on 24h of swarming (5g/l agar) and colony expansion at 3d (34g/l agar)
  ifelse(agarConc == "5 g/l", 
         
         {agar <- "5gAgar"
         tp <- "24h"
         d <- d[d$tp == "24h", ]},
         
         {agar <- "34gAgar"
         tp <- "3d"
         d <- d[d$tp == "3d", ]})
  
  # Export summary
  fileName <- paste("statsOut/dispersal_ANOVA-Tukey_agar", agar,"_", toupper(lifeCycle),"_pop", pop,".txt", sep ="")
  
  
  # Fit model and check for normality of residuals
  m <- lm(Length ~ morph.simple, data = d)
  
  # Save results to file
  sink(fileName)
  print("ANOVA")
  print(anova(m))
  
  print("TUKEY")
  marginal <- emmeans(m, ~ morph.simple)
  print(as.data.frame(pairs(marginal)))
  
  print("TUKEY LETTERS")
  print(cld(glht(m, mcp(morph.simple="Tukey"))))
  sink()
  
  
  return("Within pop. comparisons saved to file.")
  
  }

# 7.5. Cycle through tests
for (a in 1:length(agars)){

  for (t in 1:length(treats)){
   runLMM_Morphology(data, agars[a], treats[t])
    
     for(p in pops){
       #runLM_Tukey_WithinPops(data, agars[a], treats[t], pops[p])
       #runTukeyPostHoc_Morphology_WithinPops(runANOVA_Morphology_WithinPops(data, agars[a], treats[t], pops[p]), agars[a], treats[t], pops[p])
}}}