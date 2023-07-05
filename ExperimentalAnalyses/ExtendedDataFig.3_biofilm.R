################################
# Morphotype biofilms in vitro #
################################

#### = 1. Set dependencies = ####
# Set working directory to source location
currentDirectory <- dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(currentDirectory)

# Libraries & source files
library(plyr); library(ggplot2); library(gridExtra); library(reshape);library(DescTools); library(gtools); 
library(emmeans); library(lawstat); library(lme4); library(lmerTest); library(multcompView); library(multcomp); library(rcompanion)
source("my_theme.R")

# Global parameters
nRepPops <- 6
colsTreat <- c(rgb(0, 158,115, maxColorValue = 255), rgb(230, 159,0, maxColorValue = 255))
colsTreat <- c("#8d163cff", "#2e4057ff")

#### = 2. Import data = ####
data <- read.csv("morphotypes_biofilm.csv")

#### = 3. Plotting data = ####
# Relevant data
pops <- sort(unique(data$pop))
treats <- c("mono", "bi")

plotBiofilm <- function(data, lifeCycle){
  # Subset data from respective life cycle
  ifelse(lifeCycle == "mono", 
         { lifeC <- "bi"
           colsTreat <- colsTreat[2]
           ymax <- 0.5}, 
         { lifeC <- "mono"
           colsTreat <- colsTreat[1]
           ymax <- 0.85 }  )
  d <- data[data$treatE != lifeC, ]  
  
  # Reverse MT order to ensure consistent shape for ancestral MT48
  revMorphotypes <- rev(mixedsort(as.character(unique(d$morphotype))))
  d$morphotype <- factor(d$morphotype, levels = revMorphotypes)
  
  # Prepare plot
  plots <- list(length = length(unique(unique(data$tp))))
  counter <- 1
  
  for (p in 1:length(pops)){
    
    dPlot <- d[d$pop == pops[p],] 

    # Collect ancestral mean and SE for plottin as dashed line + exclude from data to be plotted
    meanANC <- mean(dPlot[dPlot$treatE == "anc", 'od_cor'])
    seANC <- sd(dPlot[dPlot$treatE == "anc", 'od_cor'])/sqrt(length(dPlot[dPlot$treatE == "anc", 'od_cor']))
    
    dPlot <- dPlot[dPlot$treatE != "anc", ]
    
    # Set population label within plot      
    ifelse(lifeCycle == "bi", 
           labelPop <- paste("Pop. B", dPlot$pop[1], sep =""), 
           labelPop <- paste("Pop. M", dPlot$pop[1], sep =""))
    
    plots[[counter]] <-ggplot(dPlot, aes(x = morph.simple, y = od_cor, col = treatE)) + 
      geom_rect(xmin = 0 , xmax = 6, ymin = meanANC - seANC, ymax = meanANC + seANC, col = NA, fill = "lightgrey", alpha = 0.5)+
      geom_hline(yintercept = meanANC, lty = 2, size = 1, alpha = 0.5)+
      geom_boxplot()+
      stat_summary(data = dPlot,
                   fun = mean,
                   fun.min = function(x) mean(x) - sd(x), 
                   fun.max = function(x) mean(x) + sd(x), 
                   geom = "pointrange",alpha = 0.5, size = 1, 
                   aes(group = morphotype),
                   position = position_jitter(height=0, width = 0.1))+
      scale_y_continuous(limits=c(-0.01, 0.85))+
      labs(x = "Morphology", y = "OD(600nm)")+
      scale_color_manual(values = "black")+
      my_theme(dPlot)+
      theme(
        legend.position = "none", 
        axis.text.x = element_text(hjust = 1, angle = 45))
    
    # Advance counter to cycle through time points (tps)
    counter <- counter + 1
    
  }
  
  # Save plot to file
  fileExtensions <- c(".png", ".svg")
  for(i in 1:length(fileExtensions)){
    ggsave(paste("biofilm~morphologySIMPLE_2022_",toupper(lifeCycle), fileExtensions[i], sep =""),
           do.call("grid.arrange", c(plots, nrow=6)), 
           dpi = 300, height = 20, width = 2) }
  
  return(print(plots))
  
}

allPlots <- c(plotBiofilm(data, treats[2]), plotBiofilm(data, treats[1]))

# Save plot to file
fileExtensions <- c(".png", ".pdf")
for(i in 1:length(fileExtensions)){
  ggsave(paste("biofilm~morphologySIMPLE_all_2022_simple", fileExtensions[i], sep =""),
         do.call("grid.arrange", 
                 c(allPlots[7], allPlots[8],allPlots[9], allPlots[10],allPlots[11], allPlots[12], 
                   allPlots[1], allPlots[2], allPlots[3], allPlots[4], allPlots[5], allPlots[6],
                   ncol=6)), 
         dpi = 300, height = 5, width = 9.5)}

# Clean up
rm(colsTreat, my_theme, plotBiofilm)

#### = 4. Statistical analysis = ####
# Check for differences across populations
runLMM_Morphology <- function(data, lifeCycle){
  # Subset life cycle
  ifelse(lifeCycle == "mono", 
         lifeC <- "bi",
         lifeC <- "mono")
  
  d <- data[data$treatE != lifeC, ]  
  
  # Export summary
  fileName <- paste("stats/biofilm_LMM.effectAcrossPops-withTukey",toupper(lifeCycle),".txt", sep ="")
  
  # Fit model and check for normality of residuals
  m <- lmer(od_cor ~ morph.simple + (1|pop), data = d)
  
  marginal <- emmeans(m, ~ morph.simple)
  
  # Save results to file
  sink(fileName)
  print(anova(m))
  print("")
  print(as.data.frame(pairs(marginal)))
  sink()
  
  # Return model for Tukey post-hoc contrasts
  return(m)
}

# Posthoc contrasts of morphologies across populations
runTukeyPostHoc_Morphology_AcrossPops <- function(m, lifeCycle){
  marginal <- emmeans(m, ~ morph.simple)
  
  # Export summary
  fileName <- paste("stats/biofilm_Tukey_effectAcrossPops",toupper(lifeCycle),".txt", sep ="")
  
  # Save to file
  write.csv(pairs(marginal), fileName, row.names = F)
  sink(fileName)
  print(as.data.frame(pairs(marginal)))
  sink()
  
  return("Tukey saved to file.")
}

# Check for differences within populations
runLM_Tukey_WithinPops <- function(data, lifeCycle, numPop){
  # Subset life cycle
  ifelse(lifeCycle == "mono", 
         lifeC <- "bi",
         lifeC <- "mono")
  
  d <- data[data$treatE != lifeC & data$pop == numPop, ]
  
  # Export summary
  fileName <- paste("stats/MT simple categories/biofilm_ANOVA_Tukey.WithinPops",toupper(lifeCycle),"_", numPop, ".txt", sep ="")
  
  # Fit model and check for normality of residuals
  m <- lm(od_cor ~ morph.simple, data = d)
  
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
  
  return("Stats saved to file.")
}


# 7.5. Cycle through tests
treats <- c("mono", "bi")
pops <- 1:6

for (t in 1:length(treats)){
  runLMM_Morphology(data, treats[t])
    
  for(p in pops){
    runLM_Tukey_WithinPops(data, treats[t], pops[p])
    }}