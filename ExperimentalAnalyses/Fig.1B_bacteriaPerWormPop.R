###################################
#                                 #
#            Figure 1B            # 
#      Bacterial colonization     #
#                                 #
###################################

#### 1. Set dependencies ####
# Set working directory to source location
currentDirectory <- dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(currentDirectory)
source("my_theme.R")
library(plyr); library(ggplot2); library(reshape); library(stringr); library(lme4); library(multcomp);
library(dplyr)

#### 2. Read data ####
dataWorm <- read.csv("colWormPop_relData.csv")

#### 3. Plot relative load (evo/anc) ####
plotRelLoad <- function(data, nameSubset){
  # Prepare plot
  ggplot(data, aes(x = toupper(treatE), y = rel.cfu.perWormInPop)) +
    geom_hline(yintercept = 1, lty=2, col = "grey50", size = 1)+
    geom_hline(yintercept =max(data$rel.cfu.perWormInPop)*1.05, alpha = 0)+ # keep space for significance
    geom_boxplot(aes(fill = treatE), alpha = 0.85,outlier.shape = NA)+
    geom_jitter(size = 1, alpha = 0.7, width = 0.1)+
    scale_fill_manual(values = colsTreat[2:3]) +
    labs(x = "Treatment", y = "Evo/anc")+
    my_theme(data)
  
  # Save to file
  fileExtensions <- c(".png", ".svg")
  for(i in 1:length(fileExtensions)){
    ggsave(paste("rel.cfuWormPop~treatE_10_", nameSubset, "_MYb11", fileExtensions[i], sep=""), height=3.2, width=2.5, dpi=300) 
  }
  
  return("Relative load plotted to file.")
}

plotRelLoad(dataWorm, "perWorm")

#### 4. Statistics: relative load (evo/anc) ####

testDifferencesRel <- function(data, subsetName){
  # Prepare file name
  dirStatsOut <- "./stats out"
  fileName <- paste(dirStatsOut, "/colPop_rel.tests_tTest-ANOVA_C10_",subsetName, "_MYb11.txt", sep="")
  
  # Run tests and save to file alongside
  sink(fileName)
  
  # 9.1. Does either treatment differ from the ancestor?
  print("Does either treatment differ from the ancestor?")
  print("MONO vs ANC")
  print(t.test(log10(data[data$treatE == "mono", 'diff_cfu']), mu = 0, alternative = "two.sided"))
  
  print("BI vs ANC")
  print(t.test(log10(data[data$treatE == "bi", 'diff_cfu']), mu = 0, alternative = "two.sided"))
  
  print("corrected pvalues (FDR) - mono, bi")
  print(p.adjust(c(t.test(log10(data[data$treatE == "bi", 'diff_cfu']), mu = 0, alternative = "two.sided")$p.val,
                   t.test(log10(data[data$treatE == "mono", 'diff_cfu']), mu = 0, alternative = "two.sided")$p.val),
                 method = "fdr"))
  
  # 9.2. Do treatments differ in their evolutionary response?
  print("Do treatments differ in their evolutionary response?")
  m <- lm(log10(diff_cfu) ~ treatE, data = data)
  print(anova(m))
  print("")
  print(summary(glht(m, mcp(treatE="Tukey"))))
  
  # Close file
  sink()
  
  return("Test outcomes saved to file.")
}

testDifferencesRel(dataWorm, "perWorm")

# Outcome: Only biphasic treatment differs from the ancestor.
# Outcome: Yes, There is a significant difference between biphasically and monophasically evolved populations.