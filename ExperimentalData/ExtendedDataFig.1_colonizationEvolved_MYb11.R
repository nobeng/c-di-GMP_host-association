###################################
#                                 #
# Quantify evolved colonization   #
#       in worm population        #
#                                 #
###################################

#### 1. Set dependencies ####
# Set working directory to source location
currentDirectory <- dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(currentDirectory)
source("my_theme.R")
library(plyr); library(ggplot2); library(reshape); library(stringr); library(lme4); library(multcomp);
library(dplyr)

colorsMYb11 <- c("#9E163C", "#CF1332","#2e4057","#4A678C","#66947C")

#### 2. Read data ####
cfuWorm <- read.csv("cfuPerWormInPop_MYb11.csv")

#### 3. Plot abundance of CFUs/worm in population samples #####

# Color scheme
colsTreat <- c("grey25", rgb(0, 158,115, maxColorValue = 255), rgb(230, 159,0, maxColorValue = 255))

# Proper axis labels
fancy_scientific <- function(l) {
  # turn in to character string in scientific notation
  l <- format(l, scientific = TRUE)
  # quote the part before the exponent to keep all the digits
  l <- gsub("e\\+","e",l)
  # turn the 'e+' into plotmath format
  l <- gsub("e", "%*%10^", l)
  # return this as an expression
  parse(text=l)}

# Ensure universal scale labels
scale_log10_labels <- function(l) {
  # turn in to character string in scientific notation
  l <- format(l, scientific = F)
  # return this as an expression
  parse(text=l)}

# Plot Bacterial population size/ worm (within pop)
ggplot(cfuWorm, aes(x = toupper(treatE), y = cfuWorm))+
  geom_boxplot(aes(fill = treatE), alpha = 0.85)+
  geom_jitter(size = 1, alpha = 0.7, width = 0.1)+
  geom_hline(yintercept = max(cfuWorm$cfuWorm)*2, alpha = 0)+ # place holder for space of significance indication
  scale_y_log10(labels = fancy_scientific)+
  scale_fill_manual(values = colsTreat) +
  labs(x = "Evo. life cycle", y = "Mean CFU/ worm") +
  my_theme(cfuWorm)

# Save plot to file
fileExtensions <- c(".svg", ".png")
for(i in 1:length(fileExtensions)){
  ggsave(paste("cfuWormPop~treatE_perWorm_MYb11", fileExtensions[i], sep=""), dpi = 300, height = 2.5, width = 3)  
}

#### 4. Analyze statistics ####
testDifferencesAbs <- function(data, subsetName){
  # Prepare file name
  dirStatsOut <- "./stats out"
  fileName <- paste(dirStatsOut, "/colPop_abs.tests_GLM_C10_",subsetName, "_MYb11.txt", sep="")
  
  # Run tests and save to file alongside
  sink(fileName)
  
  # 6.1. Are evolved populations from ancestral?
  print("6.1. Are evolved populations from ancestral?")
  
  m <- glm(log10(cfu) ~ treatE, family = "Gamma", data = data)
  print(summary(m))
  
  print("ANOVA")
  print(anova(m, test = "Chisq"))
  
  print("post-hoc: TUKEY")
  print(summary(glht(m, mcp(treatE="Tukey"))))
  
  print("post-hoc: DUNNET")
  print(summary(glht(m, mcp(treatE="Dunnet"))))
  
  # Close file
  sink()
  
  return("Test outcomes saved to file.")
}
testDifferencesAbs(dataWorm, "perWorm")

# Outcome:  Biphasic differ from ancestral and monophasic.
#           Amount of bacteria is rather low given that we are looking at a whole population.