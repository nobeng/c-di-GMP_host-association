################################
#   Morphotype colonization    #
################################

#### = 1. Set dependencies + open data ####
# Set working directory to source location
currentDirectory <- dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(currentDirectory)

# Specify libraries & source files
library(plyr); library(ggplot2); library(gridExtra); library(lawstat); library(multcomp); library(multcompView)
library(lme4); library(lmerTest); library(emmeans); library(gtools)
library(finalfit)
source("C:/Users/Nancy/Documents/phd/data/experimental evolution/G. pseudomonas diversification project/B) Student sub-projects/data/my_theme.R")

# General variables
colsTreat <- c("#8d163cff", "#2e4057ff")
colsTyps <- c()
pops <- 1:6 # replicate populations assayed

#### = 2. Import data = ####
earlyColonization <- read.csv("earlyColonization_acrossMTs.csv")
persistAndRelease <- read.csv("persistAndRelease_acrossMTs.csv")

#### = 3. Plot data ####
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

# Functions to plot CFUs
plotCol <- function(data, assay, lifeCycle){
  data <- as.data.frame(data)
  # Subset data to test according to assay 
  ifelse(assay == "earlyCol", 
         d <- data,
         ifelse(assay == "shortPersist", 
                { d <- subset(data, select = c(pop, Replicate, morphotype, treatE, morph.simple, CFUs_Col))
                colnames(d)[colnames(d) == "CFUs_Col"] <- "CFUs" }, 
                ifelse(assay == "release", 
                       { d <- subset(data, select = c(pop, Replicate, morphotype, treatE, morph.simple, CFUs_Rel))
                       colnames(d)[colnames(d) == "CFUs_Rel"] <- "CFUs" },
                       "Assay not recognized.")))
  # Subset data from respective life cycle
  ifelse(lifeCycle == "mono", 
         { lifeC <- "bi"
         colsTreat <- colsTreat[2]}, 
         { lifeC <- "mono"
         colsTreat <- colsTreat[1]}  )
  
  d <- d[d$treatE != lifeC, ]  
  
  # Prepare data to plot on log10-scale
  d$CFUs[d$CFUs == 0]<-1
  d$CFUs[d$CFUs < 0]<-NA
  d <- d[!is.na(d$CFUs),]
  
  # Reverse MT order to ensure consistent shape for ancestral MT48
  revMorphotypes <- rev(mixedsort(as.character(unique(d$morphotype))))
  d$morphotype <- factor(d$morphotype, levels = revMorphotypes)
  
  # Prepare plot
  plots <- list(length = length(unique(unique(data$pop))))
  counter <- 1
  
  for (p in 1:length(pops)){
    dPlot <- d[d$pop == pops[p],] 
    spaceOnTop <- max(dPlot$CFUs)*3.5
    
    # Collect ancestral mean and SE for plottin as dashed line + exclude from data to be plotted
    meanANC <- mean(dPlot[dPlot$treatE == "anc", 'CFUs'])
    seANC <- sd(dPlot[dPlot$treatE == "anc", 'CFUs'])/sqrt(length(dPlot[dPlot$treatE == "anc", 'CFUs']))
    
    dPlot <- dPlot[dPlot$treatE != "anc", ]
    
    # Set population label within plot      
    ifelse(lifeCycle == "bi", 
           labelPop <- paste("Pop. T", dPlot$pop[1], sep =""), 
           labelPop <- paste("Pop. C", dPlot$pop[1], sep =""))
    
    dPlot$morphotype <- as.factor(dPlot$morphotype)
    
    # Use indirect log10 scale for y axis, due to problems with geom_rect and scale_y_log10
    plots[[counter]] <-
      ggplot(dPlot, aes(x = morph.simple, y = log10(CFUs), col = treatE)) + 
        geom_rect(xmin = 0 , xmax = length(unique(dPlot$morph.simple))+1,
                  ymin = log10(meanANC - seANC),
                  ymax = log10(meanANC + seANC), col = NA, fill = "lightgrey", alpha = 0.5)+
      geom_point(alpha = 0)+ # place holder plotting to allow putting annotate layer on bottom w/o errors
      geom_hline(yintercept = log10(meanANC), lty = 2, size = 1, alpha = 0.5)+
      geom_boxplot(outliers.shape = NA)+
      stat_summary(data = dPlot,
                   fun = mean,
                   fun.min = function(x) mean(x) - sd(x), 
                   fun.max = function(x) mean(x) + sd(x), 
                   geom = "pointrange",alpha = 0.5, size = 1, 
                   aes(group = morphotype),
                   position = position_jitter(height=0, width = 0.1))+
      #   geom_point(size = 4, alpha = 0.5,
      #            position = position_jitterdodge(jitter.width = 1, dodge.width = 0.5)) +
      scale_y_continuous(limits = c(1,6), labels = fancy_scientific)+
              labs(x = "Morphology", y = "CFU/worm")+
      #ggtitle(labelPop)+
      scale_color_manual(values = "black")+
      my_theme(dPlot)+
      theme(legend.position = "none", 
        axis.text.x = element_text(hjust = 1, angle = 45))
    
    # Advance counter to cycle through time points (tps)
    counter <- counter + 1
  }
  
  # Select directory for plots
  dirPlots <- paste(currentDirectory, "/plots/", sep = "")
  
  # Save plot to file
  ggsave(paste(dirPlots, assay, "~morphology_",toupper(lifeCycle),".png",sep =""),
         do.call("grid.arrange", c(plots, nrow=6)), 
         dpi = 300, height = 20, width = 2)
  
  ggsave(paste(dirPlots, assay, "~morphology_",toupper(lifeCycle),".svg",sep =""),
         do.call("grid.arrange", c(plots, nrow=6)), 
         dpi = 300, height = 20, width = 2)
  
  return(print(plots)) }

# Run plotting functions
datasets <- list(earlyColonization, persistAndRelease, persistAndRelease)
assays <- c("earlyCol", "shortPersist", "release")
treats <- c("mono", "bi")

persistAndRelease <- rbind(persistAndRelease, 
                           data.frame(pop = 3, Replicate = 1, morphotype = "fuzzy", treatE = "bi", 
                                      morphology = "fuzzy", morph.simple = "fuzzy", 
                                      CFUs_Col = NA, CFUs_Rel = 0))
for(a in 1:length(assays)){
  allPlots <- data.frame()
  for (i in 1:length(treats)){
      allPlots <- c(allPlots, plotCol(datasets[a], assays[a], treats[i])) }
  
  # Save plot to file
  fileExtensions <- c(".png", ".pdf")
  for(i in 1:length(fileExtensions)){
    ggsave(paste("colonization~morphologySIMPLE_",assays[a] ,"_all_simple", fileExtensions[i], sep =""),
           do.call("grid.arrange", 
                   c(allPlots[7], allPlots[8],allPlots[9], allPlots[10],allPlots[11], allPlots[12], 
                     allPlots[1], allPlots[2], allPlots[3], allPlots[4], allPlots[5], allPlots[6],
                     ncol=6)), 
           dpi = 300, height = 5, width = 12) }}

# Clean up 
rm(plotCol, colsTreat, colsTyps, dirUnblind, my_theme)

#### = 4. Statistical analysis #####
runLMM_morphology <- function(data, assay, lifeCycle){
  data <- as.data.frame(data)
  # Subset data to test according to assay 
  ifelse(assay == "earlyCol", 
         d <- data,
         ifelse(assay == "shortPersist", 
                { d <- subset(data, select = c(pop, Replicate, morphotype, treatE, morph.simple, CFUs_Col))
                colnames(d)[colnames(d) == "CFUs_Col"] <- "CFUs" }, 
                ifelse(assay == "release", 
                       { d <- subset(data, select = c(pop, Replicate, morphotype, treatE, morph.simple, CFUs_Rel))
                       colnames(d)[colnames(d) == "CFUs_Rel"] <- "CFUs" },
                       "Assay not recognized.")))
  # Work on log10-scale
  d$CFUs <- d$CFUs+1
  d$CFUs[d$CFUs < 0] <- NA
  d <- d[!is.na(d$CFUs),]
  
  # Subset life cycle
  ifelse(lifeCycle == "mono", 
         lifeC <- "bi",
         lifeC <- "mono")
  
  d <- d[d$treatE != lifeC, ]
  
  # Set file name for export of summary
  fileName <- paste("stats/", assay, "_LMM-Tukey.effectAcrossPops_",toupper(lifeCycle), ".txt", sep ="")
  
  # Fit model and check for normality of residuals
  m <- lmer(log10(CFUs) ~ morph.simple + (1|pop), data = d)
  
  marginal <- emmeans(m, ~ morph.simple)
  
  # Save results to file
  sink(fileName)
  print(anova(m))
  print("")
  print(as.data.frame(pairs(marginal)))
  sink()
  
  return("Effect of morphology across pops: GLMM + Tukey saved to file.")
}
runPairwise_withinPops <- function(data, assay, lifeCycle, pop){
  data <- as.data.frame(data)
  # Subset data to test according to assay 
  ifelse(assay == "earlyCol", 
         d <- data,
         ifelse(assay == "shortPersist", 
                { d <- subset(data, select = c(pop, Replicate, morphotype, treatE, morph.simple, CFUs_Col))
                colnames(d)[colnames(d) == "CFUs_Col"] <- "CFUs" }, 
                ifelse(assay == "release", 
                       { d <- subset(data, select = c(pop, Replicate, morphotype, treatE, morph.simple, CFUs_Rel))
                       colnames(d)[colnames(d) == "CFUs_Rel"] <- "CFUs" },
                       "Assay not recognized.")))
  
  # Subset life cycle and agar concentration
  ifelse(lifeCycle == "mono", 
         lifeC <- "bi",
         lifeC <- "mono")
  
  d <- d[d$treatE != lifeC & d$pop == pop, ]  
  
  # Set file name for export of summary
  fileName <- paste("stats/", assay, "_pairwiseWilcoxWithinPops",toupper(lifeCycle),"_pop", pop,".txt", sep ="")
  fileName2 <- paste("stats/", assay, "_AOV-WithinPops",toupper(lifeCycle),"_pop", pop,".txt", sep ="")
  
  # Work on log10-scale
  d$CFUs <- d$CFUs+1
  d$CFUs[d$CFUs < 0] <- NA
  d <- d[!is.na(d$CFUs),]
  
  # Run (G)LM
  # Fit model and check for normality of residuals
  d$morph.simple <- as.factor(d$morph.simple)
  m <- lm(log10(CFUs) ~ morph.simple, data =d)# family = "Gamma", data = d)
  #m <- glm(log10(CFUs) ~ morph.simple, family = "Gamma", data = d)
  
  # Tukey post-hoc contrasts
  marginal <- emmeans(m, ~ morph.simple)
  
  sink(fileName2)
  print("ANOVA")
  print(anova(m))
  
  print("TUKEY")
  print(as.data.frame(pairs(marginal)))
  
  print("TUKEY LETTERS")
  print(cld(glht(m, mcp(morph.simple="Tukey"))))
  sink()
  
  
  return("Pairwise within pops saved to file.")
}

for(a in 1:length(assays)){
  for (t in 1:length(treats)){
    runLMM_morphology(datasets[a], assays[a], treats[t])
    for(p in pops){
      runPairwise_withinPops(datasets[a], assays[a], treats[t], pops[p]) }}}

#### Final clean up ####
rm(runPairwise_withinPops, runLMM_morphology, datasets, earlyColonization, persistAndRelease, assays, currentDirectory, i, pops)