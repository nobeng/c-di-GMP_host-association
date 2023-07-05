########################################
# Analysis of colonization experiments #
#    4. ANAYLSIS: RELATIVE FITNESS     #
#        EVOLVED VS. ANCESTRAL         #
########################################

#### == 1. Set dependencies == ####
# Set working directory to source location
currentDirectory <- dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(currentDirectory)

# Specify libraries
library(plyr); library(ggplot2); library(lawstat); library(emmeans)

# Source files
source("my_theme.R")

# Variables
cycs <- c(4, 10)
treats <- c("bi", "mono")
colsMYb <- c("#FC2F62", "grey27", "#189E00", "blue")
colsTreat <- c("grey25", rgb(0, 158, 115, maxColorValue = 255), rgb(230, 159, 0, maxColorValue = 255))

#### == 2. Import data == ####
# Import early colonization data
dataLoad <- read.csv("col.2h_relData.csv", header = T)
# Import persistence and release data
dataReleas <- read.csv("releas_relData.csv", header = T)

#### == 3. Function definitions == ####

# Ensure universal scale labels
scale_log10_labels <- function(l) {
  # turn in to character string in scientific notation
  l <- format(l, scientific = F)
  # return this as an expression
  parse(text=l)}

# Plot relative bacterial load per worm (evo/anc)
plotRelLoad <- function(data, assayName, cycle, cols, dimW, dimH){
  # Set cycle and subset respective cycle and lawn
  ifelse(cycle == 10, cyc <- 4, cyc <- 10)
  data <- data[data$cycle != cyc, ]
  
  # Leave 20% space on top of max. data point for significance stars later
  top <- max(data$diff_cfuWorm)*1.2
  
  # Print plot
  p <- ggplot(data, aes(x = toupper(treatE), y = diff_cfuWorm)) +
    geom_hline(yintercept = 1, lty=2, col = "grey50", size = 1)+
    geom_boxplot(aes(fill = treatE), alpha = 0.85, outlier.color = NA)+
    geom_jitter(size = 1, alpha = 0.7, width = 0.1)+
    geom_hline(yintercept =max(data$diff_cfuWorm)*1.05, alpha = 0)+ # keep space for significance
    scale_fill_manual(values = colsTreat[2:3])+
    scale_y_log10(labels=scale_log10_labels)+
    labs(x = "Treatment", y = expression("Evo/anc"))+
    my_theme(data)

  # Save plot to directory with rel. loads
  fileExtensions <- c(".png", ".svg")
  for(i in 1:length(fileExtensions)){
    
    # Define file name & save to file
    fileName <- paste(assayName, "_rel.load~treatE_C", cycle, "_MYb11", fileExtensions[i], sep = "")
    ggsave(fileName, p, height=dimH, width=dimW,dpi=300)
  }
  
  return(p)}
# Plot relative bacteria released per worm (evo/anc)
plotRelReleas <- function(data, assayName, cycle, cols, dimW, dimH){
  # Set cycle and subset respective cycle and lawn
  ifelse(cycle == 10, cyc <- 4, cyc <- 10)
  data <- data[data$cycle != cyc, ]
  
  # Leave 20% space on top of max. data point for significance stars later
  top <- max(data$diff_cfuReleas)*1.2
  
  # Print plot
  p <- ggplot(data, aes(x = toupper(treatE), y = diff_cfuReleas)) +
    geom_hline(yintercept = 1, lty=2, col = "grey50", size = 1)+
    geom_boxplot(aes(fill = treatE), alpha = 0.85, outlier.color = NA)+
    geom_jitter(size = 1, alpha = 0.7, width = 0.1)+
    geom_hline(yintercept =max(data$diff_cfuReleas)*1.05, alpha = 0)+ # keep space for significance
    scale_fill_manual(values = colsTreat[2:3])+
    scale_y_log10(labels=scale_log10_labels)+
    labs(x = "Treatment", y = "Evo/anc")+
    my_theme(data)
  
  # Save plot to directory with rel. loads
  fileExtensions <- c(".png", ".svg")
  for(i in 1:length(fileExtensions)){
    # Define file name
    fileName <- paste(assayName, "_rel.releas~treatE_C", cycle, "_MYb11", fileExtensions[i], sep = "")
  
    # Save plot to directory with rel. loads
    ggsave(fileName, p, height=dimH, width=dimW,dpi=300)
  }
  
  return(p)}

# Run one-sample t-tests
runOneSample.tTest.LOAD <- function(data,assayName, c, t){
  
  # Subset data
  d <- with(data, data[cycle == cycs[c] & treatE == treats[t], 'diff_cfuWorm'])
  
  # Get log10 and exclude inf.
  d <- log10(d)[is.finite(log10(d))]
  
  # Run t-test
  ttestOUT <- t.test(d, mu = 0, var.equal = T, alternative = "two.sided")
  
  # Return data frame with test details and outcomes
  return(data.frame(assay = assayName,
                    cycle = cycs[c],
                    treatE = treats[t],
                    phenotype = "diff_cfuWorm",
                    t = signif(as.numeric(ttestOUT$statistic),3),
                    df = as.numeric(ttestOUT$parameter),
                    pval = ttestOUT$p.value))
}
runOneSample.tTest.RELEASE <- function(data,assayName, c, t){
  
  # Subset data
  d <- with(data, data[cycle == cycs[c] & treatE == treats[t], 'diff_cfuReleas'])
  
  # Get log10 and exclude inf.
  d <- log10(d)[is.finite(log10(d))]
  
  # Run t-test
  ttestOUT <- t.test(d, mu = 0, var.equal = T, alternative = "two.sided")
  
  # Return data frame with test details and outcomes
  return(data.frame(assay = paste(assayName, "releas", sep="_"),
                    cycle = cycs[c],
                    treatE = treats[t],
                    phenotype = "diff_cfuReleas",
                    t = signif(as.numeric(ttestOUT$statistic),3),
                    df = as.numeric(ttestOUT$parameter),
                    pval = ttestOUT$p.value))
}

# Subset data from a specific evol. cycle
subsetCycle <- function(data, cycle){
  ifelse(cycle == 10,
         cyc <- 4, 
         cyc <- 10)
  
  return(data[data$cycle != cyc, ])
}

#### == 4. Plot rel. fitness == ####
for (i in 1:length(cycs)){ # Loop through evol. cycles
  plotRelLoad(dataLoad, "earlyCol", cycs[i], colsMYb, "MYb11", 2.5, 3.2)
  plotRelLoad(dataRelease, "persistence", cycs[i], colsMYb, "MYb11", 2.5, 3.2)
   
  # Also plot CFU released, when appropriate
  plotRelReleas(dataReleas, "release", cycs[i], colsMYb, "MYb11", 2.5, 3.2)}

#### == 5. Calculate one-sample t-tests == ####
# Pre-allocate data frame
resultsTtest <- data.frame()

for(c in 1:length(cycs)){ # Evol. cycles
    for(t in 1:length(treats)){ # Evol. treatments
        resultsTtest <- rbind(resultsTtest, runOneSample.tTest.LOAD(dataLoad, "earlyCol", c, t))
        resultsTtest <- rbind(resultsTtest, runOneSample.tTest.LOAD(dataLoad, "persistence", c, t))
        resultsTtest <- rbind(resultsTtest, runOneSample.tTest.RELEASE(dataReleas, "release", c, t)) }}

# Correct p-values for multiple testing
resultsTtest$pval.FDR <- rep(100, length(resultsTtest$assay))

for(c in 1:length(cycs)){ # Evol. cycles
  for(t in 1:length(treats)){ # Evol. treatments
        
        # Subset pvalues to be corrected for multiple comparisons
        pvalUncor <- with(resultsTtest, resultsTtest[cycle == cycs[c], 'pval'])
        
        # Correct using false-discovery rate
        pvalCor <- p.adjust(pvalUncor, method = "fdr")
        
        # Attach to overall data frame
        resultsTtest[resultsTtest$cycle == cycs[c], 'pval.FDR'] <- pvalCor
      }}

# Save test results to file
# Check presence of directory for stats results
statsOut <- paste(currentDirectory, "/stats out", sep="")
ifelse(dir.exists(statsOut) == F,
       dir.create(statsOut), "Directory for statistics exists.")

# Define file name
fileName <- paste(statsOut, "/", assayName, "_relativeFitness_oneSample.tTests_MYb11.csv", sep = "")

# Write to file
write.csv(resultsTtest, fileName, row.names = F)