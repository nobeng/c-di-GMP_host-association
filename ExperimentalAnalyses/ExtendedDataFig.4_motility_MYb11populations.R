###########################################
#             Bacterial motility          #
#               Sept. 2019                #
#              Nancy Obeng                #
###########################################

#### = 1. Dependencies, directories and variables = ####
library(ggplot2);library(plyr); library(gtools); library(lawstat); library(emmeans)

# Set working directory to source location
currentDirectory <- dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(currentDirectory)

# Source files
source("my_theme.R")

# Variables
colsTreat <- c("grey25", rgb(0, 158,115, maxColorValue = 255), rgb(230, 159,0, maxColorValue = 255))
cycs <- c(4, 10)
assayName <- "motility"

#### Import data ####
dataRel <- read.csv("motility_MYb11populations.csv", header = T)

#### === Plot: Evol. of max. colony diameter === ####
colsMYb <- c("#FC2F62", "grey27", "#189E00")

# Ensure universal scale labels
scale_log10_labels <- function(l) {
  # turn in to character string in scientific notation
  l <- format(l, scientific = F)
  # return this as an expression
  parse(text=l)}

plotEvo <- function(data, strains, cagar, cyc,tp, h, w){
  ifelse(cyc == 4, excyc <- "10", excyc <- "4")
  data <- data[data$cycle.x != excyc &
         data$agar == cagar &
         data$strain%in%strains &
         data$tp == tp, ]
  
  p <- ggplot(data,aes(x = toupper(treatE.x), y = relTrans)) +
    geom_hline(yintercept = 1, lty=2, col = "grey50", size = 1)+
    geom_point(aes(x = "BI", y= max(data$relTrans)*1.05), alpha = 0)+ # keep space for significance
    geom_boxplot(aes(fill = treatE.x), alpha = 0.85)+
    geom_jitter(size = 1, alpha = 0.7, width = 0.1)+
    scale_fill_manual(values = colsTreat[2:3]) +
    xlab("Life cycle") + ylab(bquote("Colony diameter"[(evo/anc)]))+
    my_theme(data)
  
  # Define file name
  fileExtensions <- c(".png", ".svg")
  for(i in 1:length(fileExtensions)){
    # Check presence of directory
    relPlots <- paste(currentDirectory, "/relative phenotype plots", sep="")
    ifelse(dir.exists(relPlots) == F, 
          dir.create(relPlots), "Directory for relative data already exists.")
  
    ggsave(paste(relPlots, "/", assayName,"~treatE_", cyc,"_",strains,"_",cagar,"_",tp,fileExtensions[i], sep=""),
           p, height=h, width=w, dpi=300)
  }
  }

# Plotting
for (c in c(4,10)){
  plotEvo(dataRel,"MYb11",0.5,c,"24h", 3.2,2.5) # Plot 0.5% agar
  
  for(t in unique(dataRel$tp)){
    plotEvo(dataRel,"MYb11",3.4,c,t, 3.2,2.5) # Plot 3.4% agar time points
    }}

# Clean up
rm(plotEvo,scale_log10_labels, my_theme)

#### === Statistics: One-sided t-tests of relative motility === ####
dataRel <- dataRel[dataRel$strain == "MYb11", ]

runOneSample.tTest <- function(data, cyc, cagar, strains, tp, treatE){
  # Subset data
  ifelse(cyc == 4, excyc <- "10", excyc <- "4")
  
  d <- data[data$cycle.x != excyc &
                 data$agar == cagar &
                 data$strain%in%strains &
                 data$tp == tp &
                 data$treatE.x == treatE, 
            'relTrans']
  
  # Get mean for possible direction of evolutionary change
  evolChange <- "two.sided"
  
  # Run t-test
  ttestOUT <- t.test(d, mu = 1, var.equal = T, alternative = evolChange)
  
  # Return data frame with test details and outcomes
  return(data.frame(assay = "motility",
                    agar = cagar,
                    strain = strains,
                    cycle = cyc,
                    treatE = treatE,
                    tp = tp,
                    evolChange = evolChange,
                    t = signif(as.numeric(ttestOUT$statistic),3),
                    df = as.numeric(ttestOUT$parameter),
                    pval = ttestOUT$p.value))
  }

cycles <- c(4, 10)
treats <- c("bi", "mono")
tps <- unique(dataRel$tp)

resultsTtest <- data.frame()

for (l in 1:length(treats)){
  for (c in 1:length(cycles)){

  resultsTtest <- rbind(resultsTtest, runOneSample.tTest(dataRel, cycles[c], 0.5, "MYb11", "24h", treats[l])) # test 0.5% agar
  
  for(t in 1:length(tps)){
    resultsTtest <- rbind(resultsTtest, runOneSample.tTest(dataRel, cycles[c], 3.4, "MYb11", tps[t], treats[l])) # test 3.4% agar
  }}}


# Correct for multiple testing
resultsTtest$pval.FDR <- rep(100, length(resultsTtest$assay))

for (l in 1:length(treats)){
  for (c in 1:length(cycles)){
    
    # 0.5 % agar
    # Subset pvalues to be corrected for multiple comparisons
    pvalUncor <- with(resultsTtest, resultsTtest[cycle == cycles[c] & agar == 0.5 & tp == "24h", 'pval'])
    # Correct using false-discovery rate
    pvalCor <- p.adjust(pvalUncor, method = "fdr")
    
    # Attach to overall data frame
    resultsTtest[resultsTtest$cycle == cycles[c] & resultsTtest$agar == 0.5 & resultsTtest$tp == "24h", 'pval.FDR'] <- pvalCor
    
    for(t in 1:length(tps)){
      # 3.4 % agar
      # Subset pvalues to be corrected for multiple comparisons
      pvalUncor <- resultsTtest[resultsTtest$cycle == cycles[c] & resultsTtest$agar == 3.4 & resultsTtest$tp == tps[t], 'pval']
      # Correct using false-discovery rate
      pvalCor <- p.adjust(pvalUncor, method = "fdr")
      
      # Attach to overall data frame
      resultsTtest[resultsTtest$cycle == cycles[c] & resultsTtest$agar == 3.4 & resultsTtest$tp == tps[t], 'pval.FDR'] <- pvalCor
    }}}

# Save test results to file
# Check presence of directory for stats results
statsOut <- paste(currentDirectory, "/stats out", sep="")
ifelse(dir.exists(statsOut) == F,
       dir.create(statsOut), "Directory for statistics exists.")

# Define file name
fileName <- paste(statsOut, "/", assayName, "_relativeFitness_oneSample.tTests_MYb11.csv", sep = "")

# Write to file
write.csv(resultsTtest, fileName, row.names = F)

##### Clean up ####
rm(statsOut, dataRel, resultsTtest, assayName, colsMYb, colsTreat, currentDirectory, cycles, cycs, tps, treats, runOneSample.tTest)