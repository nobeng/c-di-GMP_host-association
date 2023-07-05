#################################
#                               #
#     MY316 sizes @ L4 + 3d     #
#         (evo/anc.)            #
#                               #
#################################

#### = 1. Dependencies, directories and variables = ####
library(plyr); library(ggplot2); library(gridExtra); library(lawstat); library(emmeans)

# Set working directory to source location
currentDirectory <- dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(currentDirectory)

# Source files
source("my_theme.R")

# Variables
colsMYb <- c("#FC2F62", "grey27", "#189E00", "blue")
colsTreat <- c("grey25", rgb(0, 158, 115, maxColorValue = 255), rgb(230, 159, 0, maxColorValue = 255))
assayName <- "wormSize"
cycs <- c(4, 10)

#### = 2. Import worm data = ####
data <- read.csv("wormPopSize_MYb11_relData.csv", sep=",", header = T)

#### = 2. Function definitions = ####

# Function to plot relative phenotypes
plotRelDim <- function(data, cycle, cols, strains, dims,labelY, dimW, dimH){
  # Set cycle and subset respective cycle and lawn
  ifelse(cycle == 10, cyc <- 4, cyc <- 10)
  data <- data[data$cycle != cyc, ]
  
  # Leave 20% space on top of max. data point for significance stars later
  top <- max(data[, dims])*1.05
  
  # Print plot
  p <- ggplot(data, aes(x = toupper(treatE), y = data[, dims])) +
    geom_hline(yintercept = 1, lty=2, col = "grey50", size = 1)+
    geom_hline(yintercept =max(data[, dims])*1.05, alpha = 0)+ # keep space for significance
    geom_boxplot(aes(fill = treatE), alpha = 0.85)+
    geom_jitter(size = 1, alpha = 0.7, width = 0.1)+
    scale_fill_manual(values = colsTreat[2:3]) +
    xlab("Life cycle") + ylab(bquote(.(labelY)*.[(evo/anc)]))+
    my_theme(d)
  
  # Check presence of directory
  relPlots <- paste(currentDirectory, "/relative phenotype plots", sep="")
  ifelse(dir.exists(relPlots) == F,
         dir.create(relPlots), "Directory for relative data already exists.")

  fileExtensions <- c(".png", ".svg")
  for(i in 1:length(fileExtensions)){
    # Define file name
    fileName <- paste(relPlots, "/", assayName, "_", dims,"~treatE_C", cycle, "_", strains, fileExtensions[i], sep = "")
  
    # Save plot to directory
    ggsave(fileName, p, height=dimH, width=dimW,dpi=300) }
  
  return(p)}

# Function to run one-sample t-tests
runOneSample.tTest <- function(data, dims, s, c, t){
  
  # Subset data
  d <- with(data, data[strain == "MYb11" & cycle == cycs[c] & treatE == treats[t], ])
  d <- d[, dims]

  # Run t-test
  ttestOUT <- t.test(d, mu = 1, var.equal = T, alternative = "two.sided")
  
  # Return data frame with test details and outcomes
  return(data.frame(assay = assayName,
                    strain = "MYb11",
                    cycle = cycs[c],
                    treatE = treats[t],
                    phenotype = dims,
                    t = signif(as.numeric(ttestOUT$statistic),3),
                    df = as.numeric(ttestOUT$parameter),
                    pval = ttestOUT$p.value))
}

#### = 8. Run plotting & statistics = ####
treats <- unique(data$treatE)

# Pre-allocate data frame
resultsTtest <- data.frame()

for(c in 1:length(cycs)){ # Evol. cycles
    for(t in 1:length(treats)){ # Evol. treatments
        resultsTtest <- rbind(resultsTtest, runOneSample.tTest(data, "relnWorms", "MYb11", c, t))
        plotRelDim(data[data$strain == "MYb11",], cycs[c], colsMYb, "MYb11", "relnWorms", "# worms",2.5, 3.2)
        }}

# Correct for multiple testing of single ancestor
resultsTtest$pval.FDR <- rep(1, length(resultsTtest$assay))
for(c in 1:length(cycs)){ # Evol. cycles
    # Subset pvalues to be corrected for multiple comparisons
      pvalUncor <- with(resultsTtest, resultsTtest[phenotype == "relnWorms" & cycle == cycs[c] & strain == "MYb11", 'pval'])
      
      # Correct using false-discovery rate
      pvalCor <- p.adjust(pvalUncor, method = "fdr")
      
      # Attach to overall data frame
      resultsTtest[resultsTtest$phenotype == "relnWorms" & resultsTtest$cycle == cycs[c] & resultsTtest$strain == "MYb11", 'pval.FDR'] <- pvalCor
}

# Save test results to file
# Check presence of directory for stats results
statsOut <- paste(currentDirectory, "/stats out", sep="")
ifelse(dir.exists(statsOut) == F,
       dir.create(statsOut), "Directory for statistics exists.")

# Define file name
fileName <- paste(statsOut, "/", assayName, "_relWormSizes_oneSample.tTests", "_MYb11.csv", sep = "")

# Write to file
write.csv(resultsTtest, fileName, row.names = F)

#### = Clean up = ####
rm(data, resultsTtest, assayName, c, colsMYb, colsTreat, currentDirectory, 
   cycs, fileName, statsOut, t, treats, my_theme, plotRelDim, runOneSample.tTest, 
   ttestOUT, dims, evolChange,pvalCor, pvalUncor)
