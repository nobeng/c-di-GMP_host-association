###################################
#                                 #
#   Morphotype props. over time   #
#       = growth on agar =        #
#                                 #
###################################

#### 1. Set dependencies ####
# Set working directory to source location
currentDirectory <- dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(currentDirectory)
dirStats <- "./stats out"

source("my_theme.R")
library(plyr); library(ggplot2); library(reshape); library(stringr); library(lme4);
library(tseries); library(multcomp); library(gamlss); library(emmeans)

#### 2. Collect data ####
data <- rbind(read.csv("growth_agar_morphotypes_t0.csv", header = T),
              read.csv("growth_agar_morphotypes_t24h.csv", header = T),
              read.csv("growth_agar_morphotypes_thalfC.csv", header = T),
              read.csv("growth_agar_morphotypes_tfullC.csv", header = T))

# Add zero counts for morphotypes at tps, at which they are not detected (only if they have been observed at at least one time point)
X <- 12 # numberNewEntries
addZeros <- data.frame(pop = c(3, 3, 5, 5, 5, 6, 6, 6, 2, 4, rep(5, 2)), 
                       sample = c(rep(2, 5), 12, 12, 12, 6, 25, rep(28,2)), 
                       assay = c(rep("growth_agar", X)),
                       tp = c(rep(c("t0", "24h"), 2), "fullC", "t0", "24h", "halfC", "24h","halfC","t0","fullC"), 
                       morphotype = c(rep("fuzzy", 5), rep("smooth", 3), rep("fuzzy",4)),
                       cfu = c(rep(0, X)), 
                       treatE = c(rep("bi", 8), rep("mono", 4)), 
                       strain = c(rep("MYb11", X)), 
                       cycle = c(rep(10, X)))

data <- rbind(data, addZeros)
rm(addZeros)

# Sum up morpotypes to get proportions per sample
dataSums <- ddply(data, .(pop, sample, assay, tp, treatE, strain, cycle), summarise, cfuTotal = sum(cfu))
data <- merge(data, dataSums, by = c("pop", "sample", "assay", "tp", "treatE", "strain", "cycle"))
data$prop <- data$cfu/data$cfuTotal
rm(dataSums)

#### 3. Format data frame for plotting ####

# Call replicate populations by # and treatE
popNames <- paste(toupper(substring(data$treatE,1,1)), data$pop, sep="")
data$pop <- popNames

# Add time scale in hours for proper x-axis
data$hours <- rep(1, length(data$tp))

tps <- unique(data$tp)
hours <- c(24, 168, 72, 0)

for(i in 1:length(tps)){
  data[data$tp%in%tps[i], 'hours'] <- hours[i] }

# Create population ID across samples, to connect lines
data$id <- with(data, paste(pop, sample, morphotype))
data$pop <- as.factor(data$pop)

#### 4. Print plots ####
# Color scheme
colsTreat <- c("grey25", rgb(0, 158,115, maxColorValue = 255), rgb(230, 159,0, maxColorValue = 255))
colsMT <- c("#FF4F76", "grey40", "#36A5A5")
colsMT <- c("#006111ff", "#7f7f7fff", "#ee9b33ff")

#### Plot temporal dynamics of morphotypes ####
plotDynamics <- function(data, treatE){
  # Subset data of selected life cycle treatment
  data <- data[data$treatE == treatE, ]
  
  # Prepare plot
  p <- ggplot(data, aes(x = hours, y = prop, col = morphotype, shape = pop, fill = morphotype))+
    geom_line(aes(group = id), lty = 5, alpha = 0.8)+
    geom_point(size = 3, fill = "white")+
    geom_point(size = 3)+
    scale_x_continuous(breaks=c(0, 24, 72, 168), labels=c("0", "24", "72", "168"))+
    scale_y_continuous(labels=scales::percent, limits = c(0, 1))+
    labs(x = "Time (h)", y = "Morph. composition", col = "Morphology", shape = "Replicate\npopulation")+
    scale_shape_manual(values =c(7, 21:25))+ #c(22, 21, 24, 15:20))+
    scale_color_manual(values = colsMT)+
    scale_fill_manual(values = colsMT)+
    my_theme(data)+
    theme(legend.position = "right")
  
  # Save to file
  fileExtensions <- c(".png", ".svg")
  for(i in 1:length(fileExtensions)){
    ggsave(paste("morphotypeProps_agar~time_", treatE,"2022", fileExtensions[i], sep = ""),
           p, dpi = 300, height = 5, width = 4.5) }
  
  return(p)
  }

plotDynamics(data, "bi")
plotDynamics(data, "mono")

#### Plot temporal dynamics of whole populations ####
dataTotal <- ddply(data, .(pop, tp, treatE, strain, cycle), summarize, cfu = sum(cfu))
dataTotal$pop <- substring(dataTotal$pop,2,2)
dataTotal$id <- paste(dataTotal$pop, dataTotal$treatE, sep="_")

# Add time scale in hours for proper x-axis
dataTotal$hours <- rep(1, length(dataTotal$tp))

tps <- unique(dataTotal$tp)
hours <- c(0, 24, 72, 168)

for(i in 1:length(tps)){
  dataTotal[dataTotal$tp%in%tps[i], 'hours'] <- hours[i] }

# Save to file
absDir <- "C:/Users/Nancy/Documents/phd/data/experimental evolution/C. post-ee experiments part II/E. data/analysis/absolute phenotype data"
#write.csv(dataTotal, paste(absDir,"/", "growthAgar_absData_MYb11_overtime.csv",sep=""), row.names = F)

dataTotal <- read.csv(paste(absDir,"/", "growthAgar_absData_popsMYb11_overtime.csv",sep=""))

ggplot(dataTotal, aes(x=hours, y=cfu, col=treatE, fill=treatE))+
  geom_point(size = 4, alpha = 0.8)+
  geom_line(aes(group=dataTotal$id), lty=2, size = 1)+
  scale_y_log10(labels=fancy_scientific)+
  scale_x_continuous(breaks=c(0, 24, 72, 168), labels=c("0", "24", "72", "168"))+
  labs(x = "Time (hours)", y = "CFU/population",col = "Evo. life cycle")+
  scale_color_manual(values=c("black", "#8d163cff", "#2e4057ff"))+
  scale_fill_manual(values=c("black", "#8d163cff", "#2e4057ff"))+
  facet_wrap(~toupper(treatE))+
  my_theme(dataTotal)
  

# Save to file
fileExtensions <- c(".png", ".svg")
for(i in 1:length(fileExtensions)){
  ggsave(paste("growth_agar~time_2022colors", fileExtensions[i], sep = ""),
         dpi = 300, height = 4, width = 9) }

#### 6. Test for statistical differences over time ####

# Note, the inoculum (i.e.t0) of a plate is connected to the subsequent time points. The others are not a classic time series, as I sampled repl. plates.

##### 6.1. Do the proportions of the morphotypes change over time? (beta regression) ####
# Note: require a model with a closed interval [0,1], as I have both zeros and 1s in the dataset
#       -> zero-inflated beta regression

betaReg_propsOverTimePerTreatE <- function(data, treatE){
  # Subset data
  d <- data[data$treatE == treatE, ]
  
  # Set up zero- and one-inflated beta-regression
  m <- gamlss(prop ~ hours * morphotype, family = BEINF, data = d)
  
  # Define fileName for saving data
  fileName <- paste("propMT~time_betaRegression_", toupper(treatE), sep="")
  
  # Save model summary to file
  sink(paste(dirStats, "/", fileName,"_betaRegression_0+1-inflated.txt", sep=""))
  print(summary(m))
  sink()
  
  # Save diagnostic plot summary to file
  png(paste(dirStats, "/", fileName,"_betaRegression_0+1-inflated_diagnostic.png", sep=""), 
      width = 800, height = 800, units = "px")
  plot(m)
  dev.off()
  
  # Save predicted values to file as plot
  # Predicted responses
  means_m <- lpred(m, type='response', what='mu', se.fit=T)
  df_fit <- data.frame(HOURS = d$hours, M = means_m$fit, SE = means_m$se.fit)
  
  # Plot predicted responses
  plotPredict <- ggplot(df_fit, aes(HOURS, M, col = d$morphotype)) + 
    geom_pointrange(aes(ymin=M-SE, ymax=M+SE), size = 0.6) + 
    geom_line(lty=2, size = 0.8, alpha = 0.8)+
    labs(x="Time (hours)",y="Morph. composition (predicted)", col = "Morphology") + 
    scale_y_continuous(labels=scales::percent, limits = c(0, 1))+
    scale_color_manual(values = colsMT)+
    my_theme()+
    theme(legend.position = "right");plotPredict
  
  fileExtensions <- c(".svg", ".png")
  for(i in 1:length(fileExtensions)){
    ggsave(paste("predicted_", fileName,"_2022", fileExtensions[i], sep=""), plotPredict, dpi = 300, height = 5, width = 4) }
  
  return(plotPredict)
}

betaReg_propsOverTimePerTreatE(data, "bi")
betaReg_propsOverTimePerTreatE(data, "mono")

betaReg_propsOverTimePerMT <- function(data, treatE, MT){
  # Subset data
  d <- data[data$treatE == treatE & data$morphotype == MT, ]
  
  # Set up zero- and one-inflated beta-regression
  m <- gamlss(prop ~ hours, family = BEINF, data = d)
  
  # Extract p-value of hours
  sum.m <- summary(m, save = T)
  pval <- sum.m$coef.table[2,4]
  pval <- data.frame(treatE = treatE, 
                     morphotype = MT, 
                     pval = pval)
  
  # Define fileName for saving data
  fileName <- paste("propMT~time_betaRegression_", toupper(treatE),"_",MT,sep="")
  
  # Save model summary to file
  sink(paste(dirStats, "/", fileName,"_betaRegression_0+1-inflated.txt", sep=""))
  print(summary(m))
  sink()
  
  # Save diagnostic plot summary to file
  png(paste(dirStats, "/", fileName,"_betaRegression_0+1-inflated_diagnostic.png", sep=""))
  plot(m)
  dev.off()
  
  return(pval)
}

treats <- c(rep("bi", 3), rep("mono", 2))
morphotypes <- c("wrinkly", rep(c("smooth", "fuzzy"), 2))

pvalues <- data.frame()
for(t in 1:length(treats)){
  pvalues <- rbind(pvalues, betaReg_propsOverTimePerMT(data, treats[t], morphotypes[t]))}

# Correct for multiple testing
treats <- unique(treats)
pvalues$pval.FDR <- rep(100, length(pvalues))
for (i in 1:length(treats)){
  pvalues[pvalues$treatE == treats[i], 'pval.FDR'] <- p.adjust(pvalues[pvalues$treatE == treats[i], 'pval'], method = "fdr")
}

# Save to file
fileName <- paste("propMT~time_betaRegression_", toupper(treatE),sep="")
write.csv(pvalues, paste(dirStats, "/", fileName,"_betaRegression_0+1-inflated_FDRcorrectedPvalues.csv", sep=""))
