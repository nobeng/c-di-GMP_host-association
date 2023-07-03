################################
# Morphotype biofilms in vitro #
#                              #
#    Open, plot and analyze    #
#        biofilm data          #
#                              #
#     January - March 2019     #
################################

#### = 1. Set dependencies = ####
# Set working directory to source location
currentDirectory <- dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(currentDirectory)

# Libraries & source files
library(ggplot2);library(plyr);  library(ggbiplot); library(lattice); library(car); library(RColorBrewer); library(GGally)
library(vegan)
source("my_theme.R")

# Global parameters
nRepPops <- 6
colsTreat <- c("#02734A", "#F28907")

#### = 2. Import & concatenate data = ####
biofilms <- read.csv("morphotypes_biofilm.csv", header = T)
biofilms <- ddply(biofilms, .(morphotype, pop, treatE, morphology),
                  summarize, 
                  sdOD = sd(od), od = mean(od), 
                  sdOD_cor = sd(od_cor), od_cor = mean(od_cor))

dispersal <- read.csv("morphotypes_dispersal.csv", header = T)
dispersal <- subset(dispersal, select = -c(morphology, treatE, rep))
dispersal <- ddply(dispersal, .(pop, morphotype), 
                   summarize, 
                   sdSwarm24h = sd(swarm24h), swarm24h = mean(swarm24h), 
                   sdExpand3d = sd(expand3d), expand3d = mean(expand3d))

colonization <- read.csv("morphotypes_colonization.csv", header = T)
colonization <- colonization[complete.cases(colonization),]

colonization <- subset(colonization, select = -c(treatE, rep))
colonization <- ddply(colonization, .(pop, morphotype), 
                      summarize, 
                      sdCFUs_earlyCol = sd(CFUs_earlyCol), CFUs_earlyCol = mean(CFUs_earlyCol), 
                      sdCFUs_shortPersist = sd(CFUs_shortPersist), CFUs_shortPersist = mean(CFUs_shortPersist),
                      sdCFUs_release = sd(CFUs_release), CFUs_release = mean(CFUs_release))

# Work with log10-values for CFUs
colonization$CFUs_earlyCol <- log10(colonization$CFUs_earlyCol)
colonization$CFUs_shortPersist <- log10(colonization$CFUs_shortPersist)
colonization$CFUs_release <- log10(colonization$CFUs_release)

data <- merge(biofilms, dispersal, by = c("morphotype", "pop"))
data <- merge(data, colonization, by = c("morphotype", "pop"))

#### = 3. Calculate evo/anc quotient per phenotype = ####
phenotypes <- c("od_cor", "swarm24h", "expand3d", "CFUs_earlyCol", "CFUs_shortPersist", "CFUs_release")

calculateEvoAncQuotient <- function(data, phenotype){
  # Subset phenotype + treatment info from overall data table
  dataQuotient <- cbind(subset(data, select = c(morphotype, pop, treatE, morphology)),
                        phenotype = data[, phenotype])
  
  # Collect ancestral data points
  dataAnc <- dataQuotient[dataQuotient$treatE == "anc",]
  dataAnc <- subset(dataAnc, select = c(pop, phenotype))
  
  # Merge data frames and calculate the quotient of evolved/ancestral phenotype
  dataQuotient <- merge(dataQuotient, dataAnc, by = c("pop"))
  dataQuotient$quotient <- dataQuotient$phenotype.x/dataQuotient$phenotype.y
  
  # Data frame to keep only quotients
  dataQuotient <- subset(dataQuotient[dataQuotient$treatE != "anc", ], select = -c(phenotype.x, phenotype.y))
  
  # Name column with quotient by phenotype
  colnames(dataQuotient)[colnames(dataQuotient) == "quotient"] <- phenotype
  
  return(dataQuotient)  
}

dataQuotients <- data.frame()
# Collect evo/anc-quotients across phenotypes
for (p in 1:length(phenotypes)){
  ifelse(p == 1, 
         dataQuotients <- rbind(dataQuotients, calculateEvoAncQuotient(data, phenotypes[p])), 
         
         dataQuotients <- merge(dataQuotients, calculateEvoAncQuotient(data, phenotypes[p]),
                                by = c("pop", "morphotype","treatE", "morphology")))
}

# Rename data columns for better graph readibility
newNames <- c("Biofilm formation", "Swarming", "Colony expansion","Early colonization", "Short-term persistence", "Release")
colnames(data)[colnames(data)%in%phenotypes] <- newNames
colnames(dataQuotients)[colnames(dataQuotients)%in%phenotypes] <- newNames

# Clean up
rm(biofilms, dispersal, colonization)

treats <- c("mono", "bi")
treats <- c("bi")

# Use simplified annotation scheme of morphotypes (smooth, wrinkly, fuzzy)
legend <- read.csv("decode.treatID.csv", header = T, sep="\t" )
legend <- legend[legend$X == "",]
legend <- subset(legend, select = c("morphotype", "morph.simple", "pop"))

data <- merge(data, legend, by = c("morphotype", "pop"))

data$morph.simple <- as.character(data$morph.simple)
data[data$morph.simple == "anc", 'morph.simple'] <- "ancestor"
data$morph.simple <- as.factor(data$morph.simple)

#### = 4. Focus on wrinkly MT12, MT14, MT22 compared to MYb11_anc(MT48) = ####
data <- data[data$morphotype%in%c("MT48", "MT12", "MT14", "MT22"),]

data <- data[data$morph.simple%in%c("ancestor", "wrinkly"),]


# post-hoc for permanova (see: install_github("pmartinezarbizu/pairwiseAdonis/pairwiseAdonis"))
pairwise.adonis <- function(x,factors, sim.function = 'vegdist', sim.method = 'bray', p.adjust.m ='bonferroni',
                            reduce=NULL,perm=1000)
{  co <- combn(unique(as.character(factors)),2)
pairs <- c()
Df <- c()
SumsOfSqs <- c()
F.Model <- c()
R2 <- c()
p.value <- c()


for(elem in 1:ncol(co)){
  if(inherits(x, 'dist')){
    x1=as.matrix(x)[factors %in% c(as.character(co[1,elem]),as.character(co[2,elem])),
                    factors %in% c(as.character(co[1,elem]),as.character(co[2,elem]))]
  }
  
  else  (
    if (sim.function == 'daisy'){
      x1 = daisy(x[factors %in% c(co[1,elem],co[2,elem]),],metric=sim.method)
    } 
    else{x1 = vegdist(x[factors %in% c(co[1,elem],co[2,elem]),],method=sim.method)}
  )
  
  ad <- adonis(x1 ~ factors[factors %in% c(co[1,elem],co[2,elem])],
               permutations = perm);
  pairs <- c(pairs,paste(co[1,elem],'vs',co[2,elem]));
  Df <- c(Df,ad$aov.tab[1,1])
  SumsOfSqs <- c(SumsOfSqs, ad$aov.tab[1,2])
  F.Model <- c(F.Model,ad$aov.tab[1,4]);
  R2 <- c(R2,ad$aov.tab[1,5]);
  p.value <- c(p.value,ad$aov.tab[1,6])
}
p.adjusted <- p.adjust(p.value,method=p.adjust.m)

sig = c(rep('',length(p.adjusted)))
sig[p.adjusted <= 0.05] <-'.'
sig[p.adjusted <= 0.01] <-'*'
sig[p.adjusted <= 0.001] <-'**'
sig[p.adjusted <= 0.0001] <-'***'
pairw.res <- data.frame(pairs,Df,SumsOfSqs,F.Model,R2,p.value,p.adjusted,sig)

if(!is.null(reduce)){
  pairw.res <- subset (pairw.res, grepl(reduce,pairs))
  pairw.res$p.adjusted <- p.adjust(pairw.res$p.value,method=p.adjust.m)
  
  sig = c(rep('',length(pairw.res$p.adjusted)))
  sig[pairw.res$p.adjusted <= 0.1] <-'.'
  sig[pairw.res$p.adjusted <= 0.05] <-'*'
  sig[pairw.res$p.adjusted <= 0.01] <-'**'
  sig[pairw.res$p.adjusted <= 0.001] <-'***'
  pairw.res <- data.frame(pairw.res[,1:7],sig)
}
class(pairw.res) <- c("pwadonis", "data.frame")
return(pairw.res)
} 
summary.pwadonis = function(object, ...) {
  cat("Result of pairwise.adonis:\n")
  cat("\n")
  print(object, ...)
  cat("\n")
  cat("Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1\n")
}

colsMT <- c("black", "#006111ff", "#7f7f7fff", "#ee9b33ff")
colsMT <- c(colsMT[1], "#f83636ff",colsMT[4])

plotPCA <- function(data, lifeCycle){
  # Subset life cycle evolutionary treatment of interest
  ifelse(lifeCycle == "bi", 
         {d <- data[data$treatE != "mono", ]
         backgroundCol <- colsTreat[1]},
         {d <- data[data$treatE != "bi", ]
         backgroundCol <- colsTreat[2]})
  
  # Collect matrix of response variable data
  d <- subset(d,
              select = c("pop", "morph.simple","Biofilm formation", "Swarming", "Colony expansion", "Short-term persistence", "Early colonization", "Release"))
  responseVars <- subset(d, select = -c(morph.simple, pop))
  responseVars <- responseVars[complete.cases(responseVars),]
  d <- d[complete.cases(d),]
  d$pop <- as.factor(d$pop)
  
  # Set up PCA
  pcaOut <- prcomp(responseVars, center = T, scale. = T)
  
  # Test for differences between groups using the PERMANOVA
  m <- adonis(responseVars~ morph.simple, data = d); m
  sink(paste("PERMANOVA_morphotypes_", lifeCycle,"_2022-wrinkly.txt", sep=""))
  print("PERMANOVA")
  print(m)
  sink()
  
  pairOut <- pairwise.adonis(responseVars,factors=d$morph.simple)
  
  # Focus on ancestor vs. evolved comparison -> add adjusted p-values for those comparisons 
  pairOut <- pairOut[grepl("ancestor", pairOut$pairs), ]
  pairOut$p.adjusted <- p.adjust(pairOut$p.value, method = "fdr")
  
  # Save to file
  write.csv(pairOut, paste("pairwise_PERMANOVA_phenotypes_MT_", lifeCycle, "_2022-wrinkly.csv", sep=""))
  
  # Check out PCA summary & save to file
  sink(paste("morphotypes_PCA_", lifeCycle, "_2022-wrinkly.txt", sep=""))
  print(pcaOut)
  print(summary(pcaOut))
  sink()
  
  # Print PCA plot to file
  fontSize <- 25
  
  pca.plot <-  my_ggbiplot(pcaOut, ellipse = T, groups = d$morph.simple, varname.size = 8, varname.adjust = 0, alpha = 0)+
    geom_point(aes(color=d$morph.simple, fill = d$morph.simple),size = 4, alpha = 0.9) +
    scale_color_manual(values = colsMT)+
    scale_fill_manual(values = colsMT)+
    scale_shape_manual(values = c(7, 21:25))+
    labs(color = "Morphology", fill = "Morphology")+
    theme_bw()+
    theme(legend.direction ="vertical", 
          legend.position = "right", 
          legend.title = element_text(size = fontSize),
          legend.text = element_text(size = fontSize),
          axis.title = element_text(size = fontSize), 
          axis.text = element_text(size = fontSize)); pca.plot
  
  # Save to file
  fileExtension <- c(".png", ".svg")
  for(i in 1:length(fileExtension)){
    ggsave(paste("morphotypes_phenotypic PCAs_", lifeCycle, "_threeTypes_2022-wrinkly", fileExtension[i], sep = ""), 
           plot = pca.plot, height = 7, width = 9, dpi = 300) }
  
  return(pca.plot)
}

for (t in 1:length(treats)){ plotPCA(data, treats[t]) }

# Clean up
rm(data, newNames, t, treats, plotPCA, dataQuotients, colsTreat, currentDirectory, dirTheme, fontSize, nRepPops, p, phenotypes, 
   pointShapes, calculateEvoAncQuotient, my_ggbiplot, my_theme)