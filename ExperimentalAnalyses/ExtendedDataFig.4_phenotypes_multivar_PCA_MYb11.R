########################################
# Analysis of colonization experiments #
#    6. ANAYLSIS: MULTIVARIATE         #
#        EVOLVED VS. ANCESTRAL         #
########################################

#### = Set dependencies + open data ####
# Set working directory to source location
currentDirectory <- dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(currentDirectory)

# Specify libraries
library(ggplot2);library(plyr);  library(ggbiplot); library(RColorBrewer); library(GGally)
library(missMDA)# see: https://francoishusson.wordpress.com/2017/08/05/multiple-imputation-for-continuous-and-categorical-data-2/
library(vegan)

# Source files
source("my_theme.R")

#### = Open absolute data and concatenate from different assays = ####
importData <- function(assayName){
  data <- subset(read.csv(paste(assayName, "_absData.csv", sep = "")))
  
  # Streamline data frame
  ifelse("assay"%in%colnames(data) & !data$assay%in%c("motility"), 
         {
           # Average technical replicates
           data <- ddply(data, .(assay, treatE, cycle, strain, repO), summarize, 
                         cfu_w = mean(cfu_w), cfu_releas = mean(cfu_releas))
     
           # Note assay name for colonization assays in cfu read out & kick out redundant columns
           ifelse(data$assay == "releas",
                  {colnames(data)[colnames(data) == "cfu_w"] <- paste(unique(data$assay), "cfu_w", sep="_")
                  colnames(data)[colnames(data) == "cfu_releas"] <- paste(unique(data$assay), "cfu_releas", sep="_")
                  data <- subset(data, select = -c(assay))}, 
                  
                  {colnames(data)[colnames(data) == "cfu_w"] <- paste(unique(data$assay), "cfu_w", sep="_")
                  data <- subset(data, select = -c(assay, cfu_releas))})
          }, 
         data <- subset(data, select = -assay))
  
  return(data)
}

assayNames <- c("col.2h", "releas", "motility")

allAbsData <- data.frame()
for (i in 1:length(assayNames)){
  newData <- importData(assayNames[i])
  newData <- newData[newData$strain == "MYb11" & newData$cycle != 4, ]
  
  ifelse(i == 1, 
         allAbsData <- newData, 
         allAbsData <- merge(allAbsData, importData(assayNames[i]), by = c("treatE", "cycle", "strain", "repO"), all.x=T))
}

allAbsData <- subset(allAbsData, select = -c(motil_3.4agar_7d, motil_3.4agar_24h))

colnames(allAbsData)[colnames(allAbsData)== "repO"] <- "pop"

# Note, for motility, all populations were assayed at once with one ancestor (with techn. reps.), 
# thus add ancestral values to other ancestral phenotypes
allAbsData[allAbsData$treatE == "anc", 'motil_0.5agar_24h'] <- rep(allAbsData[allAbsData$treatE == "anc" & allAbsData$pop == 1,
                                                                          'motil_0.5agar_24h'], length(unique(allAbsData$pop)))

allAbsData[allAbsData$treatE == "anc", 'motil_3.4agar_3d'] <- rep(allAbsData[allAbsData$treatE == "anc" & allAbsData$pop == 1,
                                                                              'motil_3.4agar_3d'], length(unique(allAbsData$pop)))
allAbsData$col.2h_cfu_w <- log10(allAbsData$col.2h_cfu_w)
allAbsData$releas_cfu_w <- log10(allAbsData$releas_cfu_w)
allAbsData$releas_cfu_releas <- log10(allAbsData$releas_cfu_releas)

# Rename phenotypes for better graph readibility
oldNames <- c("col.2h_cfu_w","releas_cfu_w", "releas_cfu_releas", "motil_0.5agar_24h", "motil_3.4agar_3d")
 
newNames <- c("Early colonization", "Short-term persistence", "Release", "Swarming", "Colony expansion")

for(i in 1:length(oldNames)){
  colnames(allAbsData)[colnames(allAbsData) == oldNames[i]] <-newNames[i]}

# Clean up
rm(assayNames, i, newNames, oldNames, importData, newData)

#### = Compute PCAs = ####

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

compute_and_plot_PCA <- function(d, dataSubset){
  # Subset data: focus MYb11 at cycle 10
  d <- d[d$strain == "MYb11" & d$cycle != 4, ]
  
  # Impute ancestral NAs
  impute <- imputePCA(subset(d, select=-c(treatE, cycle, strain, pop)))
  impute <- impute$completeObs
  d <- cbind(subset(d, select=c(treatE, cycle, strain, pop)),impute)

  # Collect matrix of response variable data
  d <- subset(d, select =-c(cycle, strain))
  responseVars <- subset(d, select = -c(treatE, pop))
  d$pop <- as.factor(d$pop)
  
  # Set up PCA
  pcaOut <- prcomp(responseVars, center = T, scale = T)

  # Test for differences between groups using the PERMANOVA
  m <- adonis(responseVars~ treatE, data = d); m
  sink(paste("PCAs/PERMANOVA_popsMYb11_", dataSubset,"_simple.txt", sep=""))
  print("PERMANOVA")
  print(m)
  sink()

  pairOut <- pairwise.adonis(responseVars,factors=d$treatE)

  # Save to file
  dirStats <- paste(currentDirectory, "/stats out", sep="")
  write.csv(pairOut, paste(dirStats, "/pairwise_PERMANOVA_phenotypes_popsMYb11_", dataSubset, "_simple.csv", sep=""))

  # Check out PCA summary & save to file
  sink(paste(dirStats, "/StagesLifeCycle_PCA_MYb11_", dataSubset,"_simple.txt", sep=""))
  print("PCA")
  print(pcaOut)
  print("IMPORTANCE OF COMPONENTS")
  print(summary(pcaOut))
  sink()

  # Plot PCA
  # Function based on ggbiplot (edited arrow thickeness, arrow and label color, label position and axis label)
  my_ggbiplot <- function (pcobj, choices = 1:2, scale = 1, pc.biplot = TRUE, 
                           obs.scale = 1 - scale, var.scale = scale, groups = NULL, 
                           ellipse = FALSE, ellipse.prob = 0.68, labels = NULL, labels.size = 3, 
                           alpha = 1, var.axes = TRUE, circle = FALSE, circle.prob = 0.69, 
                           varname.size = 3, varname.adjust = 1.5, varname.abbrev = FALSE, 
                           ...) {
    library(ggplot2)
    library(plyr)
    library(scales)
    library(grid)
    stopifnot(length(choices) == 2)
    if (inherits(pcobj, "prcomp")) {
      nobs.factor <- sqrt(nrow(pcobj$x) - 1)
      d <- pcobj$sdev
      u <- sweep(pcobj$x, 2, 1/(d * nobs.factor), FUN = "*")
      v <- pcobj$rotation
    }
    else if (inherits(pcobj, "princomp")) {
      nobs.factor <- sqrt(pcobj$n.obs)
      d <- pcobj$sdev
      u <- sweep(pcobj$scores, 2, 1/(d * nobs.factor), FUN = "*")
      v <- pcobj$loadings
    }
    else if (inherits(pcobj, "PCA")) {
      nobs.factor <- sqrt(nrow(pcobj$call$X))
      d <- unlist(sqrt(pcobj$eig)[1])
      u <- sweep(pcobj$ind$coord, 2, 1/(d * nobs.factor), FUN = "*")
      v <- sweep(pcobj$var$coord, 2, sqrt(pcobj$eig[1:ncol(pcobj$var$coord), 
                                                    1]), FUN = "/")
    }
    else if (inherits(pcobj, "lda")) {
      nobs.factor <- sqrt(pcobj$N)
      d <- pcobj$svd
      u <- predict(pcobj)$x/nobs.factor
      v <- pcobj$scaling
      d.total <- sum(d^2)
    }
    else {
      stop("Expected a object of class prcomp, princomp, PCA, or lda")
    }
    choices <- pmin(choices, ncol(u))
    df.u <- as.data.frame(sweep(u[, choices], 2, d[choices]^obs.scale, 
                                FUN = "*"))
    v <- sweep(v, 2, d^var.scale, FUN = "*")
    df.v <- as.data.frame(v[, choices])
    names(df.u) <- c("xvar", "yvar")
    names(df.v) <- names(df.u)
    if (pc.biplot) {
      df.u <- df.u * nobs.factor
    }
    r <- sqrt(qchisq(circle.prob, df = 2)) * prod(colMeans(df.u^2))^(1/4)
    v.scale <- rowSums(v^2)
    df.v <- r * df.v/sqrt(max(v.scale))
    if (obs.scale == 0) {
      u.axis.labs <- paste("PC", choices, 
                           sep = "")
    }
    else {
      u.axis.labs <- paste("PC", choices, sep = "")
    }
    u.axis.labs <- paste(u.axis.labs, sprintf("(%0.1f%%)", 
                                              100 * pcobj$sdev[choices]^2/sum(pcobj$sdev^2)))
    if (!is.null(labels)) {
      df.u$labels <- labels
    }
    if (!is.null(groups)) {
      df.u$groups <- groups
    }
    if (varname.abbrev) {
      df.v$varname <- abbreviate(rownames(v))
    }
    else {
      df.v$varname <- rownames(v)
    }
    df.v$angle <- with(df.v, (180/pi) * atan(yvar/xvar))
    df.v$hjust = with(df.v, (1 - varname.adjust * sign(xvar))/2)
    g <- ggplot(data = df.u, aes(x = xvar, y = yvar)) + xlab(u.axis.labs[1]) + 
      ylab(u.axis.labs[2]) + coord_equal()
    g <- g + geom_hline(yintercept = 0, lty= 2) +
      geom_vline(xintercept = 0, lty= 2)
    
    if (var.axes) {
      if (circle) {
        theta <- c(seq(-pi, pi, length = 50), seq(pi, -pi, 
                                                  length = 50))
        circle <- data.frame(xvar = r * cos(theta), yvar = r * 
                               sin(theta))
        g <- g + geom_path(data = circle, color = muted("white"), 
                           size = 1/2, alpha = 1/3)
      }
      g <- g + geom_segment(data = df.v, aes(x = 0, y = 0, 
                                             xend = xvar, yend = yvar), arrow = arrow(length = unit(1/2, 
                                                                                                    "picas")),
                            color = "black", size = 1.5)
    }
    if (!is.null(df.u$labels)) {
      if (!is.null(df.u$groups)) {
        g <- g + geom_text(aes(label = labels, color = groups), 
                           size = labels.size)
      }
      else {
        g <- g + geom_text(aes(label = labels), size = labels.size)
      }
    }
    else {
      if (!is.null(df.u$groups)) {
        g <- g + geom_point(aes(color = groups), alpha = alpha)
      }
      else {
        g <- g + geom_point(alpha = alpha)
      }
    }
    if (!is.null(df.u$groups) && ellipse) {
      theta <- c(seq(-pi, pi, length = 50), seq(pi, -pi, length = 50))
      circle <- cbind(cos(theta), sin(theta))
      ell <- ddply(df.u, "groups", function(x) {
        if (nrow(x) <= 2) {
          return(NULL)
        }
        sigma <- var(cbind(x$xvar, x$yvar))
        mu <- c(mean(x$xvar), mean(x$yvar))
        ed <- sqrt(qchisq(ellipse.prob, df = 2))
        data.frame(sweep(circle %*% chol(sigma) * ed, 2, 
                         mu, FUN = "+"), groups = x$groups[1])
      })
      names(ell)[1:2] <- c("xvar", "yvar")
      g <- g + geom_path(data = ell, aes(color = groups, group = groups), size = 1)
    }
    if (var.axes) {
      g <- g + geom_text(data = df.v, aes(label = varname, 
                                          x = xvar, y = yvar, angle = 0, hjust = hjust), 
                         color = "black", size = varname.size)
    }
    return(g)
  }
  
  ifelse(dataSubset == "abs", 
         colsTreat <- c("black", "#8d163cff","#666666ff"), 
         colsTreat <- c(rgb(0, 158,115, maxColorValue = 255), rgb(230, 159,0, maxColorValue = 255)))
  
  fontSize <- 25
  
  pca.plot <- my_ggbiplot(pcaOut, ellipse = T, groups = toupper(d$treatE), varname.size = 5, varname.adjust = 0, alpha = 0)+
    geom_jitter(aes(color=toupper(d$treatE), fill = toupper(d$treatE)), size = 4, alpha = 0.9, width = 0.05)+ #add tiny jitter to allow seeing overlapping ancestral samples
    scale_color_manual(values = colsTreat)+
    scale_fill_manual(values = colsTreat)+
    scale_shape_manual(values = c(7, 21:25))+
    labs(color = "Evo. life cycle", fill = "Evo. life cycle", shape = "Replicate\npopulation")+
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
    ggsave(paste("treatE_phenotypic PCAs_MYb11_simple_",dataSubset, fileExtension[i], sep = ""), 
           plot = pca.plot, height = 8, width = 10, dpi = 300) }
  
 return("PCA saved to file.") 
}

set.seed(20200714)
compute_and_plot_PCA(allAbsData, "abs")

# Clean up
rm(allAbsData, currentDirectory, compute_and_plot_PCA, my_theme)