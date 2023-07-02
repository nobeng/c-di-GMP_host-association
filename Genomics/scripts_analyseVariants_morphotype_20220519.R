#######################################
#                                     #
# Script to analyse genomic variants  #
#    (samtools, VarScan2, snpEff)     #
#                                     #
#         May 2022, Nancy Obeng       #
#                                     #
#######################################

#### Set dependencies ####
# Set working directory to source location
currentDirectory <- dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(currentDirectory)

# Load packages
library(gtools); library(data.table); library(vcfR); library(stringr); library(plyr); library(dplyr)
library(ggplot2); library(svglite); library(RColorBrewer); library(tidyverse); library(VennDiagram); library(rlist)
library(gridExtra)
#source("C:/Users/Lenovo/Documents/work/phd/data/experimental evolution/G. pseudomonas diversification project/B) Student sub-projects/data/my_theme.R")
source("C:/Users/Nancy/Documents/phd/data/experimental evolution/G. pseudomonas diversification project/B) Student sub-projects/data/my_theme.R")

# Isolate-specific locations of genomic references
refDir <-  "C:/Users/Lenovo/Documents/work/phd/methods/genomics"
references <- c("myb11_ancestralREFSEQ", 
                "myb71REFSEQ_GCF_002975205.1_ASM297520v1_genomic", 
                "op50_GCF_009496595.1_ASM949659v1_genomic")

gffFiles <- paste(refDir,"/", references,".gff", sep="")
dnaFiles <- paste(refDir, "/", references, ".fasta", sep="")
chromNames <- c("myb11_anc")
strainOfAnalysis <- 1 # this is MYb11, see line 21 ff.

# #### Preparation to open data ####
# 
# # Filter ancestral variants from evolved
# ancestralSample <- "MT48_S16"
# allSamples <- list.files(pattern = ".vcf")
# evolvedSamples <- allSamples[!allSamples%like%ancestralSample] # exclude ancestral sample
# evolvedSamples <- unique(gsub("_var.*", "" , evolvedSamples)) # keep only sample identifier
# 
# # Leave out MTs identified to be OP50
# evolvedSamples <- evolvedSamples[!evolvedSamples%in%c("MT15_S31", "MT18_S6")]
# 
# # Note: Only needed initially, then muted and saved data loaded to save computational effort & time
# 
# #Collect variant effects in single row per effect
# collectVariantINFO <- function(my_sample, entry, fileName){
#   # Subset entry of focus
#   varX <- my_sample[entry,]
#   
#   # Collect variant info, including read depth info (DP)
#   variantInfo <- str_split(as.character(my_sample$INFO[entry]), pattern = ";ANN=")[[1]][1]
#   variantInfo <- strsplit(variantInfo, ";")[[1]]
#   
#   # Collect annotion INFO
#   annotation <- str_split(as.character(my_sample$INFO[entry]), pattern = ";ANN=")[[1]][2]
#   annotation <- strsplit(annotation, ",")[[1]]
#   
#   # Focus only on best call
#   annotation <- annotation[1]
#   
#   # Split annotation string based on pipe divisions, collect gene name containing entries and select the first (all should be the same)
#   # Attach to entry
#   varX$effectSO <- strsplit(as.character(annotation), "|", fixed = T)[[1]][2]
#   varX$gene_name <- substring(strsplit(as.character(annotation), "|", fixed = T)[[1]][4], 1, 13) # select substring to get single flanking gene of intergenic call
#   varX$geneID <-substring(strsplit(as.character(annotation), "|", fixed = T)[[1]][5], 1, 13)# select substring to get single flanking gene of intergenic call
#   varX$nucleotideChange <- strsplit(as.character(annotation), "|", fixed = T)[[1]][10]
#   varX$aaChange <- strsplit(as.character(annotation), "|", fixed = T)[[1]][11]
#   varX$impact <- strsplit(as.character(annotation), "|", fixed = T)[[1]][3]
#   varX$mutFreq <- strsplit(as.character(my_sample$inBAMfile[[1]]),":", fixed = T)[[1]][7]
#   
#   ifelse(fileName == "filtered", 
#          varX$coverage <- NA,
#          varX$coverage <- as.numeric(substring(variantInfo[grepl("ADP", variantInfo) == T], 5, 9)))
#   
#   return(varX) } # is called in function compileAnnotations()
# 
# compileAnnotations <- function(sampleName, fileName){
#   
#   # Check if variants were called by checking for >24 rows in file
#   vcfFile <- file(paste(sampleName, "_var_ann_", fileName,".vcf", sep="")) 
#   checkVCFContent <- length(readLines(vcfFile))
#   on.exit(close(vcfFile)) # Close file again
#   
#   variants <- data.frame()
#   if(checkVCFContent > 24) 
#          {
#            ### === 1. Filter variants with a coverage (read depth) >= 10 === #
#            # 1.1 Open in vcfR-package
#            # --> Combine vcf-file, annotation gff and fasta of the reference
#            vcf <- read.vcfR(paste(sampleName, "_var_ann_", fileName,".vcf", sep=""), verbose = F)
#            gff <- read.table(gffFiles[strainOfAnalysis], sep="\t", quote="")
#            dna <- ape::read.dna(dnaFiles[strainOfAnalysis], format = "fasta")
#            
#            # Create chrom object (can be plotted using plot)
#            chrom <- create.chromR(name = chromNames[strainOfAnalysis], vcf = vcf, seq=dna, ann = gff)
#            
#            # Filter data (min. quality and min.read depth)
#            chrom <- masker(chrom, min_QUAL = 1, min_DP = 10)
#            
#            # Extract positions and coverage of variants left
#            postions_to_keep <- variant.table(chrom)$POS
#            
#            # 1.2. Open annotated vcf and subset filtered variants
#            my_sample <-  read.table(paste(sampleName, "_var_ann_", fileName, ".vcf", sep=""))
#            colnames(my_sample) <- c("CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", "inBAMfile")
#            my_sample <- my_sample[my_sample$POS%in%postions_to_keep, ]
#            my_sample <- cbind(sample = rep(sampleName, length(my_sample$POS)), my_sample)
#            
#            ### === 2. Deparse annotation INFO === #
#            # 2.1. Collect INFO column and names of affected genes
#            variants <- data.frame()
#            for(entry in 1:length(my_sample$POS)){
#              variants <- rbind(variants, collectVariantINFO(my_sample, entry, fileName))  }
#            
#            # Check for NCBI gene name and product of the gene in the annotation database
#            variants$product <- rep("ph", length(variants$sample))# Attach product column
#            
#            # Annotate full table
#            for(entry in 1:length(variants$POS)){
#              # Check in annotation for gene product
#              gene_of_interest <- variants[entry, 'geneID']
#              annotation_of_interest <- as.character(chrom@ann[grep(gene_of_interest, chrom@ann$V9, fixed = T), 'V9'][[2]])
#              
#              # Enter into product column
#              variants[entry, 'product'] <- gsub(";.*", "", gsub(".*product=", "", annotation_of_interest)) }}
#   return(variants)
# }
# 
# # Function to get variants found by both samtools and VarScan2 (filtered = samtools, SNPs and Indels = VarScan2)
# get_wellSupportedVariants <- function(sampleName){
#   indels <- compileAnnotations(sampleName, "Indels") %>%
#     filter(POS %in% compileAnnotations(sampleName, "filtered")$POS)
#   indels$variantType <- rep("indel", length(indels$POS))
#   
#   SNPs <- compileAnnotations(sampleName, "SNPs") %>%
#     filter(POS %in% compileAnnotations(sampleName, "filtered")$POS)
#   SNPs$variantType <- rep("SNP", length(SNPs$POS))
#   
#   variants <- rbind(indels, SNPs)
#   
#   # add unique identifier for variants
#   variants$variantID <- with(variants, paste(POS, REF, ALT, sep="_"))
#   
#   return(variants) 
# }
# 
# # 1. Collect ancestral variants
# # ancestralVariants <- get_wellSupportedVariants(ancestralSample)
# # 
# # # 2. Collect evolved variants
# # evolvedVariants <- data.frame()
# # for(i in 1:length(evolvedSamples)){
# #   
# #   # Collect detected variants
# #   vars <- get_wellSupportedVariants(evolvedSamples[i])
# #   
# #   # Add column to note whether a variant appeared or disappeared during the EE (default: appeared)
# #   vars$change <- rep("appeared", length(vars$sample))
# #   
# #   # Add and higlight variants that were present in the ancestor but reverted to the reference
# #   varsUniqueToAncestor <- ancestralVariants %>%
# #     filter(!variantID %in% vars$variantID)
# #   varsUniqueToAncestor$sample <- rep(evolvedSamples[i], length(varsUniqueToAncestor$variantID))
# #   
# #   varsUniqueToAncestor$change <- rep("disappeared", length(varsUniqueToAncestor$sample))
# #   
# #   # Deduct variants present in the ancestor
# #   vars <- vars %>%
# #     filter(!variantID %in% ancestralVariants$variantID)
# #   
# #   vars <- rbind(vars, varsUniqueToAncestor)
# #   rm(varsUniqueToAncestor)
# #   
# #   # If variants are present, attach to overall list of evolved variants
# #   evolvedVariants <- rbind(evolvedVariants, vars)
# #   
# #   # Clean up
# #   rm(vars)
# #   
# # }  
# # 
# 
# # 5. Attach treatment information and organize table for better readibility
# # Morphotype legend
# dirUnblind <- "C:/Users/Lenovo/Documents/work/phd/data/experimental evolution/G. pseudomonas diversification project/B) Student sub-projects/data/dispersal"
# legend <- read.csv(paste(dirUnblind,"decode.treatID.csv", sep="/"), header = T)
# colnames(legend)[1]<- "sample"
# legend <- legend[legend$X != "2. Umlauf",]
# legend <- subset(legend, select = -c(sample, X))
# legend$morphotype <- as.character(legend$morphotype)
# legend <- subset(legend, select = c(morphotype, pop, treatE, morphology, morph.simple))
# legend <- legend[!legend$morphotype%in%c("OP50", "MYb11-dTomato", "MT15", "MT18"), ]
# 
# oldNames <- c("ancestral", "monophasic", "biphasic")
# newNames <- c("anc", "mono", "bi")
# 
# legend$treatE <- as.character(legend$treatE)
# for(i in 1:length(oldNames)){
#   legend[legend$treatE == oldNames[i], 'treatE'] <- newNames[i]
# }
# # 
# # evolvedVariants$morphotype <- vapply(strsplit(as.character(evolvedVariants$sample), "_", fixed = T), '[', 1, FUN.VALUE=character(1))
# # evolvedVariants <- merge(legend, evolvedVariants, by = "morphotype", all.x = F)
# # 
# # # Save to file
# # write.csv(evolvedVariants, file = "evolvedVariants_MT.csv", row.names = F)
# # write.csv(ancestralVariants, file = "ancestralVariants_MT.csv", row.names = F)
# # 
# # # Clean up
# # rm(ancestralVariants)
# 
# # 
#### Plotting ####
evolvedVariants <- read.csv("evolvedVariants_MT_withGO_morphology-sorted.csv", header = T)

# Add combined product and gene ID
evolvedVariants$productAndID <- with(evolvedVariants, paste(product, " (", gene_nameManualAdditions, ")", sep=""))

# Sort genes by morphotype, then functional (GO) annotation
evolvedVariants <- evolvedVariants %>% arrange(morph.simple, gene_nameManualAdditions)#, position_in_graph, geneID)
evolvedVariants$productAndID <- factor(evolvedVariants$productAndID, levels = rev(unique(evolvedVariants$productAndID)))

evolvedVariants$morph.simple <- factor(evolvedVariants$morph.simple, levels = c("wrinkly","smooth","fuzzy"))

evolvedVariants$gene_nameManualAdditions <- factor(evolvedVariants$gene_nameManualAdditions,
                                                   levels = rev(c("wspE", "wspF", "rph","CLM75_RS05480", 
                                                              "panB", "CLM75_RS14725", "CLM75_RS13935", 
                                                              "CLM75_RS24045", "CLM75_RS02445","CLM75_RS20620", "CLM75_RS16500", 
                                                              "CLM75_RS06470","metH",
                                                              "minC","minD", "rpoB")))
                  
evolvedVariants$gene <- factor(evolvedVariants$morph.simple, levels = c("wrinkly","smooth","fuzzy"))
evolvedVariants$pop <- as.factor(evolvedVariants$pop)

# Color schemes
colsTreat <- c(rgb(0, 158,115, maxColorValue = 255), rgb(230, 159,0, maxColorValue = 255))
colsTreat <- c(rgb(141,22,60, maxColorValue = 255), rgb(46, 64,87, maxColorValue = 255)) # 2022 update

# Main figure graph
# Plot by morphology across evo. treatment (corresponds to Venn Diagram)
colsMT <- c("#7f7f7fff", "#ee9b33ff", "#006111ff")

p <- ggplot(evolvedVariants, 
            aes(x = morph.simple, y = gene_nameManualAdditions))+
  geom_jitter(size = 4, height = 0, width=0.2)+
  facet_grid(~treatE)+
  labs(x = "Morphology", y = "Gene affected")+
  #scale_color_manual(values = colsMT)+
  #scale_fill_manual(values = colsMT)+
  my_theme()+
  theme(axis.text.x = element_text(hjust = 1, angle = 45));p

# Save to file
fileExtensions <- c(".svg", ".png")
for(i in 1:length(fileExtensions)){
  ggsave(paste("variantsMTs_byMT_2022_simple",fileExtensions[i], sep = ""),
         dpi = 300, height = 7, width = 7)}


# Supplementary reference data: counts & types of mutations
# Plot number of avg. mutations per isolate
numPerType <- legend %>% select(treatE, morph.simple) %>% table() %>% data.frame()
mutPerIsol <- evolvedVariants %>% select(treatE, morph.simple) %>% table() %>% data.frame()

mutPerIsol <- merge(numPerType, mutPerIsol, by = c("treatE", "morph.simple"))
mutPerIsol$avgMutPerIsolate <- mutPerIsol$Freq.y/ mutPerIsol$Freq.x

mutPerIsol$morph.simple <- factor(mutPerIsol$morph.simple, levels = c("smooth", "wrinkly", "fuzzy"))
ggplot(mutPerIsol, aes(x = toupper(treatE), y = avgMutPerIsolate, fill = treatE)) + 
  geom_bar(stat="identity",position="dodge")+
  labs(x = "Evo. life cycle", y = "Avg. number of mutations/ isolate")+
  scale_fill_manual(values = colsTreat)+
  facet_wrap(~morph.simple)+
  my_theme()

# Save to file
fileExtensions <- c(".svg", ".png")
for(i in 1:length(fileExtensions)){
  ggsave(paste("variantsMTs_SUPPL_avgMutPerIsolate",fileExtensions[i], sep = ""), dpi = 300, height = 4, width = 5)
}

# Plot mutation types by treatE
int_breaks <- function(x, n = 5) pretty(x, n)[pretty(x, n) %% 1 == 0] 
ggplot(evolvedVariants, aes(x = effectSO, fill = treatE)) + 
  geom_histogram(stat="count")+
  labs(x = "Type of mutational change", y = "Total number of changes")+
  scale_y_continuous(breaks = int_breaks)+
  facet_grid(morph.simple~toupper(treatE))+
  scale_fill_manual(values = colsTreat)+
  my_theme()+
  theme(axis.text.x = element_text(hjust = 1, angle = 50))


# Select colors based on mutational changes present to keep consistent legend
colsMut <- c("#173F5F", "#20639B", "#3CAEA3", "#F6D55C", "#ED553B", "grey80")
typesMut <- as.factor(c("missense_variant", "frameshift_variant","stop_gained","conservative_inframe_deletion",
                        "disruptive_inframe_deletion","upstream_gene_variant"))

biSamples <- evolvedVariants[evolvedVariants$treatE == "bi", ]
biSamples$effectSO <- factor(biSamples$effectSO, typesMut)
colsMutBi <- colsMut[which(typesMut%in%biSamples$effectSO == TRUE)]

monoSamples <- evolvedVariants[evolvedVariants$treatE == "mono", ]
monoSamples$effectSO <- factor(monoSamples$effectSO, typesMut)
colsMutMono <- colsMut[which(typesMut%in%monoSamples$effectSO == TRUE)]


bi <- ggplot(biSamples, aes(x = morph.simple, fill = effectSO)) + 
  geom_histogram(stat="count")+
  labs(x = "Morphology", fill = "Type of mutational change", y = "# of mutational changes")+
  scale_y_continuous(breaks = int_breaks)+
  facet_wrap(~toupper(treatE))+
  scale_fill_manual(values = colsMutBi)+
  my_theme()+
  theme(axis.text.x = element_text(hjust = 1, angle = 50), 
        legend.position = "right")

mono <- ggplot(monoSamples, aes(x = morph.simple, fill = effectSO)) + 
  geom_histogram(stat="count")+
  labs(x = "Morphology", fill = "Type of mutational change", y = "# of mutational changes")+
  scale_y_continuous(breaks = int_breaks)+
  facet_wrap(~toupper(treatE))+
  scale_fill_manual(values =colsMutMono)+
  my_theme()+
  theme(axis.text.x = element_text(hjust = 1, angle = 50), 
        legend.position = "right")

# Save to file
fileExtensions <- c(".svg", ".png")
for(i in 1:length(fileExtensions)){
  ggsave(paste("variantsMTs_SUPPL_typeMutChangePerTreatAndType",fileExtensions[i], sep = ""), 
         grid.arrange(bi, mono, nrow=2), dpi = 300, height = 6.5, width = 7)
}

# Plot variants by evo. life cycle
ggplot(evolvedVariants, aes(x = toupper(treatE), y = productAndID,
                            fill = treatE,  col = treatE))+
  geom_jitter(size = 4, alpha = 0.7, width = 0.1, height = 0.1)+
  labs(x = "Evo. life cycle", y = "Gene affected", col = "Evo. life cycle", shape ="Replicate\npopulation")+
  scale_color_manual(values = colsTreat)+
  scale_fill_manual(values = colsTreat)+
  scale_shape_manual(values = c(7, 21:25))+
  my_theme()+
  theme(legend.position = "none")

# Save to file
fileExtensions <- c(".svg", ".png")
for(i in 1:length(fileExtensions)){
  ggsave(paste("variantsMTs_byTreatE_2022",fileExtensions[i], sep = ""), dpi = 300, height = 6, width = 13)
}


  

# Plot by morphology within evo. treatment (corresponds to Venn Diagram)
plotByMorphology <- function(data, lifeC, plotH, plotW){
  ifelse(lifeC == "bi", 
         
         {colsMT <- c("grey50", "#36A5A5", "#FF4F76")
         selectShapes <- c(7, 21:25)[-2]},
         
         {colsMT <- c("grey50","#FF4F76")
         selectShapes <- c(7, 21:25)[1:3]})
  
  p <- ggplot(data[data$treatE == lifeC,], 
         aes(x = morph.simple, y = productAndID, fill = morph.simple, col = morph.simple, shape = pop))+
    geom_jitter(size = 4, alpha = 0.7, width = 0.1, height = 0.1)+
    labs(x = "Morphology", y = "Gene affected")+
    scale_color_manual(values = colsMT)+
    scale_fill_manual(values = colsMT)+
    scale_shape_manual(values = selectShapes)+
    my_theme()+
    theme(axis.text.x = element_text(hjust = 1, angle = 45))
  
  # Save to file
  fileExtensions <- c(".svg", ".png")
  for(i in 1:length(fileExtensions)){
    ggsave(paste("variantsMTs_byMT_",lifeC,fileExtensions[i], sep = ""), dpi = 300, height = plotH, width = plotW)
  }
  
  return(p)
}

plotByMorphology(evolvedVariants, "mono", 3, 9)
plotByMorphology(evolvedVariants, "bi", 5.5, 11)

#### Venn diagram of morphological specialities ####

# 1. Plot rough Venn diagrams by evo. life cycle treatment to count overlaps
plotRawVenn <- function(variantData, lifeC, types, fileExt){
  ifelse(lifeC =="biphasic" & types == 1, 
         {dataForVD <- list(smooth = variantData %>% filter(treatE == lifeC & morph.simple == "smooth") %>% unique() %>% select(geneID) %>% unlist(),
                            wrinkly = variantData %>% filter(treatE == lifeC & morph.simple == "wrinkly") %>% unique() %>% select(geneID) %>% unlist())
         colsMT <- c("grey50", "#36A5A5")
         fillMT <- c(alpha("grey50", 0.4), alpha("#36A5A5", 0.4))}, 
         
         ifelse(lifeC =="biphasic" & types == 2, 
                {dataForVD <- list(wrinkly = variantData %>% filter(treatE == lifeC & morph.simple == "wrinkly") %>% unique() %>% select(geneID) %>% unlist(), 
                                   fuzzy = variantData %>% filter(treatE == lifeC & morph.simple == "fuzzy") %>% unique() %>% select(geneID) %>% unlist())
                colsMT <- c("#36A5A5", "#FF4F76")
                fillMT <- c(alpha("#36A5A5", 0.4), alpha("#FF4F76", 0.4))}, 
         
         {dataForVD <- list(smooth = variantData %>% filter(treatE == lifeC & morph.simple == "smooth") %>% unique() %>% select(geneID) %>% unlist(),
                           fuzzy = variantData %>% filter(treatE == lifeC & morph.simple == "fuzzy") %>% unique() %>% select(geneID) %>% unlist())
         colsMT <- c("grey50", "#FF4F76")
         fillMT <- c(alpha("grey50", 0.4), alpha("#FF4F76", 0.4))}))

  venn.diagram(
    # Input data
    x = dataForVD,
    
    # Save to file
    filename = paste("rawVenn_", as.character(lifeC),"_",types, ".", fileExt, sep=""),
    imagetype = fileExt,
    height = 480 ,
    width = 480 ,
    resolution = 300,

    # Circles
    lwd = 2,
    fill = fillMT,
    col = colsMT,

    # Numbers
    cex = .5,
    fontfamily = "sans",

    # Set names
    cat.cex = 0.6,
    cat.default.pos = "outer",
    cat.fontfamily = "sans")
}

typeSubsets <- c(1, 2, 1)
treats <- c(rep("biphasic",2), "monophasic")
extensions <- c("svg", "png")
for(i in 1:length(treats)){
  for(j in 1:length(extensions)){
    plotRawVenn(evolvedVariants, treats[i], typeSubsets[i], extensions[j]) }}

# Plot overview of all three biphasic to get numbers
lifeC <- "biphasic"
variantsByMorph <- list(smooth = as.character(evolvedVariants %>% filter(treatE == lifeC & morph.simple == "smooth") %>% select(geneID) %>% unique() %>% unlist()),
                        wrinkly = as.character(evolvedVariants %>% filter(treatE == lifeC & morph.simple == "wrinkly") %>% select(geneID)%>% unique() %>% unlist()), 
                        fuzzy = as.character(evolvedVariants %>% filter(treatE == lifeC & morph.simple == "fuzzy") %>% select(geneID)%>% unique() %>% unlist()))

venn.diagram(x = variantsByMorph, filename = "rawVenn_biphasic_all.png")

# Collect geneIDs of Venn diagram
variantsByMorph <- list.append(variantsByMorph,
                          uniquelyS = variantsByMorph$smooth[!variantsByMorph$smooth%in%c(variantsByMorph$wrinkly, variantsByMorph$fuzzy)], 
                          uniquelyW = variantsByMorph$wrinkly[!variantsByMorph$wrinkly%in%c(variantsByMorph$smooth, variantsByMorph$fuzzy)], 
                          uniquelyF = variantsByMorph$fuzzy[!variantsByMorph$fuzzy%in%c(variantsByMorph$smooth, variantsByMorph$wrinkly)],
                          all = intersect(intersect(variantsByMorph$smooth, variantsByMorph$wrinkly), variantsByMorph$fuzzy))

# Step 1 in getting intersections
variantsByMorph <- list.append(variantsByMorph, 
                               SW = intersect(variantsByMorph$smooth, variantsByMorph$wrinkly),
                               SF = intersect(variantsByMorph$smooth, variantsByMorph$fuzzy), 
                               WF = intersect(variantsByMorph$fuzzy, variantsByMorph$wrinkly)) 
variantsByMorph$SW <- variantsByMorph$SW[!variantsByMorph$SW%in%c(variantsByMorph$all)]
variantsByMorph$SF <- variantsByMorph$SF[!variantsByMorph$SF%in%c(variantsByMorph$all)]
variantsByMorph$WF <- variantsByMorph$WF[!variantsByMorph$WF%in%c(variantsByMorph$all)]


lifeC <- "monophasic"
variantsByMorphMONO <- list(smooth = as.character(evolvedVariants %>% filter(treatE == lifeC & morph.simple == "smooth") %>% select(geneID) %>% unique() %>% unlist()),
                            fuzzy = as.character(evolvedVariants %>% filter(treatE == lifeC & morph.simple == "fuzzy") %>% select(geneID)%>% unique() %>% unlist()))

variantsByMorphMONO <- list.append(variantsByMorphMONO,
                               uniquelyS = variantsByMorphMONO$smooth[!variantsByMorphMONO$smooth%in%variantsByMorphMONO$fuzzy], 
                               uniquelyF = variantsByMorphMONO$fuzzy[!variantsByMorphMONO$fuzzy%in%variantsByMorphMONO$smooth],
                               all = intersect(variantsByMorphMONO$smooth, variantsByMorphMONO$fuzzy))

# Clean up
rm(variantsByMorph, variantsByMorphMONO)
