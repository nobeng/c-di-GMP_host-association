############################
# Setting up R environment #
############################

# Set working directory (with folder 'target sequences' containing gff files & 'complete sequences.seq' + 'classification.csv' in parental directory)
wd <- "\PATH_TO_THE_WORKING_DIRECTORY"
setwd(wd)

# Load required packages
library(rentrez)
library(rvest)
library(annotate)
library(seqinr)
library(ape)
library(rBLAST)


#############################################################
# Download genomes from NCBI for identified gene candidates #
#############################################################

setwd(paste(wd, "target sequences", sep = ""))

gene.names            <- c("rph", "wspA", "wspB", "wspC", "wspD", "wspE", "wspF", "wspR")
reference.seqs        <- vector(mode = "list", length = length(gene.names))
names(reference.seqs) <- gene.names

for (i in 1:length(gene.names)) {
  # Read file with reference genome IDs & features
  features                  <- read.gff(list.files()[i])
  features                  <- features[which(features$type == "CDS"),]
  
  # Initialize data frame containing reference gene info
  reference.seqs[[i]]           <- as.data.frame(matrix(nrow = nrow(features), ncol = 6))
  colnames(reference.seqs[[i]]) <- c("NCBI ID", "gene name", "gene start", "gene end", "strand", "gene sequence")
  reference.seqs[[i]][, 1]      <- features$seqid
  reference.seqs[[i]][, 2]      <- rep(gene.names[i], nrow(features))
  reference.seqs[[i]][, 3]      <- features$start
  reference.seqs[[i]][, 4]      <- features$end
  reference.seqs[[i]][, 5]      <- features$strand
  
  # Download genome containing reference gene sequence from NCBI & extract gene sequence (+ strand coding)
  for (j in 1:nrow(features)) {
    sequence                    <- entrez_fetch(db = "nucleotide", id = features$seqid[j], rettype = "fasta")
    sequence                    <- unname(unlist(strsplit(unlist(read.fasta(textConnection(sequence), as.string = T)), split = "")))
    sequence                    <- sequence[features$start[j]:features$end[j]]
    if (features$strand[j] == "-") {
      sequence                  <- rev(comp(sequence))
    }
    reference.seqs[[i]][j, 6]   <- paste(toupper(sequence), collapse = "")
    cat(c(i, j, "\n"))
  }
}

saveRDS(reference.seqs, "reference seqs.RDS")


##################################
# Complete genome identification #
##################################

setwd(wd)
sequences <- read.table("complete sequences.seq")

no.seqs <- NULL

for (i in 1:nrow(sequences)) {
  fasta <- entrez_fetch(db = "nucleotide", id = sequences[i, 1], rettype = "fasta")
  if (fasta == "\n") {
    no.seqs <- c(no.seqs, i)
  }
  if (i %% 10 == 0) {
    print(i)
  }
}

sequences <- sequences[-no.seqs,]
write.csv(sequences, "validated whole genome sequences.csv", row.names=F)


##############################################
# Collect meta info on 2279 complete genomes #
##############################################

# Data on organism, strain, length, biosample ID & NCBI reference ID
genome.list         <- read.csv("validated whole genome sequences.csv")[,1]
sequences           <- as.data.frame(matrix(nrow = length(genome.list), ncol = 5))
colnames(sequences) <- c("organism", "strain", "length", "biosample nr", "NCBI reference nr")

for (i in 1:length(genome.list)) {
  summary         <- entrez_summary(db = "nucleotide", id = genome.list[i])
  sequences[i, 1] <- ifelse(is.null(summary$organism), NA, summary$organism)
  sequences[i, 2] <- ifelse(is.null(summary$strain), NA, summary$strain)
  sequences[i, 3] <- ifelse(is.null(summary$slen), NA, summary$slen)
  sequences[i, 4] <- ifelse(is.null(summary$biosample), NA, summary$biosample)
  sequences[i, 5] <- genome.list[i]
  if (i %% 10 == 0) {
    cat(paste(i, "/", length(genome.list), "\n"))
  }
}

write.csv(sequences, "validated whole genome sequences.csv", row.names = F)


# Crawl info on host, host disease, isolation source & sample type from NCBI's biosample database
setwd(wd)
sequences <- read.csv("validated whole genome sequences.csv")

host                                                           <- NULL
host.disease                                                   <- NULL
isolation.source                                               <- NULL
sample.type                                                    <- NULL

for (i in 1:nrow(sequences)) {
  url                                                          <- paste("https://www.ncbi.nlm.nih.gov/biosample/?term=", sequences$biosample.nr[i], sep="")
  page                                                         <- read_html(url)
  listings                                                     <- html_nodes(page, css = ".rprt")
  table                                                        <- as.data.frame(html_table(listings))
  if (nrow(table) > 0) {
    table                                                       <- table[which(table[,1] %in% c("host", "host disease", "isolation source", "sample type")),]
    table[table[,2] %in% c("unknown", "missing", "Unknown"), 2] <- NA
    
    host[i]                                                     <- ifelse(length(table[table[,1] == "host", 2]) == 0, NA, table[table[,1] == "host", 2])
    host.disease[i]                                             <- ifelse(length(table[table[,1] == "host disease", 2]) == 0, NA, table[table[,1] == "host disease", 2])
    isolation.source[i]                                         <- ifelse(length(table[table[,1] == "isolation source", 2]) == 0, NA, table[table[,1] == "isolation source", 2])
    sample.type[i]                                              <- ifelse(length(table[table[,1] == "sample type", 2]) == 0, NA, table[table[,1] == "sample type", 2])
  }
}

sequences$host                                                   <- host
sequences$host.disease                                           <- host.disease
sequences$isolation.source                                       <- isolation.source
sequences$sample.type                                            <- sample.type
write.csv(sequences, "validated whole genome sequences.csv", row.names=F)


##########################################################
# Perform BLAST with reference genes as target sequences #
##########################################################

setwd(wd)
sequences            <- read.csv("validated whole genome sequences.csv")

setwd(paste(wd, "target sequences", sep = ""))
reference.seqs       <- readRDS("reference seqs.RDS")

setwd(paste(wd, "blast", sep = ""))

blast.res            <- NULL

for (i in 1:nrow(sequences)) {
  seq                <- entrez_fetch(db = "nucleotide", id = sequences$NCBI.reference.nr[i], rettype = "fasta")
  write(seq, "blast.fasta")
  makeblastdb("blast.fasta", dbtype = "nucl")
  db                 <- blast("blast.fasta")
  for (j in 1:length(reference.seqs)) {
    blast            <- predict(db, DNAStringSet(reference.seqs[[j]]$`gene sequence`))
    if (nrow(blast) > 0) {
      bases          <- unname(unlist(strsplit(unlist(read.fasta(textConnection(seq), as.string = T)), split = "")))
      gene.seq       <- NULL
      for (k in 1:nrow(blast)) {
        geneseq      <- bases[blast$S.start[k]:blast$S.end[k]]
        if (blast$S.start[k] > blast$S.end[k]) {
          geneseq    <- comp(geneseq)
        }
        gene.seq[k]  <- paste(toupper(geneseq), collapse="")
      }
      blast$Sequence <- gene.seq
      blast$QueryID  <- reference.seqs[[j]]$`NCBI ID`[as.integer(sapply(strsplit(blast$QueryID, split = "Query_"), "[[", 2))]
      blast$Gene     <- rep(names(reference.seqs)[j], nrow(blast))
      blast.res      <- rbind.data.frame(blast.res, blast)
    }
  }
  unlink(getwd(), recursive=T)
  cat(paste("Run ", i, ": ", nrow(blast.res), " matches in total\n", sep = ""))
  if (i %% 10 == 0) {
    write.csv(blast.res, paste(wd, "blast_res.csv", sep = ""), row.names = F)
  }
}
write.csv(blast.res, paste(wd, "blast_res.csv", sep = ""), row.names = F)


#############################################################
# Exclude genome sequence duplicates (e.g. NC_xxx & NZ_xxx) #
#############################################################

setwd(wd)

sequences     <- read.csv("validated whole genome sequences.csv")
blasted       <- read.csv("blast_res.csv")
behavior      <- read.csv("classification.csv")

gene.names    <- c("rph", "wspA", "wspB", "wspC", "wspD", "wspE", "wspF", "wspR")

# identify NZ_ or NC_ duplicates
duplicate <- NULL
for (i in 1:nrow(sequences)) {
  id <- unlist(strsplit(sequences$NCBI.reference.nr[i], split = "_"))
  if (length(id) == 2) {
    if (length(which(sequences$NCBI.reference.nr == id[2])) > 0) {
      duplicate <- c(duplicate, sequences$NCBI.reference.nr[i])
    }
  }
}
blasted <- blasted[-which(blasted$SubjectID %in% duplicate),]
sequences <- sequences[-which(sequences$NCBI.reference.nr %in% duplicate),]


##########################
# Common filter criteria #
##########################

# Create plots to identify optimal filter criteria
library(lattice)
library(gridExtra)

col = colorRampPalette(c("white", "mediumblue", "midnightblue"))(100)

i = "rph"
subset     <- blasted[blasted$Gene == i,]
seqlengths <- abs(subset$S.end - subset$S.start) + 1
identity   <- subset$Perc.Ident
data       <- table(data.frame(x = cut(identity, breaks = seq(70, 100, length = 11)), y = cut(seqlengths, breaks = seq(min(seqlengths), max(seqlengths), length = 11))))
data       <- data/sum(data)
plot1      <- levelplot(data, xlab = "Percent Identity", ylab = "Sequenz-Länge", main = i, col.regions = col, scales = list(x = list(rot = 90)), panel = function(...) {
  panel.levelplot(...)
  panel.lines(x=c(3.5,10.5), y=c(9.5,9.5), lwd = 2, col = "red")
  panel.lines(x=c(3.5,10.5), y=c(10.5,10.5), lwd = 2, col = "red")
  panel.lines(x=c(3.5,3.5), y=c(9.5,10.5), lwd = 2, col = "red")
  panel.lines(x=c(10.5,10.5), y=c(9.5,10.5), lwd = 2, col = "red")
})

i = "wspA"
subset     <- blasted[blasted$Gene == i,]
seqlengths <- abs(subset$S.end - subset$S.start) + 1
identity   <- subset$Perc.Ident
data       <- table(data.frame(x = cut(identity, breaks = seq(70, 100, length = 11)), y = cut(seqlengths, breaks = seq(min(seqlengths), max(seqlengths), length = 11))))
data       <- data/sum(data)
plot2      <- levelplot(data, xlab = "Percent Identity", ylab = "Sequenz-Länge", main = i, col.regions = col, scales = list(x = list(rot = 90)), panel = function(...) {
  panel.levelplot(...)
  panel.lines(x=c(4.5,10.5), y=c(5.5,5.5), lwd = 2, col = "red")
  panel.lines(x=c(4.5,10.5), y=c(10.5,10.5), lwd = 2, col = "red")
  panel.lines(x=c(4.5,4.5), y=c(5.5,10.5), lwd = 2, col = "red")
  panel.lines(x=c(10.5,10.5), y=c(5.5,10.5), lwd = 2, col = "red")
})

i = "wspB"
subset     <- blasted[blasted$Gene == i,]
seqlengths <- abs(subset$S.end - subset$S.start) + 1
identity   <- subset$Perc.Ident
data       <- table(data.frame(x = cut(identity, breaks = seq(70, 100, length = 11)), y = cut(seqlengths, breaks = seq(min(seqlengths), max(seqlengths), length = 11))))
data       <- data/sum(data)
plot3      <- levelplot(data, xlab = "Percent Identity", ylab = "Sequenz-Länge", main = i, col.regions = col, scales = list(x = list(rot = 90)), panel = function(...) {
  panel.levelplot(...)
  panel.lines(x=c(2.5,10.5), y=c(9.5,9.5), lwd = 2, col = "red")
  panel.lines(x=c(2.5,10.5), y=c(10.5,10.5), lwd = 2, col = "red")
  panel.lines(x=c(2.5,2.5), y=c(9.5,10.5), lwd = 2, col = "red")
  panel.lines(x=c(10.5,10.5), y=c(9.5,10.5), lwd = 2, col = "red")
})

i = "wspC"
subset     <- blasted[blasted$Gene == i,]
seqlengths <- abs(subset$S.end - subset$S.start) + 1
identity   <- subset$Perc.Ident
data       <- table(data.frame(x = cut(identity, breaks = seq(70, 100, length = 11)), y = cut(seqlengths, breaks = seq(min(seqlengths), max(seqlengths), length = 11))))
data       <- data/sum(data)
plot4      <- levelplot(data, xlab = "Percent Identity", ylab = "Sequenz-Länge", main = i, col.regions = col, scales = list(x = list(rot = 90)), panel = function(...) {
  panel.levelplot(...)
  panel.lines(x=c(1.5,10.5), y=c(5.5,5.5), lwd = 2, col = "red")
  panel.lines(x=c(1.5,10.5), y=c(10.5,10.5), lwd = 2, col = "red")
  panel.lines(x=c(1.5,1.5), y=c(5.5,10.5), lwd = 2, col = "red")
  panel.lines(x=c(10.5,10.5), y=c(5.5,10.5), lwd = 2, col = "red")
})

i = "wspD"
subset     <- blasted[blasted$Gene == i,]
seqlengths <- abs(subset$S.end - subset$S.start) + 1
identity   <- subset$Perc.Ident
data       <- table(data.frame(x = cut(identity, breaks = seq(70, 100, length = 11)), y = cut(seqlengths, breaks = seq(min(seqlengths), max(seqlengths), length = 11))))
data       <- data/sum(data)
plot5      <- levelplot(data, xlab = "Percent Identity", ylab = "Sequenz-Länge", main = i, col.regions = col, scales = list(x = list(rot = 90)), panel = function(...) {
  panel.levelplot(...)
  panel.lines(x=c(0.5,10.5), y=c(9.5,9.5), lwd = 2, col = "red")
  panel.lines(x=c(0.5,10.5), y=c(10.5,10.5), lwd = 2, col = "red")
  panel.lines(x=c(0.5,0.5), y=c(9.5,10.5), lwd = 2, col = "red")
  panel.lines(x=c(10.5,10.5), y=c(9.5,10.5), lwd = 2, col = "red")
})

i = "wspE"
subset     <- blasted[blasted$Gene == i,]
seqlengths <- abs(subset$S.end - subset$S.start) + 1
identity   <- subset$Perc.Ident
data       <- table(data.frame(x = cut(identity, breaks = seq(70, 100, length = 11)), y = cut(seqlengths, breaks = seq(min(seqlengths), max(seqlengths), length = 11))))
data       <- data/sum(data)
plot6      <- levelplot(data, xlab = "Percent Identity", ylab = "Sequenz-Länge", main = i, col.regions = col, scales = list(x = list(rot = 90)), panel = function(...) {
  panel.levelplot(...)
  panel.lines(x=c(2.5,10.5), y=c(6.5,6.5), lwd = 2, col = "red")
  panel.lines(x=c(2.5,10.5), y=c(10.5,10.5), lwd = 2, col = "red")
  panel.lines(x=c(2.5,2.5), y=c(6.5,10.5), lwd = 2, col = "red")
  panel.lines(x=c(10.5,10.5), y=c(6.5,10.5), lwd = 2, col = "red")
})

i = "wspF"
subset     <- blasted[blasted$Gene == i,]
seqlengths <- abs(subset$S.end - subset$S.start) + 1
identity   <- subset$Perc.Ident
data       <- table(data.frame(x = cut(identity, breaks = seq(70, 100, length = 11)), y = cut(seqlengths, breaks = seq(min(seqlengths), max(seqlengths), length = 11))))
data       <- data/sum(data)
plot7      <- levelplot(data, xlab = "Percent Identity", ylab = "Sequenz-Länge", main = i, col.regions = col, scales = list(x = list(rot = 90)), panel = function(...) {
  panel.levelplot(...)
  panel.lines(x=c(1.5,10.5), y=c(9.5,9.5), lwd = 2, col = "red")
  panel.lines(x=c(1.5,10.5), y=c(10.5,10.5), lwd = 2, col = "red")
  panel.lines(x=c(1.5,1.5), y=c(9.5,10.5), lwd = 2, col = "red")
  panel.lines(x=c(10.5,10.5), y=c(9.5,10.5), lwd = 2, col = "red")
})

i = "wspR"
subset     <- blasted[blasted$Gene == i,]
seqlengths <- abs(subset$S.end - subset$S.start) + 1
identity   <- subset$Perc.Ident
data       <- table(data.frame(x = cut(identity, breaks = seq(70, 100, length = 11)), y = cut(seqlengths, breaks = seq(min(seqlengths), max(seqlengths), length = 11))))
data       <- data/sum(data)
plot8      <- levelplot(data, xlab = "Percent Identity", ylab = "Sequenz-Länge", main = i, col.regions = col, scales = list(x = list(rot = 90)), panel = function(...) {
  panel.levelplot(...)
  panel.lines(x=c(3.5,10.5), y=c(7.5,7.5), lwd = 2, col = "red")
  panel.lines(x=c(3.5,10.5), y=c(10.5,10.5), lwd = 2, col = "red")
  panel.lines(x=c(3.5,3.5), y=c(7.5,10.5), lwd = 2, col = "red")
  panel.lines(x=c(10.5,10.5), y=c(7.5,10.5), lwd = 2, col = "red")
})

# Apply filters
# Filtering BLAST results
rph.blast  <- blasted[blasted$Gene == "rph",]
rph.length <- abs(rph.blast$S.end - rph.blast$S.start) + 1
rph.pres   <- ifelse(rph.blast$Perc.Ident >= 79 & rph.blast$Perc.Ident <= 100 & rph.length >=683 & rph.length <= max(rph.length), 1, 0)

wspA.blast  <- blasted[blasted$Gene == "wspA",]
wspA.length <- abs(wspA.blast$S.end - wspA.blast$S.start) + 1
wspA.pres   <- ifelse(wspA.blast$Perc.Ident >= 82 & wspA.blast$Perc.Ident <= 100 & wspA.length >= 850 & wspA.length <= max(wspA.length), 1, 0)

wspB.blast  <- blasted[blasted$Gene == "wspB",]
wspB.length <- abs(wspB.blast$S.end - wspB.blast$S.start) + 1
wspB.pres   <- ifelse(wspB.blast$Perc.Ident >= 76 & wspB.blast$Perc.Ident <= 100 & wspB.length >= 488 & wspB.length <= max(wspB.length), 1, 0)

wspC.blast  <- blasted[blasted$Gene == "wspC",]
wspC.length <- abs(wspC.blast$S.end - wspC.blast$S.start) + 1
wspC.pres   <- ifelse(wspC.blast$Perc.Ident >= 73 & wspC.blast$Perc.Ident <= 100 & wspC.length >= 680 & wspC.length <= max(wspC.length), 1, 0)

wspD.blast  <- blasted[blasted$Gene == "wspD",]
wspD.length <- abs(wspD.blast$S.end - wspD.blast$S.start) + 1
wspD.pres   <- ifelse(wspD.blast$Perc.Ident >= 70 & wspD.blast$Perc.Ident <= 100 & wspD.length >= 637 & wspD.length <= max(wspD.length), 1, 0)

wspE.blast  <- blasted[blasted$Gene == "wspE",]
wspE.length <- abs(wspE.blast$S.end - wspE.blast$S.start) + 1
wspE.pres   <- ifelse(wspE.blast$Perc.Ident >= 76 & wspE.blast$Perc.Ident <= 100 & wspE.length >= 1440 & wspE.length <= max(wspE.length), 1, 0)

wspF.blast  <- blasted[blasted$Gene == "wspF",]
wspF.length <- abs(wspF.blast$S.end - wspF.blast$S.start) + 1
wspF.pres   <- ifelse(wspF.blast$Perc.Ident >= 73 & wspF.blast$Perc.Ident <= 100 & wspF.length >= 936 & wspF.length <= max(wspF.length), 1, 0)

wspR.blast  <- blasted[blasted$Gene == "wspR",]
wspR.length <- abs(wspR.blast$S.end - wspR.blast$S.start) + 1
wspR.pres   <- ifelse(wspR.blast$Perc.Ident >= 79 & wspR.blast$Perc.Ident <= 100 & wspR.length >= 742 & wspR.length <= max(wspR.length), 1, 0)

# Identify gene presence
sequences$rph  <- rep(NA, nrow(sequences))
sequences$wspA <- rep(NA, nrow(sequences))
sequences$wspB <- rep(NA, nrow(sequences))
sequences$wspC <- rep(NA, nrow(sequences))
sequences$wspD <- rep(NA, nrow(sequences))
sequences$wspE <- rep(NA, nrow(sequences))
sequences$wspF <- rep(NA, nrow(sequences))
sequences$wspR <- rep(NA, nrow(sequences))

for (i in 1:nrow(sequences)) {
  sequences$rph[i]  <- ifelse(sum(rph.pres[which(rph.blast$SubjectID %in% sequences$NCBI.reference.nr[i])]) > 0, 1, 0)
  sequences$wspA[i] <- ifelse(sum(wspA.pres[which(wspA.blast$SubjectID %in% sequences$NCBI.reference.nr[i])]) > 0, 1, 0)
  sequences$wspB[i] <- ifelse(sum(wspB.pres[which(wspB.blast$SubjectID %in% sequences$NCBI.reference.nr[i])]) > 0, 1, 0)
  sequences$wspC[i] <- ifelse(sum(wspC.pres[which(wspC.blast$SubjectID %in% sequences$NCBI.reference.nr[i])]) > 0, 1, 0)
  sequences$wspD[i] <- ifelse(sum(wspD.pres[which(wspD.blast$SubjectID %in% sequences$NCBI.reference.nr[i])]) > 0, 1, 0)
  sequences$wspE[i] <- ifelse(sum(wspE.pres[which(wspE.blast$SubjectID %in% sequences$NCBI.reference.nr[i])]) > 0, 1, 0)
  sequences$wspF[i] <- ifelse(sum(wspF.pres[which(wspF.blast$SubjectID %in% sequences$NCBI.reference.nr[i])]) > 0, 1, 0)
  sequences$wspR[i] <- ifelse(sum(wspR.pres[which(wspR.blast$SubjectID %in% sequences$NCBI.reference.nr[i])]) > 0, 1, 0)
}


################################
# Add lifestyle classification #
################################

behavior <- read.csv("classification.csv")
behavior <- behavior[-which(behavior$NCBI.reference.nr %in% duplicate),]

sequences$free.living        <- behavior$free.living
sequences$disease.associated <- behavior$disease.associated

write.csv(sequences, "final summary.csv", row.names = F)

final <- read.csv("final summary.csv")

final$class <- ifelse(final$free.living == "free.living", "free-living", ifelse(final$free.living == "unknown", "unknown", ifelse(final$disease.associated == "disease.associated", "in host - diseased", "in host - no/unknown disease")))

summary <- as.data.frame(matrix(nrow = 48, ncol = 4))
colnames(summary) <- c("gene", "gene presence", "class", "number")
summary[,1] <- c(rep("rph", 6), rep("wspA", 6), rep("wspB", 6), rep("wspC", 6), rep("wspD", 6), rep("wspE", 6), rep("wspF", 6), rep("wspR", 6))
summary[,2] <- rep(c(0,0,0,1,1,1), 8)
summary[,3] <- c("free-living", "in host - diseased", "in host - no/unknown disease")
for (i in 1:nrow(summary)) {
  summary[i, 4] <- length(which(final[,which(colnames(final) == summary$gene[i])] == summary$`gene presence`[i] & final$class == summary$class[i]))
}


############
# Analyses #
############

# Lifestyles in absence / presence of individual rph & wsp genes
summary$rel.num <- rep(NA, nrow(summary))

for (i in unique(summary$gene)) {
  vals <- summary$number[which(summary$gene == i & summary$`gene presence` == 0)]
  summary$rel.num[which(summary$gene == i & summary$`gene presence` == 0)] <- vals/sum(vals)
  vals <- summary$number[which(summary$gene == i & summary$`gene presence` == 1)]
  summary$rel.num[which(summary$gene == i & summary$`gene presence` == 1)] <- vals/sum(vals)
}

summary$`gene presence` <- as.factor(summary$`gene presence`)
summary$`gene presence` <- as.factor(ifelse(summary$`gene presence` == 0, "absent", "present"))

# Plots
library(ggplot2)
ggplot(aes(x = 1, y = rel.num, fill = class), data = summary) +
  geom_bar(stat = "identity") +
  facet_grid(gene ~ `gene presence`) +
  labs(x = "Gene presence", y = "Gene") +
  scale_fill_manual(values = c("#2e4057ff", "#8d163cff", "#ba6a6a")) +
  theme(axis.ticks.x = element_blank(), axis.text.x = element_blank(), axis.ticks.y = element_blank(), axis.text.y = element_blank(), strip.text.y = element_text(face = "italic"))

# Chi square tests
final <- read.csv("final summary.csv")

chi.tests <- as.data.frame(matrix(nrow = 8, ncol = 6))
colnames(chi.tests) <- c("Gen", "N(absent)", "N(present)", "Chi^2", "df", "p-Wert")
for (i in 1:length(unique(summary$gene))) {
  subset <- summary[which(summary$gene == unique(summary$gene)[i]),]
  chisq  <- chisq.test(x = subset$number[1:3], p = subset$number[4:6]/sum(subset$number[4:6]))
  chi.tests[i, 1] <- unique(summary$gene)[i]
  chi.tests[i, 2] <- sum(subset$number[subset$`gene presence` == "absent"])
  chi.tests[i, 3] <- sum(subset$number[subset$`gene presence` == "present"])
  chi.tests[i, 4] <- chisq$statistic
  chi.tests[i, 5] <- chisq$parameter
  chi.tests[i, 6] <- chisq$p.value
}


# Lifestyles of complete / incomplete wsp Operon
final <- read.csv("final summary.csv")

final$class <- ifelse(final$free.living == "free.living" | final$free.living == "unknown", "free-living/unknown", ifelse(final$disease.associated == "disease.associated", "in host - diseased", "in host - no/unknown disease"))

allwsps <- ifelse(rowSums(final[,11:17]) == 7, 1, 0)
summary <- cbind.data.frame(wsp.operon = c("incomplete", "incomplete", "incomplete", "complete", "complete", "complete"), class = summary$class[1:6])
summary$number <- rep(NA, 6)
summary$number[1:3] <- table(final$class[allwsps==0])
summary$number[4:6] <- table(final$class[allwsps==1])
summary$rel.num <- c(summary$number[1:3]/sum(summary$number[1:3]), summary$number[4:6]/sum(summary$number[4:6]))

# Plot
library(ggplot2)
ggplot(aes(x = 1, y = rel.num, fill = class), data = summary) +
  geom_bar(stat = "identity") +
  facet_grid(. ~ wsp.operon) +
  labs(x = "Gene presence", y = "Gene") +
  scale_fill_manual(values = c("#2e4057ff", "#8d163cff", "#ba6a6a")) +
  theme(axis.ticks.x = element_blank(), axis.text.x = element_blank(), axis.ticks.y = element_blank(), axis.text.y = element_blank(), strip.text.y = element_text(face = "italic"))

# Chi square tests
test <- chisq.test(x = summary$number[1:3], p = summary$number[4:6] / sum(summary$number[4:6]))
chi.tests <- rbind.data.frame(chi.tests, c("wsp operon", sum(summary$number[1:3]), sum(summary$number[4:6]), test$statistic, test$parameter, test$p.value))