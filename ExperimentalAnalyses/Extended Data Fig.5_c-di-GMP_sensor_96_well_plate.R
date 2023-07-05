#----------------------------------
# 2023
# Check c-di-GMP concentrations in different MYb11 morphotypes with plasmid sensor in 96-well-plates
#----------------------------------

#install.packages('openxlsx')
#install.packages("rstatix")

library(ggplot2)
library(openxlsx)
library(reshape2)
library(rstatix)
library(lme4)
library(multcomp)

# Setting the directory for reading in data
currentDirectory <- dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(currentDirectory)

#--- (1) Collect data -------------------------------

data1<-read.xlsx("Raw_Sensor_96well_1.xlsx")
data2<-read.xlsx("Raw_Sensor_96well_2.xlsx")
data3<-read.xlsx("Raw_Sensor_96well_3.xlsx")
data4<-read.xlsx("Raw_Sensor_96well_4.xlsx")

# Subset data of interest
#Bottomread
SubdataB<-function(data)
{
dataAmCyan <- data[46:54, 1:11]
dataRFP <- data[72:80, 1:11]

# Create 96-well-plate layout
dataCy <- as.matrix(dataAmCyan[2:9,2:11])
rownames(dataCy) <- dataAmCyan[2:9,1]
colnames(dataCy) <- dataAmCyan[1,2:11]

dataR <- as.matrix(dataRFP[2:9,2:11])
rownames(dataR) <- dataRFP[2:9,1]
colnames(dataR) <- dataRFP[1,2:11]

# Transform matrix into column format
dataCy<- melt(dataCy)
well<-paste(dataCy$Var1,dataCy$Var2,sep='')
tmp<-data.frame(well,valueC=dataCy$value);tmp

dataR<- melt(dataR)
well<-paste(dataR$Var1,dataR$Var2,sep='')
tmpR<-data.frame(well,valueR=dataR$value);tmp

dataRep <-(merge(tmp,tmpR, by = "well"))
return(dataRep)
}

SubdataB2<-function(data)
{
  dataAmCyan <- data[46:54, 1:5]
  dataRFP <- data[72:80, 1:5]
  
  # Create 96-well-plate layout
  dataCy <- as.matrix(dataAmCyan[2:9,2:5])
  rownames(dataCy) <- dataAmCyan[2:9,1]
  colnames(dataCy) <- dataAmCyan[1,2:5]
  
  dataR <- as.matrix(dataRFP[2:9,2:5])
  rownames(dataR) <- dataRFP[2:9,1]
  colnames(dataR) <- dataRFP[1,2:5]
  
  # Transform matrix into column format
  dataCy<- melt(dataCy)
  well<-paste(dataCy$Var1,dataCy$Var2,sep='')
  tmp<-data.frame(well,valueC=dataCy$value);tmp
  
  dataR<- melt(dataR)
  well<-paste(dataR$Var1,dataR$Var2,sep='')
  tmpR<-data.frame(well,valueR=dataR$value);tmp
  
  dataRep <-(merge(tmp,tmpR, by = "well"))
  return(dataRep)
}

dataRep1<- SubdataB(data1)
dataRep2<- SubdataB2(data2)
dataRep3<- SubdataB2(data3)
dataRep4<- SubdataB2(data4)

# Add data legend
addLegend<-function(data)
{
data <- merge(data, legend, by = "well")
data <- merge(data, legend2, by = "id")
data$id <- as.factor(data$id)
data$strain <- as.factor(data$strain)
data$valueR <- as.numeric(data$valueR)
data$valueC <- as.numeric(data$valueC)
return(data)
}

legend <- read.csv("layout_1.csv")
legend2 <- read.csv("strainIDs_1.csv", header = T)
dataRep1<-addLegend(dataRep1)
rm(legend)

legend <- read.csv("layout_2.csv")
legend2 <- read.csv("strainIDs_2.csv", header = T)
dataRep2<-addLegend(dataRep2)
rm(legend)

legend <- read.csv("layout_3.csv")
legend2 <- read.csv("strainIDs_3.csv", header = T)
dataRep3<-addLegend(dataRep3)
rm(legend)

legend <- read.csv("layout_4.csv")
legend2 <- read.csv("strainIDs_4.csv", header = T)
dataRep4<-addLegend(dataRep4)
rm(legend)

#--- (2) Determine RFI data -------------------------------

RFI<-function(data)
{
#blank raw data
datactr <- data[data$strain == 'PBS',]
datactr$C <- mean(datactr$valueC)
datactr$R <- mean(datactr$valueR)
datactr<-subset(datactr, select = c(R,C),)
datactr<-datactr[!duplicated(datactr),]
data<-cbind(data,datactr)
data$valueC<-data$valueC-data$C
data$valueR<-data$valueR-data$R

#RFI
data$RFI<-data$valueR/data$valueC

#mean RFI
data$meanRFI <-ave(data$RFI, data$strain,data$rep, FUN = median)

#relative RFI (Mutant/WT)
datactr<-data[data$strain == 'WT',]
datactr$meanRFIWT <-datactr$meanRFI
datactr<-subset(datactr, select = c(meanRFIWT,rep))
datactr<-datactr[!duplicated(datactr),]
data <- merge(data, datactr, by = "rep")
data$Rel<-data$meanRFI/data$meanRFIWT
return(data)
}

dataRep1 <- RFI(dataRep1)
dataRep1$rep<-as.factor(dataRep1$rep)

dataRep2 <- RFI(dataRep2)
dataRep2$rep<-as.factor(dataRep2$rep)

dataRep3 <- RFI(dataRep3)
dataRep3$rep<-as.factor(dataRep3$rep)

dataRep4 <- RFI(dataRep4)
dataRep4$rep<-as.factor(dataRep4$rep)

dataRep5 <-rbind(dataRep1,dataRep2,dataRep3,dataRep4)

#remove replicate 1 -> the controls did not work on this run
dataRep5<-dataRep5[dataRep5$rep !="1",]
dataRep5<-dataRep5[dataRep5$strain !='MT262',]
dataRep5<-dataRep5[dataRep5$strain !='PBS',]

#--- (3) Plot data -------------------------------
####9. Boxplot ####

my_theme <- function(data){
  
  # Set color palettes
  library(RColorBrewer)
  
  palette <- brewer.pal("Greys", n=9)
  col.back <-"grey75"
  col.grid <- palette[3]
  col.text <- palette[7]
  
  # Set plotting theme
  theme_bw(base_size=9) +
    
    # Set the entire chart region to a light gray color
    #theme(panel.background=element_rect(fill=col.back, color=col.back)) +
    #theme(plot.background=element_rect(fill=col.back, color=col.back)) +
    theme(panel.border=element_rect(color=col.back)) +
    theme(strip.background = element_rect(color="white",fill=col.back))+
    
    # Format the grid
    theme(panel.grid.major=element_line(color=col.grid,size=.25)) +
    theme(panel.grid.minor=element_blank()) +
    theme(axis.ticks=element_blank()) +
    
    # legend
    theme(legend.position="right") +
    #theme(legend.background = element_rect(fill=color.background)) +
    theme(legend.text = element_text(size=12,color=col.text)) +
    
    # Set title and axis labels, and format these and tick marks
    theme(plot.title=element_text(color=col.text,size=10,vjust=1.25,hjust=0.5, face = "italic")) +
    theme(axis.text.x=element_text(size=14,color=col.text,angle=45,hjust=1)) +
    theme(axis.text.y=element_text(size=14,color=col.text)) +
    theme(axis.title.x=element_text(size=14,color=col.text, vjust=0)) +
    theme(axis.title.y=element_text(size=14,color=col.text, vjust=1.25)) +
    theme(strip.text=element_text(size=14,color="grey30"))  }

# RFIs
scientific_10 <- function(x) {
  parse(text=gsub("e\\+*", " %*% 10^", scales::scientific_format()(x)))}

PlotRFI<- function(Data){
  Data$strain <- factor(Data$strain,levels = c('WT','MT12','MT21','MT25','MT26','MT13','MT33','MT11'),ordered = TRUE)
  Data<-Data[!duplicated(Data$meanRFI),]
  Pop_Plot <- ggplot(Data,aes(x=strain,y=meanRFI))+
    geom_boxplot(aes())+
    geom_point(aes(shape = rep), position= position_jitterdodge(jitter.width = 0.1, jitter.height = 0, dodge.width = 0.1,seed=NULL))+
    guides( shape = guide_legend("Replicate"))+
    scale_y_continuous(limits=c(0.5, 3),labels = scientific_10 )+
    #scale_x_discrete(labels=c("MT12" = "MYb11"~italic("wspF*")^"EVO", "MT14" = "MYb11"~italic("wspE*")^"EVO","MT22" = "MYb11"~italic("rph*")^"EVO","WT" = "MYb11"))+
    my_theme()+
    labs(x="Bacteria")+
    labs(y="RFI (TurboRFP / AmCyan)")
  
  return(Pop_Plot)
}

PlotRFI(dataRep5)

#ggsave("RFI_RAW.pdf", height=5, width=9,dpi=300)

PlotRel<- function(Data){
  Data$strain <- factor(Data$strain,levels = c('WT','MT12','MT21','MT25','MT26','MT13','MT33','MT11'),ordered = TRUE)
  Data<-Data[!duplicated(Data$Rel),]
  Pop_Plot <- ggplot(Data,aes(x=strain,y=Rel))+
    geom_boxplot(aes())+
    geom_point(aes(fill = "black"), position= position_jitterdodge(jitter.width = 0.1, jitter.height = 0, dodge.width = 0,seed=NULL))+
    guides( shape = guide_legend("Replicate"))+
    scale_y_continuous(limits=c(0.5, 3),labels = scientific_10 )+
    #scale_x_discrete(labels=c("MT12" = "MYb11"~italic("wspF*")^"EVO", "MT14" = "MYb11"~italic("wspE*")^"EVO","MT22" = "MYb11"~italic("rph*")^"EVO","WT" = "MYb11"))+
    my_theme()+
    labs(x="Bacteria")+
    labs(y="RFI (TurboRFP / AmCyan)")
  
  return(Pop_Plot)
}

PlotRel(dataRep5)

#ggsave("RFI_Rel.pdf", height=5, width=9,dpi=300)

######Stats
# prepare data

data<-dataRep5
data<- subset(data, select = c(rep,strain,RFI))
data$strain <- factor(data$strain,levels = c('WT','MT11','MT12','MT13','MT21','MT25','MT26','MT33'),ordered = TRUE)
data %>%
  group_by(strain) %>%
  shapiro_test(RFI)
levene_test(data, RFI ~ strain)

#RFI data

TestSensorRFI <- function(data,Set){
  
  # Prepare output file name
  fileName <- paste("stats out/Sensor/",Set,"_Sensor",".txt",sep="")
  
  # Run nested aov
  m <-aov(RFI ~ strain/rep,data=data)
  
  sink(fileName)
  print(summary(m)) # aov output
  print("")
  print(ml<-glht(m, mcp(strain="Dunnett"))) # Post-doc output: Dunnett
  print(summary(ml,adjusted("fdr")))
  sink()
  
  return("Regression analysis done and saved.")
}

TestSensorRFI(data,Set="RFI")

#normalized data

data2<-dataRep5
data2<-data2[!duplicated(data2$Rel),]
data2<- subset(data2, select = c(rep,strain,Rel))
test<-read.csv("Rel.csv")
data2<-rbind(data2,test)
data2$strain <- factor(data2$strain,levels = c('WT','MT11','MT12','MT13','MT21','MT25','MT26','MT33'),ordered = TRUE)
data3<-data2[data2$strain!='WT',]
data3 %>%
  group_by(strain) %>%
  shapiro_test(Rel)
levene_test(data2, Rel ~ strain)

TestSensorRel <- function(data,Set){
  
  # Prepare output file name
  fileName <- paste("stats out/Sensor/",Set,"_Sensor",".txt",sep="")
  
  # Run aov
  m <-aov(Rel ~ strain,data=data)
  
  sink(fileName)
  print(summary(m)) # aov output
  print("")
  print(ml<-glht(m, mcp(strain="Dunnett"))) # Post-doc output: Dunnett
  print(summary(ml,adjusted("fdr")))
  sink()
  
  return("Regression analysis done and saved.")
}

TestSensorRel(data2,Set="Rel")
