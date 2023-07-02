#########################################
# C-di-GMP concentration via fluoresces # 
#########################################

####1. get packages ####

#install.packages("rstudioapi")
#install.packages("ggplot2")
#install.packages("ggnewscale")
#install.packages("car")
#install.packages("plyr")
#install.packages("MASS")
#install.packages("lme4")
#install.packages("multcomp")
#install.packages("scales")
#install.packages("ggpubr")
#install.packages("lmtest")
#install.packages("installr")
#install.packages("rstatix")

####2. load packages####

library(rstudioapi)
library(ggplot2)
library(ggnewscale)
library(car)
library(plyr)
library(MASS)
library(lme4)
library(multcomp)
library(scales)
library(ggpubr)
library(stats)
library(rstatix)
library(installr)

####3. get name of the directory this file is in####

currentDirectory <- dirname(rstudioapi::getActiveDocumentContext()$path)

####4. set working directory to currentDirectory####

setwd(currentDirectory)

####5. get raw data ####

cdiGMP2 <-read.csv("RawData_Rep2.csv",TRUE,",")
cdiGMP3 <-read.csv("RawData_Rep3.csv",TRUE,",")
cdiGMP4 <-read.csv("RawData_Rep4.csv",TRUE,",")
cdiGMP5 <-read.csv("RawData_Rep5.csv",TRUE,",")
cdiGMP6 <-read.csv("RawData_Rep6.csv",TRUE,",")
MS<-read.csv("MS.csv",TRUE,",")


####6. function for raw data ####

Rawdata <- function(Data){
  # rename column
  names(Data)
  names(Data)[1]<-"MT"

  # edit data
  Data$MT<-as.factor(Data$MT)
  Data$Channel<-as.factor(Data$Channel)
  Data$Label<-as.factor(Data$Label)
  Data$TechRep<-as.factor(Data$TechRep)
  Data$Biorep<-as.factor(Data$Biorep)
  Data$Day<-as.factor(Data$Day)
  return(Data)
}

####6.1 apply function####

cdiGMP2<-Rawdata(cdiGMP2)
cdiGMP3<-Rawdata(cdiGMP3)
cdiGMP4<-Rawdata(cdiGMP4)
cdiGMP5<-Rawdata(cdiGMP5)
cdiGMP6<-Rawdata(cdiGMP6)

####7. function for CTCF ####

CTCF<-function(Data){
  
  # Subset two data frames by MT
  WT <- Data[Data$MT == "WT", ]
  WTBG <- Data[Data$MT == "WTBG", ]
  MT12 <- Data[Data$MT == "MT12", ]
  MT12BG <- Data[Data$MT == "MT12BG", ]
  MT14 <- Data[Data$MT == "MT14", ]
  MT14BG <- Data[Data$MT == "MT14BG", ]
  MT22 <- Data[Data$MT == "MT22", ]
  MT22BG <- Data[Data$MT == "MT22BG", ]
  
  # Mean
  WTBG$meanBG <- ave(WTBG$Mean, WTBG$Channel,WTBG$TechRep,WTBG$Day, FUN=mean)
  MT12BG$meanBG <- ave(MT12BG$Mean, MT12BG$Channel, MT12BG$TechRep, MT12BG$Day, FUN=mean)
  MT14BG$meanBG <- ave(MT14BG$Mean, MT14BG$Channel, MT14BG$TechRep, MT14BG$Day, FUN=mean)
  MT22BG$meanBG <- ave(MT22BG$Mean, MT22BG$Channel, MT22BG$TechRep, MT22BG$Day, FUN=mean)
  
  # Merge data frames
  
  MT12BG<-subset (MT12BG, select = -c(MT,Label,Area, Mean, Min, Max, IntDen, RawIntDen))
  MT12BG<-MT12BG[!duplicated(MT12BG$meanBG),]
  MT12<- merge(MT12, MT12BG, by = c("Channel","TechRep", "Day", "Biorep"),all =T)

  MT14BG<-subset (MT14BG, select = -c(MT,Label,Area, Mean, Min, Max, IntDen, RawIntDen))
  MT14BG<-MT14BG[!duplicated(MT14BG$meanBG),]
  MT14<- merge(MT14, MT14BG, by = c("Channel","TechRep", "Day", "Biorep"),all =T)
  
  MT22BG<-subset (MT22BG, select = -c(MT,Label,Area, Mean, Min, Max, IntDen, RawIntDen))
  MT22BG<-MT22BG[!duplicated(MT22BG$meanBG),]
  MT22<- merge(MT22, MT22BG, by = c("Channel","TechRep", "Day", "Biorep"),all =T)
  
  WTBG<-subset (WTBG, select = -c(MT,Label,Area, Mean, Min, Max, IntDen, RawIntDen))
  WTBG<-WTBG[!duplicated(WTBG$meanBG),]
  WT<- merge(WT, WTBG, by = c("Channel","TechRep", "Day", "Biorep"),all =T)
  
  CdiGMP<- rbind(MT12,MT14,MT22,WT)
  
  #calculate CTCF
  
  CdiGMP$CTCF<-(CdiGMP$IntDen-(CdiGMP$Area * CdiGMP$meanBG))
  
  return(CdiGMP)
}

CdiGMP2<-CTCF(cdiGMP2)
CdiGMP3<-CTCF(cdiGMP3)
CdiGMP4<-CTCF(cdiGMP4)
CdiGMP5<-CTCF(cdiGMP5)
CdiGMP6<-CTCF(cdiGMP6)

####8. Function RFI####

RFI<-function(Data){
  
  #Subset Data by Channel
  Data<-subset (Data, select = -c(Area, Mean, Min, Max, IntDen, RawIntDen,meanBG))
  Ch1 <- Data[Data$Channel == "1", ]
  Ch0 <- Data[Data$Channel == "0", ]
  
  # Merge data frames
  
  Data<- merge(Ch1, Ch0, by = c("MT","Label","TechRep", "Day", "Biorep"),all =T)
  
  #RFI
  Data$RFI<-(Data$CTCF.y/Data$CTCF.x)
  
  return(Data)
}

RFI2<-RFI(CdiGMP2)
RFI3<-RFI(CdiGMP3)
RFI4<-RFI(CdiGMP4)
RFI5<-RFI(CdiGMP5)
RFI6<-RFI(CdiGMP6)

##Randomise + mean RFI

Ran2<-function(Data){
  
D3 <- Data[Data$Day == "3", ]
D5 <- Data[Data$Day == "5", ]
set.seed(35)  
RWT<-D3[ sample( which(D3$MT=='WT'), 45 ), ]
RWT5<-D5[ sample( which(D5$MT=='WT'), 45 ), ]
R12<-D3[ sample( which(D3$MT=='MT12'), 45 ), ]
R125<-D5[ sample( which(D5$MT=='MT12'), 45 ), ]
R14<-D3[ sample( which(D3$MT=='MT14'), 45 ), ]
R145<-D5[ sample( which(D5$MT=='MT14'), 45 ), ]
R22<-D3[ sample( which(D3$MT=='MT22'), 45 ), ]
R225<-D5[ sample( which(D5$MT=='MT22'), 45 ), ]
Data<-(rbind(RWT,R12,R14,R22))

Data$meanRFI <- ave(Data$RFI, Data$MT, Data$Day, FUN=median)

return(Data)
}

Random2<-Ran2(RFI2)
Random3<-Ran2(RFI3)
Random4<-Ran2(RFI4)
Random5<-Ran2(RFI5)
Random6<-Ran2(RFI6)

Over<- (rbind(Random2, Random3, Random4, Random5,Random6))
Over<-Over[!duplicated(Over$meanRFI),]
Over$meanBiorep<- ave(Over$meanRFI, Over$MT, Over$Day, FUN=mean)
Over2<-Over[!duplicated(Over$meanBiorep),]

    Over2<-subset (Over2, select = -c(Label,TechRep,Biorep,Channel.x,Channel.y,CTCF.x,CTCF.y,RFI,meanRFI))
    WT <- Over2[Over2$MT == "WT", ]
    WT3 <- WT[WT$Day == "3", ]
    WT5 <- WT[WT$Day == "5", ]
    Over5<- Over2[Over2$Day == "5", ]
    Over3 <- Over2[Over2$Day == "3", ]
    Over3 <- cbind(Over3,WT3)
    Over5 <- cbind(Over5,WT5)
    Over2 <-rbind(Over3,Over5)
    names(Over2)[3]<-"meanBiorep.x"
    names(Over2)[4]<-"MTY"
    names(Over2)[5]<-"Dayy"
    Over2$'Rel' <- (Over2$meanBiorep.x  / Over2$meanBiorep)
    Over2$Rel<-round(Over2$Rel,digits = 1)
    
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
    
    # Hide the legend
    theme(legend.position="none") +
    #theme(legend.background = element_rect(fill=color.background)) +
    #theme(legend.text = element_text(size=12,color=col.text)) +
    
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

PlotCFUs<- function(Data){
  Data$MT <- factor(Data$MT,levels = c('WT','MT12','MT14','MT22'),ordered = TRUE)
  Pop_Plot <- ggplot(Data,aes(x=MT,y=meanRFI))+
    geom_boxplot(aes())+
    geom_point(aes(fill="black"), position= position_jitterdodge(jitter.width = 0.2, jitter.height = 0, dodge.width = 0.75,seed=NULL))+
    guides( shape = guide_legend("Replicate"))+
    scale_x_discrete(labels=c("MT12" = "MYb11"~italic("wspF*")^"EVO", "MT14" = "MYb11"~italic("wspE*")^"EVO","MT22" = "MYb11"~italic("rph*")^"EVO","WT" = "MYb11"))+
    my_theme()+
    labs(x="Bacteria")+
    labs(y="RFI (TurboRFP / AmCyan)")
  
  return(Pop_Plot)
}

PlotCFUs(Over)

PlotMS<- function(Data){
  Data$MT <- factor(Data$MT,levels = c('WT','MT12','MT14','MT22'),ordered = TRUE)
  Pop_Plot <- ggplot(Data,aes(x=MT,y=cdiGMPnom))+
    geom_boxplot(aes())+
    geom_point(aes(fill="black"), position= position_jitterdodge(jitter.width = 0.2, jitter.height = 0, dodge.width = 0,seed=NULL))+
    guides(fill=FALSE)+
    guides( shape = guide_legend("Replicate"))+
    scale_x_discrete(labels=c("MT12" = "MYb11"~italic("wspF*")^"EVO", "MT14" = "MYb11"~italic("wspE*")^"EVO","MT22" = "MYb11"~italic("rph*")^"EVO","WT" = "MYb11"))+
    my_theme()+
    labs(x="Bacteria")+
    labs(y=" Normalisied c-di-GMP (fmol c-di-GMP / ?g protein)")
  
  return(Pop_Plot)
}

PlotMS2<- function(Data){
  Data$MT <- factor(Data$MT,levels = c('WT','MT12','MT14','MT22'),ordered = TRUE)
  Pop_Plot <- ggplot(Data,aes(x=MT,y=totalcdiGMP))+
    geom_boxplot(aes())+
    geom_point(aes(fill="black"), position= position_jitterdodge(jitter.width = 0.3, jitter.height = 0, dodge.width = 0.3,seed=NULL))+
    scale_fill_manual (values = c("#ca0020","#0571b0","#018571","#a6611a"))+
    guides(fill=FALSE)+
    guides( shape = guide_legend("Replicate"))+
    scale_x_discrete(labels=c("MT12" = "MYb11"~italic("wspF*")^"EVO", "MT14" = "MYb11"~italic("wspE*")^"EVO","MT22" = "MYb11"~italic("rph*")^"EVO","WT" = "MYb11"))+
    my_theme()+
    labs(x="Bacteria")+
    labs(y=" Total c-di-GMP (pmol c-di-GMP)")
  
  return(Pop_Plot)
}

##### 10. preparation heatmap ####

MS2<-MS[!duplicated(MS$meannormalizedcdiGMP),]

names(MS)
names(MS)[1]<-"Biorep"
names(MS2)[1]<-"Biorep"

# edit data
MS$MT<-as.factor(MS$MT)
MS$Biorep<-as.factor(MS$Biorep)

PlotMS(MS)
PlotMS2(MS)

my_themeHM<- function(data){
  
  # Set color palettes
  library(RColorBrewer)
  
  palette <- brewer.pal("Greys", n=9)
  col.back <-"grey75"
  col.grid <- palette[5]
  col.text <- palette[9]
  
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
    
    # Hide the legend
    #theme(legend.position="none") +
    #theme(legend.background = element_rect(fill=color.background)) +
    #theme(legend.text = element_text(size=12,color=col.text)) +
    
    # Set title and axis labels, and format these and tick marks
    theme(plot.title=element_text(color=col.text,size=10,vjust=1.25,hjust=0.5, face = "italic")) +
    theme(axis.text.x=element_text(size=14,color=col.text,angle=45,hjust=1)) +
    theme(axis.text.y=element_blank()) +
    theme(axis.title.x=element_text(size=14,color=col.text,vjust=0)) +
    theme(axis.title.y =element_blank()) +
    theme(strip.text=element_text(size=14,color="grey30"))  }

Heat<-function(Data){
  Data$MT <- factor(Data$MT,levels = c('WT','MT12','MT14','MT22'),ordered = TRUE)
  Heat<-ggplot(Data,aes(MT,MTY, fill=Rel))+
    geom_tile(color="black")+
    scale_fill_gradient(name = "c-di-GMP relative
    to wildtype",low="white",high="red",breaks=c(1,2,3,4,5),limits=c(floor(1),ceiling(5)))+
    geom_text(aes(label = Rel), color = "black") +
    coord_fixed()+
    my_themeHM()+
    guides(fill = guide_colourbar(barwidth = 1,barheight = 2.5))+
    labs(x="Bacteria")+
    scale_x_discrete(labels=c("MT12" = "MYb11"~italic("wspF*")^"EVO", "MT14" = "MYb11"~italic("wspE*")^"EVO","MT22" = "MYb11"~italic("rph*")^"EVO","WT" = "MYb11"))
  
  return(Heat) }

Heat(Over2)

HeatMS<-function(Data){
  Data$MT <- factor(Data$MT,levels = c('WT','MT12','MT14','MT22'),ordered = TRUE)
  HeatMS<-ggplot(Data,aes(MT,Biorep, fill=cdiGMPRel))+
    geom_tile(color="black")+
    scale_fill_gradient(name = "c-di-GMP relative
    to wildtype",low="white",high="red",breaks=c(1,2,3,4,5),limits=c(floor(1),ceiling(5)))+
    geom_text(aes(label = cdiGMPRel), color = "black") +
    coord_fixed()+
    my_themeHM()+
    labs(x="Bacteria")+
    guides(fill = guide_colourbar(barwidth = 1,barheight = 2.5))+
    scale_x_discrete(labels=c("MT12" = "MYb11"~italic("wspF*")^"EVO", "MT14" = "MYb11"~italic("wspE*")^"EVO","MT22" = "MYb11"~italic("rph*")^"EVO","WT" = "MYb11"))
  
  return(HeatMS) }

HeatMS(MS2)

#### 11.Statistic ####

TestMS <- function(data,Set){
  # Subset data
  
  # Prepare output file name
  fileName <- paste("stats out/MS/",Set,"_Norm",".txt",sep="")
  
 
  m <-oneway.test(cdiGMPnom ~ MT, data = data, var.equal = FALSE)
  
  
  sink(fileName)
  print(m) # lmer output
  print("")
  print(games_howell_test(data, cdiGMPnom ~ MT, conf.level=0.95, detailed =FALSE))
  sink()
  
  return("Regression analysis done and saved.")
}


TestMS(MS,Set="c-di-GMP")

TestMS2 <- function(data,Set){
  # Subset data
  
  # Prepare output file name
  fileName <- paste("stats out/MS/",Set,"_Total",".txt",sep="")
  
  
  m <-oneway.test(totalcdiGMP ~ MT, data = data, var.equal = FALSE)
  
  
  sink(fileName)
  print(m) # output
  print("")
  print(games_howell_test(data, totalcdiGMP ~ MT, conf.level=0.95, detailed =FALSE))
  sink()
  
  return("Regression analysis done and saved.")
}


TestMS2(MS,Set="c-di-GMP")

TestSensor <- function(data,Set){
  # Subset data
  
  # Prepare output file name
  fileName <- paste("stats out/Sensor/",Set,"_Sensor",".txt",sep="")
  
  
  m <-oneway.test(meanRFI ~ MT, data = data, var.equal = FALSE)
  
  
  sink(fileName)
  print(m) # output
  print("")
  print(games_howell_test(data, meanRFI ~ MT, conf.level=0.95, detailed =TRUE))
  sink()
  
  return("Regression analysis done and saved.")
}


TestSensor(Over,Set="c-di-GMP")
