#########################################
# Short-term persistence competition    # 
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

####3. get name of the directory this file is in####

currentDirectory <- dirname(rstudioapi::getActiveDocumentContext()$path)

####4. set working directory to currentDirectory####

setwd(currentDirectory)

####5. get raw data ####

Competition <-read.csv("STP_RAW_1.csv",TRUE,";")
Competition2 <-read.csv("STP_RAW_2.csv",TRUE,";")
Competition3 <-read.csv("STP_RAW_3.csv",TRUE,";")
Competition4 <-read.csv("STP_RAW_4.csv",TRUE,";")
Competition5 <-read.csv("STP_RAW_5.csv",TRUE,";")
Competition6 <-read.csv("STP_RAW_6.csv",TRUE,",")


####6. function for raw data ####

Rawdata <- function(Data){
# rename column
names(Data)
names(Data)[1]<-"Replicate"
names(Data)[5]<-"1"
names(Data)[6]<-"2"
names(Data)[7]<-"3"
# edit data
Data$MT<-as.factor(Data$MT)
Data$Replicate<-as.factor(Data$Replicate)
Data$Set<-as.factor(Data$Set)
return(Data)
}

####6.1 apply function####

Competition<-Rawdata(Competition)
Competition2<-Rawdata(Competition2)
Competition3<-Rawdata(Competition3)
Competition4<-Rawdata(Competition4)
Competition5<-Rawdata(Competition5)
Competition6<-Rawdata(Competition6)

####7. function for CFUs/Worm ####

CFUSCompCol<-function(Data){
  
  # Subset two data frames, one with only worm data and one with only supernatant
  onlyWorms <- Data[Data$Typ == "W", ]
  onlySup <- Data[Data$Typ == "S", ]
  
  # Merge worm and sup data frames
  combiWormSup <- merge(onlyWorms, onlySup, by = c("Replicate", "MT","Set","Worms","Sample","Col"), all.x = T)
  
  # Mean
  combiWormSup$MeanW <- rowMeans(combiWormSup[(8:10)],TRUE)
  combiWormSup$MeanS <- rowMeans(combiWormSup[(13:15)],TRUE)
  
  # Mean/Worm
  combiWormSup$'MeanW/Worm' <- (combiWormSup$'MeanW'/combiWormSup$'Worms')
  combiWormSup$'MeanS/Worm' <- (combiWormSup$'MeanS'/combiWormSup$'Worms')
  
  # CFU/5?l/Worm
  combiWormSup$'MeanW/Worm/25?l' <- (combiWormSup$'MeanW/Worm'* 10^combiWormSup$'Dilution.x')
  combiWormSup$'MeanS/Worm/25?l' <- (combiWormSup$'MeanS/Worm'* 10^combiWormSup$'Dilution.y')
  
  # CFU/100?l/Worm
  combiWormSup$'MeanW/Worm/100?l' <- (combiWormSup$'MeanW/Worm/25?l'*4)
  combiWormSup$'MeanS/Worm/100?l' <- (combiWormSup$'MeanS/Worm/25?l'*4)
 
  # CFUs
  combiWormSup$'CFUs' <- (combiWormSup$'MeanW/Worm/100?l'-combiWormSup$'MeanS/Worm/100?l')
 
  #final Tabel
  CFUs <- combiWormSup[,c(1:5,25)]
  CFUs2<-log10(CFUs[6])
  CFUs<-cbind(CFUs,CFUs2)
  names(CFUs)[7]<-"LOG_CFU"
  return(CFUs)
}

CFUSCompCol2<-function(Data){
  
  # Subset two data frames, one with only worm data and one with only supernatants
  onlyWorms <- Data[Data$Typ == "W", ]
  onlySup <- Data[Data$Typ == "S", ]
  
  # Merge worm and sup data frames
  combiWormSup <- merge(onlyWorms, onlySup, by = c("Replicate", "MT","Set","Worms","Sample","Col"),
                        all.x = T)
  # all.x = T behält die komplette erste Tabelle bei
  
  # Mean
  combiWormSup$MeanW <- rowMeans(combiWormSup[(8:10)],TRUE)
  combiWormSup$MeanS <- rowMeans(combiWormSup[(13:15)],TRUE)
  
  # Mean/Worm
  combiWormSup$'MeanW/Worm' <- (combiWormSup$'MeanW'/combiWormSup$'Worms')
  combiWormSup$'MeanS/Worm' <- (combiWormSup$'MeanS'/combiWormSup$'Worms')
  
  # CFU/15µl/Worm
  combiWormSup$'MeanW/Worm/15µl' <- (combiWormSup$'MeanW/Worm'* 10^combiWormSup$'Dilution.x')
  combiWormSup$'MeanS/Worm/15µl' <- (combiWormSup$'MeanS/Worm'* 10^combiWormSup$'Dilution.y')
  
  # CFU/100µl/Worm
  combiWormSup$'MeanW/Worm/100µl' <- (combiWormSup$'MeanW/Worm/15µl'*6.666666666666667)
  combiWormSup$'MeanS/Worm/100µl' <- (combiWormSup$'MeanS/Worm/15µl'*6.666666666666667)
  
  # CFUs
  combiWormSup$'CFUs' <- (combiWormSup$'MeanW/Worm/100µl'-combiWormSup$'MeanS/Worm/100µl')
  
  #final Tabel
  CFUs <- combiWormSup[,c(1:5,25)]
  CFUs2<-log10(CFUs[6])
  CFUs<-cbind(CFUs,CFUs2)
  names(CFUs)[7]<-"LOG_CFU"
  return(CFUs)
}

####7.1 apply function####

CFUs_Comp <- CFUSCompCol(Competition)
CFUs_Comp2 <- CFUSCompCol(Competition2)
CFUs_Comp3 <- CFUSCompCol(Competition3)
CFUs_Comp4 <- CFUSCompCol(Competition4)
CFUs_Comp5 <- CFUSCompCol(Competition5)
CFUs_Comp6 <- CFUSCompCol2(Competition6)

####7.2 Prepare data for Box-Plot####

CFUs_Comp$CFUs[CFUs_Comp$CFUs< 0]<-NA
CFUs_Comp<- CFUs_Comp[!is.na(CFUs_Comp$CFUs),]
CFUs_Comp$MT<-as.factor(CFUs_Comp$MT)

CFUs_Comp2$CFUs[CFUs_Comp2$CFUs< 0]<-NA
CFUs_Comp2<- CFUs_Comp2[!is.na(CFUs_Comp2$CFUs),]
CFUs_Comp2$MT<-as.factor(CFUs_Comp2$MT)

CFUs_Comp3$CFUs[CFUs_Comp3$CFUs< 0]<-NA
CFUs_Comp3<- CFUs_Comp3[!is.na(CFUs_Comp3$CFUs),]
CFUs_Comp3$MT<-as.factor(CFUs_Comp3$MT)

CFUs_Comp4$CFUs[CFUs_Comp4$CFUs< 0]<-NA
CFUs_Comp4<- CFUs_Comp4[!is.na(CFUs_Comp4$CFUs),]
CFUs_Comp4$MT<-as.factor(CFUs_Comp4$MT)

CFUs_Comp5$CFUs[CFUs_Comp5$CFUs< 0]<-NA
CFUs_Comp5<- CFUs_Comp5[!is.na(CFUs_Comp5$CFUs),]
CFUs_Comp5$MT<-as.factor(CFUs_Comp5$MT)

CFUs_Comp6$CFUs[CFUs_Comp6$CFUs< 0]<-NA
CFUs_Comp6<- CFUs_Comp6[!is.na(CFUs_Comp6$CFUs),]
CFUs_Comp6$MT<-as.factor(CFUs_Comp6$MT)

####8. Boxplot ####

my_theme <- function(data){
  
  # Set color palettes
  library(RColorBrewer)
  
  palette <- brewer.pal("Greys", n=9)
  col.back <-"grey75"
  col.grid <- palette[6]
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


# CFUs (Quartet)
scientific_10 <- function(x) {
  parse(text=gsub("e\\+*", " %*% 10^", scales::scientific_format()(x)))}

PlotCFUs<- function(Data,Data2){
  Data$Set <- factor(Data$Set,levels = c('MTdT','MT11','MT13','MT12','MT14','MT11+MT12+MT13+MTdT','MT11+MT13+MT14+MTdT'),ordered = TRUE)
  Data$MT <- factor(Data$MT,levels = c('MTdT','MT11','MT13','MT12','MT14'),ordered = TRUE)
  Data2$Set <- factor(Data2$Set,levels = c('MTdT','MT11','MT13','MT12','MT14', 'MT11+MT12+MT13+MTdT','MT11+MT13+MT14+MTdT'),ordered = TRUE)
  Data2$MT <- factor(Data2$MT,levels = c('MTdT','MT11','MT13','MT12','MT14'),ordered = TRUE)
  Pop_Plot <- ggplot(Data,aes(x=Set,y=CFUs))+
    geom_boxplot(data=Data,aes(fill=MT))+
    geom_boxplot(data=Data2,aes(fill=MT))+
    scale_fill_manual (values = c("white","white","white","white","white"))+
    geom_point(data=Data,aes(fill=MT),size=2, position= position_jitterdodge(jitter.width = 0.2, jitter.height = 0, dodge.width = 0.75,seed=NA))+
    scale_y_log10(limits=c(1, 100000),labels = scientific_10 )+
    guides(fill=FALSE)+
    my_theme()+
    labs(y="Short-term persistence [CFU/worm]",x="Bacteria")+
    scale_x_discrete(labels=c(expression("MTdT" = "MYb11 dT","MT12" = "MYb11"~italic("wspF*")^"EVO","MT14" = "MYb11"~italic("wspE*")^"EVO",
                                         "MT11"="MYb11 fuzzy"^"EVO","MT13"="MYb11 smooth"^"EVO",
                                         "MT11+MT12+MT13+MTdT" = "MYb11 dT"+"MYb11 fuzzy"^"EVO"+"MYb11 smooth"^"EVO"+"MYb11"~italic("wspF*")^"EVO",
                                         "MT11+MT13+MT14+MTdT"="MYb11 dT"+"MYb11 fuzzy"^"EVO"+"MYb11 smooth"^"EVO"+"MYb11"~italic("wspE*")^"EVO")))
  return(Pop_Plot)
}

#### 8.1 apply function ####

PlotCFUs(CFUs_Comp[(CFUs_Comp$Set %in% c("MT11+MT12+MT13+MTdT","MT11+MT13+MT14+MTdT")),],CFUs_Comp[(CFUs_Comp$Set %in% c("MT11+MT12+MT13+MTdT","MT11+MT13+MT14+MTdT")),])

#### 8.2 plots CFU relative #####

CFUs_Comp4<- CFUs_Comp4[!is.infinite(CFUs_Comp4$LOG_CFU),]
CFUs_Comp3<- CFUs_Comp3[!is.infinite(CFUs_Comp3$LOG_CFU),]
CFUs_Comp5<- CFUs_Comp5[!is.infinite(CFUs_Comp5$LOG_CFU),]
CFUs_Comp6<- CFUs_Comp6[!is.infinite(CFUs_Comp6$LOG_CFU),]

#subset data relative competition
Rel1<- CFUs_Comp[CFUs_Comp$Set %in% c("MT12+MTdT","MT14+MTdT","MT48+MTdT"),]
Rel2<-CFUs_Comp2[CFUs_Comp2$Set %in% c("MT12*+MTdT","MT22+MTdT","SBW25+SBW25F*"),]
Rel3<-CFUs_Comp3[CFUs_Comp3$Set %in% c("MT14*+MTdT","193*+193","MT22*+MTdT"),]
Rel4<- CFUs_Comp4[CFUs_Comp4$Set %in% c("MT48+MT14S","MT187F+MT187","MT48+MTdT"),]
Rel5<- CFUs_Comp5
Rel6<-CFUs_Comp6
Rel<-rbind(Rel1,Rel2,Rel3)

REL1<-CFUs_Comp2[CFUs_Comp2$Set %in% c("MT12S+MT48","MT48+MTdT"),]
REL2<-CFUs_Comp3[CFUs_Comp3$Set %in% c("MT22S+MTdT"),]
REL3<-CFUs_Comp4[CFUs_Comp4$Set %in% c("MT48+MT14S"),]
REL<-rbind(REL1,REL2,REL3)

MTs <- Rel[!Rel$MT == "MTdT", ]
MTs <- MTs[!MTs$MT == "SBW25", ]
MTs <- MTs[!MTs$MT == "193", ]
MTs4 <- Rel4[!Rel4$MT == "MTdT", ]
MTs4 <- MTs4[!MTs4$MT == "MT14S", ]
MTs4 <- MTs4[!MTs4$MT == "MT187", ]
MTs5 <- Rel5[!Rel5$MT == "MTdT", ]
MTs6 <- CFUs_Comp6[CFUs_Comp6$MT %in% c("12PDE", "14PDE", "22PDE", "GCN4", "12PDEi", "14PDEi", "22PDEi", "GCN4i"),]

MT1<- REL[!REL$MT == "MTdT", ]
MT1<- MT1[!MT1$MT == "MT12S", ]
MT1<- MT1[!MT1$MT == "MT14S", ]

dT<-Rel[Rel$MT== "MTdT",]
dT2<-Rel[Rel$MT== "SBW25",]
dT3<-Rel[Rel$MT== "193",]
dT4<-Rel4[Rel4$MT== "MTdT",]
dT5<-Rel4[Rel4$MT== "MT14S",]
dT6<-Rel4[Rel4$MT== "MT187",]
dT<-rbind(dT,dT2,dT3)
dT7<-rbind(dT4,dT5,dT6)
dT8<-Rel5[Rel5$MT== "MTdT",]
dT9<-CFUs_Comp6[!CFUs_Comp6$MT %in% c("12PDE", "14PDE", "22PDE", "GCN4", "12PDEi", "14PDEi", "22PDEi", "GCN4i"),]

DT<-REL[REL$MT== "MTdT",]
DT2<-REL[REL$MT== "MT12S",]
DT3<-REL[REL$MT== "MT14S",]
DT<-rbind(DT,DT2,DT3)

Rel <- merge(MTs, dT, by = c("Replicate","Set","Worms","Sample"),all.x = T)
Rel4 <- merge(MTs4, dT7, by = c("Replicate","Set","Worms","Sample"),all.x = T)
REL <- merge(MT1, DT, by = c("Replicate","Set","Worms","Sample"),all.x = T)
Rel5 <- merge(MTs5, dT8, by = c("Replicate","Set","Worms","Sample"),all.x = T)
Rel6 <- merge(MTs6, dT9, by = c("Replicate","Set","Worms","Sample"),all.x = T)

Rel$'Rel' <- (Rel$'CFUs.x'/Rel$'CFUs.y')
Rel4$'Rel' <- (Rel4$'CFUs.x'/Rel4$'CFUs.y')
REL$'Rel' <- (REL$'CFUs.x'/REL$'CFUs.y')
Rel5$'Rel' <- (Rel5$'CFUs.x'/Rel5$'CFUs.y')
Rel6$'Rel' <- (Rel6$'CFUs.x'/Rel6$'CFUs.y')

PlotCFUs2<- function(Data){
  Data$Set <- factor(Data$Set,levels = c('MT48+MTdT','MT11+MTdT','MT12+MTdT','MT13+MTdT','MT14+MTdT','MT12S+MT48','MT48+MT14S','MT22S+MTdT','MT12*+MTdT','MT22+MTdT','SBW25+SBW25F*','193*+193','MT14*+MTdT','MT22*+MTdT','MT187F+MT187','12R+MTdT','14R+MTdT','22R+MTdT', 'GCN4+pJStrep', '12PDE+12pJN105', '14PDE+14pJN105','22PDE+22pJN105'),ordered = TRUE)
  Pop_Plot <- ggplot(Data,aes(x=Set,y=Rel))+
    geom_hline(yintercept = 1,linetype = "dashed",size=1)+
    geom_boxplot()+
    geom_point(aes(fill="black"),size=2, position= position_jitterdodge(jitter.width = 0, jitter.height = 0, dodge.width = 0,seed=NA))+
    scale_y_log10(limits=c(0.1, 1000),labels = scientific_10)+
    my_theme()+
    labs(y="Relative short-term persistence competition")+
    scale_x_discrete(labels=c(expression("MT48+MTdT" = "MYb11 dT vs. MYb11", "MT12+MTdT" ="MYb11 dT vs. MYb11"~italic("wspF*")^"EVO",
                                         "MT14+MTdT" = "MYb11 dT vs. MYb11"~italic("wspE*")^"EVO","SBW25+SBW25F*" = "SBW25 vs. SBW25"~italic(Delta*wspF)*"",
                                         "MT22+MTdT" = "MYb11 dT vs. MYb11"~italic("rph*")^"EVO","MT12*+MTdT"="MYb11 dT vs. MYb11"~italic("wspF*")*"",
                                         "MT14*+MTdT"="MYb11 dT vs. MYb11"~italic("wspE*")*"","MT22*+MTdT"="MYb11 dT vs. MYb11"~italic("rph*")*"",
                                         "193*+193"="MYb193 vs. MYb193"~italic(Delta*wspF)*"", "MT187F"="MYb187"~italic(Delta*wspF)*"","MT187F+MT187"="MYb187 vs. MYb187"~italic(Delta*wspF)*"",
                                         "MT12S+MT48" ="MYb11 dT vs. MYb11"~italic("wspF")*"","MT22S+MTdT" ="MYb11 dT vs. MYb11"~italic("rph")*"","MT48+MT14S" ="MYb11 dT vs. MYb11"~italic("wspE")*"")))
  return(Pop_Plot)
}


PlotCFUs3<- function(Data){
  Data$Set <- factor(Data$Set,levels = c('MT48+MTdT','MT11+MTdT','MT12+MTdT','MT13+MTdT','MT14+MTdT','MT12S+MT48','MT48+MT14S','MT22S+MTdT','MT12*+MTdT','MT22+MTdT','SBW25+SBW25F*','193*+193','MT14*+MTdT','MT22*+MTdT','MT187F+MT187','12R+MTdT','14R+MTdT','22R+MTdT', '12PDE+12pJN105', '14PDE+14pJN105','22PDE+22pJN105', 'GCN4+pJStrep'),ordered = TRUE)
  Pop_Plot <- ggplot(Data,aes(x=Set,y=Rel))+
    geom_hline(yintercept = 1,linetype = "dashed",size=1)+
    geom_boxplot()+
    geom_point(aes(fill="black"),size=2, position= position_jitterdodge(jitter.width = 0, jitter.height = 0, dodge.width = 0,seed=NA))+
    scale_y_log10(limits=c(0.01, 100),labels = scientific_10)+
    my_theme()+
    labs(y="Relative short-term persistence competition")+
    scale_x_discrete(labels=c(expression("12R+MTdT" = "MYb11 dT vs. MYb11"~italic("wspF*")~italic(Delta*wspR)*"","14R+MTdT" = "MYb11 dT vs. MYb11"~italic("wspE*")~italic(Delta*wspR)*"","22R+MTdT" = "MYb11 dT vs. MYb11"~italic("rph*")~italic(Delta*wspR)*"",
                                         "GCN4+pJStrep" = "MYb11 pJStrep vs. MYb11 GCN4wspR","12PDE+12pJN105" = "MYb11"~italic("wspF*")*"pJN105 vs. MYb11"~italic("wspF*")*"pPA2133",
                                         "14PDE+14pJN105" = "MYb11"~italic("wspE*")*"pJN105 vs. MYb11"~italic("wspE*")*"pPA2133","22PDE+22pJN105" = "MYb11"~italic("rph*")*"pJN105 vs. MYb11"~italic("rph*")*"pPA2133")))
    return(Pop_Plot)
}

PlotCFUs2(Rel[Rel$Set %in% c("MT12+MTdT","MT14+MTdT","MT22+MTdT"),])
PlotCFUs2(Rel[Rel$Set %in% c("MT12*+MTdT","MT14*+MTdT","MT22*+MTdT"),])
PlotCFUs2(Rel[Rel$Set %in% c("SBW25+SBW25F*"),])
PlotCFUs2(Rel[Rel$Set %in% c("193*+193"),])
PlotCFUs2(Rel4[Rel4$Set %in% c("MT187F+MT187"),])
PlotCFUs2(REL[REL$Set %in% c("MT12S+MT48","MT48+MT14S","MT22S+MTdT"),])
PlotCFUs3(Rel5[Rel5$Set %in% c("12R+MTdT","14R+MTdT","22R+MTdT"),])
A<-PlotCFUs3(Rel6[Rel6$Set %in% c("12PDE+12pJN105"),])
B<-PlotCFUs3(Rel6[Rel6$Set %in% c("14PDE+14pJN105"),])
C<-PlotCFUs3(Rel6[Rel6$Set %in% c("22PDE+22pJN105"),])
D<-PlotCFUs3(Rel6[Rel6$Set %in% c("GCN4+pJStrep"),])
PlotCFUs3(Rel6)

ggarrange(A,B,C,D, ncol=4,heights=c(1),widths = c(1,1,1,1),align= "hv")

#####9. Statistic ######

#Quartet competition

Test1 <- function(data,Set){
  # Subset data
  d <- data
  d <- subset(d, select = c(Replicate, Set,MT, LOG_CFU))
  d$MT <- droplevels(d$MT)
  d$MT <- factor(d$MT, levels = c("MTdT", levels(d$MT)[levels(d$MT)!= "MTdT"]))
  
  # Prepare output file name
  fileName <- paste("stats out/Duo_Quartet/",Set,"_Worm",".txt",sep="")
  
  # Run GLM
  m <-lmer(LOG_CFU ~ MT + (1|Replicate),data=d)
  
  sink(fileName)
  print(as.data.frame(coef(summary(m)))) # lmer output
  print("")
  print(summary(glht(m, mcp(MT="Tukey")))) # Post-hoc output: Tukey
  print("")
  print(ml<-glht(m, mcp(MT="Dunnett"))) # Post-doc output: Dunnett
  print(summary(ml,adjusted("fdr")))
  sink()
  
  return("Regression analysis done and saved.")
}

Test1(CFUs_Comp[CFUs_Comp$Set =="MT11+MT13+MT14+MTdT",],Set="MT11+MT13+MT14+MTdT" )
Test1(CFUs_Comp[CFUs_Comp$Set =="MT11+MT12+MT13+MTdT",],Set="MT11+MT12+MT13+MTdT" )



# Relative Competition#

Test2 <- function(data,Set){
  # Subset data
  d <- data
  d1<-read.csv("test.csv",TRUE,",") #Adds 1 (anc/anc) as a relative competition index for comparisons of mutants with the ancestor 
  names(d1)[1]<-"Replicate"
  d<-rbind(d,d1)
  
  d1<-log10(d[11])
  d<-cbind(d,d1)
  names(d)[12]<-"LOG_Rel"
  
  d<- d[!is.na(d$Rel),]
  d <- subset(d, select = c(Replicate, Set, LOG_Rel))
  d$Set <- droplevels(d$Set)
  d$Set <- factor(d$Set, levels = c("MTdT", levels(d$Set)[levels(d$Set)!= "MTdT"]))
  
  # Prepare output file name
  fileName <- paste("stats out/Rel Comp/",Set,"_Worm",".txt",sep="")
  
  # Run GLM
  m <-lmer(LOG_Rel ~ Set + (1|Replicate),data=d)
  
  sink(fileName)
  print(as.data.frame(coef(summary(m)))) # lmer output
  print("")
  print(summary(glht(m, mcp(Set="Tukey")))) # Post-hoc output: Tukey
  print("")
  print(ml<-glht(m, mcp(Set="Dunnett")))# Post-doc output: Dunnett
  print(summary(ml,adjusted("fdr")))
  sink()
  
  return("Regression analysis done and saved.")
}

Test2(Rel[Rel$Set %in% c("MT48+MTdT","MT12+MTdT","MT14+MTdT","MT22+MTdT"),],Set="Evo")
Test2(Rel[Rel$Set %in% c("MT48+MTdT","MT12*+MTdT","MT14*+MTdT","MT22*+MTdT"),],Set="Con")
Test2(Rel[Rel$Set %in% c("SBW25+SBW25F*"),],Set="SBW25")
Test2(Rel4[Rel4$Set %in% c("MT187F+MT187"),],Set="MT187")
Test2(Rel[Rel$Set %in% c("193*+193"),],Set="MT193")
Test2(REL,Set="Rescue")
Test2(Rel5[Rel5$Set %in% c("12R+MTdT","14R+MTdT","22R+MTdT"),],Set="wspR")
Test2(Rel6,Set="PDE")

