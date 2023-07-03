### ===    EE MYB bacteria    === ###
###          Nancy Obeng          ###
#####################################

#### Set dependencies ####
# Set working directory to source location
currentDirectory <- dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(currentDirectory)

library(plyr); library(ggplot2); library(gridExtra)
source("my_theme.R")

#### === Load data === ####
evoBS <- read.csv("expEvo_MYb11_withinWorms.csv")
evoMS <- read.csv("expEvo_MYb11_onPlates.csv")

#### === Plotting == ####
# Proper axis labels
fancy_scientific <- function(l) {
  # turn in to character string in scientific notation
  l <- format(l, scientific = TRUE)
  # quote the part before the exponent to keep all the digits
  l <- gsub("e\\+","e",l)
  # turn the 'e+' into plotmath format
  l <- gsub("e", "%*%10^", l)
  # return this as an expression
  parse(text=l)}

colsTreat <- c(rgb(0, 158,115, maxColorValue = 255), rgb(230, 159,0, maxColorValue = 255))

mM <- ddply(evoMS, .(strain, cycle, treatment), summarize, mcfu = mean(cfu))
mB <- ddply(evoBS, .(strain, cycle, treatment), summarize, mcfu = mean(cfu))

a <- ggplot(evoBS, aes(x = cycle, y = cfu, group = popID, col = strain))+
  geom_line(alpha = 0.5, lty =2, size =1.5)+
  geom_line(data = mB, aes(x = cycle, y = mcfu, group = strain), size = 1.5)+
  scale_y_log10(labels = fancy_scientific)+
  xlab("Cycle") + ylab("CFU/ worm pop. (3.5 d)") +
  ggtitle("A")+
  scale_color_manual(values = colsTreat[1])+
  my_theme(evoBS)+
  theme(plot.title = element_text(hjust = 0, size = 25)); a

b <- ggplot(evoMS, aes(x = cycle, y = cfu, group = popID, col = strain))+
  geom_line(alpha = 0.5, lty =2, size =1.5)+
  geom_line(data = mM, aes(x = cycle, y = mcfu, group = strain), size = 1.5)+
  scale_y_log10(labels = fancy_scientific)+
  xlab("Cycle") + ylab("CFU/ plate (7 d)") +
  ggtitle("B")+
  scale_color_manual(values = colsTreat[2])+
  my_theme(evoMS)+
  theme(plot.title = element_text(hjust = 0, size = 25)); b

#Save plot to file
fileExtensions <- c(".svg", ".png")
for(i in 1:length(fileExtensions))
ggsave(paste("Fitness during EE_rand.selected pops_MYb11", fileExtensions[i],sep =""),
       grid.arrange(a, b, nrow=1), 
       dpi = 300, height = 4, width = 9)
