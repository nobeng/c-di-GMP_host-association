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
  theme(axis.text.x=element_text(size=14,color=col.text)) +
  theme(axis.text.y=element_text(size=14,color=col.text)) +
  theme(axis.title.x=element_text(size=14,color=col.text, vjust=0)) +
  theme(axis.title.y=element_text(size=14,color=col.text, vjust=1.25)) +
  theme(strip.text=element_text(size=14,color="grey30"))  }