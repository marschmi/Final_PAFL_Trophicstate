###  This is a file of functions
####  Functions in this file include:
# 1. multiplot:  plot multiple ggplots in one 
# 2. makeCategories:  add categorical metadata to my dataframes
# 3. makeCategories_dups:  add categorical metadata to my dataframe of MERGED DUPLICATES
# 4. deSEQ:  Calculate the log2-fold ratio 
# 5. plot_deSEQ - plot the OTU level log2-fold ratio
# 6. plot_phylum_deSEQ - plot the PHYLUM level log2-fold ratio
# 7. scale_reads - McMurdie & Holmes scaling to dataset for uneven sequencing depth across samples
# 8. nacols - look for NAs in the dataframe
# 9. merge_samples_mean - mean the counts between two samples 



#################################################################################### 1
#################################################################################### 1
# Multiple plot function is from the R-Cookbook website 
# http://www.cookbook-r.com/Graphs/Multiple_graphs_on_one_page_(ggplot2)/
#
# ggplot objects can be passed in ..., or to plotlist (as a list of ggplot objects)
# - cols:   Number of columns in layout
# - layout: A matrix specifying the layout. If present, 'cols' is ignored.
#
# If the layout is something like matrix(c(1,2,3,3), nrow=2, byrow=TRUE),
# then plot 1 will go in the upper left, 2 will go in the upper right, and
# 3 will go all the way across the bottom.
#
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  require(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}
#################################################################################### 1
#################################################################################### 1




#################################################################################### 2 
#################################################################################### 2
## This function adds the categorical metadata to a dataframe based on the sample name
      # IMPORTANT!!!  The dataframe MUST have column named "names"
makeCategories <- function(dataframe){ 
  dataframe$lakenames<-substr(dataframe$names,1,3) # Create a new row called "lakenames" with first 3 letters of string
  dataframe$limnion <- substr(dataframe$names, 4, 4) # Create a column called limnon with top or bottom 
  dataframe$filter <- substr(dataframe$names, 6, 8) # Create a column called "filter" with pore size
  dataframe$limnion <- as.character(dataframe$limnion)  
  dataframe$limnion[dataframe$limnion == "E"] <- "Top"
  dataframe$limnion[dataframe$limnion == "H"] <- "Bottom"
  dataframe$filter <- as.character(dataframe$filter)
  dataframe$filter[dataframe$filter == "3um"] <- "Particle"
  dataframe$filter[dataframe$filter == ""] <- "Free"
  dataframe$lakenames <- as.character(dataframe$lakenames)
  dataframe$lakenames[dataframe$lakenames == "WIN"] <- "Wintergreen"
  dataframe$lakenames[dataframe$lakenames == "SIX"] <- "Sixteen"
  dataframe$lakenames[dataframe$lakenames == "SHE"] <- "Sherman"
  dataframe$lakenames[dataframe$lakenames == "PAY"] <- "Payne"
  dataframe$lakenames[dataframe$lakenames == "LON"] <- "LittleLong"
  dataframe$lakenames[dataframe$lakenames == "LEE"] <- "Lee"
  dataframe$lakenames[dataframe$lakenames == "GUL"] <- "Gull"
  dataframe$lakenames[dataframe$lakenames == "BRI"] <- "Bristol"
  dataframe$lakenames[dataframe$lakenames == "BAK"] <- "Baker"
  dataframe$lakenames[dataframe$lakenames == "BAS"] <- "Baseline"
  dataframe$lakenames[dataframe$lakenames == "BST"] <- "Bassett"
  ## Fix Sherman Lake
  dataframe$limnion[dataframe$lakenames == "Sherman" & dataframe$limnion == "Epilimnion"] <- "Mixed"
  dataframe$limnion[dataframe$lakenames == "Sherman" & dataframe$limnion == "Hypolimnion"] <- "Mixed"
  ## TO MAKE QUADRANT
  for(i in 1:length(dataframe$limnion)){
    dataframe$quadrant[i]<-paste(as.character(dataframe$filter[i]),as.character(dataframe$limnion[i]))}
  #Add Trophic State column by using the name of the lake
  dataframe <- data.table(dataframe)
  library(data.table)
  dataframe[, trophicstate := ifelse(lakenames %in% c("Wintergreen", "Baker", "Baseline"), "Eutrophic",
                                     ifelse(lakenames %in% c("Bassett", "Bristol", "Payne", "Sherman"), "Mesotrophic",
                                            #ifelse(lakenames %in% c("Sherman"), "Mixed",
                                                   ifelse(lakenames %in% c("Gull", "Sixteen", "LittleLong", "Lee"), "Oligotrophic", NA)))]
  dataframe$ProdLevel <- dataframe$trophicstate
  dataframe$ProdLevel <- as.character(dataframe$ProdLevel)
  dataframe$ProdLevel[dataframe$ProdLevel == "Eutrophic"] <- "Productive"
  dataframe$ProdLevel[dataframe$ProdLevel == "Mesotrophic"] <- "Productive"
  dataframe$ProdLevel[dataframe$ProdLevel == "Oligotrophic"] <- "Unproductive"
  dataframe$ProdLevel[dataframe$ProdLevel == "Mixed"] <- "Mixed"
  dataframe$ProdLevel <- as.factor(dataframe$ProdLevel)
  dataframe$ProdLevel <-factor(dataframe$ProdLevel,levels=c("Productive","Unproductive","Mixed"))
}
#################################################################################### 2
#################################################################################### 2




#################################################################################### 3
#################################################################################### 3
####  This function applies the same info as the above function -> but with combined duplicates!
######
makeCategories_dups <- function(dataframe){ # IMPORTANT!!!  dataframe MUST have column named "names"
  dataframe$lakenames<-substr(dataframe$names,1,3) # Create a new row called "lakenames" with first 3 letters of string
  dataframe$limnion <- substr(dataframe$names, 4, 4) # Create a column called limnon with top or bottom 
  dataframe$filter <- substr(dataframe$names, 5, 7) # Create a column called "filter" with pore size
  dataframe$limnion <- as.character(dataframe$limnion)  
  dataframe$limnion[dataframe$limnion == "E"] <- "Epilimnion"
  dataframe$limnion[dataframe$limnion == "H"] <- "Hypolimnion"
  dataframe$filter <- as.character(dataframe$filter)
  dataframe$filter[dataframe$filter == "3um"] <- "Particle"
  dataframe$filter[dataframe$filter == ""] <- "Free"
  dataframe$lakenames <- as.character(dataframe$lakenames)
  dataframe$lakenames[dataframe$lakenames == "WIN"] <- "Wintergreen"
  dataframe$lakenames[dataframe$lakenames == "SIX"] <- "Sixteen"
  dataframe$lakenames[dataframe$lakenames == "SHE"] <- "Sherman"
  dataframe$lakenames[dataframe$lakenames == "PAY"] <- "Payne"
  dataframe$lakenames[dataframe$lakenames == "LON"] <- "LittleLong"
  dataframe$lakenames[dataframe$lakenames == "LEE"] <- "Lee"
  dataframe$lakenames[dataframe$lakenames == "GUL"] <- "Gull"
  dataframe$lakenames[dataframe$lakenames == "BRI"] <- "Bristol"
  dataframe$lakenames[dataframe$lakenames == "BAK"] <- "Baker"
  dataframe$lakenames[dataframe$lakenames == "BAS"] <- "Baseline"
  dataframe$lakenames[dataframe$lakenames == "BST"] <- "Bassett"
  ## Fix Sherman Lake
  dataframe$limnion[dataframe$lakenames == "Sherman" & dataframe$limnion == "Epilimnion"] <- "Mixed"
  dataframe$limnion[dataframe$lakenames == "Sherman" & dataframe$limnion == "Hypolimnion"] <- "Mixed"
  ## TO MAKE QUADRANT
  for(i in 1:length(dataframe$limnion)){
    dataframe$quadrant[i]<-paste(as.character(dataframe$filter[i]),as.character(dataframe$limnion[i]))}
  #Add Trophic State column by using the name of the lake
  dataframe <- data.table(dataframe)
  library(data.table)
  dataframe[, trophicstate := ifelse(lakenames %in% c("Wintergreen", "Baker", "Baseline"), "Eutrophic",
                                     ifelse(lakenames %in% c("Bassett", "Bristol", "Payne", "Sherman"), "Mesotrophic",
                                            #ifelse(lakenames %in% c("Sherman"), "Mixed",
                                                   ifelse(lakenames %in% c("Gull", "Sixteen", "LittleLong", "Lee"), "Oligotrophic", NA)))]
  #dataframe$ProdLevel <- as.character(dataframe$trophicstate)
  #dataframe$ProdLevel[dataframe$trophicstate == "Eutrophic"] <- "Productive"
  #dataframe$ProdLevel[dataframe$trophicstate == "Mesotrophic"] <- "Productive"
  #dataframe$ProdLevel[dataframe$trophicstate == "Oligotrophic"] <- "Unproductive"
  #dataframe$ProdLevel[dataframe$ProdLevel == "Mixed"] <- "Productive"
  #dataframe$ProdLevel <- as.factor(dataframe$ProdLevel)
  #dataframe$ProdLevel <-factor(dataframe$ProdLevel,levels=c("Productive","Unproductive"))
}
#################################################################################### 3
#################################################################################### 3




#################################################################################### 4 + 5
#################################################################################### 4 + 5
#OTU OTU OTU OTU OTU OTU OTU OTU OTU OTU OTU OTU OTU OTU OTU OTU OTU OTU OTU OTU OTU 
## I wrote these 2 functions off of the following tutorial from the Phyloseq GitHub page:
#http://joey711.github.io/phyloseq-extensions/DESeq2.html

deSEQ <- function(data, valuetest){
  data_pruned=prune_taxa(taxa_sums(data)>0,data)
  de_data = phyloseq_to_deseq2(data_pruned, valuetest)
  de_data2 = DESeq(de_data, test="Wald", fitType="parametric")
  res_data = results(de_data2, cooksCutoff = FALSE)
  alpha = 0.01
  sig_data = res_data[which(res_data$padj < alpha), ]
  sigtab_sherm = cbind(as(sig_data, "data.frame"), as(tax_table(data_pruned)[rownames(sig_data), ], "matrix"))
} 


plot_deSEQ <- function(deSEQdata, title){
  y = tapply(deSEQdata$log2FoldChange, deSEQdata$Genus, function(x) max(x))
  y = sort(y, TRUE)
  deSEQdata$Genus = factor(as.character(deSEQdata$Genus), levels=names(y))
  ggplot(deSEQdata, aes(x=Genus, y=log2FoldChange, color=Phylum)) + 
    geom_point(size=6) +  ggtitle(title) +  
    scale_color_manual(values = phylum.colors,name="Phylum") +
    scale_y_continuous(breaks=seq(-15,15,1)) + 
    theme(axis.title.x = element_text(face="bold", size=16),
          axis.text.x = element_text(angle=30, colour = "black", vjust=1, hjust = 1, size=12),
          axis.text.y = element_text(colour = "black", size=16),
          axis.title.y = element_text(face="bold", size=16),
          plot.title = element_text(face="bold", size = 20),
          legend.title = element_text(size=12, face="bold"),
          legend.text = element_text(size = 12),
          strip.background = element_rect(colour="black"),
          legend.position="right")
}

#################################################################################### 4 + 5
#################################################################################### 4 + 5




#################################################################################### 7
#################################################################################### 7
#PHYLUM PHYLUM PHYLUM PHYLUM PHYLUM PHYLUM PHYLUM PHYLUM PHYLUM PHYLUM PHYLUM PHYLUM 
plot_phylum_deSEQ <- function(deSEQdata, title){
  x = tapply(deSEQdata$log2FoldChange, deSEQdata$Phylum, function(x) max(x))
  x = sort(x, TRUE)
  deSEQdata$Phylum = factor(as.character(deSEQdata$Phylum), levels=names(x))
  ggplot(deSEQdata, aes(x=Phylum, y=log2FoldChange, color=Phylum)) + 
    geom_point(size=6) +  ggtitle(title) +  
    scale_y_continuous(breaks=seq(-15,15,1)) + 
    scale_color_manual(values = phylum.colors,name="Phylum") +
    theme(axis.title.x = element_text(face="bold", size=16),
          axis.text.x = element_text(angle=30, colour = "black", vjust=1, hjust = 1, size=12),
          axis.text.y = element_text(colour = "black", size=16),
          axis.title.y = element_text(face="bold", size=16),
          plot.title = element_text(face="bold", size = 20),
          legend.title = element_text(size=12, face="bold"),
          legend.text = element_text(size = 12),
          strip.background = element_rect(colour="black"),
          legend.position="right")
}
#################################################################################### 7
#################################################################################### 7



#################################################################################### 8
#################################################################################### 8
## This function was written by Michelle Berry and will scale the reads to the minimum number of sample reads in the dataset
#scale_reads_floor <- function(physeq,n){  #n is the size you want to scale to. Usually, do min(sample_sums(physeq))
#  physeq.scale <- transform_sample_counts(physeq, function(x) {(n*x/sum(x))})
#  otu_table(physeq.scale) = floor(otu_table(physeq.scale))
#  physeq.scale = prune_taxa(taxa_sums(physeq.scale)>0,physeq.scale)
#  return(physeq.scale)
#}



#scale_reads_round <- function(physeq,n){  #n is the size you want to scale to. Usually, do min(sample_sums(physeq))
#  physeq.scale <- transform_sample_counts(physeq, function(x) {(n*x/sum(x))})
#  otu_table(physeq.scale) = round(otu_table(physeq.scale))
#  physeq.scale = prune_taxa(taxa_sums(physeq.scale)>0,physeq.scale)
#  return(physeq.scale)
#}
#################################################################################### 8
#################################################################################### 8




#################################################################################### 9
#################################################################################### 9
#Function to find NAs in a data frame
nacols <- function(df) {
  colnames(df)[unlist(lapply(df, function(x) any(is.na(x))))]
}
#################################################################################### 9
#################################################################################### 9


# This function was written by Michelle Berry 
# Better rounding function than R's base round
matround <- function(x){trunc(x+0.5)}

#################################################################################### 10
#################################################################################### 10
# This function was written by Michelle Berry and edited by Marian Schmidt and can be found at https://github.com/joey711/phyloseq/issues/465
merge_samples_mean <- function(physeq, group, round = "none"){
  # Calculate the number of samples in each group
  group_sums <- as.matrix(table(sample_data(physeq)[ ,group]))[,1]
  
  # Merge samples by summing
  merged <- merge_samples(physeq, group)
  
  # Divide summed OTU counts by number of samples in each group to get mean
  # Calculation is done while taxa are columns, but then transposed at the end
  x <- as.matrix(otu_table(merged))
  if(taxa_are_rows(merged)){ x<-t(x) }
  
  print(round)
  
  # Pick the rounding function
  if (round == "none"){
    out <- t(x/group_sums)
  } else if (round == "round"){
    out <- t(round(x/group_sums))
  } else if (round == "matround"){
    out <- t(matround(x/group_sums))
  } else if (round == "floor"){
    out <- t(floor(x/group_sums))
  }
  
  # Return new phyloseq object with taxa as rows
  out <- otu_table(out, taxa_are_rows = TRUE)
  otu_table(merged) <- out
  return(merged)
}
#################################################################################### 11
#################################################################################### 11


# Scales reads by 
# 1) taking proportions
# 2) multiplying by a given library size of n
# 3) rounding 
# Default for n is the minimum sample size in your library
# Default for round is floor
scale_reads <- function(physeq, n = min(sample_sums(physeq)), round = "floor") {
  
  # transform counts to n
  physeq.scale <- transform_sample_counts(physeq, 
                                          function(x) {(n * x/sum(x))}
  )
  
  # Pick the rounding functions
  if (round == "floor"){
    otu_table(physeq.scale) <- floor(otu_table(physeq.scale))
  } else if (round == "round"){
    otu_table(physeq.scale) <- round(otu_table(physeq.scale))
  } else if (round == "matround"){
    otu_table(physeq.scale) <- matround(otu_table(physeq.scale))
  }
  
  # Prune taxa and return new phyloseq object
  physeq.scale <- prune_taxa(taxa_sums(physeq.scale) > 0, physeq.scale)
  return(physeq.scale)
}




