#Load libraries
library(phyloseq)
library(ggplot2)
library(ape)
library(vegan)
library(plyr)
library(scales)
library(grid)
library(reshape2)
library(data.table)
library(sciplot)
library("DESeq2")
library(tidyr)
library(dplyr)


### Set the Working Directory
setwd("~/Final_PAFL_Trophicstate")


### Source written functions in the file Functions.PAFL.R
source("Functions_PAFL.R")


### Set our theme for our plots down the road
set_ggtheme <- theme_set(theme_bw() + 
      theme(plot.title = element_text(face="bold", size = 20),  #Set the plot title
        strip.text.x = element_text(size=12, face="bold"),  #Set the facet titles on x-axis 
        strip.text.y = element_text(size=12, face="bold"),  #Set the facet titles on x-axis 
        strip.background = element_blank(),  #Set the facet background to no background
        axis.title.x = element_text(face="bold", size=16),  #Set the x-axis title
        axis.title.y = element_text(face="bold", size=16),  #Set the y-axis title
        axis.text.x = element_text(colour = "black", size=14),  #Set the x-axis labels
        axis.text.y = element_text(colour = "black", size=14),  #Set the y-axis labels
        legend.title = element_text(size=14, face="bold"),  #Set the legend title 
        legend.text = element_text(size = 14),  #Set the legend text
        legend.position="right"))  #Default the legend position to the right


### Data import
#sharedfile = "mothur.normalized.shared"
sharedfile = "~/Final_PAFL_Trophicstate/raw_data/stability.trim.contigs.good.unique.good.filter.precluster.pick.pick.pick.an.unique_list.shared"
taxfile = "~/Final_PAFL_Trophicstate/raw_data/stability.trim.contigs.good.unique.good.filter.precluster.pick.pick.pick.an.unique_list.0.03.cons.taxonomy"
mothurdata = import_mothur(mothur_shared_file = sharedfile ,mothur_constaxonomy_file = taxfile )


# We need to change the taxonomy names
tax_table(mothurdata)<-cbind(tax_table(mothurdata),row.names(tax_table(mothurdata)))
colnames(tax_table(mothurdata)) <- c("Kingdom","Phylum","Class","Order","Family","Genus","Species")

#We need to change some of the filenames to match our metadata
samplenames <- sample_names(mothurdata)
names2 <- gsub("GULL", "GUL", samplenames)
names3 <- gsub("LONG", "LON", names2)
sample_names(mothurdata) <- names3


# Import metadata file and merge with mothurdata object
mapfile <- "~/Final_PAFL_Trophicstate/raw_data/metadata"
map <- import_qiime_sample_data(mapfile)
merge <- merge_phyloseq(mothurdata,map)

###### Analysis ##########
# lets look at only samples (removing blanks and mock and samples that didn't amplify)
pruned <- prune_taxa(taxa_sums(merge) > 0, merge)
bact_samples <-subset_taxa(pruned, Kingdom == "Bacteria")

#Summary of my good_samples object
bact_samples

## Raw Sample read counts
# Histogram of RAW sample read counts
ggplot(data.frame(sum=sample_sums(bact_samples)),aes(sum)) + ylab("Number of Sequences per Sample") +
  geom_histogram(colour="black",fill="dodgerblue2")  + xlab("Total Sequences") + ggtitle("Raw Sample Read Counts")

# mean, max and min of sample read counts
mins<-min(sample_sums(bact_samples))
paste(c("The minimum sample read count is",mins))
means<-mean(sample_sums(bact_samples))
paste(c("The mean sample read count is",means))
maxs<-max(sample_sums(bact_samples))
paste(c("The max sample read count is",maxs))

### Now let's scale our read counts
good_samples <- scale_reads(physeq = bact_samples, n = min(sample_sums(bact_samples)))

# Histogram of SCALED sample read counts
ggplot(data.frame(sum=sample_sums(good_samples)),aes(sum, fill = s)) + ylab("Number of Sequences per Sample") +
  geom_histogram(colour="black",fill="mediumorchid3") + ggtitle("McMurdie & Holme's Scaled Reads") + xlab("Total Sequences")

# mean, max and min of sample read counts
mins<-min(sample_sums(good_samples))
paste(c("The minimum sample read count is",mins))
means<-mean(sample_sums(good_samples))
paste(c("The mean sample read count is",means))
maxs<-max(sample_sums(good_samples))
paste(c("The max sample read count is",maxs))





