#Load libraries
#library(phyloseq)  #must load the devtools version
library("devtools")
install_github("phyloseq", "joey711")
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

# lets look at only samples (removing blanks and mock and samples that didn't amplify)
pruned <- prune_taxa(taxa_sums(merge) > 0, merge)
bact_samples <-subset_taxa(pruned, Kingdom == "Bacteria")
#nrow(otu_table(pruned)) - nrow(otu_table(bact_samples))  ## Check how many rows were not bacteria
bact_samples2 <-subset_taxa(bact_samples, Class != "Chloroplast")
#nrow(otu_table(bact_samples)) - nrow(otu_table(bact_samples2))  # Check how many rows were chloroplasts

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


## Not all samples have a duplicate:  Let's see which samples do/don't
data <- data.frame(sample_data(good_samples))
nrow(subset(data, lakenames == "Baker")) # BAKER has all replicates!
nrow(subset(data, lakenames == "Bristol")) # BRISTOL has all replicates!
nrow(subset(data, lakenames == "Payne")) # PAYNE has all replicates!
nrow(subset(data, lakenames == "Sherman")) # SHERMAN has all replicates!
nrow(subset(data, lakenames == "Wintergreen")) # WINTERGREEN has all replicates!

# MISSING REPLICATES
# nrow(subset(data, lakenames == "Baseline")) # MISSING A REPLICATE
# row.names(subset(data, lakenames == "Baseline"))  # We are missing "BASE23um"
# nrow(subset(data, lakenames == "Bassett")) 
# row.names(subset(data, lakenames == "Bassett"))  # We are missing "BSTE2"
# nrow(subset(data, lakenames == "Gull"))
# row.names(subset(data, lakenames == "Gull"))  # We are missing "GULH1"
# nrow(subset(data, lakenames == "Lee"))
# row.names(subset(data, lakenames == "Lee"))  # We are missing "LEEE1"
# nrow(subset(data, lakenames == "LittleLong")) #MISSING 3 samples!!!
# row.names(subset(data, lakenames == "LittleLong"))  # We are missing  "LONE13um","LONE23um", "LONH23um"
# nrow(subset(data, lakenames == "Sixteen")) # Missing 2 samples!!
# row.names(subset(data, lakenames == "Sixteen"))  # We are missing  "SIXE1","SIXH2"

missing <- c("BASE23um", "BSTE2", "GULH1", "LEEE1", "LONE13um","LONE23um", "LONH23um", "SIXE1","SIXH2")

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



### Time to merge the replicate samples!!!  
##  In our metadata there's a column called dups (aka duplicates) that we can use to merge our samples
merged_samps <- merge_samples(good_samples, "dups", fun = "mean") 
merged_samps <- prune_taxa(taxa_sums(merged_samps) > 0, good_merged)

ggplot(data.frame(sum=sample_sums(good_merged)),aes(sum, fill = s)) + ylab("Number of Sequences per Sample") +
  geom_histogram(colour="black",fill="red") + ggtitle("Merged Duplicates:  McMurdie & Holme's Scaled Reads") + xlab("Total Sequences")

# mean, max and min of sample read counts
mins<-min(sample_sums(good_merged))
paste(c("The minimum sample read count is",mins))
means<-mean(sample_sums(good_merged))
paste(c("The mean sample read count is",means))
maxs<-max(sample_sums(good_merged))
paste(c("The max sample read count is",maxs))

missing_dups <- c("BASE3um", "BSTE", "GULH", "LEEE", "LONE3um","LONE3um", "LONH3um", "SIXE","SIXH")









dups_metadata <- sample_data(good_merged)
dups_metadata$names <- row.names(dups_metadata)
dups_metadata$dups <- row.names(dups_metadata)
good_metadata <-makeCategories_dups(dups_metadata)







