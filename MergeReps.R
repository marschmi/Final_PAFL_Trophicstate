#Load libraries
#library(phyloseq)  #must load the devtools version
library("devtools")
install_github("phyloseq", "joey711")
library(ggplot2)
library(ape)
library(vegan)


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
sharedfile <- "~/Final_PAFL_Trophicstate/raw_data/stability.trim.contigs.good.unique.good.filter.precluster.pick.pick.pick.an.unique_list.shared"
taxfile <- "~/Final_PAFL_Trophicstate/raw_data/stability.trim.contigs.good.unique.good.filter.precluster.pick.pick.pick.an.unique_list.0.03.cons.taxonomy"
mothurdata <- import_mothur(mothur_shared_file = sharedfile ,mothur_constaxonomy_file = taxfile )


# We need to change the taxonomy names
tax_table(mothurdata) <- cbind(tax_table(mothurdata),row.names(tax_table(mothurdata)))
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

### We will scale our read counts, combine (sum) our reads, take the average, and then re-scale.
good_samples <- scale_reads(physeq = bact_samples, n = min(sample_sums(bact_samples)))

## Not all samples have a duplicate:  Let's see which samples do/don't
#data <- data.frame(sample_data(good_samples))
#nrow(subset(data, lakenames == "Baker")) # BAKER has all replicates!
#nrow(subset(data, lakenames == "Bristol")) # BRISTOL has all replicates!
#nrow(subset(data, lakenames == "Payne")) # PAYNE has all replicates!
#nrow(subset(data, lakenames == "Sherman")) # SHERMAN has all replicates!
#nrow(subset(data, lakenames == "Wintergreen")) # WINTERGREEN has all replicates!

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
noreps <- c("BASE13um", "BSTE1", "GULH2", "LEEE2", "LONH13um", "SIXE2","SIXH1")

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


### Time to merge the replicate samples
##  In our metadata there's a column called dups (aka duplicates) that we can use to merge our samples
good_merged <- merge_samples(good_samples, "dups", fun = "mean") #  
good_merged <- prune_taxa(taxa_sums(good_merged) > 0, good_merged)

ggplot(data.frame(sum=sample_sums(good_merged)),aes(sum)) + ylab("Number of Sequences per Sample") +
  geom_histogram(colour="black",fill="red") + ggtitle("Merged Duplicates: \n McMurdie & Holme's Scaled Reads") + 
  xlab("Total Sequences")
## As we can see, the fun = "mean" DOESN'T WORK.  Instead, merge_samples is only summing the sample counts!

# mean, max and min of sample read counts
mins<-min(sample_sums(good_merged))
paste(c("The minimum sample read count is",mins))
means<-mean(sample_sums(good_merged))
paste(c("The mean sample read count is",means))
maxs<-max(sample_sums(good_merged))
paste(c("The max sample read count is",maxs))


### AVERAGING AND SCALING THE READ COUNTS 
otus <- data.frame(otu_table(merged_samps))
missing_dups <- c("BASE3um", "BSTE", "GULH", "LEEE", "LONE3um","LONE3um", "LONH3um", "SIXE","SIXH")
otus$names <- row.names(otus)
nodups <- otus[otus$names %in% missing_dups, ] # collect only the samples that DO NOT have duplicates
dups_only <- otus[!otus$names %in% missing_dups, ] # collect samples that DO have duplicates
dups_only$names = NULL
nodups$names = NULL
rowSums(dups_only)
dups_half <- round(dups_only/2)  #Here we're taking the average.
norm_all <- rbind(dups_half, nodups)
sums <- data.frame(rowSums(norm_all))
colnames(sums) <- c("Totals")

sums_range <- max(sums) -min(sums)
paste(c("The range of sample read counts when scaled, merged (summed), averaged and then scaled is",sums_range))

ggplot(sums, aes(Totals)) + ylab("Number of Sequences per Sample") +
  geom_histogram(colour="black",fill="gold") + xlab("Total Sequences") +
  ggtitle(expression(atop("Averaged Merged Duplicates:\n  McMurdie & Holme's Scaled Reads", atop("Range = 508"), ""))) 
 

data_dup <- read.table("~/Final_PAFL_Trophicstate/raw_data/duplicate_metadata", header = TRUE, sep = "\t")
data_dup <- sample_data(data_dup)
tnorm_all <- t(norm_all)
otu_manual <- otu_table(tnorm_all, taxa_are_rows = TRUE)
taxo <- tax_table(good_merged)

manual_merged <- merge_phyloseq(taxo, otu_manual, data_dup)
manual_merged <- prune_taxa(taxa_sums(manual_merged) > 0, manual_merged)



################################################################
### COMBINING BEFORE MCMURDIE HOLMES SCALING WITH OUR RAW READS
noscale_merge <- merge_samples(bact_samples, "dups", fun = "mean") 
noscale_merge <- prune_taxa(taxa_sums(noscale_merge) > 0, noscale_merge)
range_nonscale <- max(sample_sums(noscale_merge)) - min(sample_sums(noscale_merge))

##  Sample read counts with merging samples WITHOUT SCALING!!!! 
# Histogram of RAW sample read counts
ggplot(data.frame(sum=sample_sums(noscale_merge)),aes(sum)) + ylab("Number of Sequences per Sample") +
  geom_histogram(colour="black",fill="violetred")  + xlab("Total Sequences") + 
  ggtitle(expression(atop("Merged Raw Reads", atop("Range = 54,691"), ""))) 

# mean, max and min of sample read counts
mins<-min(sample_sums(noscale_merge))
paste(c("The minimum sample read count is",mins))
means<-mean(sample_sums(noscale_merge))
paste(c("The mean sample read count is",means))
maxs<-max(sample_sums(noscale_merge))
paste(c("The max sample read count is",maxs))

### Now let's scale our read counts from our MERGED RAW data
scaled_merged <- scale_reads(physeq = noscale_merge, n = min(sample_sums(noscale_merge)))

##  Sample read counts with merging samples WITH SCALING!!!! 
# Histogram of RAW sample read counts
jpeg(filename="raw_merged_scaled.jpeg", width= 20, height=15, units= "cm", pointsize= 14, res=500)
ggplot(data.frame(sum=sample_sums(scaled_merged)),aes(sum)) + ylab("Number of Sequences per Sample") +
  geom_histogram(colour="black",fill="magenta", binwidth = 18)  + xlab("Total Sequences") + 
  scale_x_continuous(seq(14600, 15200, by =100), labels = seq(14600, 15100, by =100)) +
  ggtitle(expression(atop("Raw Merged then Scaled Reads", atop("Range = 397"), ""))) 
dev.off()


# mean, max and min of sample read counts
mins<-min(sample_sums(scaled_merged))
paste(c("The minimum sample read count is",mins))
means<-mean(sample_sums(scaled_merged))
paste(c("The mean sample read count is",means))
maxs<-max(sample_sums(scaled_merged))
paste(c("The max sample read count is",maxs))

range_scaledmerged <- max(sample_sums(scaled_merged)) - min(sample_sums(scaled_merged))
paste(c("The range of sample read counts when merged (summed) and then scaled is",range_scaledmerged))





####  So now we have 2 types of merging:
### 1. Where we scaled our read counts, summed our reads between replicates, take the average + rounded, and then re-scale.
manual_merged

# CLUSTERING ANALYSIS
manual_otu <- otu_table(manual_merged) 
norm_manual_bc <- vegdist(manual_otu, method = "bray")  # calculates the Bray-Curtis Distances

### 2. Where we merged our samples (doubled the reads with samples that had replicates) and then scaled the read counts.
scaled_merged
data_dup <- sample_data(data_dup)
sc_tax <- otu_table(scaled_merged) 
sc_tax <- tax_table(scaled_merged)
merged_final <- merge_phyloseq(sc_tax, sc_otu, data_dup)
scaled_otu <- otu_table(merged_final) 
norm_scaled_bc <- vegdist(scaled_otu, method = "bray")  # calculates the Bray-Curtis Distances

### Comparing the two methods with a clustering analysis 
#jpeg(filename="clustering_merged_Comparison.jpeg", width= 45, height=32, units= "cm", pointsize= 14, res=500)
par(mfrow = c(2,1))
plot(hclust(norm_manual_bc), main = "Bray-Curtis Distance: Manual")
plot(hclust(norm_scaled_bc), main = "Bray-Curtis Distance: Scaled")
#dev.off()

nowin_manual <- subset_samples(manual_merged, names != "WINH" & names != "WINH3um")
ordu <- ordinate(nowin_manual, "NMDS", "bray")
plot_ordination(nowin_manual, ordu, shape = "trophicstate", color = "quadrant") + geom_point(size = 6) + ggtitle("Manual Merged")

nowin_scaled <- subset_samples(merged_final, names != "WINH" & names != "WINH3um")
ordu2 <- ordinate(nowin_scaled, "NMDS", "bray")
plot_ordination(nowin_scaled, ordu2, shape = "trophicstate", color = "quadrant") + geom_point(size = 6) + ggtitle("Scaled Merged")







