###  PLEASE CAREFULLY READ LINES 20-34 TO RUN THE FOLLOWING CODE.
#Load libraries within R version 3.2.0 (2015-04-16) -- "Full of Ingredients"
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
library(scales)
library(pgirmess)
library(multcompView)

### Set the Working Directory  
###  Within the "Final_PAFL_Trophicstate" directory there should be the following 2 sub-directories:
        #   1.  "raw_data":  includes 5 files:  taxonomy table, OTU (shared) table, "metadata", "duplicate_metadata.txt", "AllLakes_depthprofile.csv" 
        #   2.  "alpha_data":  includes 2 files:  "InvSimpson100" and "ObservedRichness100"  (These files are made between lines 261-295 and take ~20 minutes to calculate.  These files are added for reviewer convenience.)
setwd("~/Final_PAFL_Trophicstate")

### Source written functions in the file Functions_PAFL.R that is housed within the "Final_PAFL_Trophicstate"
source("Functions_PAFL.R")

## Figure 1: See line 631                Figure S1: See line 236
## Figure 2: See line 827                Figure S2: See line 1838
## Figure 3: See line 1051               Figure S3: See line 2816  
## Figure 4: See line 1987               Figure S4: See line 2917 
## Figure 5: See line 2283               Figure S5: See line 2496  
## Figure 6: See line 2659


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

#set phylum plotting colors
phylum.colors <- c(Acidobacteria = "grey26", Actinobacteria = "royalblue", Alphaproteobacteria = "plum2", Armatimonadetes = "red", Bacteroidetes = "darkorange",
                   "BD1-5" = "chartreuse", Betaproteobacteria = "slateblue2", Caldiserica = "black","Candidate_division_BRC1" = "violetred4",
                   "Candidate_division_JS1" = "aquamarine1",
                   "Candidate_division_OD1" = "#6DDE88", "Candidate_division_OP3" = "hotpink", "Candidate_division_OP8" = "goldenrod1", "Candidate_division_OP11" = "chocolate4",
                   "Candidate_division_SR1" = "tan3", "Candidate_division_TM7" = "skyblue1", "Candidate_division_WS3" = "magenta",
                   Chlamydiae="violet", Chlorobi="cyan2", Chloroflexi="darkgreen", Cyanobacteria = "chartreuse3", 
                   Deferribacteres = "slateblue3", "Deinococcus-Thermus" = "violetred", Dictyoglomi = "cornsilk4", Deltaproteobacteria = "deepskyblue", 
                   Elusimicrobia = "yellow3", Epsilonproteobacteria = "lightskyblue", Fibrobacteres = "darkred", Firmicutes = "blue4", FGL7S = "palevioletred1",
                   Fusobacteria = "slateblue1", Gammaproteobacteria = "steelblue4", Gemmatimonadetes="black", GOUTA4 = "plum1", "Hyd24-12" = "sienna2", JTB23 = "seashell2",
                   Lentisphaerae = "yellow1", "NPL-UPA2"="#652926", OC31 = "mediumpurple4", Planctomycetes = "mediumorchid3", Proteobacteria = "deepskyblue",
                   "SHA-109" = "lightsalmon3", SM2F11 = "lightskyblue2", SPOTSOCT00m83 = "orangered",
                   Spirochaetae = "gold3", Tenericutes="pink", Thermotogae = "chocolate1", TA06 = "lightslateblue",TA18 = "rosybrown3", TM6 = "olivedrab",
                   unclassified = "grey", Verrucomicrobia = "purple4", "WCHB1-60" = "palegreen")


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
  geom_histogram(colour="black",fill="dodgerblue2", binwidth = 750)  + xlab("Total Sequences") + ggtitle("Raw Sample Read Counts")

# mean, max and min of sample read counts
mins<-min(sample_sums(bact_samples))
paste(c("The minimum sample read count is",mins))
means<-mean(sample_sums(bact_samples))
paste(c("The mean sample read count is",means))
maxs<-max(sample_sums(bact_samples))
paste(c("The max sample read count is",maxs))


################################################################
### COMBINING BEFORE MCMURDIE HOLMES SCALING WITH OUR RAW READS
noscale_merge <- merge_samples(bact_samples, "dups", fun = "mean") #  THE MEAN DOES NOT WORK
noscale_merge <- prune_taxa(taxa_sums(noscale_merge) > 0, noscale_merge)
range_nonscale <- max(sample_sums(noscale_merge)) - min(sample_sums(noscale_merge)); range_nonscale
#######  Fixing the Metadata!
data_dup <- read.table("~/Final_PAFL_Trophicstate/raw_data/duplicate_metadata.txt", header = TRUE, sep = "\t")
row.names(data_dup) <- data_dup$names
data_dup$quadrant <- factor(data_dup$quadrant,levels = c("Free Epilimnion", "Free Mixed",  "Free Hypolimnion", "Particle Epilimnion", "Particle Mixed", "Particle Hypolimnion"))
data_dup <- sample_data(data_dup)  # Metadata in our phyloseq object
scale_otu <- otu_table(noscale_merge)  # OTU table in our phyloseq object
scale_tax <- tax_table(noscale_merge)  # Taxonomy Table in our phyloseq object
raw_merged <- merge_phyloseq(scale_otu, scale_tax, data_dup) #####  THIS IS OUR RAW MERGED SAMPLES PHYLOSEQ OBJECT.


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
#jpeg(filename="raw_merged_scaled.jpeg", width= 20, height=15, units= "cm", pointsize= 14, res=500)
ggplot(data.frame(sum=sample_sums(scaled_merged)),aes(sum)) + ylab("Number of Sequences per Sample") +
  geom_histogram(colour="black",fill="magenta", binwidth = 15)  + xlab("Total Sequences") + 
  scale_x_continuous(seq(14600, 15200, by =100), labels = seq(14600, 15100, by =100)) +
  ggtitle(expression(atop("Raw Merged then Scaled Reads", atop("Range = 397"), ""))) 
#dev.off()


# mean, max and min of sample read counts
mins<-min(sample_sums(scaled_merged))
paste(c("The minimum sample read count is",mins))
means<-mean(sample_sums(scaled_merged))
paste(c("The mean sample read count is",means))
maxs<-max(sample_sums(scaled_merged))
paste(c("The max sample read count is",maxs))

range_scaledmerged <- max(sample_sums(scaled_merged)) - min(sample_sums(scaled_merged))
paste(c("The range of sample read counts when merged (summed) and then scaled is",range_scaledmerged))


###  Time to create our phyloseq object with the merged and scaled sample reads.
data_dup <- read.table("~/Final_PAFL_Trophicstate/raw_data/duplicate_metadata.txt", header = TRUE, sep = "\t")
row.names(data_dup) <- data_dup$names
data_dup$quadrant <- factor(data_dup$quadrant,levels = c("Free Epilimnion", "Free Mixed",  "Free Hypolimnion", "Particle Epilimnion", "Particle Mixed", "Particle Hypolimnion"))
data_dup <- sample_data(data_dup)  # Metadata in our phyloseq object
sc_otu <- otu_table(scaled_merged)  # OTU table in our phyloseq object
sc_tax <- tax_table(scaled_merged)  # Taxonomy Table in our phyloseq object
merged_final <- merge_phyloseq(sc_tax, sc_otu, data_dup) 

########## ADD THE PROTEOBACTERIA TO THE PHYLA
phy <- data.frame(tax_table(merged_final))
Phylum <- as.character(phy$Phylum)
Class <- as.character(phy$Class)

for  (i in 1:length(Phylum)){ 
  if (Phylum[i] == "Proteobacteria"){
    Phylum[i] <- Class[i]
  } 
}

phy$Phylum <- Phylum
t <- tax_table(as.matrix(phy))

tax_table(merged_final) <- t

# Our phyloseq object!
merged_final 



####################################################  DEPTH PROFILE  ####################################################  Working up to Figure S1
####################################################  DEPTH PROFILE  ####################################################  Working up to Figure S1
####################################################  DEPTH PROFILE  ####################################################  Working up to Figure S1

profile <- read.csv(file = "~/Final_PAFL_Trophicstate/raw_data/AllLakes_depthprofile.csv", header = TRUE); 
profile <- subset(profile, select = c('lakename',"trophicstate", "depth", "temp", "DO", "pH", "SpC"))
#profile$depth <- as.numeric(profile$depth)
#profile$depth <- as.factor(profile$depth)
profile$DO. = NULL

profile$trophicstate <- as.character(profile$trophicstate)
profile$trophicstate[profile$trophicstate == "Eutrophic"] <- "Productive"
profile$trophicstate[profile$trophicstate == "Mesotrophic"] <- "Productive"
profile$trophicstate[profile$trophicstate == "Oligotrophic"] <- "Unproductive"
profile$trophicstate[profile$trophicstate == "Mixed"] <- "Mixed"

# ALL LAKES 
profile_all <- profile %>%  gather(Variable, Value, 4:7)
profile_all$lakename <- as.character(profile_all$lakename)
profile_all$lakename[profile_all$lakename == "WIN"] <- "Wintergreen"
profile_all$lakename[profile_all$lakename == "SIX"] <- "Sixteen"
profile_all$lakename[profile_all$lakename == "SHE"] <- "Sherman"
profile_all$lakename[profile_all$lakename == "PAY"] <- "Payne"
profile_all$lakename[profile_all$lakename == "LON"] <- "Little Long"
profile_all$lakename[profile_all$lakename == "LEE"] <- "Lee"
profile_all$lakename[profile_all$lakename == "GUL"] <- "Gull"
profile_all$lakename[profile_all$lakename == "BRI"] <- "Bristol"
profile_all$lakename[profile_all$lakename == "BAK"] <- "Baker"
profile_all$lakename[profile_all$lakename == "BAS"] <- "Baseline"
profile_all$lakename[profile_all$lakename == "BST"] <- "Bassett"
profile_all$lakename[profile_all$lakename == "BST"] <- "Bassett"
profile_all$trophicstate[profile_all$trophicstate == "Productive"] <- "High-Nutrient"
profile_all$trophicstate[profile_all$trophicstate == "Unproductive"] <- "Low-Nutrient"
profile_all$trophicstate[profile_all$trophicstate == "Mixed"] <- "High-Nutrient"


profile_all <- subset(profile_all, Variable != "SpC")

#####  Plotting FIGURE S1  #####  Plotting FIGURE S1  #####  Plotting FIGURE S1  #####  Plotting FIGURE S1  #####  Plotting FIGURE S1
## All LAKES
#tiff(filename="~/Final_PAFL_Trophicstate/Final_Figures/Fig.S1_AllLake_profiles.tiff", width= 30, height=40, units= "cm",pointsize= 18, res=200)
ggplot(profile_all, aes(x=Value, y = depth, color = lakename)) +   
  geom_path(size=2, alpha = 0.8) + ylab("Depth (m)") + xlab("") + 
  theme_bw() +  geom_point(size=4, alpha = 0.8) + geom_hline(h=0) +
  scale_y_reverse(breaks=seq(0, 30, 4), expand = c(0, 0)) + #breaks=seq(0, 30, 5), lim = c(30,0), expand = c(0, 0)
  scale_color_manual(name = "Lake Name", 
                     labels = c("Baker", "Baseline", "Bassett", "Bristol", "Gull", "Lee", "Little Long", "Payne", "Sherman", "Sixteen", "Wintergreen"), 
                     values = c("blue", "red", "black", "forestgreen","violet","limegreen", "purple", "orange", "maroon1", "turquoise3", "thistle4")) +
  facet_grid(trophicstate ~ Variable, scales = "free", labeller = profile_labeller, space = "free_y") +  
  theme(axis.title.x = element_text(face="bold", size=14),
        axis.text.x = element_text(colour = "black",size=14),
        axis.text.y = element_text(colour = "black", size=14),
        axis.title.y = element_text(face="bold", size=16),
        legend.title = element_text(size=12, face="bold"),
        legend.text = element_text(size = 12),
        strip.text = element_text(size = 16, face = "bold", colour = "black"),
        strip.background = element_blank(),
        legend.position=c(0.925, 0.16));
#dev.off()

####################################################  WITHIN-SAMPLE (ALPHA) DIVERSITY  ####################################################  Working up to Figure 1
####################################################  WITHIN-SAMPLE (ALPHA) DIVERSITY  ####################################################  Working up to Figure 1
####################################################  WITHIN-SAMPLE (ALPHA) DIVERSITY  ####################################################  Working up to Figure 1
# For alpha diversity we will RAREFY our phyloseq object.  Our raw-data phyloseq object was:
#raw_merged  # Raw Merged Samples, we have 14,378 OTUs
#raw_nowin <- subset_samples(raw_merged, names != "WINH" & names != "WINH3um")
#raw_nowin <- prune_taxa(taxa_sums(raw_nowin) > 0, raw_nowin)

#Data Import for rarefied data 
#raw_nowin  
#min(sample_sums(raw_nowin))

# Following code from Michelle's Butterflygut website: http://rstudio-pubs-static.s3.amazonaws.com/25722_9fc9bdc48f0d4f13b05fa61baeda57a0.html#alpha-diversity
# Rarefy to 14937 reads with replacement 100 times to estimate species richness
# Since we are rarefying to 14937 reads we need to remove the two samples with less than 1000 reads
#raredata14936 <- prune_samples(sample_sums(raw_nowin) > 14936, raw_nowin)

# Initialize matrices to store richness and evenness estimates
#richness <- matrix(nrow = 41,ncol = 100)  # Store richness:  We have 41 samples 
#row.names(richness) <- sample_names(raredata14936)
#evenness <- matrix(nrow = 41,ncol = 100)  #Store evenness
#row.names(evenness) <- sample_names(raredata14936)

# We want to be reproducible - so let's set the seed.
#set.seed(3)

# For 100 replications, rarefy the OTU table to 14936 reads and store the richness and evennes estimates in our 2 matrices we just created.
#The default for the rarefy_even_depth command is to pick with replacement so I set it to false. Picking without replacement is more computationally intensive 
#for (i in 1:100) {
#  r <- rarefy_even_depth(raredata14936, sample.size = 14936, verbose = FALSE, replace = FALSE)
#  rich <- as.numeric(as.matrix(estimate_richness(r, measures="Observed")))
#  richness[,i] <- rich
#  even <- as.numeric(as.matrix(estimate_richness(r, measures="InvSimpson")))
#  evenness[,i] <- even
#}

#write.table(richness, "~/Final_PAFL_Trophicstate/alpha_data/ObservedRichness100", row.names = TRUE)
#write.table(evenness, "~/Final_PAFL_Trophicstate/alpha_data/InvSimpson100", row.names = TRUE)

richness <- read.table("~/Final_PAFL_Trophicstate/alpha_data/ObservedRichness100", header = TRUE)
evenness <- read.table("~/Final_PAFL_Trophicstate/alpha_data/InvSimpson100", header = TRUE)

# Create a new matrix to hold the means and standard deviations of all the richness estimates
rich_stats = matrix(nrow= nrow(richness), ncol=1)
rich_stats[,1] = apply(richness, 1, mean)
#rich_stats[,2] = apply(richness, 1, sd)
rich_stats = data.frame(row.names(richness), rich_stats)
colnames(rich_stats) = c("samples","mean")
rich_stats$Test <- "Observed Richness"

# Create a new matrix to hold the means and standard deviations of the evenness estimates
even_stats = matrix(nrow = nrow(evenness), ncol = 1)
even_stats[,1] = apply(evenness, 1, mean)
#even_stats[,2] = apply(evenness, 1, sd)
even_stats = data.frame(row.names(evenness), even_stats)
colnames(even_stats) = c("samples","mean")
even_stats$Test <- "Inverse Simpson"


###  Add SIMPSON'S EVENNESS!
#simps_even = as.data.frame(matrix(nrow = nrow(evenness), ncol = 2))
simps_even <- data.frame(matrix(nrow = nrow(evenness), ncol = 3))
colnames(simps_even) = c("samples","mean", "Test")
simps_even$samples <- even_stats$samples
simps_even$mean <-  even_stats$mean/rich_stats$mean
simps_even$Test <- "Simpson's Evenness"

#Combine the rich.stats and the even.stats dataframes
alpha_stats <- rbind(rich_stats, even_stats, simps_even)
alpha_stats$names <- alpha_stats$samples
alpha_stats <- makeCategories_dups(alpha_stats)
alpha_stats$ProdLevel <- as.character(alpha_stats$trophicstate)
alpha_stats$ProdLevel[alpha_stats$trophicstate == "Eutrophic"] <- "Productive"
alpha_stats$ProdLevel[alpha_stats$trophicstate == "Mesotrophic"] <- "Productive"
alpha_stats$ProdLevel[alpha_stats$trophicstate == "Oligotrophic"] <- "Unproductive"
alpha_stats$ProdLevel[alpha_stats$lakenames == "Sherman"] <- "Mixed"
alpha_stats$trophicstate[alpha_stats$lakenames == "Sherman"] <- "Mixed"


####  Add all information together!
for(i in 1:length(alpha_stats$limnion)){
  alpha_stats$troph_lim[i]<-paste(as.character(alpha_stats$trophicstate[i]),
                                  as.character(alpha_stats$limnion[i]),
                                  as.character(alpha_stats$filter[i]))}


### Observed Richness and InvSimpson for Troph_lim
Meantroph_lim <- ddply(alpha_stats, ~Test+troph_lim, function(x){data.frame(Meantroph_lim = mean(x$mean))})
alpha_stats <- join(alpha_stats,Meantroph_lim) 
SEtroph_lim <- ddply(alpha_stats, ~Test+troph_lim, function(x){data.frame(SEtroph_lim = se(x$mean))})
alpha_stats <- join(alpha_stats,SEtroph_lim)
SDtroph_lim <- ddply(alpha_stats, ~Test+troph_lim, function(x){data.frame(SDtroph_lim = sd(x$mean))})
alpha_stats <- join(alpha_stats,SDtroph_lim)


###  OVERALL FILTER ALPHA STATS
Meanfilter_troph <- ddply(alpha_stats, ~Test+filter+trophicstate, function(x){data.frame(Meanfilter_troph = mean(x$mean))})
alpha_stats <- join(alpha_stats,Meanfilter_troph) 
SEfilter_troph <- ddply(alpha_stats, ~Test+filter+trophicstate, function(x){data.frame(SEfilter_troph = se(x$mean))})
alpha_stats <- join(alpha_stats,SEfilter_troph)
SDfilter_troph <- ddply(alpha_stats, ~Test+filter+trophicstate, function(x){data.frame(SDfilter_troph = sd(x$mean))})
alpha_stats <- join(alpha_stats,SDfilter_troph)

## Make sure our factors are in the right order!
alpha_stats$trophicstate <-factor(alpha_stats$trophicstate,levels=c("Eutrophic", "Mesotrophic", "Oligotrophic", "Mixed"))
alpha_stats$troph_lim <-factor(alpha_stats$troph_lim,levels=c("Eutrophic Epilimnion Particle", "Eutrophic Epilimnion Free", "Eutrophic Hypolimnion Particle", "Eutrophic Hypolimnion Free",
                                                              "Mesotrophic Epilimnion Particle", "Mesotrophic Epilimnion Free", "Mesotrophic Hypolimnion Particle", "Mesotrophic Hypolimnion Free",
                                                              "Oligotrophic Epilimnion Particle", "Oligotrophic Epilimnion Free", "Oligotrophic Hypolimnion Particle", "Oligotrophic Hypolimnion Free",
                                                              "Mixed Mixed Particle", "Mixed Mixed Free"))


alpha_filt <- ggplot(alpha_stats, aes(x = filter, y = Meanfilter_troph, color = filter)) + geom_point(size = 5, alpha = 0.7) +
  facet_grid(Test ~trophicstate, scales="free", space="free_x") + ggtitle("Rarefied: SD") +
  geom_errorbar(aes(ymin=Meanfilter_troph-SDfilter_troph, ymax=Meanfilter_troph+SDfilter_troph), width=.2, position=position_dodge(.9)) +
  xlab("Habitat") + ylab("Within Sample Diversity") + theme_bw() + 
  theme(axis.title.x = element_text(face="bold", size=16),
        axis.text.x = element_text(angle=30, colour = "black", vjust=1, hjust = 1, size=14),
        axis.text.y = element_text(colour = "black", size=14),
        axis.title.y = element_text(face="bold", size=16),
        plot.title = element_text(face="bold", size = 20),
        legend.title = element_text(size=14, face="bold"),
        legend.text = element_text(size = 14),
        strip.background = element_blank(), strip.text = element_text(size=14, face="bold"),
        legend.position="none"); alpha_filt

#jpeg(filename="~/Final_PAFL_Trophicstate/Figures/Fig.3b_alpha_TROPH_SD.jpeg", width= 28, height=28, units= "cm", pointsize= 14, res=500)
facet_alpha_stats <- ggplot(alpha_stats, aes(x = troph_lim, y = Meantroph_lim, color = troph_lim)) + geom_point(size = 5, alpha = 0.7) +
  facet_grid(Test ~ trophicstate, scales="free", space="free_x") + 
  geom_errorbar(aes(ymin=Meantroph_lim-SDtroph_lim, ymax=Meantroph_lim+SDtroph_lim), width=.2, position=position_dodge(.9)) +
  ggtitle("Alpha Diversity: Standard Deviation") + theme_bw() +   xlab("Habitat") + ylab("") + 
  scale_color_manual(name = "", limits=c("Eutrophic Epilimnion Particle", "Eutrophic Epilimnion Free", "Eutrophic Hypolimnion Particle", "Eutrophic Hypolimnion Free",
                                         "Mesotrophic Epilimnion Particle", "Mesotrophic Epilimnion Free", "Mesotrophic Hypolimnion Particle", "Mesotrophic Hypolimnion Free",
                                         "Oligotrophic Epilimnion Particle", "Oligotrophic Epilimnion Free", "Oligotrophic Hypolimnion Particle", "Oligotrophic Hypolimnion Free",
                                         "Mixed Mixed Particle", "Mixed Mixed Free"), 
                     values = c("deeppink", "deeppink", "deeppink", "deeppink","orange","orange","orange","orange", 
                                "turquoise3","turquoise3","turquoise3","turquoise3", "blue", "blue"))+
  scale_x_discrete(breaks=c("Eutrophic Epilimnion Particle", "Eutrophic Epilimnion Free", "Eutrophic Hypolimnion Particle", "Eutrophic Hypolimnion Free",
                            "Mesotrophic Epilimnion Particle", "Mesotrophic Epilimnion Free", "Mesotrophic Hypolimnion Particle", "Mesotrophic Hypolimnion Free",
                            "Oligotrophic Epilimnion Particle", "Oligotrophic Epilimnion Free", "Oligotrophic Hypolimnion Particle", "Oligotrophic Hypolimnion Free",
                            "Mixed Mixed Particle", "Mixed Mixed Free"), 
                   labels=c("Epilimnion Particle", "Epilimnion Free", "Hypolimnion Particle", "Hypolimnion Free",
                            "Epilimnion Particle", "Epilimnion Free", "Hypolimnion Particle", "Hypolimnion Free",
                            "Epilimnion Particle", "Epilimnion Free", "Hypolimnion Particle", "Hypolimnion Free",
                            "Particle", "Free")) + 
  theme(axis.title.x = element_text(face="bold", size=16),
        axis.text.x = element_text(angle=30, colour = "black", vjust=1, hjust = 1, size=14),
        axis.text.y = element_text(colour = "black", size=14),
        axis.title.y = element_text(face="bold", size=16),
        plot.title = element_text(face="bold", size = 20),
        strip.text.x = element_text(size=12, face="bold"),
        strip.text.y = element_text(size=14, face="bold"),
        legend.title = element_text(size=14, face="bold"),
        legend.text = element_text(size = 14),
        strip.background = element_blank(),
        legend.position="none"); facet_alpha_stats
#dev.off()



####   ALPHA DIVERSITY:  COMBINING EUTROPHIC AND MESOTROPHIC into High-Nutrient
####   ALPHA DIVERSITY:  COMBINING EUTROPHIC AND MESOTROPHIC into High-Nutrient 
####   ALPHA DIVERSITY:  COMBINING EUTROPHIC AND MESOTROPHIC Into High-Nutrient
prodalpha_stats <- subset(alpha_stats, select= c(1:10))
prodalpha_stats$trophicstate <- as.character(prodalpha_stats$trophicstate)
prodalpha_stats$trophicstate[prodalpha_stats$trophicstate == "Mesotrophic"] <- "Eutrophic"
prodalpha_stats$trophicstate[prodalpha_stats$trophicstate == "Eutrophic"] <- "Productive"
prodalpha_stats$trophicstate[prodalpha_stats$trophicstate == "Oligotrophic"] <- "Unproductive"
unique(prodalpha_stats$trophicstate)

for(i in 1:length(prodalpha_stats$limnion)){
  prodalpha_stats$troph_lim[i]<-paste(as.character(prodalpha_stats$trophicstate[i]),
                                      as.character(prodalpha_stats$limnion[i]),
                                      as.character(prodalpha_stats$filter[i]))}

### Observed Richness and InvSimpson for Troph_lim
Meantroph_lim <- ddply(prodalpha_stats, ~Test+troph_lim, function(x){data.frame(Meantroph_lim = mean(x$mean))})
prodalpha_stats <- join(prodalpha_stats,Meantroph_lim) 
SEtroph_lim <- ddply(prodalpha_stats, ~Test+troph_lim, function(x){data.frame(SEtroph_lim = se(x$mean))})
prodalpha_stats <- join(prodalpha_stats,SEtroph_lim)
SDtroph_lim <- ddply(prodalpha_stats, ~Test+troph_lim, function(x){data.frame(SDtroph_lim = sd(x$mean))})
prodalpha_stats <- join(prodalpha_stats,SDtroph_lim)


prodalpha_stats$trophicstate <-factor(prodalpha_stats$trophicstate,levels=c("Productive", "Unproductive", "Mixed"))
prodalpha_stats$troph_lim <-factor(prodalpha_stats$troph_lim,levels=c("Productive Epilimnion Particle", "Productive Epilimnion Free", "Productive Hypolimnion Particle", "Productive Hypolimnion Free",
                                                                      "Unproductive Epilimnion Particle", "Unproductive Epilimnion Free", "Unproductive Hypolimnion Particle", "Unproductive Hypolimnion Free",
                                                                      "Mixed Mixed Particle", "Mixed Mixed Free"))

#jpeg(filename="~/Final_PAFL_Trophicstate/Figures/Fig.3a_alpha_PROD_SD.jpeg", width= 25, height=28, units= "cm", pointsize= 14, res=500)
facet_prod_stats <- ggplot(prodalpha_stats, aes(x = troph_lim, y = Meantroph_lim, color = troph_lim)) + geom_point(size = 5, alpha = 0.7) +
  facet_grid(Test ~ trophicstate, scales="free", space="free_x") + 
  geom_errorbar(aes(ymin=Meantroph_lim-SDtroph_lim, ymax=Meantroph_lim+SDtroph_lim), width=.2, position=position_dodge(.9)) +
  ggtitle("Alpha Diversity: Standard Deviation") + theme_bw() +   xlab("Habitat") + ylab("") + 
  scale_color_manual(name = "", limits=c("Productive Epilimnion Particle", "Productive Epilimnion Free", "Productive Hypolimnion Particle", "Productive Hypolimnion Free",
                                         "Unproductive Epilimnion Particle", "Unproductive Epilimnion Free", "Unproductive Hypolimnion Particle", "Unproductive Hypolimnion Free",
                                         "Mixed Mixed Particle", "Mixed Mixed Free"), 
                     values = c("deeppink", "deeppink", "deeppink", "deeppink",
                                "turquoise3","turquoise3","turquoise3","turquoise3", "blue", "blue"))+
  scale_x_discrete(breaks=c("Productive Epilimnion Particle", "Productive Epilimnion Free", "Productive Hypolimnion Particle", "Productive Hypolimnion Free",
                            "Unproductive Epilimnion Particle", "Unproductive Epilimnion Free", "Unproductive Hypolimnion Particle", "Unproductive Hypolimnion Free",
                            "Mixed Mixed Particle", "Mixed Mixed Free"), 
                   labels=c("Epilimnion \nParticle-Associated", "Epilimnion \nFree-Living", "Hypolimnion \nParticle-Associated", "Hypolimnion \nFree-Living",
                            "Epilimnion \nParticle-Associated", "Epilimnion \nFree-Living", "Hypolimnion \nParticle-Associated", "Hypolimnion \nFree-Living",
                            "Particle", "Free")) + 
  theme(axis.title.x = element_text(face="bold", size=16),
        axis.text.x = element_text(angle=30, colour = "black", vjust=1, hjust = 1, size=14),
        axis.text.y = element_text(colour = "black", size=14),
        axis.title.y = element_text(face="bold", size=16),
        plot.title = element_text(face="bold", size = 20),
        strip.text.x = element_text(size=12, face="bold"),
        strip.text.y = element_text(size=14, face="bold"),
        legend.title = element_text(size=14, face="bold"),
        legend.text = element_text(size = 14),
        strip.background = element_blank(),
        legend.position="none"); facet_prod_stats
#dev.off()

prod_even <- subset(prodalpha_stats, Test == "Simpson's Evenness")
prod_richobs <-subset(prodalpha_stats, Test == "Observed Richness") 

prod_even$trophicstate <- as.character(prod_even$trophicstate)
prod_even$trophicstate[prod_even$trophicstate == "Productive"] <-"High-Nutrient"
prod_even$trophicstate[prod_even$trophicstate == "Unproductive"] <-"Low-Nutrient"

prod_even$trophicstate <-factor(prod_even$trophicstate,levels=c("High-Nutrient", "Low-Nutrient", "Mixed"))


####################################################################  STATS TIME!  Simpson's Meausre of Evenness
####################################################################  STATS TIME!  Simpson's Meausre of Evenness
####################################################################  STATS TIME!  Simpson's Meausre of Evenness
####  Run the test on ALL the data that goes into the mean!  #### Check it out here:  https://aquaticr.wordpress.com/2012/12/18/multiple-comparison-test-for-non-parametric-data/
hist(prod_even$mean, breaks = 40)  # Not normally distributed!!!
prod_even$mean <- as.numeric(prod_even$mean)
prod_even$troph_lim1 <- as.factor(prod_even$troph_lim)
## Do the KW test
even_prod_KW <- kruskal.test(prod_even$mean ~ prod_even$troph_lim) # Kruskal Wallis test on evensen!
print(even_prod_KW)  # show Kruskal Wallis result
### Which samples are significantly different from each other?  Significant???  YES! WOOOHOOOOOOO!
even_prod_KW_MC <- kruskalmc(prod_even$mean ~ prod_even$troph_lim)  ## Defaults to P < 0.05
print(even_prod_KW_MC)
### Time to figure out letters to represent significance in a plot!
even_test <- even_prod_KW_MC$dif.com$difference # select logical vector
names(even_test) <- row.names(even_prod_KW_MC$dif.com) # add comparison names
# create a list with "homogenous groups" coded by letter
even_letters <- multcompLetters(even_test, compare="<", threshold=0.05, Letters=c(letters, LETTERS, "."), reversed = FALSE)#['Letters']
###  Let's extract the values from the multcompLetters object
even_sigs_dataframe <-  data.frame(as.vector(names(even_letters$Letters)), as.vector(even_letters$Letters))
colnames(even_sigs_dataframe) <- c("troph_lim", "siglabel")
even_try <- merge(prod_even, even_sigs_dataframe, by = "troph_lim")

#geom_text(aes(label = siglabel, x = troph_lim1, y = ((mean+sd) + 0.03)), size = 5) +


plot_even_sigs <- ggplot(even_try, aes(x = troph_lim, y = Meantroph_lim, color = troph_lim)) + geom_point(size = 5, alpha = 0.7) +
  facet_grid(Test ~ trophicstate, scales="free", space="free_x") + 
  geom_text(aes(label = siglabel, x = troph_lim, y = ((Meantroph_lim+SDtroph_lim) + 0.006)), size = 5) +
  geom_errorbar(aes(ymin=Meantroph_lim-SDtroph_lim, ymax=Meantroph_lim+SDtroph_lim), width=.2, position=position_dodge(.9)) +
  theme_bw() +   xlab("Habitat") + ylab("Simpson's Evenness") + #ggtitle("Within-Sample Diversity") +
  scale_color_manual(name = "", limits=c("Productive Epilimnion Particle", "Productive Epilimnion Free", "Productive Hypolimnion Particle", "Productive Hypolimnion Free",
                                         "Unproductive Epilimnion Particle", "Unproductive Epilimnion Free", "Unproductive Hypolimnion Particle", "Unproductive Hypolimnion Free",
                                         "Mixed Mixed Particle", "Mixed Mixed Free"), 
                     values = c("black", "black", "black", "black",
                                "gray48","gray48","gray48","gray48", "gray63", "gray63"))+
  scale_x_discrete(breaks=c("Productive Epilimnion Particle", "Productive Epilimnion Free", "Productive Hypolimnion Particle", "Productive Hypolimnion Free",
                            "Unproductive Epilimnion Particle", "Unproductive Epilimnion Free", "Unproductive Hypolimnion Particle", "Unproductive Hypolimnion Free",
                            "Mixed Mixed Particle", "Mixed Mixed Free"), 
                   labels=c("Epilimnion \nParticle-Associated", "Epilimnion \nFree-Living", "Hypolimnion \nParticle-Associated", "Hypolimnion \nFree-Living",
                            "Epilimnion \nParticle-Associated", "Epilimnion \nFree-Living", "Hypolimnion \nParticle-Associated", "Hypolimnion \nFree-Living",
                            "Particle-Associated", "Free-Living")) + 
  theme(axis.title.x = element_text(face="bold", size=16),
        axis.text.x = element_text(angle=30, colour = "black", vjust=1, hjust = 1, size=14),
        axis.text.y = element_text(colour = "black", size=14),
        axis.title.y = element_text(face="bold", size=16),
        plot.title = element_text(face="bold", size = 20),
        strip.text.x = element_text(size=16, face="bold"),
        strip.text.y = element_blank(),
        legend.title = element_text(size=14, face="bold"),
        legend.text = element_text(size = 14),
        strip.background = element_blank(),
        legend.position="none"); plot_even_sigs


prod_richobs <-subset(prodalpha_stats, Test == "Observed Richness") 

prod_richobs$trophicstate <- as.character(prod_richobs$trophicstate)
prod_richobs$trophicstate[prod_richobs$trophicstate == "Productive"] <-"High-Nutrient"
prod_richobs$trophicstate[prod_richobs$trophicstate == "Unproductive"] <-"Low-Nutrient"

prod_richobs$trophicstate <-factor(prod_richobs$trophicstate,levels=c("High-Nutrient", "Low-Nutrient", "Mixed"))


####################################################################  STATS TIME!  Richness
####################################################################  STATS TIME!  Richness
####################################################################  STATS TIME!  Richness
####  Run the test on ALL the data that goes into the mean!  #### Check it out here:  https://aquaticr.wordpress.com/2012/12/18/multiple-comparison-test-for-non-parametric-data/
hist(prod_richobs$mean, breaks = 40)  # Not normally distributed!!!
prod_richobs$mean <- as.numeric(prod_richobs$mean)
prod_richobs$troph_lim1 <- as.factor(prod_richobs$troph_lim)
## Do the KW test
richobs_prod_KW <- kruskal.test(prod_richobs$mean ~ prod_richobs$troph_lim) # Kruskal Wallis test on Observed Richness!
print(richobs_prod_KW)  # show Kruskal Wallis result
### Which samples are significantly different from each other?  Significant???  YES! WOOOHOOOOOOO!
richobs_prod_KW_MC <- kruskalmc(prod_richobs$mean ~ prod_richobs$troph_lim)  ## Defaults to P < 0.05
print(richobs_prod_KW_MC)
### Time to figure out letters to represent significance in a plot!
richobs_test <- richobs_prod_KW_MC$dif.com$difference # select logical vector
names(richobs_test) <- row.names(richobs_prod_KW_MC$dif.com) # add comparison names
# create a list with "homogenous groups" coded by letter
richobs_letters <- multcompLetters(richobs_test, compare="<", threshold=0.05, Letters=c(letters, LETTERS, "."), reversed = FALSE)#['Letters']
###  Let's extract the values from the multcompLetters object
richobs_sigs_dataframe <-  data.frame(as.vector(names(richobs_letters$Letters)), as.vector(richobs_letters$Letters))
colnames(richobs_sigs_dataframe) <- c("troph_lim", "siglabel")
richobs_try <- merge(prod_richobs, richobs_sigs_dataframe, by = "troph_lim")


plot_richobs_sigs <- ggplot(richobs_try, aes(x = troph_lim, y = Meantroph_lim, color = troph_lim)) + geom_point(size = 5, alpha = 0.7) +
  facet_grid(Test ~ trophicstate, scales="free", space="free_x") + 
  geom_text(aes(label = siglabel, x = troph_lim, y = ((Meantroph_lim+SDtroph_lim) + 75)), size = 5) +
  geom_errorbar(aes(ymin=Meantroph_lim-SDtroph_lim, ymax=Meantroph_lim+SDtroph_lim), width=.2, position=position_dodge(.9)) +
  theme_bw() +   xlab("Habitat") + ylab("Observed Richness") + 
  scale_color_manual(name = "", limits=c("Productive Epilimnion Particle", "Productive Epilimnion Free", "Productive Hypolimnion Particle", "Productive Hypolimnion Free",
                                         "Unproductive Epilimnion Particle", "Unproductive Epilimnion Free", "Unproductive Hypolimnion Particle", "Unproductive Hypolimnion Free",
                                         "Mixed Mixed Particle", "Mixed Mixed Free"), 
                     values = c("black", "black", "black", "black",
                                "gray48","gray48","gray48","gray48", "gray63", "gray63"))+
  scale_x_discrete(breaks=c("Productive Epilimnion Particle", "Productive Epilimnion Free", "Productive Hypolimnion Particle", "Productive Hypolimnion Free",
                            "Unproductive Epilimnion Particle", "Unproductive Epilimnion Free", "Unproductive Hypolimnion Particle", "Unproductive Hypolimnion Free",
                            "Mixed Mixed Particle", "Mixed Mixed Free"), 
                   labels=c("Epilimnion \nParticle-Associated", "Epilimnion \nFree-Living", "Hypolimnion \nParticle-Associated", "Hypolimnion \nFree-Living",
                            "Epilimnion \nParticle-Associated", "Epilimnion \nFree-Living", "Hypolimnion \nParticle-Associated", "Hypolimnion \nFree-Living",
                            "Particle-\nAssociated", "Free-\nLiving")) + 
  theme(axis.title.x = element_text(face="bold", size=16),
        axis.text.x = element_text(angle=30, colour = "black", vjust=1, hjust = 1, size=14),
        axis.text.y = element_text(colour = "black", size=14),
        axis.title.y = element_text(face="bold", size=16),
        plot.title = element_text(face="bold", size = 20),
        strip.text.x = element_text(size=16, face="bold"),
        strip.text.y = element_blank(),
        legend.title = element_text(size=14, face="bold"),
        legend.text = element_text(size = 14),
        strip.background = element_blank(),
        legend.position="none"); plot_richobs_sigs




########## Almost final PLOT!  (Bray-Curtis on top)
plot_even_sigs2 <- plot_even_sigs + scale_y_continuous(breaks=seq(0.01, 0.09, 0.02), lim = c(0.01,0.093)) +
  theme(axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank(), 
        plot.margin = unit(c(0.1, 0.1, 0.1, 0.23), "cm")) #top, right, bottom, left)

plot_richobs_sigs2 <- plot_richobs_sigs + scale_y_continuous(breaks=seq(400, 1200, 200), lim = c(300,1300)) +
  theme(strip.text.x = element_blank(), plot.margin = unit(c(-0.8, 0.1, 0.2, 0.1), "cm")) #top, right, bottom, left)

#tiff(filename="~/Final_PAFL_Trophicstate/Final_Figures/Fig.1_alpha_SIGS.tiff", width= 30, height=22, units= "cm", pointsize= 14, res=200)
#grid.newpage()
#pushViewport(viewport(layout=grid.layout(2,1,height=c(0.45,0.55))))
#print(plot_even_sigs2, vp=viewport(layout.pos.row=1,layout.pos.col=1))
#print(plot_richobs_sigs2, vp=viewport(layout.pos.row=2,layout.pos.col=1))
#dev.off()




########## FLIP THE FINAL PLOT  --> Sorensen on top, Bray-Curtis on the bottom
plot_richobs_sigs_FLIP <- plot_richobs_sigs + scale_y_continuous(breaks=seq(400, 1200, 200), lim = c(300,1300)) +
  theme(axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank(), 
        plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), "cm")) #top, right, bottom, left)


plot_even_sigs_FLIP <- plot_even_sigs + scale_y_continuous(breaks=seq(0.01, 0.09, 0.02), lim = c(0.01,0.093)) +
  theme(strip.text.x = element_blank(), plot.margin = unit(c(-0.8, 0.1, 0.2, 0.23), "cm")) #top, right, bottom, left)

#####  Plotting FIGURE 1  #####  Plotting FIGURE 1  #####  Plotting FIGURE 1  #####  Plotting FIGURE 1  #####  Plotting FIGURE 1
#tiff(filename="~/Final_PAFL_Trophicstate/Final_Figures/Fig.1_alpha_SIGS.tiff", width= 30, height=22, units= "cm", pointsize= 14, res=200)
pdf(file="~/Final_PAFL_Trophicstate/Final_Figures/Fig.1_alpha_SIGS.pdf", width= 4, height=3)
grid.newpage()
pushViewport(viewport(layout=grid.layout(2,1,height=c(0.45,0.55))))
print(plot_richobs_sigs_FLIP, vp=viewport(layout.pos.row=1,layout.pos.col=1))
print(plot_even_sigs_FLIP, vp=viewport(layout.pos.row=2,layout.pos.col=1))
dev.off()







####################################################################################  OTU OVERLAP #############  Working up to Figure 2
####################################################################################  OTU OVERLAP #############  Working up to Figure 2
####################################################################################  OTU OVERLAP #############  Working up to Figure 2
# For alpha diversity we will RAREFY our phyloseq object.  Our raw-data phyloseq object was:
raw_merged  # Raw Merged Samples, we have 14,378 OTUs
raw_nowin <- subset_samples(raw_merged, names != "WINH" & names != "WINH3um")
raw_nosherwin <- subset_samples(raw_nowin, lakenames != "Sherman")
raw_nosherwin <- prune_taxa(taxa_sums(raw_nosherwin) > 0, raw_nosherwin)

#Data Import for rarefied data 
raw_nosherwin  
min(sample_sums(raw_nowin))

# Following code from Michelle's Butterflygut website: http://rstudio-pubs-static.s3.amazonaws.com/25722_9fc9bdc48f0d4f13b05fa61baeda57a0.html#alpha-diversity
# Rarefy to 14937 reads with replacement 100 times to estimate species richness
# Since we are rarefying to 14937 reads we need to remove the two samples with less than 1000 reads
nosherwin_raredata14936 <- prune_samples(sample_sums(raw_nosherwin) > 14936, raw_nosherwin)
ggplot(data.frame(sum=sample_sums(nosherwin_raredata14936)),aes(sum)) + ylab("Number of Sequences per Sample") +
  geom_histogram(colour="black",fill="violetred")  + xlab("Total Sequences") 


nosherwin_raredata14936 <- rarefy_even_depth(nosherwin_raredata14936, sample.size = 14936, rngseed =3)
ggplot(data.frame(sum=sample_sums(nosherwin_raredata14936)),aes(sum)) + ylab("Number of Sequences per Sample") +
  geom_histogram(colour="black",fill="deepskyblue")  + xlab("Total Sequences") 

nosherwin_raredata14936 = prune_taxa(taxa_sums(nosherwin_raredata14936) > 0, nosherwin_raredata14936)
good_OTU <- nosherwin_raredata14936

#################  FREE LIVING VS PARTICLE ASSOCIATED 
#Create FREE LIVING DF 
freeOTU <- subset_samples(good_OTU, filter == "Free")
freeOTU2 <- prune_taxa(taxa_sums(freeOTU) > 0, freeOTU)
freeOTU_melt <- psmelt(freeOTU2)
freeOTU_melt <- subset(freeOTU_melt, select = c("Sample", "filter", "Species"))
dim(freeOTU_melt)
# Creating the df 
free_vector <- as.vector(unique(freeOTU_melt$Species))
erp <- rep("Free", len = length(free_vector))
unique_free <- as.data.frame(cbind(erp, free_vector))
colnames(unique_free) <- c("Preference", "OTU")
unique_free$OTU <- as.character(unique_free$OTU)

## CREATE PARTICLE ASSOCIATED DF 
partOTU <- subset_samples(good_OTU, filter == "Particle")
partOTU2 <- prune_taxa(taxa_sums(partOTU) > 0, partOTU)
partOTU_melt <- psmelt(partOTU2)
partOTU_melt <- subset(partOTU_melt, select = c("Sample", "filter", "Species"))
dim(partOTU_melt)
#creating the df 
part_vector <- as.vector(unique(partOTU_melt$Species))
erppart <- rep("Particle", len = length(part_vector))
unique_part <- as.data.frame(cbind(erppart, part_vector))
colnames(unique_part) <- c("Preference", "OTU")
unique_part$OTU <- as.character(unique_part$OTU)

overlap <- length(intersect(unique_part$OTU, unique_free$OTU))
free_only <- length(setdiff(unique_free$OTU, unique_part$OTU))
part_only <- length(setdiff(unique_part$OTU, unique_free$OTU))
free_total <- length(free_vector)
part_total <- length(part_vector)


# Create Data Frame with NA's
PAFL_otu <- data.frame(matrix(NA, nrow=6, ncol=3))
# Confirm Size of Data Frame
dim(PAFL_otu)
# Change Variable Names
names(PAFL_otu) <- c("SampleType","Type", "NumOTUs")
PAFL_otu$SampleType <- c("Particle-Associated", "Particle-Associated", "Particle-Associated", "Free-Living", "Free-Living", "Free-Living")
PAFL_otu$Type <- c("Total", "Particle-Associated Only", "Both", "Total", "Free-Living Only", "Both")
PAFL_otu$NumOTUs <- c(part_total, part_only, overlap, free_total, free_only, overlap)
PAFL_otu2 <- subset(PAFL_otu, Type != "Total")
#Fisher's test = INSIGNIFICANT

#################  TOP VS BOTTOM 
#Create EPILIMNION  DF 
topOTU <- subset_samples(good_OTU, limnion == "Epilimnion")
topOTU2 <- prune_taxa(taxa_sums(topOTU) > 0, topOTU)
topOTU_melt <- psmelt(topOTU2)
topOTU_melt <- subset(topOTU_melt, select = c("Sample", "limnion", "Species"))
dim(topOTU_melt)
# Creating the df 
top_vector <- as.vector(unique(topOTU_melt$Species))
toperp <- rep("Epilimnion", len = length(top_vector))
unique_top <- as.data.frame(cbind(toperp, top_vector))
colnames(unique_top) <- c("Preference", "OTU")
unique_top$OTU <- as.character(unique_top$OTU)

## CREATE PARTICLE ASSOCIATED DF 
bottomOTU <- subset_samples(good_OTU, limnion == "Hypolimnion")
bottomOTU2 <- prune_taxa(taxa_sums(bottomOTU) > 0, bottomOTU)
bottomOTU_melt <- psmelt(bottomOTU2)
bottomOTU_melt <- subset(bottomOTU_melt, select = c("Sample", "limnion", "Species"))
dim(bottomOTU_melt)
#creating the df 
bot_vector <- as.vector(unique(bottomOTU_melt$Species))
boterp <- rep("Hypolimnion", len = length(bot_vector))
unique_bottom <- as.data.frame(cbind(boterp, bot_vector))
colnames(unique_bottom) <- c("Preference", "OTU")
unique_bottom$OTU <- as.character(unique_bottom$OTU)


library(dplyr)
TB_overlap <- length(intersect(unique_bottom$OTU, unique_top$OTU))
top_only <- length(setdiff(unique_top$OTU, unique_bottom$OTU))
bottom_only <- length(setdiff(unique_bottom$OTU, unique_top$OTU))
top_total <- length(top_vector)
bottom_total <- length(bot_vector)

# Create Data Frame with NA's
TB_otu <- data.frame(matrix(NA, nrow=6, ncol=3))
# Confirm Size of Data Frame
dim(TB_otu)
# Change Variable Names
names(TB_otu) <- c("SampleType","Type", "NumOTUs")
TB_otu$SampleType <- c("Epilimnion", "Epilimnion", "Epilimnion", "Hypolimnion", "Hypolimnion", "Hypolimnion")
TB_otu$Type <- c("Total", "Epilimnion Only", "Both", "Total", "Hypolimnion Only", "Both")
TB_otu$NumOTUs <- c(top_total, top_only, TB_overlap, bottom_total, bottom_only, TB_overlap)
TB_otu2 <- subset(TB_otu, Type != "Total")
# Fisher's test:  p = 0.0001

#################  EU/MESO/OLIGO + PRODUCTIVE VS UNPRODUCTIVE 
#Create EPILIMNION  DF 
prodOTU <- subset_samples(good_OTU, ProdLevel == "Productive")
prodOTU2 <- prune_taxa(taxa_sums(prodOTU) > 0, prodOTU)
prodOTU_melt <- psmelt(prodOTU2)
prodOTU_melt <- subset(prodOTU_melt, select = c("Sample", "ProdLevel", "Species"))
dim(prodOTU_melt)
# Creating the df 
prod_vector <- as.vector(unique(prodOTU_melt$Species))
proderp <- rep("Productve", len = length(prod_vector))
unique_prod <- as.data.frame(cbind(proderp, prod_vector))
colnames(unique_prod) <- c("Preference", "OTU")
unique_prod$OTU <- as.character(unique_prod$OTU)

## CREATE PARTICLE ASSOCIATED DF 
unprodOTU <- subset_samples(good_OTU, ProdLevel == "Unproductive")
unprodOTU2 <- prune_taxa(taxa_sums(unprodOTU) > 0, unprodOTU)
unprodOTU_melt <- psmelt(unprodOTU2)
unprodOTU_melt <- subset(unprodOTU_melt, select = c("Sample", "ProdLevel", "Species"))
dim(unprodOTU_melt)
#creating the df 
unprod_vector <- as.vector(unique(unprodOTU_melt$Species))
unproderp <- rep("Unproductive", len = length(unprod_vector))
unique_unprod <- as.data.frame(cbind(unproderp, unprod_vector))
colnames(unique_unprod) <- c("Preference", "OTU")
unique_unprod$OTU <- as.character(unique_unprod$OTU)

Prod_overlap <- length(intersect(unique_unprod$OTU, unique_prod$OTU))
prod_only <- length(setdiff(unique_prod$OTU, unique_unprod$OTU))
unprod_only <- length(setdiff(unique_unprod$OTU, unique_prod$OTU))
prod_total <- length(prod_vector)
unprod_total <- length(unprod_vector)

# Create Data Frame with NA's
prod_otu <- data.frame(matrix(NA, nrow=6, ncol=3))
# Confirm Size of Data Frame
dim(prod_otu)
# Change Variable Names
names(prod_otu) <- c("SampleType","Type", "NumOTUs")
prod_otu$SampleType <- c("High-Nutrient", "High-Nutrient", "High-Nutrient", "Low-Nutrient", "Low-Nutrient", "Low-Nutrient")
prod_otu$Type <- c("Total", "High-Nutrient Only", "Both", "Total", "Low-Nutrient Only", "Both")
prod_otu$NumOTUs <- c(prod_total, prod_only, Prod_overlap, unprod_total, unprod_only, Prod_overlap)
prod_otu2 <- subset(prod_otu, Type != "Total")
#Fisher's test:  p = 0.0001


###  PUT IT ALL TOGETHER! 
PAFL_otu2$Comparison <- "Particle-Associated \n vs. \nFree-Living"
TB_otu2$Comparison <- "Epilimnion \n vs. \nHypolimnion"
prod_otu2$Comparison <- "High-Nutrient\n  vs. \nLow-Nutrient"
otusums <- rbind(PAFL_otu2, TB_otu2, prod_otu2)

otusums$Comparison <- factor(otusums$Comparison,levels = c("Particle-Associated \n vs. \nFree-Living", "High-Nutrient\n  vs. \nLow-Nutrient", "Epilimnion \n vs. \nHypolimnion"))
otusums$SampleType <- factor(otusums$SampleType,levels = c("Free-Living", "Particle-Associated", "Low-Nutrient", "High-Nutrient", "Epilimnion", "Hypolimnion"))
otusums$Type <- factor(otusums$Type, levels = c("Both","Free-Living Only", "Particle-Associated Only", "Low-Nutrient Only", "High-Nutrient Only", "Epilimnion Only", "Hypolimnion Only"))
#otusums$Type <- factor(otusums$Type, levels = rev(levels(otusums$Type)))


label.df <- data.frame(Group = c("Productive \n vs. Unproductive", "Epilimnion \n vs. Hypolimnion"),
                       Value = c(8600, 8000))

#####  Plotting FIGURE 2  #####  Plotting FIGURE 2  #####  Plotting FIGURE 2  #####  Plotting FIGURE 2  #####  Plotting FIGURE 2
#tiff(filename="~/Final_PAFL_Trophicstate/Final_Figures/Fig.2_DetectedOTUs_rarefied.tiff", width= 25, height=19, units= "cm", pointsize= 10, res=200)
ggplot(otusums, aes(y=NumOTUs , x=SampleType, fill=Type, order=Type)) +
  facet_grid(. ~ Comparison,scales = "free") + #geom_text(x = 2, y = 8750, label = "***") +
  xlab("Sample Type ") + ylab("Number of Deteceted UniqueOTUs") + 
  geom_bar(stat="identity") +   geom_bar(stat="identity", colour="black", show_guide=FALSE) +  
  scale_y_continuous(expand = c(0,0), breaks=seq(0, 6000, 1000), lim = c(0, 6000)) + 
  guides(fill = guide_legend(keywidth = 1.5, keyheight = 1.5)) + 
  scale_fill_manual(limits = c("Free-Living Only", "Particle-Associated Only", "Low-Nutrient Only", "High-Nutrient Only", "Epilimnion Only", "Hypolimnion Only", "Both"),
                    breaks = c("Free-Living Only", "Particle-Associated Only", "Low-Nutrient Only", "High-Nutrient Only", "Epilimnion Only", "Hypolimnion Only", "Both"),
                    values = c( "goldenrod1", "firebrick1", "turquoise3", "green4","palevioletred1","cornflowerblue", "gray39")) +
  #values = c("orange", "red", "darkgreen", "limegreen", "deepskyblue", "blue4", "gray39")) +
  theme_classic() + theme(axis.title.x = element_text(face="bold", size=16),
                          axis.text.x = element_text(angle=20, vjust = 1, hjust = 1, colour = "black", size=14),
                          axis.text.y = element_text(colour = "black", size=14),
                          axis.title.y = element_text(face="bold", size=16),
                          plot.title = element_text(face="bold", size = 20),
                          strip.text.x = element_text(size=14, face="bold"),
                          legend.title = element_blank(),
                          legend.text = element_text(size = 12),
                          strip.background = element_blank(),
                          legend.position="right")
#dev.off()
#

####  Test for Significance!
PAFL_chisq <-matrix(c(2734,2620,2263,2263),nrow=2)
PAFL_chisq_res <- chisq.test(PAFL_chisq, correct = FALSE)
prod_chisq <-matrix(c(4289,1502,1826,1502),nrow=2)
prod_chisq_res <- chisq.test(prod_chisq, correct = FALSE)
limnion_chisq <-matrix(c(2347,1540,3730,1540),nrow=2)
limnion_chisq_res <- chisq.test(limnion_chisq, correct = FALSE)
shared_chisq <- matrix(c(2263,1502,1540),nrow=1)
shared_chisq_res <- chisq.test(shared_chisq, correct = FALSE)







####################################################  ORDINATIONS  ####################################################  Working up to Figure 3
####################################################  ORDINATIONS  ####################################################  Working up to FIgure 3
####################################################  ORDINATIONS  ####################################################  Working up to Figure 3
### Clustering 
### Get rid of the Wintergreen HYPOLIMNION Samples
nowin_merged <- subset_samples(merged_final, names != "WINH" & names != "WINH3um")
nowin_merged <- prune_taxa(taxa_sums(nowin_merged) > 0, nowin_merged)
nowinOTU <- otu_table(nowin_merged)
#weighted
norm_bray <- vegdist(nowinOTU, method = "bray", binary = FALSE)  # calculates the Bray-Curtis Distances
nowinOTU_df <- data.frame(otu_table(nowin_merged))  
norm_soren <- vegdist(nowinOTU_df, method = "bray", binary = TRUE)  # SORENSEN INDEX


### Clustering based on 
#jpeg(filename="clustering_bray+soren.jpeg", width= 45, height=32, units= "cm", pointsize= 14, res=500)
par(mfrow = c(2,1))
plot(hclust(norm_bray), main = "Bray-Curtis Distance")
plot(hclust(norm_soren), main = "Sorensen Distance")
#dev.off()

# vector of colors
mypal = c("black", "red", "blue", "purple", "green3", "orange", "maroon1", "gold")
# cutting dendrogram in 5 clusters
hc <- hclust(otu_bray)
clus5 = cutree(hc, 5)
par(mfrow = c(1,1))
plot(as.phylo(hc), tip.color = mypal[clus5],  main = "Bray Curtis Dissimilarity") 



########  NMDS + PCoA Plots
#########  PLOT ORDINATIONS FOR BOTH SCALED AND MANUAL
bray_pcoa <- pcoa(norm_bray)
bray_pcoa2 <- bray_pcoa$vectors
bray_pcoa3 <- data.frame(bray_pcoa2[, 1:3])
bray_pcoa3$names <- row.names(bray_pcoa3)
bray_pcoa4 <- makeCategories_dups(bray_pcoa3)
bray_pcoa4$quadrant <- factor(bray_pcoa4$quadrant,levels = c("Free Epilimnion", "Free Mixed",  "Free Hypolimnion", "Particle Epilimnion", "Particle Mixed", "Particle Hypolimnion"))
### Let's extract the values of each axis for the pcoa
bray_values <- bray_pcoa$values
bray_rel_eigens <- bray_values$Relative_eig
bray_rel_eigen1 <- bray_rel_eigens[1]
bray_rel_eigen1_percent <- round(bray_rel_eigen1 * 100, digits = 1)
bray_rel_eigen2 <- bray_rel_eigens[2]
bray_rel_eigen2_percent <- round(bray_rel_eigen2 * 100, digits = 1)
bray_axis1 <- paste("PCoA1:",bray_rel_eigen1_percent,"%")
bray_axis2 <- paste("PCoA2:",bray_rel_eigen2_percent,"%")

pcoa_bray <- ggplot(bray_pcoa4, aes(Axis.1, Axis.2 * -1, color = quadrant, shape = trophicstate)) +
  xlab(bray_axis1) + ylab(bray_axis2) + #ggtitle("Bray-Curtis: Trophic State") +
  geom_point(size= 6, alpha=0.9) + theme_bw() + #ylim(1, -1) + xlim(1, -1) +
  scale_color_manual(name = "Habitat", breaks=c("Free Epilimnion", "Free Mixed",  "Free Hypolimnion", "Particle Epilimnion", "Particle Mixed", "Particle Hypolimnion"),
                     labels = c("Free-Living Epilimnion", "Free-Living Mixed",  "Free-Living Hypolimnion", "Particle-Associated Epilimnion", "Particle-Associated Mixed", "Particle-Associated Hypolimnion"), 
                     values = c("purple3", "mediumblue", "darkgreen", "orchid1", "deepskyblue", "green")) +
  scale_shape_manual(name = "Trophic State", breaks = c("Eutrophic", "Mesotrophic", "Oligotrophic", "Mixed"), 
                     labels = c("Eutrophic", "Mesotrophic", "Oligotrophic", "Mixed"),
                     values = c(15, 19, 17)) +
  theme(axis.text.x = element_text(colour="black", vjust=0.5, size=14), 
        axis.text.y = element_text(colour="black", vjust=0.5, size=14),
        axis.title.x = element_text(face="bold", size=16),
        axis.title.y = element_text(face="bold", size=16),
        legend.title = element_text(size=12, face="bold"),
        legend.text = element_text(size = 12),
        ###LEGEND TOP RIGHT CORNER
        legend.position = "none"); pcoa_bray

##########  SORENSEN'S DISSIMILARITY
soren_pcoa <- pcoa(norm_soren)
soren_pcoa2 <- soren_pcoa$vectors
soren_pcoa3 <- data.frame(soren_pcoa2[, 1:3])
soren_pcoa3$names <- row.names(soren_pcoa3)
soren_pcoa4 <- makeCategories_dups(soren_pcoa3)
soren_pcoa4$quadrant <- factor(soren_pcoa4$quadrant,levels = c("Free Epilimnion", "Free Mixed",  "Free Hypolimnion", "Particle Epilimnion", "Particle Mixed", "Particle Hypolimnion"))
### Let's extract the values of each axis for the pcoa
soren_values <- soren_pcoa$values
soren_rel_eigens <- soren_values$Relative_eig
soren_rel_eigen1 <- soren_rel_eigens[1]
soren_rel_eigen1_percent <- round(soren_rel_eigen1 * 100, digits = 1)
soren_rel_eigen2 <- soren_rel_eigens[2]
soren_rel_eigen2_percent <- round(soren_rel_eigen2 * 100, digits = 1)
soren_axis1 <- paste("PCoA1:",soren_rel_eigen1_percent,"%")
soren_axis2 <- paste("PCoA2:",soren_rel_eigen2_percent,"%")

pcoa_soren <- ggplot(soren_pcoa4, aes(Axis.1 * -1, Axis.2 * -1, color = quadrant, shape = trophicstate)) +
  xlab(soren_axis1) + ylab(soren_axis2) + #ggtitle("Bray-Curtis: Trophic State") +
  geom_point(size= 6, alpha=0.9) + theme_bw() + #ylim(1, -1) + xlim(1, -1) +
  scale_color_manual(name = "Habitat", breaks=c("Free Epilimnion", "Free Mixed",  "Free Hypolimnion", "Particle Epilimnion", "Particle Mixed", "Particle Hypolimnion"),
                     labels = c("Free-Living Epilimnion", "Free-Living Mixed",  "Free-Living Hypolimnion", "Particle-Associated Epilimnion", "Particle-Associated Mixed", "Particle-Associated Hypolimnion"), 
                     values = c("purple3", "mediumblue", "darkgreen", "orchid1", "deepskyblue", "green")) +
  scale_shape_manual(name = "Trophic State", breaks = c("Eutrophic", "Mesotrophic", "Oligotrophic", "Mixed"), 
                     labels = c("Eutrophic", "Mesotrophic", "Oligotrophic", "Mixed"),
                     values = c(15, 19, 17)) +
  theme(axis.text.x = element_text(colour="black", vjust=0.5, size=14), 
        axis.text.y = element_text(colour="black", vjust=0.5, size=14),
        axis.title.x = element_text(face="bold", size=16),
        axis.title.y = element_text(face="bold", size=16),
        legend.title = element_text(size=12, face="bold"),
        legend.text = element_text(size = 12),
        ###LEGEND TOP RIGHT CORNER
        legend.position = "right"); pcoa_soren

#jpeg(filename="~/Final_PAFL_Trophicstate/Figures/Fig.2_NMDS_bc+soren.jpeg", width= 45, height=18, units= "cm", pointsize= 14, res=500)
grid.newpage()
pushViewport(viewport(layout=grid.layout(1,2,width=c(0.4,0.6))))
print(pcoa_bray, vp=viewport(layout.pos.row=1,layout.pos.col=1))
print(pcoa_soren, vp=viewport(layout.pos.row=1,layout.pos.col=2))
#dev.off()


############## NMDS NMDS NMDS NMDS NDMS
set.seed(3)
nmds_bray <- metaMDS(nowinOTU, distance="bray")  #, autotransform = FALSE

nmds_bray4 <- metaMDS(nowinOTU, distance="bray", k = 4)  #, autotransform = FALSE

plot(nmds_bray4, choices = c(1,3))



nmds_bray <- data.frame(nmds_bray$points) #http://strata.uga.edu/software/pdf/mdsTutorial.pdf
nmds_bray$names<-row.names(nmds_bray) #new names column
nmds_bray <- makeCategories_dups(nmds_bray) #will add our categorical information:  lakenames, limnion, filter, quadrant and trophicstate
nmds_bray$ProdLevel <- as.character(nmds_bray$trophicstate)
nmds_bray$ProdLevel[nmds_bray$trophicstate == "Eutrophic"] <- "Productive"
nmds_bray$ProdLevel[nmds_bray$trophicstate == "Mesotrophic"] <- "Productive"
nmds_bray$ProdLevel[nmds_bray$trophicstate == "Oligotrophic"] <- "Unproductive"
nmds_bray$quadrant <- factor(nmds_bray$quadrant,levels = c("Free Epilimnion", "Free Mixed",  "Free Hypolimnion", "Particle Epilimnion", "Particle Mixed", "Particle Hypolimnion"))

nmds_bc_quad <- ggplot(nmds_bray, aes(MDS1*10, MDS2*10, color = quadrant, shape = ProdLevel)) +
  xlab("NMDS1") + ylab("NMDS2") + ggtitle("Bray-Curtis Dissimilarity") + 
  geom_point(size= 6, alpha=0.9) + theme_bw() +   geom_point(colour="white", size = 2) +
  annotate("text", label = "B", x = (min(nmds_bray$MDS1*10) + 0.05), y = (max(nmds_bray$MDS2*10) -0.02), face = "bold",size = 10, colour = "black") +
  annotate("text", label = "Stress = 0.17", x = (max(nmds_bray$MDS1*10) -0.2), y = (max(nmds_bray$MDS2*10) -0.02), size = 6, colour = "black") +
  scale_color_manual(name = "Filter Fraction and Lake Layer", breaks=c("Free Epilimnion", "Free Mixed",  "Free Hypolimnion", "Particle Epilimnion", "Particle Mixed", "Particle Hypolimnion"),
                     labels = c("Free-Living Epilimnion", "Free-Living Mixed",  "Free-Living Hypolimnion", "Particle-Associated Epilimnion", "Particle-Associated Mixed", "Particle-Associated Hypolimnion"), 
                     values = c("darkgreen", "mediumblue", "purple3", "green", "deepskyblue", "orchid1")) +
  scale_shape_manual(name = "Nutrient Level", breaks = c("Productive", "Unproductive"), 
                     labels = c("High", "Low"),
                     values = c(15, 17)) +
  guides(shape = guide_legend(order=1), color = guide_legend(order=2)) +
  scale_y_continuous(breaks=seq(-0.6, 0.6, 0.3), lim = c(-0.6,0.7)) +
  theme(axis.text.x = element_text(colour="black", vjust=0.5, size=14), 
        axis.text.y = element_text(colour="black", vjust=0.5, size=14),
        axis.title.x = element_text(face="bold", size=16),
        axis.title.y = element_text(face="bold", size=16),
        legend.title = element_text(size=12, face="bold"),
        legend.text = element_text(size = 12),
        plot.title = element_text(size = 16, face="bold"),
        ###LEGEND TOP RIGHT CORNER
        legend.position = "right");  nmds_bc_quad


# UNWEIGHTED
set.seed(3)
nmds_soren <- metaMDS(nowinOTU, distance="bray", binary = TRUE)
nmds_soren <- data.frame(nmds_soren$points) #http://strata.uga.edu/software/pdf/mdsTutorial.pdf
nmds_soren$names<-row.names(nmds_soren) #new names column
nmds_soren <- makeCategories_dups(nmds_soren) #will add our categorical information:  lakenames, limnion, filter, quadrant and trophicstate
nmds_soren$ProdLevel <- as.character(nmds_soren$trophicstate)
nmds_soren$ProdLevel[nmds_soren$trophicstate == "Eutrophic"] <- "Productive"
nmds_soren$ProdLevel[nmds_soren$trophicstate == "Mesotrophic"] <- "Productive"
nmds_soren$ProdLevel[nmds_soren$trophicstate == "Oligotrophic"] <- "Unproductive"
nmds_soren$quadrant <- factor(nmds_soren$quadrant,levels = c("Free Epilimnion", "Free Mixed",  "Free Hypolimnion", "Particle Epilimnion", "Particle Mixed", "Particle Hypolimnion"))

nmds_soren_quad <- ggplot(nmds_soren, aes(MDS1, MDS2, color = quadrant, shape = ProdLevel)) +
  xlab("NMDS1") + ylab("NMDS2") + ggtitle("Sørensen Dissimilarity") +
  geom_point(size= 6, alpha=0.9) + theme_bw() +   geom_point(colour="white", size = 2) +
  annotate("text", label = "A", x = (min(nmds_soren$MDS1) + 0.05), y = (max(nmds_soren$MDS2) -0.02), face = "bold",size = 10, colour = "black") +
  annotate("text", label = "Stress = 0.14", x = (max(nmds_soren$MDS1) -0.18), y = (max(nmds_soren$MDS2) -0.01), size = 6, colour = "black") +
  scale_color_manual(name = "Habitat", breaks=c("Free Epilimnion", "Free Mixed",  "Free Hypolimnion", "Particle Epilimnion", "Particle Mixed", "Particle Hypolimnion"),
                     labels = c("Free-Living Epilimnion", "Free-Living Mixed",  "Free-Living Hypolimnion", "Particle-Associated Epilimnion", "Particle-Associated Mixed", "Particle-Associated Hypolimnion"), 
                     values = c("darkgreen", "mediumblue", "purple3", "green", "deepskyblue", "orchid1")) +
  scale_shape_manual(name = "Nutrient Level", breaks = c("Productive", "Unproductive"), 
                     labels = c("High", "Low"),
                     values = c(15, 17)) +
  guides(shape = guide_legend(order=1), color = guide_legend(order=2)) +
  theme(axis.text.x = element_text(colour="black", vjust=0.5, size=14), 
        axis.text.y = element_text(colour="black", vjust=0.5, size=14),
        axis.title.x = element_text(face="bold", size=16),
        axis.title.y = element_text(face="bold", size=16),
        legend.title = element_text(size=12, face="bold"),
        legend.text = element_text(size = 12),
        plot.title = element_text(size = 16, face="bold"),
        ###LEGEND TOP RIGHT CORNER
        legend.position = "none");  nmds_soren_quad


#multiplot(nmds_bc_quad, nmds_soren_quad, cols = 2)

#####  Plotting FIGURE 3  #####  Plotting FIGURE 3  #####  Plotting FIGURE 3  #####  Plotting FIGURE 3  #####  Plotting FIGURE 3
#tiff(filename="~/Final_PAFL_Trophicstate/Final_Figures/Fig.3_NMDS_bc+soren_prod_small.tiff", width= 45, height=18, units= "cm", pointsize= 14, res=300, compression = "lzw")
grid.newpage()
pushViewport(viewport(layout=grid.layout(1,2,width=c(0.41,0.59))))
print(nmds_soren_quad, vp=viewport(layout.pos.row=1,layout.pos.col=1))
print(nmds_bc_quad, vp=viewport(layout.pos.row=1,layout.pos.col=2))
#dev.off()



####################################################################################  ADONIS  ####################################################################################
####################################################################################  ADONIS  ####################################################################################
####################################################################################  ADONIS  ####################################################################################
#  First we need to organize our environmental data
nowin_merged <- subset_samples(merged_final, names != "WINH" & names != "WINH3um")
nowin_merged <- prune_taxa(taxa_sums(nowin_merged) > 0, nowin_merged)
nosherwin_merged <- subset_samples(nowin_merged, lakenames != "Sherman")
nosherwin_merged <- prune_taxa(taxa_sums(nosherwin_merged) > 0, nosherwin_merged)  # 594 OTUs are unique to SHERMAN!

#nowin_merged
environ <- sample_data(nosherwin_merged)
environ <- subset(environ, select = c("names", "lakenames", "limnion", "filter", "trophicstate","ProdLevel", "totaldepth", "DO", "SpC", "temp", "pH")) # , "nitrate", "chla", "ammonia", "SRP", "TP", "Ncoord", "Wcoord"))
environ <- data.frame(environ)
## TO MAKE QUADRANT
for(i in 1:length(environ$limnion)){
  environ$quadrant[i]<-paste(as.character(environ$filter[i]),as.character(environ$limnion[i]))}
## TO MAKE HABITAT
for(i in 1:length(environ$limnion)){
  environ$habitat[i]<-paste(as.character(environ$ProdLevel[i]),as.character(environ$filter[i]),as.character(environ$limnion[i]))}


##  Calculate the BC dissimilarity and the Sorensen Dissimilarity
nosherwinOTU <- otu_table(nosherwin_merged)  # This is our OTU table that we will use for Adonis
set.seed(3)
#BCdist <- vegdist(nosherwinOTU, method = "bray", binary = FALSE)  # calculates the Bray-Curtis Distances

df_nosherwinOTU <- data.frame(nosherwinOTU)
otu_soren <- vegdist(df_nosherwinOTU, method = "bray", binary = TRUE)  ##SORENSEN DISTANCE --> Test's the presence/absence and makes 
BCdist <- otu_soren  ### This way we don't need to re-type the code for adonis, BUT - **BE CAREFUL** with this!  

##  Do you have the right object for BCdist (BC vs sorensen???)
##  Are you sure?!  

# #Run an ADONIS test!
adonis_PAFL_mult <- adonis(BCdist ~ limnion+filter +ProdLevel+ DO+temp + pH, data=environ) # R2 = 0.36245
adonis_PAFL_hab <- adonis(BCdist ~ habitat, data=environ) # R2 = 0.36245

adonis_PAFL_quad <- adonis(BCdist~quadrant,data=environ) #  R2 = 
adonis_PAFL_filt <- adonis(BCdist~filter,data=environ)  #R2 = 
adonis_PAFL_limnion <- adonis(BCdist~limnion,data=environ) #R2 = 
adonis_PAFL_DO <- adonis(BCdist~DO,data=environ) # R2 = 
adonis_PAFL_trophicstate <- adonis(BCdist~trophicstate, data = environ) # R2 = 
adonis_PAFL_prod <- adonis(BCdist~ProdLevel, data = environ) # R2 = 
adonis_PAFL_temp <- adonis(BCdist~ temp,data=environ) # R2 = 
adonis_PAFL_pH <- adonis(BCdist~ pH,data=environ) # R2 = 


# Particle-Associated!
part <- subset_samples(nosherwin_merged, filter == "Particle")
set.seed(3)
partOTU <- otu_table(part)
#part_BC <- vegdist(partOTU, method = "bray", binary = FALSE)  # BRAY CURTIS!

#df_partOTU <- data.frame(partOTU)
#part_soren <- vegdist(df_partOTU, method = "bray", binary = TRUE)  ##SORENSEN DISTANCE --> Test's the presence/absence and makes 
#part_BC <- part_soren  ### This way we don't need to re-type the code for adonis, BUT - **BE CAREFUL** with this!  

##  Do you have the right object for part_BC (BC vs sorensen???)
##  Are you sure?! 

part_env <- subset(environ, filter == "Particle")  #Load Environmental Data

part_adon_mult<- adonis(part_BC ~ limnion+ProdLevel+ DO + temp + pH, data=part_env) #R2 = 
part_adon_limn <- adonis(part_BC ~ limnion, data=part_env) # R2 = 
part_adon_troph <- adonis(part_BC ~ trophicstate, data=part_env) # R2 =  
part_adon_prod <- adonis(part_BC ~ ProdLevel, data=part_env) # R2 =  
part_adon_DO <- adonis(part_BC ~ DO, data=part_env) # R2 = 0.086
part_adon_temp <- adonis(part_BC ~ temp, data=part_env) # R2 = 0.086
part_adon_pH <- adonis(part_BC ~ pH, data=part_env) # R2 = 0.086



# FREE-LIVING!!!
free <- subset_samples(nosherwin_merged, filter == "Free")
set.seed(3)
freeOTU <- otu_table(free)
free_BC <- vegdist(freeOTU, method = "bray", binary = FALSE)  # BRAY CURTIS!

#df_freeOTU <- data.frame(freeOTU)
#free_soren <- vegdist(df_freeOTU, method = "bray", binary = TRUE)  ##SORENSEN DISTANCE --> Test's the presence/absence and makes 
#free_BC <- free_soren  ### This way we don't need to re-type the code for adonis, BUT - **BE CAREFUL** with this!  

##  Do you have the right object for free_BC (BC vs sorensen???)
##  Are you sure?! 

free_env <- subset(environ, filter == "Free")  #Load Environmental Data

free_adon_mult<- adonis(free_BC ~ limnion+ProdLevel+ DO + temp + pH, data=free_env) #R2 = 
free_adon_limn <- adonis(free_BC ~ limnion, data=free_env) # R2 = 
free_adon_troph <- adonis(free_BC ~ trophicstate, data=free_env) # R2 =  
free_adon_prod <- adonis(free_BC ~ ProdLevel, data=free_env) # R2 =  
free_adon_DO <- adonis(free_BC ~ DO, data=free_env) # R2 = 0.086
free_adon_temp <- adonis(free_BC ~ temp, data=free_env) # R2 = 0.086
free_adon_pH <- adonis(free_BC ~ pH, data=free_env) # R2 = 0.086



# Epilimnion!
epi <- subset_samples(nosherwin_merged, limnion == "Epilimnion")
set.seed(3)
epiOTU <- otu_table(epi)
epi_BC <- vegdist(epiOTU, method = "bray", binary = FALSE)  # BRAY CURTIS!

#df_epiOTU <- data.frame(epiOTU)
#epi_soren <- vegdist(df_epiOTU, method = "bray", binary = TRUE)  ##SORENSEN DISTANCE --> Test's the presence/absence and makes 
#epi_BC <- epi_soren  ### This way we don't need to re-type the code for adonis, BUT - **BE CAREFUL** with this!  

##  Do you have the right object for epi_BC (BC vs sorensen???)
##  Are you sure?! 

epi_env <- subset(environ, limnion == "Epilimnion")  #Load Environmental Data

epi_adon_mult<- adonis(epi_BC ~ filter+ProdLevel+ DO + temp + pH, data=epi_env) #R2 = 
epi_adon_filt <- adonis(epi_BC ~ filter, data=epi_env) # R2 = 
epi_adon_quad <- adonis(epi_BC ~ quadrant, data=epi_env) # R2 = 
epi_adon_troph <- adonis(epi_BC ~ trophicstate, data=epi_env) # R2 =  
epi_adon_prod <- adonis(epi_BC ~ ProdLevel, data=epi_env) # R2 =  
epi_adon_DO <- adonis(epi_BC ~ DO, data=epi_env) # R2 = 0.086
epi_adon_temp <- adonis(epi_BC ~ temp, data=epi_env) # R2 = 0.086
epi_adon_pH <- adonis(epi_BC ~ pH, data=epi_env) # R2 = 0.086


# hypolimnion!
hypo <- subset_samples(nosherwin_merged, limnion == "Hypolimnion")
hypoOTU <- otu_table(hypo)
set.seed(3)
hypo_BC <- vegdist(hypoOTU, method = "bray", binary = FALSE) # BRAY CURTIS!

#df_hypoOTU <- data.frame(hypoOTU)
#hypo_soren <- vegdist(df_hypoOTU, method = "bray", binary = TRUE)  ##SORENSEN DISTANCE --> Test's the presence/absence and makes 
#hypo_BC <- hypo_soren  ### This way we don't need to re-type the code for adonis, BUT - **BE CAREFUL** with this!  

##  Do you have the right object for hypo_BC (BC vs sorensen???)
##  Are you sure?! 

hypo_env <- subset(environ, limnion == "Hypolimnion") #Load Environmental Data

hypo_adon_mult <- adonis(hypo_BC ~ filter+ProdLevel+DO + temp + pH, data=hypo_env) # R2 = 
hypo_adon_filt <- adonis(hypo_BC ~ filter, data=hypo_env) # R2 = 
hypo_adon_quad <- adonis(hypo_BC ~ quadrant, data=hypo_env) # R2 = 
hypo_adon_troph <- adonis(hypo_BC ~ trophicstate, data=hypo_env) # R2 = 
hypo_adon_prod <- adonis(hypo_BC ~ ProdLevel, data=hypo_env) # R2 = 
hypo_adon_DO <- adonis(hypo_BC ~ DO, data=hypo_env) # R2 = 
hypo_adon_temp <- adonis(hypo_BC ~ temp, data=hypo_env) # R2 = 
hypo_adon_pH <- adonis(hypo_BC ~ pH, data=hypo_env) # R2 = 


# Oligotrophic!!
oligo <- subset_samples(nosherwin_merged, trophicstate == "Oligotrophic") 
oligoOTU <- otu_table(oligo)
set.seed(3)
#oligo_BC <- vegdist(oligoOTU, method = "bray", binary = FALSE) # BRAY CURTIS!

df_oligoOTU <- data.frame(oligoOTU)
oligo_soren <- vegdist(df_oligoOTU, method = "bray", binary = TRUE)  ##SORENSEN DISTANCE --> Test's the presence/absence and makes 
oligo_BC <- oligo_soren  ### This way we don't need to re-type the code for adonis, BUT - **BE CAREFUL** with this!  

##  Do you have the right object for oligo_BC (BC vs sorensen???)
##  Are you sure?! 

oligo_env <- subset(environ, trophicstate == "Oligotrophic")

oligo_adon_mult <- adonis(oligo_BC ~ limnion+filter+DO + temp + pH, data=oligo_env) # R2 = 
oligo_adon_filt <- adonis(oligo_BC ~ filter, data=oligo_env) # R2 = 
oligo_adon_quad <- adonis(oligo_BC ~ quadrant, data=oligo_env) # R2 = 
oligo_adon_limnion <- adonis(oligo_BC ~ limnion, data=oligo_env) # R2 = 
oligo_adon_DO <- adonis(oligo_BC ~ DO, data=oligo_env) # R2 = 
oligo_adon_temp <- adonis(oligo_BC ~ temp, data=oligo_env) # R2 = 
oligo_adon_pH <- adonis(oligo_BC ~ pH, data=oligo_env) # R2 = 


# MESO + EUTROPHIC!!
prod <- subset_samples(nosherwin_merged, ProdLevel == "Productive")
prodOTU <- otu_table(prod)
set.seed(3)
#prod_BC <- vegdist(prodOTU, method = "bray", binary = FALSE) # BRAY CURTIS DISTANCE

df_prodOTU <- data.frame(prodOTU)
prod_soren <- vegdist(df_prodOTU, method = "bray", binary = TRUE)  ##SORENSEN DISTANCE --> Test's the presence/absence and makes 
prod_BC <- prod_soren  ### This way we don't need to re-type the code for adonis, BUT - **BE CAREFUL** with this!  

##  Do you have the right object for prod_BC (BC vs sorensen???)
##  Are you sure?! 

prod_env <- subset(environ, ProdLevel == "Productive")

prod_adon_mult <- adonis(prod_BC ~ limnion+filter+DO + temp + pH, data=prod_env) # R2 = 
prod_adon_filt <- adonis(prod_BC ~ filter, data=prod_env) # R2 = 
prod_adon_quad <- adonis(prod_BC ~ quadrant, data=prod_env) # R2 = 
prod_adon_limnion <- adonis(prod_BC ~ limnion, data=prod_env) # R2 = 
prod_adon_DO <- adonis(prod_BC ~ DO, data=prod_env) # R2 = 
prod_adon_temp <- adonis(prod_BC ~ temp, data=prod_env) # R2 = 
prod_adon_pH <- adonis(prod_BC ~ pH, data=prod_env) # R2 = 


#  First we need to organize our environmental data
sherman_merged <- subset_samples(merged_final, lakenames == "Sherman")
sherman_merged <- prune_taxa(taxa_sums(sherman_merged) > 0, sherman_merged)

#nowin_merged
sherm_environ <- sample_data(sherman_merged)
sherm_environ <- subset(sherm_environ, select = c("names", "lakenames", "limnion", "filter", "trophicstate","ProdLevel", "totaldepth", "DO", "SpC", "temp", "pH")) # , "nitrate", "chla", "ammonia", "SRP", "TP", "Ncoord", "Wcoord"))
sherm_environ <- data.frame(sherm_environ)
sherm_environ$limnion <- as.character(sherm_environ$limnion)
sherm_environ$limnion[sherm_environ$limnion == "Mixed" & sherm_environ$names == "SHEE"] <- "Surface"
sherm_environ$limnion[sherm_environ$limnion == "Mixed" & sherm_environ$names == "SHEH"] <- "Bottom"
sherm_environ$limnion[sherm_environ$limnion == "Mixed" & sherm_environ$names == "SHEE3um"] <- "Surface"
sherm_environ$limnion[sherm_environ$limnion == "Mixed" & sherm_environ$names == "SHEH3um"] <- "Bottom"

## TO MAKE QUADRANT
for(i in 1:length(sherm_environ$limnion)){
  sherm_environ$quadrant[i]<-paste(as.character(sherm_environ$filter[i]),as.character(sherm_environ$limnion[i]))}


#########  PLOT ORDINATIONS FOR BOTH SCALED AND MANUAL
sherm_OTU <- otu_table(sherman_merged) # This is our OTU table that we will use for Adonis 
set.seed(3)
sherm_BCdist <- vegdist(sherm_OTU, method = "bray", binary = FALSE)

#df_sherm_OTU <- data.frame(sherm_OTU)
#sherm_soren <- vegdist(df_sherm_OTU, method = "bray", binary = TRUE)  ##SORENSEN DISTANCE --> Test's the presence/absence and makes 
#sherm_BCdist <- sherm_soren  ### This way we don't need to re-type the code for adonis, BUT - **BE CAREFUL** with this!  

##  Do you have the right object for sherm_BCdist (BC vs sorensen???)
##  Are you sure?! 

# 
# #Run an ADONIS test!
adonis_sherm_mult <- adonis(sherm_BCdist ~ filter+limnion+DO+temp+pH, data=sherm_environ, nperm = 999) # R2 = 0.36245
adonis_sherm_filt <- adonis(sherm_BCdist ~ filter, data=sherm_environ, nperm = 999) # R2 = 0.36245
adonis_sherm_limnion <- adonis(sherm_BCdist ~ limnion, data=sherm_environ, nperm = 999) # R2 = 0.36245
adonis_sherm_quad <- adonis(sherm_BCdist ~ quadrant, data=sherm_environ, nperm = 999) # R2 = 0.36245
adonis_sherm_temp <- adonis(sherm_BCdist ~ temp, data=sherm_environ, nperm = 999) # R2 = 0.36245
adonis_sherm_pH <- adonis(sherm_BCdist ~ pH, data=sherm_environ, nperm = 999) # R2 = 0.36245
adonis_sherm_DO <- adonis(sherm_BCdist ~ DO, data=sherm_environ, nperm = 999) # R2 = 0.36245





####################################################  DIFFERENCES IN COMMUNITY COMPOSITION  ####################################################  Working up to Figure 4 and Figure S2
####################################################  DIFFERENCES IN COMMUNITY COMPOSITION  ####################################################  Working up to Figure 4 and Figure S2
####################################################  DIFFERENCES IN COMMUNITY COMPOSITION  ####################################################  Working up to Figure 4 and FIgure S2

#BC Distance on Shared file = community composition
bray.dist <- norm_bray  # Bray Curtis dissimiliarity matrix for McMurdie & Holmes scaled, NO HYPO WINTERGREEN, samples

#http://stackoverflow.com/questions/23474729/convert-object-of-class-dist-into-data-frame-in-r
bray <- melt(as.matrix(bray.dist), varnames = c("samp1", "samp2"))
bray <- subset(bray, value > 0)

bray$lakenames1 <- substr(bray$samp1,1,3) # Create a new row called "lakenames" with first 3 letters of string
bray$lakenames2 <- substr(bray$samp2,1,3) # Create a new row called "lakenames" with first 3 letters of string
bray$limnion1 <- substr(bray$samp1, 4, 4) # Create a column called limnon with hypo or epi
bray$limnion2 <- substr(bray$samp2, 4, 4) # Create a column called limnon with hypo or epi
bray$filter1 <- substr(bray$samp1, 5, 7) 
bray$filter2 <- substr(bray$samp2, 5, 7) 


bray$lakenames1 <- as.character(bray$lakenames1)
bray$lakenames1[bray$lakenames1 == "WIN"] <- "Wintergreen"
bray$lakenames1[bray$lakenames1 == "SIX"] <- "Sixteen"
bray$lakenames1[bray$lakenames1 == "SHE"] <- "Sherman"
bray$lakenames1[bray$lakenames1 == "PAY"] <- "Payne"
bray$lakenames1[bray$lakenames1 == "LON"] <- "LittleLong"
bray$lakenames1[bray$lakenames1 == "LEE"] <- "Lee"
bray$lakenames1[bray$lakenames1 == "GUL"] <- "Gull"
bray$lakenames1[bray$lakenames1 == "BRI"] <- "Bristol"
bray$lakenames1[bray$lakenames1 == "BAK"] <- "Baker"
bray$lakenames1[bray$lakenames1 == "BAS"] <- "Baseline"
bray$lakenames1[bray$lakenames1 == "BST"] <- "Bassett"

bray$lakenames2 <- as.character(bray$lakenames2)
bray$lakenames2[bray$lakenames2 == "WIN"] <- "Wintergreen"
bray$lakenames2[bray$lakenames2 == "SIX"] <- "Sixteen"
bray$lakenames2[bray$lakenames2 == "SHE"] <- "Sherman"
bray$lakenames2[bray$lakenames2 == "PAY"] <- "Payne"
bray$lakenames2[bray$lakenames2 == "LON"] <- "LittleLong"
bray$lakenames2[bray$lakenames2 == "LEE"] <- "Lee"
bray$lakenames2[bray$lakenames2 == "GUL"] <- "Gull"
bray$lakenames2[bray$lakenames2 == "BRI"] <- "Bristol"
bray$lakenames2[bray$lakenames2 == "BAK"] <- "Baker"
bray$lakenames2[bray$lakenames2 == "BAS"] <- "Baseline"
bray$lakenames2[bray$lakenames2 == "BST"] <- "Bassett"



#Add Trophic State column by using the name of the lake
bray <- data.table(bray)
library(data.table)
bray[, trophicstate1 := ifelse(lakenames1 %in% c("Wintergreen", "Baker", "Baseline"), "Eutrophic",
                               ifelse(lakenames1 %in% c("Bassett", "Bristol", "Payne"), "Mesotrophic",
                                      ifelse(lakenames1 %in% c("Sherman"), "Mixed",
                                             ifelse(lakenames1 %in% c("Gull", "Sixteen", "LittleLong", "Lee"), "Oligotrophic", NA))))]
bray[, trophicstate2 := ifelse(lakenames2 %in% c("Wintergreen", "Baker", "Baseline"), "Eutrophic",
                               ifelse(lakenames2 %in% c("Bassett", "Bristol", "Payne"), "Mesotrophic",
                                      ifelse(lakenames2 %in% c("Sherman"), "Mixed",
                                             ifelse(lakenames2 %in% c("Gull", "Sixteen", "LittleLong", "Lee"), "Oligotrophic", NA))))]

bray$limnion1[bray$limnion1 == "E"] <- "Epilimnion"
bray$limnion1[bray$limnion1 == "H"] <- "Hypolimnion"
bray$limnion2[bray$limnion2 == "E"] <- "Epilimnion"
bray$limnion2[bray$limnion2 == "H"] <- "Hypolimnion"


###  Pull out the SHERMAN LAKE SAMPLES
bray$limnion1[bray$lakenames1 == "Sherman" & bray$limnion1 == "Epilimnion"] <- "Mixed"
bray$limnion1[bray$lakenames1 == "Sherman" & bray$limnion1 == "Hypolimnion"] <- "Mixed"
bray$limnion2[bray$lakenames2 == "Sherman" & bray$limnion2 == "Epilimnion"] <- "Mixed"
bray$limnion2[bray$lakenames2 == "Sherman" & bray$limnion2 == "Hypolimnion"] <- "Mixed"


bray$filter1[bray$filter1 == "3um"] <- "Particle"
bray$filter1[bray$filter1 == ""] <- "Free"
bray$filter2[bray$filter2 == "3um"] <- "Particle"
bray$filter2[bray$filter2 == ""] <- "Free"

#Add the combined lake layer
bray$combined<-rep(NA,length(bray$limnion1))
bray$combined[bray$limnion1=="Epilimnion"&bray$limnion2=="Epilimnion"]<-"Epilimnion"
bray$combined[bray$limnion1=="Hypolimnion"&bray$limnion2=="Hypolimnion"]<-"Hypolimnion"
bray$combined[bray$limnion1=="Hypolimnion"&bray$limnion2=="Epilimnion"]<-"EH"
bray$combined[bray$limnion1=="Epilimnion"&bray$limnion2=="Hypolimnion"]<-"EH"

bray$combined[bray$limnion1=="Mixed"&bray$limnion2=="Mixed"]<-"Mixed"
bray$combined[bray$limnion1=="Mixed"&bray$limnion2=="Hypolimnion"]<-"Mixed Hypolimnion"
bray$combined[bray$limnion1=="Hypolimnion"&bray$limnion2=="Mixed"]<-"Mixed Hypolimnion"
bray$combined[bray$limnion1=="Epilimnion"&bray$limnion2=="Mixed"]<-"Mixed Epilimnion"
bray$combined[bray$limnion1=="Mixed"&bray$limnion2=="Epilimnion"]<-"Mixed Epilimnion"


##  Add the commbined filter.
bray$filt_comb<-rep(NA,length(bray$filter1))
bray$filt_comb[bray$filter1=="Free"&bray$filter2=="Free"]<-"Free"
bray$filt_comb[bray$filter1=="Particle"&bray$filter2=="Particle"]<-"Particle"
bray$filt_comb[bray$filter1=="Particle"&bray$filter2=="Free"]<-"PF"
bray$filt_comb[bray$filter1=="Free"&bray$filter2=="Particle"]<-"PF"

#creating a column for each combination of top bottom particle free
for(i in 1:length(bray$limnion1)){
  bray$cat[i]<-paste(as.character(bray$filt_comb[i]),as.character(bray$combined[i]))}


#  Add for trophicstate combintions
bray$troph_comb<-rep(NA,length(bray$lakenames2))
bray$troph_comb[bray$trophicstate1=="Eutrophic" & bray$trophicstate2=="Eutrophic"]<-"Eutrophic"
bray$troph_comb[bray$trophicstate1=="Eutrophic" & bray$trophicstate2=="Mesotrophic"]<-"Eutrophic-Mesotrophic"
bray$troph_comb[bray$trophicstate1=="Mesotrophic" & bray$trophicstate2=="Eutrophic"]<-"Eutrophic-Mesotrophic"
bray$troph_comb[bray$trophicstate1=="Eutrophic" & bray$trophicstate2=="Oligotrophic"]<-"Eutrophic-Oligotrophic"
bray$troph_comb[bray$trophicstate1=="Oligotrophic" & bray$trophicstate2=="Eutrophic"]<-"Eutrophic-Oligotrophic"
bray$troph_comb[bray$trophicstate1=="Eutrophic" & bray$trophicstate2=="Mixed"]<-"Eutrophic-Mixed"
bray$troph_comb[bray$trophicstate1=="Mixed" & bray$trophicstate2=="Eutrophic"]<-"Eutrophic-Mixed"

bray$troph_comb[bray$trophicstate1=="Mesotrophic" & bray$trophicstate2=="Mesotrophic"]<-"Mesotrophic"
bray$troph_comb[bray$trophicstate1=="Mesotrophic" & bray$trophicstate2=="Oligotrophic"]<-"Mesotrophic-Oligotrophic"
bray$troph_comb[bray$trophicstate1=="Oligotrophic" & bray$trophicstate2=="Mesotrophic"]<-"Mesotrophic-Oligotrophic"
bray$troph_comb[bray$trophicstate1=="Mesotrophic" & bray$trophicstate2=="Mixed"]<-"Mesotrophic-Mixed"
bray$troph_comb[bray$trophicstate1=="Mixed" & bray$trophicstate2=="Mesotrophic"]<-"Mesotrophic-Mixed"

bray$troph_comb[bray$trophicstate1=="Oligotrophic" & bray$trophicstate2=="Oligotrophic"]<-"Oligotrophic"
bray$troph_comb[bray$trophicstate1=="Mixed" & bray$trophicstate2=="Oligotrophic"]<-"Oligotrophic-Mixed"
bray$troph_comb[bray$trophicstate1=="Oligotrophic" & bray$trophicstate2=="Mixed"]<-"Oligotrophic-Mixed"

bray$troph_comb[bray$trophicstate1=="Mixed" & bray$trophicstate2=="Mixed"]<-"Mixed"

##  Add combined trophicstate, filter, and limnion 
for(i in 1:length(bray$limnion1)){
  bray$troph_lim1[i]<-paste(as.character(bray$trophicstate1[i]),
                            as.character(bray$limnion1[i]),
                            as.character(bray$filter1[i]))}
for(i in 1:length(bray$limnion2)){
  bray$troph_lim2[i]<-paste(as.character(bray$trophicstate2[i]),
                            as.character(bray$limnion2[i]),
                            as.character(bray$filter2[i]))}

bray$samp1 <- as.character(bray$samp1)
bray$samp2 <- as.character(bray$samp2)

#creating a column for each combination of top bottom particle free
for(i in 1:length(bray$limnion1)){
  bray$cat[i]<-paste(as.character(bray$filt_comb[i]),as.character(bray$combined[i]))}

###  Subsetting out the samples that are in the same categories as each other
beta <- subset(bray, troph_lim1 == troph_lim2)

### Put everything in the order we would like.
beta$trophicstate1 <-factor(beta$trophicstate1,levels=c("Eutrophic", "Mesotrophic", "Oligotrophic", "Mixed"))
beta$trophicstate2 <-factor(beta$trophicstate2,levels=c("Eutrophic", "Mesotrophic", "Oligotrophic", "Mixed"))
beta$troph_lim2 <-factor(beta$troph_lim2,levels=c("Eutrophic Epilimnion Particle", "Eutrophic Epilimnion Free", "Eutrophic Hypolimnion Particle", "Eutrophic Hypolimnion Free",
                                                  "Mesotrophic Epilimnion Particle", "Mesotrophic Epilimnion Free", "Mesotrophic Hypolimnion Particle", "Mesotrophic Hypolimnion Free",
                                                  "Oligotrophic Epilimnion Particle", "Oligotrophic Epilimnion Free", "Oligotrophic Hypolimnion Particle", "Oligotrophic Hypolimnion Free",
                                                  "Mixed Mixed Particle", "Mixed Mixed Free"))

#### Subsetting out the mixed lake!
nomix_beta <- subset(beta, trophicstate1 != "Mixed")
nomix_beta2 <- subset(nomix_beta, trophicstate2 != "Mixed")


#######  Creating geom_point plot for Beta Diverstiy.
# STATS ON BETA
ddply_beta <- ddply(nomix_beta2, c("troph_lim1", "trophicstate1", "filter1"), summarise, 
                    N = length(value),
                    mean = mean(value),
                    sd   = sd(value),
                    se   = sd / sqrt(N))



####  Run the test on ALL the data that goes into the mean!
hist(nomix_beta2$value, breaks = 30)  # Not normally distributed!!!
nomix_beta2$value <- as.numeric(nomix_beta2$value)
nomix_beta2$troph_lim1 <- as.factor(nomix_beta2$troph_lim1)
nomix_bray_beta <- nomix_beta2
## Do the KW test
bray_nomix_KW <- kruskal.test(nomix_bray_beta$value ~ nomix_bray_beta$troph_lim1) # Kruskal Wallis test on braysen!
print(bray_nomix_KW)  # show Kruskal Wallis result
### Which samples are significantly different from each other?  Significant???  YES! WOOOHOOOOOOO!
bray_nomix_KW_MC <- kruskalmc(nomix_bray_beta$value ~ nomix_bray_beta$troph_lim1)  ## Defaults to P < 0.05
print(bray_nomix_KW_MC)
### Time to figure out letters to represent significance in a plot!
nomix_bray_test <- bray_nomix_KW_MC$dif.com$difference # select logical vector
names(nomix_bray_test) <- row.names(bray_nomix_KW_MC$dif.com) # add comparison names
# create a list with "homogenous groups" coded by letter
nomix_bray_letters <- multcompLetters(nomix_bray_test, compare="<", threshold=0.05, Letters=c(letters, LETTERS, "."), reversed = FALSE)#['Letters']
###  Let's extract the values from the multcompLetters object
nomix_bray_sigs_dataframe <-  data.frame(as.vector(names(nomix_bray_letters$Letters)), as.vector(nomix_bray_letters$Letters))
colnames(nomix_bray_sigs_dataframe) <- c("troph_lim1", "siglabel")
nomix_bray_try <- merge(ddply_beta, nomix_bray_sigs_dataframe, by = "troph_lim1")


#tiff(filename="~/Final_PAFL_Trophicstate/Final_Figures/Fig.S2_beta_TROPH_SD.tiff", width= 35, height= 20, units= "cm", pointsize= 14, res=200)
beta_plot <- ggplot(nomix_bray_try, aes(x = troph_lim1, y = mean, color = troph_lim1)) + geom_point(size = 5) +
  geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=.2, position=position_dodge(.9)) +
  geom_text(aes(label = siglabel, x = troph_lim1, y = ((mean+sd) + 0.03)), size = 5) +
  scale_color_manual(name = "", limits=c("Eutrophic Epilimnion Particle", "Eutrophic Epilimnion Free", "Eutrophic Hypolimnion Particle", "Eutrophic Hypolimnion Free",
                                         "Mesotrophic Epilimnion Particle", "Mesotrophic Epilimnion Free", "Mesotrophic Hypolimnion Particle", "Mesotrophic Hypolimnion Free",
                                         "Oligotrophic Epilimnion Particle", "Oligotrophic Epilimnion Free", "Oligotrophic Hypolimnion Particle", "Oligotrophic Hypolimnion Free"), 
                     values = c("firebrick", "firebrick", "firebrick", "firebrick","darkorange","darkorange","darkorange","darkorange", 
                                "cornflowerblue","cornflowerblue","cornflowerblue","cornflowerblue"))+
  scale_x_discrete(breaks=c("Eutrophic Epilimnion Particle", "Eutrophic Epilimnion Free", "Eutrophic Hypolimnion Particle", "Eutrophic Hypolimnion Free",
                            "Mesotrophic Epilimnion Particle", "Mesotrophic Epilimnion Free", "Mesotrophic Hypolimnion Particle", "Mesotrophic Hypolimnion Free",
                            "Oligotrophic Epilimnion Particle", "Oligotrophic Epilimnion Free", "Oligotrophic Hypolimnion Particle", "Oligotrophic Hypolimnion Free"),
                   labels=c("Epilimnion \nParticle-Associated", "Epilimnion \nFree-Living", "Hypolimnion \nParticle-Associated", "Hypolimnion \nFree-Living",
                            "Epilimnion \nParticle-Associated", "Epilimnion \nFree-Living", "Hypolimnion \nParticle-Associated", "Hypolimnion \nFree-Living",
                            "Epilimnion \nParticle-Associated", "Epilimnion \nFree-Living", "Hypolimnion \nParticle-Associated", "Hypolimnion \nFree-Living")) + 
  xlab("Habitat") + ylab("Bray-Curtis Dissimilarity") + theme_bw() + #scale_fill_brewer(palette="Paired") + 
  facet_grid(. ~ trophicstate1, scale = "free", space = "free") +
  theme(axis.title.x = element_text(face="bold", size=16),
        axis.text.x = element_text(angle=30, colour = "black", vjust=1, hjust = 1, size=14),
        axis.text.y = element_text(colour = "black", size=14),
        axis.title.y = element_text(face="bold", size=16),
        plot.title = element_text(face="bold", size = 20),
        strip.text.x = element_text(size=16, face = "bold", colour = "black"),
        strip.background = element_blank(),
        plot.margin = unit(c(0.2, 0.2, 0.2, 0.5), "cm"), #top, right, bottom, left)
        legend.position="none"); beta_plot
#dev.off()



###################################  NSF DDIG PLOT!
ddig_bray_try <- nomix_bray_try
ddig_bray_try$trophicstate1 <- as.character(ddig_bray_try$trophicstate1)
ddig_bray_try$trophicstate1[ddig_bray_try$trophicstate1 == "Eutrophic"] <- "High-Nutrient"
ddig_bray_try$trophicstate1[ddig_bray_try$trophicstate1 == "Mesotrophic"] <- "Mid-Nutrient"
ddig_bray_try$trophicstate1[ddig_bray_try$trophicstate1 == "Oligotrophic"] <- "Low-Nutrient"
ddig_bray_try$trophicstate1 <-factor(ddig_bray_try$trophicstate1,levels=c("High-Nutrient", "Mid-Nutrient", "Low-Nutrient"))


DDIG_betaplot <-  ggplot(ddig_bray_try, aes(x = troph_lim1, y = mean, color = troph_lim1)) + geom_point(size = 3) +
  geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=.2, position=position_dodge(.9)) +
  geom_text(aes(label = siglabel, x = troph_lim1, y = ((mean+sd) + 0.035)), size = 3) +
  scale_color_manual(name = "", limits=c("Eutrophic Epilimnion Particle", "Eutrophic Epilimnion Free", "Eutrophic Hypolimnion Particle", "Eutrophic Hypolimnion Free",
                                         "Mesotrophic Epilimnion Particle", "Mesotrophic Epilimnion Free", "Mesotrophic Hypolimnion Particle", "Mesotrophic Hypolimnion Free",
                                         "Oligotrophic Epilimnion Particle", "Oligotrophic Epilimnion Free", "Oligotrophic Hypolimnion Particle", "Oligotrophic Hypolimnion Free"), 
                     values = c("firebrick", "firebrick", "firebrick", "firebrick","darkorange","darkorange","darkorange","darkorange", 
                                "cornflowerblue","cornflowerblue","cornflowerblue","cornflowerblue"))+
  scale_x_discrete(breaks=c("Eutrophic Epilimnion Particle", "Eutrophic Epilimnion Free", "Eutrophic Hypolimnion Particle", "Eutrophic Hypolimnion Free",
                            "Mesotrophic Epilimnion Particle", "Mesotrophic Epilimnion Free", "Mesotrophic Hypolimnion Particle", "Mesotrophic Hypolimnion Free",
                            "Oligotrophic Epilimnion Particle", "Oligotrophic Epilimnion Free", "Oligotrophic Hypolimnion Particle", "Oligotrophic Hypolimnion Free"),
                   labels=c("Epilimnion \nParticle-Associated", "Epilimnion \nFree-Living", "Hypolimnion \nParticle-Associated", "Hypolimnion \nFree-Living",
                            "Epilimnion \nParticle-Associated", "Epilimnion \nFree-Living", "Hypolimnion \nParticle-Associated", "Hypolimnion \nFree-Living",
                            "Epilimnion \nParticle-Associated", "Epilimnion \nFree-Living", "Hypolimnion \nParticle-Associated", "Hypolimnion \nFree-Living")) + 
  xlab("Habitat") + ylab("Bray-Curtis Dissimilarity") + theme_bw() + #scale_fill_brewer(palette="Paired") + 
  facet_grid(. ~ trophicstate1, scale = "free", space = "free") +
theme(axis.title.x = element_text(face="bold", size=10),  #Set the x-axis title
      axis.title.y = element_text(face="bold", size=10, vjust=0.5),  #Set the y-axis title
      axis.text.x = element_text(colour = "black", size=8, angle = 30, hjust = 1, vjust = 1),  #Set the x-axis labels
      axis.text.y = element_text(colour = "black", size=8),  #Set the y-axis labels
      legend.title = element_text(size=7, face="bold"),  #Set the legend title 
      legend.text = element_text(size = 7),
      strip.text.x = element_text(size=8, face="bold"),  #Set the facet titles on x-axis 
      strip.text.y = element_text(size=8, face="bold"),  #Set the facet titles on x-axis 
      legend.position = "none",#"left",
      plot.margin = unit(c(0.1, 1, 0.1, 0.1), "cm"),
      strip.text.x = element_text(size=14, face = "bold", colour = "black"),
      strip.background = element_blank()); DDIG_betaplot
ggsave(DDIG_betaplot, filename = "~/Final_PAFL_Trophicstate/bray-trophic.pdf", height = 3, width = 6.5, dpi = 300)



##  MORE DDIG -- SUBSET OUT ONLY THE LOW NUTRIENT LAKES: 
oligo_beta <- subset(nomix_beta2, trophicstate1 == "Oligotrophic")
####  Run the test on ALL the data that goes into the mean!
hist(oligo_beta$value, breaks = 60)  # Not normally distributed!!!
oligo_beta$value <- as.numeric(oligo_beta$value)
oligo_beta$troph_lim1 <- as.factor(oligo_beta$troph_lim1)
oligo_beta_bray <- oligo_beta
## Do the KW test
oligo_KW <- kruskal.test(oligo_beta_bray$value ~ oligo_beta_bray$troph_lim1) # Kruskal Wallis test on braysen!
print(oligo_KW)  # show Kruskal Wallis result
### Which samples are significantly different from each other?  Significant???  YES! WOOOHOOOOOOO!
oligo_KW_MC <- kruskalmc(oligo_beta_bray$value ~ oligo_beta_bray$troph_lim1)  ## Defaults to P < 0.05
oligo_KW_MC <- subset(oligo_KW_MC$dif.com, difference != "NA")
print(oligo_KW_MC)
### Time to figure out letters to represent significance in a plot!
oligo_KW_test <- oligo_KW_MC$difference # select logical vector
names(oligo_KW_test) <- row.names(oligo_KW_MC) # add comparison names
# create a list with "homogenous groups" coded by letter
oligo_KW_letters <- multcompLetters(oligo_KW_test, compare="<", threshold=0.05, Letters=c(letters, LETTERS, "."), reversed = FALSE)#['Letters']
###  Let's extract the values from the multcompLetters object
oligo_bray_sigs_dataframe <-  data.frame(as.vector(names(oligo_KW_letters$Letters)), as.vector(oligo_KW_letters$Letters))
colnames(oligo_bray_sigs_dataframe) <- c("troph_lim1", "siglabel")
## Make oligotrophic df with ddply
oligo_ddply <- subset(ddply_beta, trophicstate1 == "Oligotrophic")
oligo_bray_try <- merge(oligo_ddply, oligo_bray_sigs_dataframe, by = "troph_lim1")


oligo_bray_try$trophicstate1 <- as.character(oligo_bray_try$trophicstate1)
oligo_bray_try$trophicstate1[oligo_bray_try$trophicstate1 == "Eutrophic"] <- "High-Nutrient"
oligo_bray_try$trophicstate1[oligo_bray_try$trophicstate1 == "Mesotrophic"] <- "Mid-Nutrient"
oligo_bray_try$trophicstate1[oligo_bray_try$trophicstate1 == "Oligotrophic"] <- "Low-Nutrient"
oligo_bray_try$trophicstate1 <-factor(oligo_bray_try$trophicstate1,levels=c("High-Nutrient", "Mid-Nutrient", "Low-Nutrient"))
nlabel <- c("n=12", "n=6", "n=12", "n=12")
oligo_bray_try$nlabel <- nlabel

oligo_beta_plot <-  ggplot(oligo_bray_try, aes(x = troph_lim1, y = mean, color = troph_lim1)) + geom_point(size = 3) +
  geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=.2, position=position_dodge(.9)) +
  geom_text(aes(label = siglabel, x = troph_lim1, y = ((mean+sd) + 0.035)), size = 3) +
  geom_text(aes(label = nlabel, x = troph_lim1, y = ((mean+sd) + 0.065)), size = 3) +
  scale_color_manual(name = "", limits=c("Eutrophic Epilimnion Particle", "Eutrophic Epilimnion Free", "Eutrophic Hypolimnion Particle", "Eutrophic Hypolimnion Free",
                                         "Mesotrophic Epilimnion Particle", "Mesotrophic Epilimnion Free", "Mesotrophic Hypolimnion Particle", "Mesotrophic Hypolimnion Free",
                                         "Oligotrophic Epilimnion Particle", "Oligotrophic Epilimnion Free", "Oligotrophic Hypolimnion Particle", "Oligotrophic Hypolimnion Free"), 
                     values = c("firebrick", "firebrick", "firebrick", "firebrick","darkorange","darkorange","darkorange","darkorange", 
                                "cornflowerblue","cornflowerblue","cornflowerblue","cornflowerblue"))+
  scale_x_discrete(breaks=c("Eutrophic Epilimnion Particle", "Eutrophic Epilimnion Free", "Eutrophic Hypolimnion Particle", "Eutrophic Hypolimnion Free",
                            "Mesotrophic Epilimnion Particle", "Mesotrophic Epilimnion Free", "Mesotrophic Hypolimnion Particle", "Mesotrophic Hypolimnion Free",
                            "Oligotrophic Epilimnion Particle", "Oligotrophic Epilimnion Free", "Oligotrophic Hypolimnion Particle", "Oligotrophic Hypolimnion Free"),
                   labels=c("Epilimnion \nParticle-Associated", "Epilimnion \nFree-Living", "Hypolimnion \nParticle-Associated", "Hypolimnion \nFree-Living",
                            "Epilimnion \nParticle-Associated", "Epilimnion \nFree-Living", "Hypolimnion \nParticle-Associated", "Hypolimnion \nFree-Living",
                            "Epilimnion \nParticle-Associated", "Epilimnion \nFree-Living", "Hypolimnion \nParticle-Associated", "Hypolimnion \nFree-Living")) + 
  xlab("Habitat") + ylab("Bray-Curtis Dissimilarity") + theme_bw() + #scale_fill_brewer(palette="Paired") + 
  facet_grid(. ~ trophicstate1, scale = "free", space = "free") +
  theme(axis.title.x = element_text(face="bold", size=10),  #Set the x-axis title
        axis.title.y = element_text(face="bold", size=10, vjust=0.5),  #Set the y-axis title
        axis.text.x = element_text(colour = "black", size=8, angle = 30, hjust = 1, vjust = 1),  #Set the x-axis labels
        axis.text.y = element_text(colour = "black", size=8),  #Set the y-axis labels
        legend.title = element_text(size=7, face="bold"),  #Set the legend title 
        legend.text = element_text(size = 7),
        strip.text.x = element_text(size=8, face="bold"),  #Set the facet titles on x-axis 
        strip.text.y = element_text(size=8, face="bold"),  #Set the facet titles on x-axis 
        legend.position = "none",#"left",
        plot.margin = unit(c(0.1, 1, 0.1, 0.1), "cm"),
        strip.text.x = element_text(size=14, face = "bold", colour = "black"),
        strip.background = element_blank()); oligo_beta_plot
ggsave(oligo_beta_plot, filename = "~/Final_PAFL_Trophicstate/bray-oligo.pdf", height = 3, width = 3, dpi = 300)










################# BRAY CURTIS DISTANCE FOR PRODUCTIVE VS UNPRODUCTIVE (High VS Low Nutrients )
################# BRAY CURTIS DISTANCE FOR PRODUCTIVE VS UNPRODUCTIVE (High VS Low Nutrients )
################# BRAY CURTIS DISTANCE FOR PRODUCTIVE VS UNPRODUCTIVE (High VS Low Nutrients )
prod_beta <- nomix_beta2
unique(prod_beta$troph_lim1)
prod_beta$troph_lim1 <- as.character(prod_beta$troph_lim1)
prod_beta$troph_lim1[prod_beta$troph_lim1 == "Eutrophic Epilimnion Free" | prod_beta$troph_lim1 =="Mesotrophic Epilimnion Free"] <- "Productive Epilimnion Free"
prod_beta$troph_lim1[prod_beta$troph_lim1 == "Eutrophic Epilimnion Particle" | prod_beta$troph_lim1 =="Mesotrophic Epilimnion Particle"] <- "Productive Epilimnion Particle"
prod_beta$troph_lim1[prod_beta$troph_lim1 == "Eutrophic Hypolimnion Free" | prod_beta$troph_lim1 =="Mesotrophic Hypolimnion Free"] <- "Productive Hypolimnion Free"
prod_beta$troph_lim1[prod_beta$troph_lim1 == "Eutrophic Hypolimnion Particle" | prod_beta$troph_lim1 =="Mesotrophic Hypolimnion Particle"] <- "Productive Hypolimnion Particle"
prod_beta$troph_lim1[prod_beta$troph_lim1 == "Oligotrophic Epilimnion Free"] <- "Unproductive Epilimnion Free"
prod_beta$troph_lim1[prod_beta$troph_lim1 == "Oligotrophic Epilimnion Particle"] <- "Unproductive Epilimnion Particle"
prod_beta$troph_lim1[prod_beta$troph_lim1 == "Oligotrophic Hypolimnion Free"] <- "Unproductive Hypolimnion Free"
prod_beta$troph_lim1[prod_beta$troph_lim1 == "Oligotrophic Hypolimnion Particle"] <- "Unproductive Hypolimnion Particle"
prod_beta$trophicstate1 <- as.character(prod_beta$trophicstate1)
prod_beta$trophicstate1[prod_beta$trophicstate1 == "Eutrophic"] <- "Productive"
prod_beta$trophicstate1[prod_beta$trophicstate1 == "Mesotrophic"] <- "Productive"
prod_beta$trophicstate1[prod_beta$trophicstate1 == "Oligotrophic"] <- "Unproductive"


####  Run the test on ALL the data that goes into the mean!
hist(prod_beta$value, breaks = 30)  # Not normally distributed!!!
prod_beta$value <- as.numeric(prod_beta$value)
prod_beta$troph_lim1 <- as.factor(prod_beta$troph_lim1)
## Significant???  YES!
prodKW <- kruskal.test(prod_beta$value ~ prod_beta$troph_lim1) 
### Which samples are significantly different from each other?
KW_bray_samps <- kruskalmc(prod_beta$value ~ prod_beta$troph_lim1)
KW_bray_samps_sigs <- subset(KW_bray_samps$dif.com, difference == TRUE)



#############  Time to calculate the mean and standard deviation for all of them!
ddply_prodbeta <- ddply(prod_beta, c("troph_lim1", "trophicstate1"), summarise, 
                        N = length(value),
                        mean = mean(value),
                        sd   = sd(value),
                        se   = sd / sqrt(N))

ddply_prodbeta$troph_lim1 <-factor(ddply_prodbeta$troph_lim1,levels=c("Productive Epilimnion Particle", "Productive Epilimnion Free", "Productive Hypolimnion Particle", "Productive Hypolimnion Free",
                                                                      "Unproductive Epilimnion Particle", "Unproductive Epilimnion Free", "Unproductive Hypolimnion Particle", "Unproductive Hypolimnion Free"))


#jpeg(filename="~/Final_PAFL_Trophicstate/Figures/Fig.3c_beta_PROD_SD.jpeg", width= 25, height=15, units= "cm", pointsize= 14, res=500)
prodbeta_plot <- ggplot(ddply_prodbeta, aes(x = troph_lim1, y = mean, color = troph_lim1)) + geom_point(size = 5) +
  geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=.2, position=position_dodge(.9)) +
  scale_color_manual(name = "", limits=c("Productive Epilimnion Particle", "Productive Epilimnion Free", "Productive Hypolimnion Particle", "Productive Hypolimnion Free",
                                         "Unproductive Epilimnion Particle", "Unproductive Epilimnion Free", "Unproductive Hypolimnion Particle", "Unproductive Hypolimnion Free"), 
                     values = c("black", "black", "black", "black",
                                "gray48","gray48","gray48","gray48"))+
  scale_x_discrete(breaks=c("Productive Epilimnion Particle", "Productive Epilimnion Free", "Productive Hypolimnion Particle", "Productive Hypolimnion Free",
                            "Unproductive Epilimnion Particle", "Unproductive Epilimnion Free", "Unproductive Hypolimnion Particle", "Unproductive Hypolimnion Free"),
                   labels=c("Epilimnion Particle", "Epilimnion Free", "Hypolimnion Particle", "Hypolimnion Free",
                            "Epilimnion Particle", "Epilimnion Free", "Hypolimnion Particle", "Hypolimnion Free")) + 
  xlab("Habitat") + ylab("Bray Curtis Dissimilarity") + theme_bw() + #scale_fill_brewer(palette="Paired") + 
  facet_grid(. ~ trophicstate1, scale = "free", space = "free") +
  theme(axis.title.x = element_text(face="bold", size=14),
        axis.text.x = element_text(angle=30, colour = "black", vjust=1, hjust = 1, size=14),
        axis.text.y = element_text(colour = "black", size=14),
        axis.title.y = element_text(face="bold", size=16),
        plot.title = element_text(face="bold", size = 20),
        strip.text.x = element_text(size=12, face = "bold", colour = "black"),
        strip.background = element_blank(),
        legend.position="none");   prodbeta_plot
#dev.off()


######  Test for significance 
#oligo_bray <- subset(ddply_prodbeta, trophicstate1 == "Unproductive")
#oligo_KW_bray <- kruskal.test(mean ~ troph_lim1, data = ddply_prodbeta)  #Still NOT significant!
prod_bray <- subset(ddply_prodbeta, trophicstate1 == "Productive")
prod_KW_bray <- kruskal.test(mean ~ troph_lim1, data = ddply_prodbeta)  #Still NOT significant!
#  Bray-Curtis 
KW_bray <- kruskal.test(mean ~ troph_lim1, data = ddply_prodbeta)

### Which samples are significantly different from each other?
KW_bray_samps <- kruskalmc(ddply_prodbeta$mean, ddply_prodbeta$troph_lim1)
KW_bray_samps_sigs <- subset(KW_bray_samps$dif.com, difference == TRUE)



####################################################  SORENSEN
####################################################  SORENSEN
####################################################  SORENSEN

#SOREN Distance on Shared file = community composition
nowin_merged <- subset_samples(merged_final, names != "WINH" & names != "WINH3um")
nowin_merged <- prune_taxa(taxa_sums(nowin_merged) > 0, nowin_merged)
nowinOTU <- otu_table(nowin_merged)
nowinOTU_df <- data.frame(otu_table(nowin_merged))  
norm_soren <- vegdist(nowinOTU_df, method = "bray", binary = TRUE)  # SORENSEN INDEX
soren.dist <- norm_soren  # soren Curtis dissimiliarity matrix for McMurdie & Holmes scaled, NO HYPO WINTERGREEN, samples

#http://stackoverflow.com/questions/23474729/convert-object-of-class-dist-into-data-frame-in-r
soren <- melt(as.matrix(soren.dist), varnames = c("samp1", "samp2"))
soren <- subset(soren, value > 0)

soren$lakenames1 <- substr(soren$samp1,1,3) # Create a new row called "lakenames" with first 3 letters of string
soren$lakenames2 <- substr(soren$samp2,1,3) # Create a new row called "lakenames" with first 3 letters of string
soren$limnion1 <- substr(soren$samp1, 4, 4) # Create a column called limnon with hypo or epi
soren$limnion2 <- substr(soren$samp2, 4, 4) # Create a column called limnon with hypo or epi
soren$filter1 <- substr(soren$samp1, 5, 7) 
soren$filter2 <- substr(soren$samp2, 5, 7) 


soren$lakenames1 <- as.character(soren$lakenames1)
soren$lakenames1[soren$lakenames1 == "WIN"] <- "Wintergreen"
soren$lakenames1[soren$lakenames1 == "SIX"] <- "Sixteen"
soren$lakenames1[soren$lakenames1 == "SHE"] <- "Sherman"
soren$lakenames1[soren$lakenames1 == "PAY"] <- "Payne"
soren$lakenames1[soren$lakenames1 == "LON"] <- "LittleLong"
soren$lakenames1[soren$lakenames1 == "LEE"] <- "Lee"
soren$lakenames1[soren$lakenames1 == "GUL"] <- "Gull"
soren$lakenames1[soren$lakenames1 == "BRI"] <- "Bristol"
soren$lakenames1[soren$lakenames1 == "BAK"] <- "Baker"
soren$lakenames1[soren$lakenames1 == "BAS"] <- "Baseline"
soren$lakenames1[soren$lakenames1 == "BST"] <- "Bassett"

soren$lakenames2 <- as.character(soren$lakenames2)
soren$lakenames2[soren$lakenames2 == "WIN"] <- "Wintergreen"
soren$lakenames2[soren$lakenames2 == "SIX"] <- "Sixteen"
soren$lakenames2[soren$lakenames2 == "SHE"] <- "Sherman"
soren$lakenames2[soren$lakenames2 == "PAY"] <- "Payne"
soren$lakenames2[soren$lakenames2 == "LON"] <- "LittleLong"
soren$lakenames2[soren$lakenames2 == "LEE"] <- "Lee"
soren$lakenames2[soren$lakenames2 == "GUL"] <- "Gull"
soren$lakenames2[soren$lakenames2 == "BRI"] <- "Bristol"
soren$lakenames2[soren$lakenames2 == "BAK"] <- "Baker"
soren$lakenames2[soren$lakenames2 == "BAS"] <- "Baseline"
soren$lakenames2[soren$lakenames2 == "BST"] <- "Bassett"



#Add Trophic State column by using the name of the lake
soren <- data.table(soren)
library(data.table)
soren[, trophicstate1 := ifelse(lakenames1 %in% c("Wintergreen", "Baker", "Baseline"), "Eutrophic",
                                ifelse(lakenames1 %in% c("Bassett", "Bristol", "Payne"), "Mesotrophic",
                                       ifelse(lakenames1 %in% c("Sherman"), "Mixed",
                                              ifelse(lakenames1 %in% c("Gull", "Sixteen", "LittleLong", "Lee"), "Oligotrophic", NA))))]
soren[, trophicstate2 := ifelse(lakenames2 %in% c("Wintergreen", "Baker", "Baseline"), "Eutrophic",
                                ifelse(lakenames2 %in% c("Bassett", "Bristol", "Payne"), "Mesotrophic",
                                       ifelse(lakenames2 %in% c("Sherman"), "Mixed",
                                              ifelse(lakenames2 %in% c("Gull", "Sixteen", "LittleLong", "Lee"), "Oligotrophic", NA))))]

soren$limnion1[soren$limnion1 == "E"] <- "Epilimnion"
soren$limnion1[soren$limnion1 == "H"] <- "Hypolimnion"
soren$limnion2[soren$limnion2 == "E"] <- "Epilimnion"
soren$limnion2[soren$limnion2 == "H"] <- "Hypolimnion"


###  Pull out the SHERMAN LAKE SAMPLES
soren$limnion1[soren$lakenames1 == "Sherman" & soren$limnion1 == "Epilimnion"] <- "Mixed"
soren$limnion1[soren$lakenames1 == "Sherman" & soren$limnion1 == "Hypolimnion"] <- "Mixed"
soren$limnion2[soren$lakenames2 == "Sherman" & soren$limnion2 == "Epilimnion"] <- "Mixed"
soren$limnion2[soren$lakenames2 == "Sherman" & soren$limnion2 == "Hypolimnion"] <- "Mixed"


soren$filter1[soren$filter1 == "3um"] <- "Particle"
soren$filter1[soren$filter1 == ""] <- "Free"
soren$filter2[soren$filter2 == "3um"] <- "Particle"
soren$filter2[soren$filter2 == ""] <- "Free"

#Add the combined lake layer
soren$combined<-rep(NA,length(soren$limnion1))
soren$combined[soren$limnion1=="Epilimnion"&soren$limnion2=="Epilimnion"]<-"Epilimnion"
soren$combined[soren$limnion1=="Hypolimnion"&soren$limnion2=="Hypolimnion"]<-"Hypolimnion"
soren$combined[soren$limnion1=="Hypolimnion"&soren$limnion2=="Epilimnion"]<-"EH"
soren$combined[soren$limnion1=="Epilimnion"&soren$limnion2=="Hypolimnion"]<-"EH"

soren$combined[soren$limnion1=="Mixed"&soren$limnion2=="Mixed"]<-"Mixed"
soren$combined[soren$limnion1=="Mixed"&soren$limnion2=="Hypolimnion"]<-"Mixed Hypolimnion"
soren$combined[soren$limnion1=="Hypolimnion"&soren$limnion2=="Mixed"]<-"Mixed Hypolimnion"
soren$combined[soren$limnion1=="Epilimnion"&soren$limnion2=="Mixed"]<-"Mixed Epilimnion"
soren$combined[soren$limnion1=="Mixed"&soren$limnion2=="Epilimnion"]<-"Mixed Epilimnion"


##  Add the commbined filter.
soren$filt_comb<-rep(NA,length(soren$filter1))
soren$filt_comb[soren$filter1=="Free"&soren$filter2=="Free"]<-"Free"
soren$filt_comb[soren$filter1=="Particle"&soren$filter2=="Particle"]<-"Particle"
soren$filt_comb[soren$filter1=="Particle"&soren$filter2=="Free"]<-"PF"
soren$filt_comb[soren$filter1=="Free"&soren$filter2=="Particle"]<-"PF"

#creating a column for each combination of top bottom particle free
for(i in 1:length(soren$limnion1)){
  soren$cat[i]<-paste(as.character(soren$filt_comb[i]),as.character(soren$combined[i]))}


#  Add for trophicstate combintions
soren$troph_comb<-rep(NA,length(soren$lakenames2))
soren$troph_comb[soren$trophicstate1=="Eutrophic" & soren$trophicstate2=="Eutrophic"]<-"Eutrophic"
soren$troph_comb[soren$trophicstate1=="Eutrophic" & soren$trophicstate2=="Mesotrophic"]<-"Eutrophic-Mesotrophic"
soren$troph_comb[soren$trophicstate1=="Mesotrophic" & soren$trophicstate2=="Eutrophic"]<-"Eutrophic-Mesotrophic"
soren$troph_comb[soren$trophicstate1=="Eutrophic" & soren$trophicstate2=="Oligotrophic"]<-"Eutrophic-Oligotrophic"
soren$troph_comb[soren$trophicstate1=="Oligotrophic" & soren$trophicstate2=="Eutrophic"]<-"Eutrophic-Oligotrophic"
soren$troph_comb[soren$trophicstate1=="Eutrophic" & soren$trophicstate2=="Mixed"]<-"Eutrophic-Mixed"
soren$troph_comb[soren$trophicstate1=="Mixed" & soren$trophicstate2=="Eutrophic"]<-"Eutrophic-Mixed"

soren$troph_comb[soren$trophicstate1=="Mesotrophic" & soren$trophicstate2=="Mesotrophic"]<-"Mesotrophic"
soren$troph_comb[soren$trophicstate1=="Mesotrophic" & soren$trophicstate2=="Oligotrophic"]<-"Mesotrophic-Oligotrophic"
soren$troph_comb[soren$trophicstate1=="Oligotrophic" & soren$trophicstate2=="Mesotrophic"]<-"Mesotrophic-Oligotrophic"
soren$troph_comb[soren$trophicstate1=="Mesotrophic" & soren$trophicstate2=="Mixed"]<-"Mesotrophic-Mixed"
soren$troph_comb[soren$trophicstate1=="Mixed" & soren$trophicstate2=="Mesotrophic"]<-"Mesotrophic-Mixed"

soren$troph_comb[soren$trophicstate1=="Oligotrophic" & soren$trophicstate2=="Oligotrophic"]<-"Oligotrophic"
soren$troph_comb[soren$trophicstate1=="Mixed" & soren$trophicstate2=="Oligotrophic"]<-"Oligotrophic-Mixed"
soren$troph_comb[soren$trophicstate1=="Oligotrophic" & soren$trophicstate2=="Mixed"]<-"Oligotrophic-Mixed"

soren$troph_comb[soren$trophicstate1=="Mixed" & soren$trophicstate2=="Mixed"]<-"Mixed"

##  Add combined trophicstate, filter, and limnion 
for(i in 1:length(soren$limnion1)){
  soren$troph_lim1[i]<-paste(as.character(soren$trophicstate1[i]),
                             as.character(soren$limnion1[i]),
                             as.character(soren$filter1[i]))}
for(i in 1:length(soren$limnion2)){
  soren$troph_lim2[i]<-paste(as.character(soren$trophicstate2[i]),
                             as.character(soren$limnion2[i]),
                             as.character(soren$filter2[i]))}

soren$samp1 <- as.character(soren$samp1)
soren$samp2 <- as.character(soren$samp2)

#creating a column for each combination of top bottom particle free
for(i in 1:length(soren$limnion1)){
  soren$cat[i]<-paste(as.character(soren$filt_comb[i]),as.character(soren$combined[i]))}

###  Subsetting out the samples that are in the same categories as each other
soren_beta <- subset(soren, troph_lim1 == troph_lim2)

### Put everything in the order we would like.
soren_beta$trophicstate1 <-factor(soren_beta$trophicstate1,levels=c("Eutrophic", "Mesotrophic", "Oligotrophic", "Mixed"))
soren_beta$trophicstate2 <-factor(soren_beta$trophicstate2,levels=c("Eutrophic", "Mesotrophic", "Oligotrophic", "Mixed"))
soren_beta$troph_lim2 <-factor(soren_beta$troph_lim2,levels=c("Eutrophic Epilimnion Particle", "Eutrophic Epilimnion Free", "Eutrophic Hypolimnion Particle", "Eutrophic Hypolimnion Free",
                                                              "Mesotrophic Epilimnion Particle", "Mesotrophic Epilimnion Free", "Mesotrophic Hypolimnion Particle", "Mesotrophic Hypolimnion Free",
                                                              "Oligotrophic Epilimnion Particle", "Oligotrophic Epilimnion Free", "Oligotrophic Hypolimnion Particle", "Oligotrophic Hypolimnion Free",
                                                              "Mixed Mixed Particle", "Mixed Mixed Free"))

#### Subsetting out the mixed lake!
nomix_soren_beta <- subset(soren_beta, trophicstate1 != "Mixed")
nomix_soren_beta2 <- subset(nomix_soren_beta, trophicstate2 != "Mixed")


#######  Creating geom_point plot for Beta Diverstiy.
# STATS ON BETA
ddply_soren_beta <- ddply(nomix_soren_beta2, c("troph_lim1", "trophicstate1", "filter1"), summarise, 
                          N = length(value),
                          mean = mean(value),
                          sd   = sd(value),
                          se   = sd / sqrt(N))


####  Run the test on ALL the data that goes into the mean!
hist(nomix_soren_beta2$value, breaks = 30)  # Not normally distributed!!!
nomix_soren_beta2$value <- as.numeric(nomix_soren_beta2$value)
nomix_soren_beta2$troph_lim1 <- as.factor(nomix_soren_beta2$troph_lim1)
nomix_soren_beta3 <- nomix_soren_beta2
## Do the KW test
soren_nomix_KW <- kruskal.test(nomix_soren_beta3$value ~ nomix_soren_beta3$troph_lim1) # Kruskal Wallis test on sorensen!
print(soren_nomix_KW)  # show Kruskal Wallis result
### Which samples are significantly different from each other?  Significant???  YES! WOOOHOOOOOOO!
soren_nomix_KW_MC <- kruskalmc(nomix_soren_beta3$value ~ nomix_soren_beta3$troph_lim1)  ## Defaults to P < 0.05
print(soren_nomix_KW_MC)
### Time to figure out letters to represent significance in a plot!
nomix_soren_test <- soren_nomix_KW_MC$dif.com$difference # select logical vector
names(nomix_soren_test) <- row.names(soren_nomix_KW_MC$dif.com) # add comparison names
# create a list with "homogenous groups" coded by letter
nomix_soren_letters <- multcompLetters(nomix_soren_test, compare="<", threshold=0.05, Letters=c(letters, LETTERS, "."), reversed = FALSE)#['Letters']
###  Let's extract the values from the multcompLetters object
nomix_soren_sigs_dataframe <-  data.frame(as.vector(names(nomix_soren_letters$Letters)), as.vector(nomix_soren_letters$Letters))
colnames(nomix_soren_sigs_dataframe) <- c("troph_lim1", "siglabel")
nomix_soren_try <- merge(ddply_soren_beta, nomix_soren_sigs_dataframe, by = "troph_lim1")



#jpeg(filename="~/Final_PAFL_Trophicstate/Figures/Fig.3c_beta_TROPH_SD.jpeg", width= 25, height=15, units= "cm", pointsize= 14, res=500)
soren_beta_plot <- ggplot(nomix_soren_try, aes(x = troph_lim1, y = mean, color = troph_lim1)) + geom_point(size = 5) +
  geom_text(aes(label = siglabel, x = troph_lim1, y = ((mean+sd) + 0.03)), size = 5) +
  geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=.2, position=position_dodge(.9)) +
  scale_color_manual(name = "", limits=c("Eutrophic Epilimnion Particle", "Eutrophic Epilimnion Free", "Eutrophic Hypolimnion Particle", "Eutrophic Hypolimnion Free",
                                         "Mesotrophic Epilimnion Particle", "Mesotrophic Epilimnion Free", "Mesotrophic Hypolimnion Particle", "Mesotrophic Hypolimnion Free",
                                         "Oligotrophic Epilimnion Particle", "Oligotrophic Epilimnion Free", "Oligotrophic Hypolimnion Particle", "Oligotrophic Hypolimnion Free"), 
                     values = c("firebrick", "firebrick", "firebrick", "firebrick","darkorange","darkorange","darkorange","darkorange", 
                                "cornflowerblue","cornflowerblue","cornflowerblue","cornflowerblue"))+
  scale_x_discrete(breaks=c("Eutrophic Epilimnion Particle", "Eutrophic Epilimnion Free", "Eutrophic Hypolimnion Particle", "Eutrophic Hypolimnion Free",
                            "Mesotrophic Epilimnion Particle", "Mesotrophic Epilimnion Free", "Mesotrophic Hypolimnion Particle", "Mesotrophic Hypolimnion Free",
                            "Oligotrophic Epilimnion Particle", "Oligotrophic Epilimnion Free", "Oligotrophic Hypolimnion Particle", "Oligotrophic Hypolimnion Free"),
                   labels=c("Epilimnion Particle", "Epilimnion Free", "Hypolimnion Particle", "Hypolimnion Free",
                            "Epilimnion Particle", "Epilimnion Free", "Hypolimnion Particle", "Hypolimnion Free",
                            "Epilimnion Particle", "Epilimnion Free", "Hypolimnion Particle", "Hypolimnion Free")) + 
  xlab("Habitat") + ylab("Sørensen\nDissimilarity") + theme_bw() + #scale_fill_brewer(palette="Paired") + 
  facet_grid(. ~ trophicstate1, scale = "free", space = "free") +
  theme(axis.title.x = element_text(face="bold", size=14),
        axis.text.x = element_text(angle=30, colour = "black", vjust=1, hjust = 1, size=14),
        axis.text.y = element_text(colour = "black", size=14),
        axis.title.y = element_text(face="bold", size=16),
        plot.title = element_text(face="bold", size = 20),
        plot.margin = unit(c(0.2, 0.2, 0.2, 0.5), "cm"), #top, right, bottom, left)
        strip.text.x = element_text(size=16, face = "bold", colour = "black"),
        strip.background = element_blank(),
        legend.position="none"); soren_beta_plot
#dev.off()



soren_beta_plot2 <- soren_beta_plot +   theme(axis.title.x = element_blank(), axis.ticks.x = element_blank(),
                                              axis.text.x = element_blank(),  plot.title = element_text(face="bold", size = 20),
                                              strip.background = element_blank(),
                                              plot.margin = unit(c(0.2, 0.2, -0.37, 0.22), "cm"), #top, right, bottom, left)
                                              legend.position="none")


beta_plot2 <- beta_plot + ylab("Bray-Curtis \nDissimilarity") + theme(strip.background = element_blank(), 
                                                                      plot.margin = unit(c(-0.37, 0.2, 0.2, 0.5), "cm"), #top, right, bottom, left)
                                                                      strip.text.x = element_blank())

#####  Plotting FIGURE S2  #####  Plotting FIGURE S2  #####  Plotting FIGURE S2  #####  Plotting FIGURE S2  #####  Plotting FIGURE S2
#tiff(filename="~/Final_PAFL_Trophicstate/Final_Figures/Fig.S2_BC+soren_beta_trophicstate_SIGS.tiff", width= 33, height=25, units= "cm", pointsize= 14, res=200)
grid.newpage()
pushViewport(viewport(layout=grid.layout(2,1,height=c(0.44,0.56))))
print(soren_beta_plot2, vp=viewport(layout.pos.row=1,layout.pos.col=1))  
print(beta_plot2, vp=viewport(layout.pos.row=2,layout.pos.col=1))  ##  See 1491 for beta_plot
#dev.off()



######  Test for significance 
#  soren-Curtis 
KW_soren_troph <- kruskal.test(mean ~ troph_lim2, data = ddply_beta)

### Which samples are significantly different from each other?
KW_soren_samps_troph <- kruskalmc(ddply_beta$mean, ddply_beta$troph_lim2)
KW_soren_samps_troph_sigs <- subset(KW_soren_samps_troph$dif.com, difference == TRUE)


prod_soren_beta <- nomix_soren_beta2
unique(prod_soren_beta$troph_lim1)
prod_soren_beta$troph_lim1[prod_soren_beta$troph_lim1 == "Eutrophic Epilimnion Free" | prod_soren_beta$troph_lim1 =="Mesotrophic Epilimnion Free"] <- "Productive Epilimnion Free"
prod_soren_beta$troph_lim1[prod_soren_beta$troph_lim1 == "Eutrophic Epilimnion Particle" | prod_soren_beta$troph_lim1 =="Mesotrophic Epilimnion Particle"] <- "Productive Epilimnion Particle"
prod_soren_beta$troph_lim1[prod_soren_beta$troph_lim1 == "Eutrophic Hypolimnion Free" | prod_soren_beta$troph_lim1 =="Mesotrophic Hypolimnion Free"] <- "Productive Hypolimnion Free"
prod_soren_beta$troph_lim1[prod_soren_beta$troph_lim1 == "Eutrophic Hypolimnion Particle" | prod_soren_beta$troph_lim1 =="Mesotrophic Hypolimnion Particle"] <- "Productive Hypolimnion Particle"
prod_soren_beta$troph_lim1[prod_soren_beta$troph_lim1 == "Oligotrophic Epilimnion Free"] <- "Unproductive Epilimnion Free"
prod_soren_beta$troph_lim1[prod_soren_beta$troph_lim1 == "Oligotrophic Epilimnion Particle"] <- "Unproductive Epilimnion Particle"
prod_soren_beta$troph_lim1[prod_soren_beta$troph_lim1 == "Oligotrophic Hypolimnion Free"] <- "Unproductive Hypolimnion Free"
prod_soren_beta$troph_lim1[prod_soren_beta$troph_lim1 == "Oligotrophic Hypolimnion Particle"] <- "Unproductive Hypolimnion Particle"
prod_soren_beta$trophicstate1 <- as.character(prod_soren_beta$trophicstate1)
prod_soren_beta$trophicstate1[prod_soren_beta$trophicstate1 == "Eutrophic"] <- "High-Nutrient"
prod_soren_beta$trophicstate1[prod_soren_beta$trophicstate1 == "Mesotrophic"] <- "High-Nutrient"
prod_soren_beta$trophicstate1[prod_soren_beta$trophicstate1 == "Oligotrophic"] <- "Low-Nutrient"



####################################################################  STATS TIME!  SORENSEN
####################################################################  STATS TIME!  SORENSEN
####################################################################  STATS TIME!  SORENSEN
####  Run the test on ALL the data that goes into the mean!  #### Check it out here:  https://aquaticr.wordpress.com/2012/12/18/multiple-comparison-test-for-non-parametric-data/
hist(prod_soren_beta$value, breaks = 30)  # Not normally distributed!!!
prod_soren_beta$value <- as.numeric(prod_soren_beta$value)
prod_soren_beta$troph_lim1 <- as.factor(prod_soren_beta$troph_lim1)
## Do the KW test
soren_prod_KW <- kruskal.test(prod_soren_beta$value ~ prod_soren_beta$troph_lim1) # Kruskal Wallis test on sorensen!
print(soren_prod_KW)  # show Kruskal Wallis result
### Which samples are significantly different from each other?  Significant???  YES! WOOOHOOOOOOO!
soren_prod_KW_MC <- kruskalmc(prod_soren_beta$value ~ prod_soren_beta$troph_lim1)  ## Defaults to P < 0.05
print(soren_prod_KW_MC)
### Time to figure out letters to represent significance in a plot!
soren_test <- soren_prod_KW_MC$dif.com$difference # select logical vector
names(soren_test) <- row.names(soren_prod_KW_MC$dif.com) # add comparison names
# create a list with "homogenous groups" coded by letter
soren_letters <- multcompLetters(soren_test, compare="<", threshold=0.05, Letters=c(letters, LETTERS, "."), reversed = FALSE)#['Letters']
###  Let's extract the values from the multcompLetters object
soren_sigs_dataframe <-  data.frame(as.vector(names(soren_letters$Letters)), as.vector(soren_letters$Letters))
colnames(soren_sigs_dataframe) <- c("troph_lim1", "siglabel")
soren_try <- merge(ddply_prodbeta_soren, soren_sigs_dataframe)


####################################################################  STATS TIME!  BRAY CURTIS
####################################################################  STATS TIME!  BRAY CURTIS
####################################################################  STATS TIME!  BRAY CURTIS
####  Run the test on ALL the data that goes into the mean!
hist(prod_beta$value, breaks = 30)  # Not normally distributed!!!
prod_beta$value <- as.numeric(prod_beta$value)
prod_beta$troph_lim1 <- as.factor(prod_beta$troph_lim1)
prod_bray_beta <- prod_beta
## Do the KW test
bray_prod_KW <- kruskal.test(prod_bray_beta$value ~ prod_bray_beta$troph_lim1) # Kruskal Wallis test on braysen!
print(bray_prod_KW)  # show Kruskal Wallis result
### Which samples are significantly different from each other?  Significant???  YES! WOOOHOOOOOOO!
bray_prod_KW_MC <- kruskalmc(prod_bray_beta$value ~ prod_bray_beta$troph_lim1)  ## Defaults to P < 0.05
print(bray_prod_KW_MC)
### Time to figure out letters to represent significance in a plot!
bray_test <- bray_prod_KW_MC$dif.com$difference # select logical vector
names(bray_test) <- row.names(bray_prod_KW_MC$dif.com) # add comparison names
# create a list with "homogenous groups" coded by letter
bray_letters <- multcompLetters(bray_test, compare="<", threshold=0.05, Letters=c(letters, LETTERS, "."), reversed = FALSE)#['Letters']
###  Let's extract the values from the multcompLetters object
bray_sigs_dataframe <-  data.frame(as.vector(names(bray_letters$Letters)), as.vector(bray_letters$Letters))
colnames(bray_sigs_dataframe) <- c("troph_lim1", "siglabel")
bray_try <- merge(ddply_prodbeta, bray_sigs_dataframe)


#############  Time to calculate the mean and standard deviation for all of them!
ddply_prodbeta_soren <- ddply(prod_soren_beta, c("troph_lim1", "trophicstate1"), summarise, 
                              N = length(value),
                              mean = mean(value),
                              sd   = sd(value),
                              se   = sd / sqrt(N))

ddply_prodbeta_soren$troph_lim1 <-factor(ddply_prodbeta_soren$troph_lim1,levels=c("Productive Epilimnion Particle", "Productive Epilimnion Free", "Productive Hypolimnion Particle", "Productive Hypolimnion Free",
                                                                                  "Unproductive Epilimnion Particle", "Unproductive Epilimnion Free", "Unproductive Hypolimnion Particle", "Unproductive Hypolimnion Free"))


#jpeg(filename="~/Final_PAFL_Trophicstate/Figures/Fig.3c_beta_PROD_SD.jpeg", width= 25, height=15, units= "cm", pointsize= 14, res=500)
soren_prodbeta_plot <- ggplot(soren_try, aes(x = troph_lim1, y = mean, color = troph_lim1)) + geom_point(size = 5) +
  geom_text(aes(label = siglabel, x = troph_lim1, y = ((mean+sd) + 0.04)), size =5) +
  geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=.2, position=position_dodge(.9)) +
  scale_color_manual(name = "", limits=c("Productive Epilimnion Particle", "Productive Epilimnion Free", "Productive Hypolimnion Particle", "Productive Hypolimnion Free",
                                         "Unproductive Epilimnion Particle", "Unproductive Epilimnion Free", "Unproductive Hypolimnion Particle", "Unproductive Hypolimnion Free"), 
                     values = c("black", "black", "black", "black",
                                "gray48","gray48","gray48","gray48"))+
  scale_x_discrete(breaks=c("Productive Epilimnion Particle", "Productive Epilimnion Free", "Productive Hypolimnion Particle", "Productive Hypolimnion Free",
                            "Unproductive Epilimnion Particle", "Unproductive Epilimnion Free", "Unproductive Hypolimnion Particle", "Unproductive Hypolimnion Free"),
                   labels=c("Epilimnion Particle", "Epilimnion Free", "Hypolimnion Particle", "Hypolimnion Free",
                            "Epilimnion Particle", "Epilimnion Free", "Hypolimnion Particle", "Hypolimnion Free")) + 
  xlab("Habitat") + ylab("Sørensen\nDissimilarity") + theme_bw() + #scale_fill_brewer(palette="Paired") + 
  scale_y_continuous(breaks=seq(0.4, 0.8, 0.1), lim = c(0.4, 0.83)) + 
  facet_grid(. ~ trophicstate1, scale = "free", space = "free") +
  theme(axis.title.x = element_blank(), axis.ticks.x = element_blank(), 
        axis.text.x = element_blank(),
        axis.text.y = element_text(colour = "black", size=14),
        axis.title.y = element_text(face="bold", size=16),
        plot.margin = unit(c(0.1, 0.1, 0.2, 0.1), "cm"), #top, right, bottom, left
        plot.title = element_text(face="bold", size = 20),
        strip.text.x = element_text(size=16, face="bold"),  #Set the facet titles on x-axis 
        strip.background = element_blank(),
        legend.position="none");   soren_prodbeta_plot
#dev.off()



prod_beta_plot <- ggplot(bray_try, aes(x = troph_lim1, y = mean, color = troph_lim1)) + geom_point(size = 5) +
  geom_text(aes(label = siglabel, x = troph_lim1, y = ((mean+sd) + 0.03)), size = 5) +
  geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=.2, position=position_dodge(.9)) + 
  scale_color_manual(name = "", limits=c("Productive Epilimnion Particle", "Productive Epilimnion Free", "Productive Hypolimnion Particle", "Productive Hypolimnion Free",
                                         "Unproductive Epilimnion Particle", "Unproductive Epilimnion Free", "Unproductive Hypolimnion Particle", "Unproductive Hypolimnion Free"), 
                     values = c("black", "black", "black", "black",
                                "gray48","gray48","gray48","gray48"))+
  scale_x_discrete(breaks=c("Productive Epilimnion Particle", "Productive Epilimnion Free", "Productive Hypolimnion Particle", "Productive Hypolimnion Free",
                            "Unproductive Epilimnion Particle", "Unproductive Epilimnion Free", "Unproductive Hypolimnion Particle", "Unproductive Hypolimnion Free"),
                   labels=c("Epilimnion \nParticle-Associated", "Epilimnion \nFree-Living", "Hypolimnion \nParticle-Associated", "Hypolimnion \nFree-Living",
                            "Epilimnion \nParticle-Associated", "Epilimnion \nFree-Living", "Hypolimnion \nParticle-Associated", "Hypolimnion \nFree-Living")) + 
  xlab("Habitat") + ylab("Bray-Curtis \nDissimilarity") + theme_bw() +  #scale_fill_brewer(palette="Paired") + 
  facet_grid(. ~ trophicstate1, scale = "free", space = "free") +
  scale_y_continuous(breaks=seq(0.4, 0.8, 0.1), lim = c(0.4, 0.84)) + 
  theme(axis.title.x = element_text(face="bold", size=16),
        axis.text.x = element_text(angle=30, colour = "black", vjust=1, hjust = 1, size=14),
        axis.text.y = element_text(colour = "black", size=14),
        axis.title.y = element_text(face="bold", size=16),
        plot.title = element_text(face="bold", size = 20),
        #axis.ticks.x = element_blank(),
        plot.margin = unit(c(-0.9, 0.1, 0.1, 0.1), "cm"), #top, right, bottom, left
        strip.background = element_blank(), strip.text = element_blank(),
        legend.position="none");prod_beta_plot


#####  Plotting FIGURE 4  #####  Plotting FIGURE 4  #####  Plotting FIGURE 4  #####  Plotting FIGURE 4  #####  Plotting FIGURE 4
#tiff(filename="~/Final_PAFL_Trophicstate/Final_Figures/Fig.4_BC+soren_beta_SIGS.tiff", width= 26, height=22, units= "cm", pointsize= 14, res=200)
grid.newpage()
pushViewport(viewport(layout=grid.layout(2,1,height=c(0.46,0.54))))
print(soren_prodbeta_plot, vp=viewport(layout.pos.row=1,layout.pos.col=1))
print(prod_beta_plot, vp=viewport(layout.pos.row=2,layout.pos.col=1))
#dev.off()






####################################################  PHYLUM LEVEL LOG-2-FOLD RATIO ANALYSIS  ####################################################  Working up to FIGURE 5
####################################################  PHYLUM LEVEL LOG-2-FOLD RATIO ANALYSIS  ####################################################  Working up to FIGURE 5
####################################################  PHYLUM LEVEL LOG-2-FOLD RATIO ANALYSIS  ####################################################  Working up to FIGURE 5
#View(data.frame(tax_table(merged_final)))  #Sanity check
good_phylum_proteo <-tax_glom(merged_final,taxrank = "Phylum")

### Subsetting Sherman lake for differences between particle and free living 
phy_sherm <- subset_samples(good_phylum_proteo, lakenames == "Sherman")  
phy_shemeister <- deSEQ(phy_sherm, ~ filter)
physherm_plot <- plot_phylum_deSEQ(phy_shemeister, "Sherman: Phylum-Level")

########## SUBSET OUT SHERMAN AND WINTERGREEN HYPOLIMNION
good_phylum_nosher <- subset_samples(good_phylum_proteo, lakenames != "Sherman")
good_phylum_nosherwin <- subset_samples(good_phylum_nosher, names != "WINH" & names != "WINH3um")
#View(data.frame(sample_data(good_phylum_nosherwin))) #Sanity check

############################
###########################################################PA VS FL
#1. Top prod 
phytopProd <- subset_samples(good_phylum_nosherwin, ProdLevel == "Productive" & limnion == "Epilimnion") 
de_phytopProd<- deSEQ(phytopProd, ~ filter)
phytopProd_plot <- plot_phylum_deSEQ(de_phytopProd, "PA vs FL:  Surface Productive (Phylum)")

#2 Top Oligo 
phytopUNProd <- subset_samples(good_phylum_nosherwin, ProdLevel == "Unproductive" & limnion == "Epilimnion") 
de_phytopUNProd <- deSEQ(phytopUNProd, ~ filter)
phytopUNProd_plot <- plot_phylum_deSEQ(de_phytopUNProd, "PA vs FL:  Surface Unproductive (Phylum)")

#3 Bottom Productive 
phybotProd <- subset_samples(good_phylum_nosherwin, ProdLevel == "Productive" & limnion == "Hypolimnion")
de_phybotProd <- deSEQ(phybotProd, ~ filter)
phybotProd_plot <- plot_phylum_deSEQ(de_phybotProd, "PA vs FL:  Bottom Productive (Phylum)")

#4 Bottom Unproductive 
phybotUNProd <- subset_samples(good_phylum_nosherwin, ProdLevel == "Unproductive" & limnion == "Hypolimnion") 
de_phybotUNProd <- deSEQ(phybotUNProd, ~ filter)
phybotUNProd_plot <- plot_phylum_deSEQ(de_phybotUNProd, "PA vs FL:  Bottom Unproductive (Phylum)")

pafl_phylum <- multiplot(phytopProd_plot, phybotProd_plot, phytopUNProd_plot, phybotUNProd_plot,cols = 2)


############################
###########################################################TOP VS BOTTOM
#1. PA prod 
phyprodPA <- subset_samples(good_phylum_nosherwin, ProdLevel == "Productive" & filter == "Particle") 
de_phyprodPA <- deSEQ(phyprodPA, ~ limnion)
phyprodPA_plot <- plot_phylum_deSEQ(de_phyprodPA, "Surface vs Bottom:  PA Productive (Phylum)")

#2 PA Oligo 
phyunprodPA <- subset_samples(good_phylum_nosherwin, ProdLevel == "Unproductive" & filter == "Particle") 
de_phyunprodPA <- deSEQ(phyunprodPA, ~ limnion)
phyunprodPA_plot <- plot_phylum_deSEQ(de_phyunprodPA, "Surface vs Bottom:  PA Unproductive (Phylum)")

#3 FL Productive 
phyprodFL <- subset_samples(good_phylum_nosherwin, ProdLevel == "Productive" & filter == "Free") 
de_phyprodFL <- deSEQ(phyprodFL, ~ limnion)
phyprodFL_plot <- plot_phylum_deSEQ(de_phyprodFL, "Surface vs Bottom:  FL Productive (Phylum)")

#4 FL Unproductive 
phyunprodFL <- subset_samples(good_phylum_nosherwin, ProdLevel == "Unproductive" & filter == "Free") 
de_phyunprodFL <- deSEQ(phyunprodFL, ~ limnion)
phyunprodFL_plot <- plot_phylum_deSEQ(de_phyunprodFL, "Surface vs Bottom:  FL Unproductive (Phylum)")

topbot_phylum <- multiplot(phyprodPA_plot, phyprodFL_plot, phyunprodPA_plot, phyunprodFL_plot,cols = 2)


############################
###########################################################PROD VS OLIGO
#1. Top PA 
phytopPA <- subset_samples(good_phylum_nosherwin, limnion == "Epilimnion" & filter == "Particle") ## DOES NOT WORK!!!!!!!!!!
de_phytopPA <- deSEQ(phytopPA, ~ ProdLevel)   ## DOES NOT WORK!!!!!!!!!!
phytopPA_plot <- plot_phylum_deSEQ(de_phytopPA, "Productive vs Unproductive:  PA Surface (Phylum)")

#2 BOTTOM PA 
phybotPA <- subset_samples(good_phylum_nosherwin, limnion == "Hypolimnion" & filter == "Particle") 
de_phybotPA <- deSEQ(phybotPA, ~ ProdLevel)
phybotPA_plot <- plot_phylum_deSEQ(de_phybotPA, "Productive vs Unproductive:  PA Bottom (Phylum)")

#3 TOP FL
phytopFL <- subset_samples(good_phylum_nosherwin, limnion == "Epilimnion" & filter == "Free") 
de_phytopFL <- deSEQ(phytopFL, ~ ProdLevel)  ## DOES NOT WORK!!!!!!!!!!
phytopFL_plot <- plot_phylum_deSEQ(de_phytopFL, "Productive vs Unproductive:  FL Surface (Phylum)")

#4 Bottom FL 
phybotFL <- subset_samples(good_phylum_nosherwin, limnion == "Hypolimnion" & filter == "Free") 
de_phybotFL <- deSEQ(phybotFL, ~ ProdLevel)
phybotFL_plot <- plot_phylum_deSEQ(de_phybotFL, "Productive vs Unproductive:  FL Bottom (Phylum)")

prodoligo_phylum <- multiplot(phybotFL_plot, phybotPA_plot, cols = 2)



#PA VS FL:  
#Surface productive
pafl_topprod <- subset(de_phytopProd, select = c(Phylum, log2FoldChange, padj))
pafl_topprod$Habitat <- "PA vs. FL: Top Productive"
#Surface unproductive
pafl_topUNprod <- subset(de_phytopUNProd, select = c(Phylum, log2FoldChange, padj))
pafl_topUNprod$Habitat <- "PA vs. FL: Top Unproductive"
#Bottom Productive
pafl_botprod <- subset(de_phybotProd, select = c(Phylum, log2FoldChange, padj))
pafl_botprod$Habitat <- "PA vs. FL: Bottom Productive"
#Bottom Unproductive
pafl_botUNprod <- subset(de_phybotUNProd, select = c(Phylum, log2FoldChange, padj))
pafl_botUNprod$Habitat <- "PA vs. FL: Bottom Unproductive"
## SHERMAN
pafl_isothermal <- subset(phy_shemeister, select = c(Phylum, log2FoldChange, padj))
pafl_isothermal$Habitat <- "PA vs. FL: Mixed"

###Top vs Bottom:  
topbot1 <- subset(de_phyprodPA, select = c(Phylum, log2FoldChange, padj))
topbot1$Habitat <- "Top vs. Bottom: PA Productive"
# Particle-Associated UNproductive
topbot2 <- subset(de_phyunprodPA, select = c(Phylum, log2FoldChange, padj))
topbot2$Habitat <- "Top vs. Bottom: PA Unproductive"
# Free-Living Productive 
topbot3 <- subset(de_phyprodFL, select = c(Phylum, log2FoldChange, padj))
topbot3$Habitat <- "Top vs. Bottom: FL Productive"
# Free-Living UNproductive
topbot4 <- subset(de_phyunprodFL, select = c(Phylum, log2FoldChange, padj))
topbot4$Habitat <- "Top vs. Bottom: FL Unproductive"


#Prod vs Oligo:  
#troph1 <- subset(de_phytopPA, select = c(Phylum, log2FoldChange, padj))  ## NO SIGNIFICANT CHANGES HERE!
#troph1$Habitat <- "Prod vs. Unprod: PA Epilimnion"
# Particle-Associated BOTTOM
troph2 <- subset(de_phybotPA, select = c(Phylum, log2FoldChange, padj))
troph2$Habitat <- "Prod vs. Unprod: PA Bottom"
# Free-Living TOP 
#troph3 <- subset(de_phytopFL, select = c(Phylum, log2FoldChange, padj))   ## NO SIGNIFICANT CHANGES HERE!
#troph3$Habitat <- "Prod vs. Unprod: FL Top"
# Free-Living BOTTOM
troph4 <- subset(de_phybotFL, select = c(Phylum, log2FoldChange, padj))
troph4$Habitat <- "Prod vs. Unprod: FL Bottom"
trophs <- rbind(troph2, troph4)
newlog2foldchange <- as.numeric(trophs$log2FoldChange * -1)  ##  To make PRODUCTIVE changes POSITIVE 
trophs$log2FoldChange <- newlog2foldchange

df_ratio <-rbind(pafl_topprod, pafl_topUNprod, pafl_botprod, pafl_botUNprod,pafl_isothermal,
                 topbot1, topbot2, topbot3, topbot4,
                 trophs)

paste(c("The range of the log2foldChange is",min(df_ratio$log2FoldChange), "to", max(df_ratio$log2FoldChange)))

#Split of the habitat column to 2 columns named comparison and habitat
split_cols <- colsplit(df_ratio$Habitat, ":", c("Comparison", "Habitat"))
df_ratio$Habitat = NULL
dfrat <- cbind(df_ratio, split_cols)

mf_labeller <- function(var, value){
  value <- as.character(value)
  if (var=="Comparison") { 
    value[value=="PA vs. FL"] <- "Particle-Associated \n vs. \nFree-Living"
    value[value=="Top vs. Bottom"]   <- "Hypolimnion \n vs. \nEpilimnion"
    value[value=="Prod vs. Unprod"]   <- "High-Nutrient \n vs.\nLow-Nutrient"
  }
  return(value)
}

uncla <- subset(dfrat, Phylum == "unclassified")

dfrat$Phylum <- factor(dfrat$Phylum,levels = c("Bacteroidetes", "Cyanobacteria", "Verrucomicrobia", "Betaproteobacteria", "Actinobacteria", "Planctomycetes", "Alphaproteobacteria", "Deltaproteobacteria",
                                               "Gammaproteobacteria", "Chloroflexi", "Lentisphaerae", "Chlorobi", "Armatimonadetes", "Firmicutes", "Acidobacteria", "Spirochaetae", "Candidate_division_OD1",
                                               "NPL-UPA2", "Deinococcus-Thermus", "Candidate_division_OP3", "TM6", "Chlamydiae", "Epsilonproteobacteria", "TA18", "Fibrobacteres", "Candidate_division_SR1", 
                                               "Candidate_division_BRC1", "BD1-5", "Gemmatimonadetes", "Fusobacteria", "Candidate_division_WS3", "Tenericutes", "Elusimicrobia", "WCHB1-60", "Deferribacteres", 
                                               "Candidate_division_TM7", "Candidate_division_OP8", "SPOTSOCT00m83", "Thermotogae", "Candidate_division_OP11", "Dictyoglomi", "unclassified"))

dfrat$Phylum = with(dfrat, factor(Phylum, levels = rev(levels(Phylum))))


#################################################################################  OVERALL ABUNDANCE PLOT
good_phylum_proteo <-tax_glom(merged_final,taxrank = "Phylum")
goodsamps_phylum <- tax_glom(good_phylum_proteo,taxrank = "Phylum")
goodsamps_phy_melt <- psmelt(goodsamps_phylum)  ##  Melt it into a dataframe 
sub_goodsamps_phy_melt <- subset(goodsamps_phy_melt, select = c("Sample", "ProdLevel", "quadrant", "Phylum","Abundance"))
TOTALS <- ddply(sub_goodsamps_phy_melt, c("Sample"), summarise, ##  Let's get the McMurdie and Holmes Scaled  Total for each sample 
                total   = sum(Abundance))   
sub_phy_melt_totals <- merge(sub_goodsamps_phy_melt, TOTALS, by = "Sample")  ### Merge phylum sums and the total numbers -> so we can calculate the Relative Abundance
sub_phy_melt_totals$RelAbundance <- sub_phy_melt_totals$Abundance/sub_phy_melt_totals$total  ## Calculate the relative abundance
sub_phy_melt_totals$PercentAbund <- sub_phy_melt_totals$RelAbundance * 100  ##  Calculate the Percent Abundance

####  Make a new dataframe with the percent abudance within the entire dataset!
phy_stats <- ddply(sub_phy_melt_totals, c("Phylum"), summarise, 
                   N = length(PercentAbund),
                   PercentPhy = mean(PercentAbund),
                   sd   = sd(PercentAbund),
                   se   = sd / sqrt(N))
abund <- subset(phy_stats,PercentPhy > 0.001)  # Only take the phyla that are more than 0.01% abundant

abund_ordered <- arrange(abund, desc(PercentPhy))

abund$Phylum <- factor(abund$Phylum,levels = c("Bacteroidetes", "Cyanobacteria", "Verrucomicrobia", "Betaproteobacteria", "Actinobacteria", "Planctomycetes",  "Alphaproteobacteria", "Deltaproteobacteria",
                                               "Gammaproteobacteria", "Chloroflexi", "Lentisphaerae", "Chlorobi", "Armatimonadetes", "Firmicutes", "Acidobacteria", "Spirochaetae", "Candidate_division_OD1",
                                               "NPL-UPA2", "Deinococcus-Thermus", "Candidate_division_OP3", "TM6", "Chlamydiae", "Epsilonproteobacteria", "TA18", "Fibrobacteres", "Candidate_division_SR1", 
                                               "Candidate_division_BRC1", "BD1-5", "Gemmatimonadetes", "Fusobacteria", "Candidate_division_WS3", "Tenericutes", "Elusimicrobia", "WCHB1-60", "Deferribacteres", 
                                               "Candidate_division_TM7", "Candidate_division_OP8", "SPOTSOCT00m83", "Thermotogae", "Candidate_division_OP11", "Dictyoglomi", "unclassified"))   

abund$Phylum = with(abund, factor(Phylum, levels = rev(levels(Phylum)))) 

#Relative abundance plot 
abund_plot <- ggplot(abund, aes(y=PercentPhy , x=Phylum))  +
  #geom_boxplot(fill = "magenta4", colour = "black") + 
  geom_bar(stat="identity", position=position_dodge(),  fill = "magenta4", colour = "black") +
  theme_bw() + ggtitle("Phyla Above 0.1% in All Samples") +
  xlab("Phylum") + ylab("Mean Relative Abundance (%)") +
  geom_errorbar(aes(ymin = PercentPhy -se, ymax = PercentPhy +se), width = 0.25) + coord_flip() +
  theme(axis.title.x = element_text(face="bold", size=16),
        axis.text.x = element_text(angle=0, colour = "black", size=14),
        axis.text.y = element_text(colour = "black", size=14),
        axis.title.y = element_text(face="bold", size=16),
        plot.title = element_text(face="bold", size = 20),
        legend.title = element_text(size=12, face="bold"),
        legend.text = element_text(size = 12),
        legend.position="none"); abund_plot


#################################################################################
#################################################################################
dfrat$Habitat <- as.character(dfrat$Habitat)
dfrat$Habitat[dfrat$Habitat == " Mixed"] <- "Mixed"
dfrat$Habitat[dfrat$Habitat == " Top Productive"] <- "Epilimnion \nHigh-Nutrient"
dfrat$Habitat[dfrat$Habitat == " Top Unproductive"] <- "Epilimnion \nLow-Nutrient"
dfrat$Habitat[dfrat$Habitat == " Bottom Productive"] <- "Hypolimnion \nHigh-Nutrient"
dfrat$Habitat[dfrat$Habitat == " Bottom Unproductive"] <- "Hypolimnion \nLow-Nutrient"
dfrat$Habitat[dfrat$Habitat == " PA Productive"] <- "High-Nutrient \nParticle-Associated"
dfrat$Habitat[dfrat$Habitat == " PA Unproductive"] <- "Low-Nutrient \nParticle-Associated"
dfrat$Habitat[dfrat$Habitat == " FL Productive"] <- "High-Nutrient \nFree-Living"
dfrat$Habitat[dfrat$Habitat == " FL Unproductive"] <- "Low-Nutrient \nFree-Living"
dfrat$Habitat[dfrat$Habitat == " PA Top"] <- "Particle-Associated \nEpilimnion"
dfrat$Habitat[dfrat$Habitat == " FL Top"] <- "Free-Living \nEpilimnion"
dfrat$Habitat[dfrat$Habitat == " FL Bottom"] <- "Free-Living \nHypolimnion"
dfrat$Habitat[dfrat$Habitat == " PA Bottom"] <- "Particle-Associated \nHypolimnion"

dfrat$Habitat <- factor(dfrat$Habitat,levels = c("Epilimnion \nHigh-Nutrient", "Hypolimnion \nHigh-Nutrient", "Mixed", "Epilimnion \nLow-Nutrient", "Hypolimnion \nLow-Nutrient",
                                                 "Particle-Associated \nHypolimnion", "Free-Living \nEpilimnion", "Free-Living \nHypolimnion", 
                                                 "High-Nutrient \nParticle-Associated", "High-Nutrient \nFree-Living", "Low-Nutrient \nParticle-Associated", "Low-Nutrient \nFree-Living"))

#jpeg(filename="~/Final_PAFL_Trophicstate/Figures/Fig.4_heat_only.jpeg", width= 30, height=35, units= "cm", pointsize= 8, res=250)
heat <- ggplot(dfrat, aes(Habitat, Phylum)) + geom_tile(aes(fill = log2FoldChange)) + 
  scale_fill_gradient2(name = "Odds\nRatio", mid = "gray", low = "darkorange", high = "blue4",  na.value = "white", guide = guide_colorbar(barwidth = 3, barheight = 18)) + #scale_y_reverse() + 
  theme_bw(base_size = 12) + scale_x_discrete(expand = c(0, 0)) + scale_y_discrete(expand = c(0, 0)) + 
  ylab("Phylum") + xlab("Habitat") + 
  #geom_text(aes(fill = splif2$Transformed, label = splif2$Transformed, size = 8)) +
  #scale_y_discrete(limits=phys) + xlab("Habitat") + ylab("Phylum") + 
  facet_grid(. ~ Comparison, scales = "free", space = "free", labeller=mf_labeller) + 
  theme(axis.text.x = element_text(colour="black", size=12, angle = 30, hjust = 1, vjust = 1), 
        axis.text.y = element_text(colour="black", vjust=0.5, size=14),
        axis.title.x = element_text(face="bold", size=16),
        legend.title = element_text(face="bold", size=14),
        legend.text = element_text(size = 14),
        legend.position = c(0.94, 0.13),#"left",
        axis.title.y = element_text(face="bold", size=16),
        plot.margin = unit(c(0.5, 1, 0.5, 0.5), "cm"),
        strip.text.x = element_text(size=14, face = "bold", colour = "black"),
        strip.background = element_blank()); 
#dev.off()

#http://stackoverflow.com/questions/4559229/drawing-pyramid-plot-using-r-and-ggplot2

setdiff(abund$Phylum, dfrat$Phylum)  # Discover the phyla that are in the abundance plot but NOT significantly differentially abundant at the phylum level
phylum_delete <- c("Candidate_division_OP11", "Candidate_division_OP8", "Candidate_division_TM7", "Deferribacteres", "Deltaproteobacteria", "Dictyoglomi", "SPOTSOCT00m83", "Thermotogae", "TM6", "unclassified", "WCHB1-60")

subset_abundPhylum <- subset(abund, !(Phylum %in% phylum_delete))
#setdiff(subset_abundPhylum$Phylum, dfrat$Phylum)  # Sanity check
#setdiff(dfrat$Phylum, subset_abundPhylum$Phylum)  # Double sanity check


relabun_plot <- ggplot(subset_abundPhylum, aes(y=PercentPhy , x=Phylum)) + #coord_cartesian(xlim = c(0, 30)) + 
  geom_bar(stat="identity", position=position_dodge(),fill = "gray", colour = "black") +
  ylab("Mean Percent \n Relative \n Abundance (%)") + coord_flip() + theme_bw() + 
  geom_errorbar(aes(ymin = PercentPhy -se, ymax = PercentPhy +se), width = 0.25, color = "black") + 
  scale_y_continuous(expand= c(0,0), limits = c(0,25)) +
  theme(axis.title.x = element_text(face="bold", size=16),
        axis.text.x = element_text(angle=0, colour = "black", size=12),
        axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        strip.background = element_rect(colour="black", fill = "black"),
        plot.margin = unit(c(2.09, 2, 1.51, -1.25), "cm"), #top, right, bottom, left    it was (2, 2, 1.65, -1.25)
        #panel.grid.minor=element_blank(), #panel.grid.major=element_blank(),
        legend.position="none"); #relabun_plot


#####  Plotting FIGURE 5  #####  Plotting FIGURE 5  #####  Plotting FIGURE 5  #####  Plotting FIGURE 5  #####  Plotting FIGURE 5
#http://stackoverflow.com/questions/20817094/how-to-control-width-of-multiple-plots-in-ggplot2
#tiff(filename="~/Final_PAFL_Trophicstate/Final_Figures/Fig.5_heat+abund.tiff", width= 40, height=35, units= "cm", pointsize= 8, res=200)
grid.newpage()
pushViewport(viewport(layout=grid.layout(1,2,width=c(0.8,0.2))))
print(heat, vp=viewport(layout.pos.row=1,layout.pos.col=1))
print(relabun_plot, vp=viewport(layout.pos.row=1,layout.pos.col=2))
#dev.off()


#### Subsetting out the top 14 phyla. 
phy14 <- c("Bacteroidetes", "Cyanobacteria", "Verrucomicrobia", "Betaproteobacteria", "Actinobacteria", "Planctomycetes",  "Alphaproteobacteria", "Deltaproteobacteria",
           "Gammaproteobacteria", "Chloroflexi", "Lentisphaerae", "Chlorobi", "Armatimonadetes", "Firmicutes")
phy14_dfrat <- subset(dfrat, Phylum %in% phy14)
phy14_dfrat <- subset(phy14_dfrat, Habitat != "Mixed")

phy14_dfrat$Phylum <- factor(phy14_dfrat$Phylum, levels = c("Bacteroidetes", "Cyanobacteria", "Verrucomicrobia", "Betaproteobacteria", "Actinobacteria", "Planctomycetes",  "Alphaproteobacteria", 
                                                            "Deltaproteobacteria", "Gammaproteobacteria", "Chloroflexi", "Lentisphaerae", "Chlorobi", "Armatimonadetes", "Firmicutes"))

phy14_dfrat$Phylum = with(phy14_dfrat, factor(Phylum, levels = rev(levels(Phylum)))) ### Reverse the order so its from the top down


phy_heat <- ggplot(phy14_dfrat, aes(Habitat, Phylum)) + geom_tile(aes(fill = log2FoldChange)) + 
  scale_fill_gradient2(name = "Odds-Ratio", mid = "gray", low = "darkorange", high = "blue4",  na.value = "white", guide = guide_colorbar(barwidth = 1, barheight = 4)) + #scale_y_reverse() + 
  theme_bw(base_size = 12) + scale_x_discrete(expand = c(0, 0)) + scale_y_discrete(expand = c(0, 0)) + 
  ylab("Phylum") + xlab("Habitat") + 
  #geom_text(aes(fill = splif2$Transformed, label = splif2$Transformed, size = 8)) +
  #scale_y_discrete(limits=phys) + xlab("Habitat") + ylab("Phylum") + 
  facet_grid(. ~ Comparison, scales = "free", space = "free", labeller=mf_labeller) + 
  theme(axis.title.x = element_text(face="bold", size=10),  #Set the x-axis title
        axis.title.y = element_text(face="bold", size=10, vjust=0.5),  #Set the y-axis title
        axis.text.x = element_text(colour = "black", size=8, angle = 30, hjust = 1, vjust = 1),  #Set the x-axis labels
        axis.text.y = element_text(colour = "black", size=8),  #Set the y-axis labels
        legend.title = element_text(size=7, face="bold"),  #Set the legend title 
        legend.text = element_text(size = 7),
        strip.text.x = element_text(size=8, face="bold"),  #Set the facet titles on x-axis 
        strip.text.y = element_text(size=8, face="bold"),  #Set the facet titles on x-axis 
        legend.position = c(0.92, 0.75),#"left",
        plot.margin = unit(c(0.5, 1, 0.5, 0.5), "cm"),
        strip.text.x = element_text(size=14, face = "bold", colour = "black"),
        strip.background = element_blank()); phy_heat
ggsave(phy_heat, filename = "~/Final_PAFL_Trophicstate/phy14_heat.pdf", height = 4, width = 6.5, dpi = 300)


#######  Making the abundance plot!
phy14_abund <- subset(abund, Phylum %in% phy14)
phy14_abund$Phylum <- factor(phy14_abund$Phylum, levels = c("Bacteroidetes", "Cyanobacteria", "Verrucomicrobia", "Betaproteobacteria", "Actinobacteria", "Planctomycetes",  "Alphaproteobacteria", 
                                                            "Deltaproteobacteria", "Gammaproteobacteria", "Chloroflexi", "Lentisphaerae", "Chlorobi", "Armatimonadetes", "Firmicutes"))

phy14_abund$Phylum = with(phy14_abund, factor(Phylum, levels = rev(levels(Phylum)))) ### Reverse the order so its from the top down


phy14_relabnd <- ggplot(phy14_abund, aes(y=PercentPhy , x=Phylum)) + #coord_cartesian(xlim = c(0, 30)) + 
  geom_bar(stat="identity", position=position_dodge(),fill = "gray", colour = "black") +
  ylab("Mean Percent \n Relative Abundance") + coord_flip() + theme_bw() + 
  geom_errorbar(aes(ymin = PercentPhy -se, ymax = PercentPhy +se), width = 0.25, color = "black") + 
  scale_y_continuous(expand= c(0,0), limits = c(0,25)) +
  theme(axis.title.x = element_text(face="bold", size=10),
        axis.text.x = element_text(angle=0, colour = "black", size=8),
        axis.text.y = element_text(colour = "black", size=8),  #Set the y-axis labels
        axis.title.y = element_text(face="bold", size=10, vjust=0.5),  #Set the y-axis title
        strip.background = element_rect(colour="black", fill = "black"),
        #plot.margin = unit(c(2.09, 2, 1.51, -1.25), "cm"), #top, right, bottom, left    it was (2, 2, 1.65, -1.25)
        #panel.grid.minor=element_blank(), #panel.grid.major=element_blank(),
        legend.position="none"); 
ggsave(phy14_relabnd, filename = "~/Final_PAFL_Trophicstate/phy14_abundance.pdf", height = 3.5, width = 3, dpi = 300)





####################################################  OTU LEVEL LOG-2-FOLD RATIO ANALYSIS  ####################################################  Working up to FIGURE 6 and Figure S5
####################################################  OTU LEVEL LOG-2-FOLD RATIO ANALYSIS  ####################################################  Working up to FIGURE 6 and Figure S5
####################################################  OTU LEVEL LOG-2-FOLD RATIO ANALYSIS  ####################################################  Working up to FIGURE 6 and Figure S5
### Subsetting Sherman lake for differences between particle and free living 
sherm <- subset_samples(merged_final, lakenames == "Sherman")  
shemeister <- deSEQ(sherm, ~ filter)
sherm_plot <- plot_deSEQ(shemeister, "Sherman: OTU-Level")
sherm_plots <- multiplot(physherm_plot,sherm_plot, cols = 2)


########## SUBSET OUT SHERMAN AND WINTERGREEN HYPOLIMNION
good_samps_nosher <- subset_samples(merged_final, lakenames != "Sherman")
good_samps_nosherwin <- subset_samples(good_samps_nosher, names != "WINH" & names != "WINH3um")
#View(data.frame(sample_data(good_samps_nosherwin))) 


############################
###########################################################PA VS FL
#1. Top prod 
topProd <- subset_samples(good_samps_nosherwin, ProdLevel == "Productive" & limnion == "Epilimnion")
de_topProd <- deSEQ(topProd, ~ filter)
topProd_plot <- plot_deSEQ(de_topProd, "PA vs FL:  Surface Productive (OTU)")

#2 Top Oligo 
topUNProd <- subset_samples(good_samps_nosherwin, ProdLevel == "Unproductive" & limnion == "Epilimnion")
de_topUNProd <- deSEQ(topUNProd, ~ filter)
topUNProd_plot <- plot_deSEQ(de_topUNProd, "PA vs FL:  Surface Unproductive (OTU)")

#3 Bottom Productive 
botProd <- subset_samples(good_samps_nosherwin, ProdLevel == "Productive" & limnion == "Hypolimnion") 
de_botProd <- deSEQ(botProd, ~ filter)
botProd_plot <- plot_deSEQ(de_botProd, "PA vs FL:  Bottom Productive (OTU)")

#4 Bottom Unproductive 
botUNProd <- subset_samples(good_samps_nosherwin, ProdLevel == "Unproductive" & limnion == "Hypolimnion") 
de_botUNProd <- deSEQ(botUNProd, ~ filter)
botUNProd_plot <- plot_deSEQ(de_botUNProd, "PA vs FL:  Bottom Unproductive (OTU)")

pafl_otu <- multiplot(topProd_plot, botProd_plot, topUNProd_plot, botUNProd_plot,cols = 2)


############################
###########################################################TOP VS BOTTOM
#1. Top prod 
prodPA <- subset_samples(good_samps_nosherwin, ProdLevel == "Productive" & filter == "Particle") 
de_prodPA <- deSEQ(prodPA, ~ limnion)
prodPA_plot <- plot_deSEQ(de_prodPA, "Surface vs Bottom:  PA Productive (OTU)")

#2 Top Oligo 
unprodPA <- subset_samples(good_samps_nosherwin, ProdLevel == "Unproductive" & filter == "Particle") 
de_unprodPA <- deSEQ(unprodPA, ~ limnion)
unprodPA_plot <- plot_deSEQ(de_unprodPA, "Surface vs Bottom:  PA Unproductive (OTU)")

#3 Bottom Productive 
prodFL <- subset_samples(good_samps_nosherwin, ProdLevel == "Productive" & filter == "Free") 
de_prodFL <- deSEQ(prodFL, ~ limnion)
prodFL_plot <- plot_deSEQ(de_prodFL, "Surface vs Bottom:  FL Productive (OTU)")

#4 Bottom Unproductive 
unprodFL <- subset_samples(good_samps_nosherwin, ProdLevel == "Unproductive" & filter == "Free") 
de_unprodFL <- deSEQ(unprodFL, ~ limnion)
unprodFL_plot <- plot_deSEQ(de_unprodFL, "Surface vs Bottom:  FL Unproductive (OTU)")

topbot_otu <- multiplot(prodPA_plot, prodFL_plot, unprodPA_plot, unprodFL_plot,cols = 2)

############################
###########################################################PROD VS OLIGO
#1. Top prod 
topPA <- subset_samples(good_samps_nosherwin, limnion == "Epilimnion" & filter == "Particle") 
de_topPA <- deSEQ(topPA, ~ ProdLevel)
topPA_plot <- plot_deSEQ(de_topPA, "Productive vs Unproductive:  PA Surface (OTU)")

#2 Top Oligo 
botPA <- subset_samples(good_samps_nosherwin, limnion == "Hypolimnion" & filter == "Particle") 
de_botPA <- deSEQ(botPA, ~ ProdLevel)
botPA_plot <- plot_deSEQ(de_botPA, "Productive vs Unproductive:  PA Bottom (OTU)")

#3 Bottom Productive 
topFL <- subset_samples(good_samps_nosherwin, limnion == "Epilimnion" & filter == "Free") 
de_topFL <- deSEQ(topFL, ~ ProdLevel)
topFL_plot <- plot_deSEQ(de_topFL, "Productive vs Unproductive:  FL Surface (OTU)")

#4 Bottom Unproductive 
botFL <- subset_samples(good_samps_nosherwin, limnion == "Hypolimnion" & filter == "Free") 
de_botFL <- deSEQ(botFL, ~ ProdLevel)
botFL_plot <- plot_deSEQ(de_botFL, "Productive vs Unproductive:  FL Bottom (OTU)")

prodoligo_otu <- multiplot(topFL_plot, botFL_plot, topPA_plot, botPA_plot, cols = 2)



#PA VS FL:  
#Surface Productive
pafl1 <- subset(de_topProd, select = c(Phylum, Genus, Species, log2FoldChange, padj))
pafl1$Habitat <- "PA vs. FL: Epilimnion \nHigh-Nutrient"
#Surface Low-Nutrient
pafl2 <- subset(de_botProd, select = c(Phylum, Genus, Species, log2FoldChange, padj))
pafl2$Habitat <- "PA vs. FL: Epilimnion \nLow-Nutrient"
#Bottom High-Nutrient
pafl3 <- subset(de_topUNProd, select = c(Phylum, Genus, Species, log2FoldChange, padj))
pafl3$Habitat <- "PA vs. FL: Hypolimnion \nHigh-Nutrient"
#Bottom Low-Nutrient
pafl4 <- subset(de_botUNProd, select = c(Phylum, Genus, Species, log2FoldChange, padj))
pafl4$Habitat <- "PA vs. FL: Hypolimnion \nLow-Nutrient"



###Top vs Bottom:  
otu_topbot1 <- subset(de_prodPA, select = c(Phylum, Genus, Species, log2FoldChange, padj))
otu_topbot1$Habitat <- "Top vs. Bottom: High-Nutrient \n Particle-Associated"
# Particle-Associated Low-Nutrient
otu_topbot2 <- subset(de_unprodPA, select = c(Phylum, Genus, Species, log2FoldChange, padj))
otu_topbot2$Habitat <- "Top vs. Bottom: Low-Nutrient \n Particle-Associated"
# Free-Living High-Nutrient 
otu_topbot3 <- subset(de_prodFL, select = c(Phylum, Genus, Species, log2FoldChange, padj))
otu_topbot3$Habitat <- "Top vs. Bottom: High-Nutrient \n Free-Living"
# Free-Living Low-Nutrient
otu_topbot4 <- subset(de_unprodFL, select = c(Phylum, Genus, Species, log2FoldChange, padj))
otu_topbot4$Habitat <- "Top vs. Bottom: Low-Nutrient \n Free-Living"



#Prod vs Oligo:  
otu_troph1 <- subset(de_topPA, select = c(Phylum, Genus, Species, log2FoldChange, padj)) 
otu_troph1$Habitat <- "Prod vs. Unprod: Particle-Associated \n Epilimnion"
# Particle-Associated BOTTOM
otu_troph2 <- subset(de_botPA, select = c(Phylum, Genus, Species, log2FoldChange, padj))
otu_troph2$Habitat <- "Prod vs. Unprod: Particle-Associated \n Hypolimnion"
# Free-Living TOP 
otu_troph3 <- subset(de_topFL, select = c(Phylum, Genus, Species, log2FoldChange, padj))
otu_troph3$Habitat <- "Prod vs. Unprod: Free-Living \nEpilimnion"
# Free-Living BOTTOM
otu_troph4 <- subset(de_botFL, select = c(Phylum, Genus, Species, log2FoldChange, padj))
otu_troph4$Habitat <- "Prod vs. Unprod: Free-Living \nHypolimnion"


######  Combing all into one dataframe named out_ratio
otu_trophs <- rbind(otu_troph1, otu_troph2, otu_troph3, otu_troph4)
newlog2foldchange <- as.numeric(otu_trophs$log2FoldChange * -1) # To correlate with Top vs bottom
otu_trophs$log2FoldChange <- newlog2foldchange

otu_ratio <-rbind(pafl1, pafl2, pafl3, pafl4,
                  otu_topbot1, otu_topbot2, otu_topbot3, otu_topbot4,
                  otu_trophs) #otu_troph1, otu_troph2, otu_troph3, otu_troph4)

#What's the min and max ratios?
paste(c("The range of the log2foldChange is", min(otu_ratio$log2FoldChange), "to",  max(otu_ratio$log2FoldChange)))


#Split of the habitat column to 2 columns named comparison and habitat
otu_cols <- colsplit(otu_ratio$Habitat, ":", c("Comparison", "Habitat"))
otu_ratio$Habitat = NULL
otu_ratios <- cbind(otu_ratio, otu_cols)

mf_labeller <- function(var, value){
  value <- as.character(value)
  if (var=="Comparison") { 
    value[value=="PA vs. FL"] <- "Particle-Associated \n vs. \n Free-Living"
    value[value=="Top vs. Bottom"]   <- "Hypolimnion \n vs. \nEpilimnion"
    value[value=="Prod vs. Unprod"]   <- "High-Nutrient \n vs. \nLow-Nutrient"
  }
  return(value)
}

#no_unclass_genus <- subset(otu_ratios, Genus != "unclassified")
phy_otu <- otu_ratios 
phy_otu$Phylum <- factor(phy_otu$Phylum,levels = c("Bacteroidetes", "Cyanobacteria", "Verrucomicrobia", "Betaproteobacteria", "Actinobacteria", "Planctomycetes", "Alphaproteobacteria", "Deltaproteobacteria",
                                               "Gammaproteobacteria", "Chloroflexi", "Lentisphaerae", "Chlorobi", "Armatimonadetes", "Firmicutes", "Acidobacteria", "Spirochaetae", "Candidate_division_OD1",
                                               "NPL-UPA2", "Deinococcus-Thermus", "Candidate_division_OP3", "TM6", "Chlamydiae", "Epsilonproteobacteria", "TA18", "Fibrobacteres", "Candidate_division_SR1", 
                                               "Candidate_division_BRC1", "BD1-5", "Gemmatimonadetes", "Fusobacteria", "Candidate_division_WS3", "Tenericutes", "Elusimicrobia", "WCHB1-60", "Deferribacteres", 
                                               "Candidate_division_TM7", "Candidate_division_OP8", "SPOTSOCT00m83", "Thermotogae", "Candidate_division_OP11", "Dictyoglomi", "unclassified"))
phy_otu$Phylum = with(phy_otu, factor(Phylum, levels = rev(levels(Phylum)))) 



########  GENUS LEVEL PLOT ORGANIZED BY PHYLUM
## GENUS LEVEL HEAT PLOT
sub_abund <- subset(abund, select = c("Phylum", "PercentPhy"))
oturats_abund <- merge(otu_ratios, sub_abund, by = "Phylum")
ordered_otu_ratios <- arrange(oturats_abund, desc(PercentPhy)) 

ordered_otu_ratios$Genus <- factor(ordered_otu_ratios$Genus,levels = c("vadinBC27_wastewater-sludge_group", "unclassified", "Pseudarcicella", "Solitalea", "Haliscomenobacter", "Ferruginibacter",
                                                   "Candidatus_Aquirestis", "Paludibacter", "Flavobacterium", "Sediminibacterium", "Fluviicola", "Pedobacter",
                                                   "Algoriphagus", "Owenweeksia", "Emticicia", "Chitinophaga", "Planktothrix","Synechococcus", "Snowella", "Anabaena",
                                                   "Pseudanabaena", "Luteolibacter","Opitutus","Candidatus_Methylacidiphilum", "LD28_freshwater_group","Iodobacter","Nitrosospira",
                                                   "Sterolibacterium","Dechloromonas","Sulfuritalea","Candidatus_Accumulibacter","PRD01a011B","Sideroxydans","MWH-UniP1_aquatic_group",
                                                   "Polynucleobacter","GKS98_freshwater_group",
                                                   "Candidatus_Branchiomonas","Candidatus_Planktophila","hgcI_clade","Candidatus_Microthrix","CL500-29_marine_group","Mycobacterium",
                                                   "Candidatus_Aquiluna","Candidatus_Limnoluna","Planctomyces","Pirellula","Candidatus_Anammoximicrobium","Candidatus_Nostocoida",
                                                   "CL500-3","Pir4_lineage","Roseospirillum","Magnetospirillum","Brevundimonas","Porphyrobacter", "Haematobacter",
                                                   "Roseomonas","Rhizomicrobium","Rickettsia","OM27_clade","Desulfobacula",
                                                   "Desulfobulbus","Syntrophus","Peredibacter","Desulfomonile","Desulfocapsa","Desulfatirhabdium","Geobacter",
                                                   "Desulfovibrio","Phaselicystis","Sorangium","Lamprocystis","Pseudospirillum",
                                                   "Thiodictyon","Thiocystis","Candidatus_Competibacter","Coxiella","Methylocaldum","Rheinheimera",
                                                   "Methylobacter","Crenothrix","Oscillochloris","Anaerolinea","Chloronema","Victivallis",
                                                   "Chlorobium","Ignavibacterium","Armatimonas","Fastidiosipila","Fonticella","Incertae_Sedis",
                                                   "Anaerofustis","marine_group","Geothrix","Treponema","Leptospira","Spirochaeta",
                                                   "Brevinema","Sulfurospirillum","Sulfurimonas","possible_genus_03","Gemmatimonas","Candidatus_Latescibacter","SC103"))
ordered_otu_ratios$Genus = with(ordered_otu_ratios, factor(Genus, levels = rev(levels(Genus)))) 

sub_ordered_oturats <- subset(ordered_otu_ratios, Genus != "unclassified")


#####  Plotting FIGURE S5  #####  Plotting FIGURE S5  #####  Plotting FIGURE S5  #####  Plotting FIGURE S5  #####  Plotting FIGURE S5
#tiff(filename="~/Final_PAFL_Trophicstate/Final_Figures/Fig.S5_genus_heat.tiff", width= 40, height=60, units= "cm", pointsize= 8, res=200)
ggplot(sub_ordered_oturats, aes(Habitat, Genus)) + geom_tile(aes(fill = log2FoldChange)) + 
  scale_fill_gradient2(name = "Odds-\nRatio", mid = "gray", low = "darkorange", high = "blue4",  na.value = "white", guide = guide_colorbar(barwidth = 3, barheight = 15)) + #scale_y_reverse() + 
  theme_bw(base_size = 12) + scale_x_discrete(expand = c(0, 0)) + scale_y_discrete(expand = c(0, 0)) + ylab(NULL) + 
  #geom_text(aes(fill = splif2$Transformed, label = splif2$Transformed, size = 8)) +
  xlab("Habitat") + ylab("Genus") + 
  facet_grid(. ~ Comparison, scales = "free", space = "free", labeller=mf_labeller) + 
  theme(axis.text.x = element_text(colour="black", size=14, angle = 30, hjust = 1, vjust = 1), 
        axis.text.y = element_text(colour="black", vjust=0.5, size=14),
        axis.title.x = element_text(face="bold", size=16),
        legend.title = element_text(face="bold", size=12),
        legend.text = element_text(size = 12),
        legend.position = c(0.96, 0.06), #c(0.1, 0.93),
        axis.title.y = element_text(face="bold", size=16),
        plot.margin = unit(c(0.5, 1, 0.5, 0.5), "cm"),
        strip.text.x = element_text(size=16, face = "bold", colour = "black"),
        strip.background = element_blank());  
#dev.off()



########   Summed-OTU Plot 
### PA (positive) vs FL (negative)
### Top (negative) vs Bottom (positive)
###Prod (negative) vs Unprod (positive)
#head(otu_ratios)
#subset out prod vs unprod
sig_oturats <- otu_ratios
pvu <- subset(sig_oturats, Comparison == "Prod vs. Unprod")
tvb <- subset(sig_oturats, Comparison == "Top vs. Bottom")
pvf <- subset(sig_oturats, Comparison == "PA vs. FL")
#plot_deSEQ(pvu, "Prod vs. Unprod")
#plot_deSEQ(tvb, "Top vs Bottom")
#plot_deSEQ(pvf, "PA vs FL")


# Let's start with PA vs FL as it has the clearest trend and the fewest Phyla
#we need a df with the sums of significant OTUs in each phylum for PA AND FL
pvf_FL <- subset(pvf, log2FoldChange < 0)  # There are 
length(unique(pvf_FL$Species))  ### 35 OTUs are significantly overrepresented in FL across all habitats
pvf_FL$Preference <- "Free-Living"
pvf_PA <- subset(pvf, log2FoldChange > 0)
length(unique(pvf_PA$Species))  ### 67 OTUs are significantly overrepresented in FL across all habitats
pvf_PA$Preference <- "Particle-Associated"
pvf2 <- rbind(pvf_FL, pvf_PA)
# Make a new data frame with each phylum and the number of sigs in each environment
phylumlist <- as.character(unique(pvf2$Phylum))
res_PA <- data.frame(matrix(NA,length(pvf2),3)) #Create a dummy table to fill in values with during for loop.
names(res_PA) <- c("Phylum", "NumSigOTUs", "Preference")
res_FL <- data.frame(matrix(NA,length(pvf2),3))
names(res_FL) <- c("Phylum", "NumSigOTUs", "Preference")
for(i in 1:length(phylumlist)){
  ## PARTICLE ASSOCIATED
  df_PA <- subset(pvf_PA, Phylum == phylumlist[i]) #particle-associated phylum i
  num_sigOTUs <- length(unique(df_PA$Species))
  res_PA[i,1] <- phylumlist[i] #put the name of the phylum in the table
  res_PA[i,2] <- num_sigOTUs #put the number of sig OTUs in the table. 
  res_PA[i,3] <- "Particle-Associated"
  ##FREE LIVING
  df_FL <- subset(pvf_FL, Phylum == phylumlist[i]) #free phylum i
  num_sigOTUs <- length(unique(df_FL$Species))
  res_FL[i,1] <- phylumlist[i] #put the name of the phylum in the table.
  res_FL[i,2] <- num_sigOTUs * -1 #put the number of sig OTUs in the table. 
  res_FL[i,3] <- "Free-Living"
  PAFL_OTUresults <- rbind(res_PA, res_FL)
  PAFL_OTUresults$Comparison <- "Particle-Associated \nvs. \nFree-Living"
}


#we need a df with the sums of significant OTUs in each phylum for PROD AND UNPROD
###  Prod (POSITIVE) vs Unprod (NEGATIVE)
#Select Unproductive
pvu_oligo <- subset(pvu, log2FoldChange < 0) 
length(unique(pvu_oligo$Species))  ### 63 OTUs are significantly overrepresented in FL across all habitats
pvu_oligo$Preference <- "Unproductive"
#Select Productive
pvu_prod <- subset(pvu, log2FoldChange > 0)
length(unique(pvu_prod$Species))  ### 67 OTUs are significantly overrepresented in FL across all habitats
pvu_prod$Preference <- "Productive"
pvu2 <- rbind(pvu_oligo, pvu_prod)
# Make a new data frame with each phylum and the number of sigs in each environment
pvu_phylumlist <- as.character(unique(pvu2$Phylum))
res_prod <- data.frame(matrix(NA,length(pvu_phylumlist),3))
names(res_prod) <- c("Phylum", "NumSigOTUs", "Preference")
res_oligo <- data.frame(matrix(NA,length(pvu_phylumlist),3))
names(res_oligo) <- c("Phylum", "NumSigOTUs", "Preference")
for(i in 1:length(pvu_phylumlist)){
  ## PRODUCTIVE
  df_prod <- subset(pvu_prod, Phylum == pvu_phylumlist[i]) #free phylum i
  num_sigOTUs <- length(unique(df_prod$Species))
  res_prod[i,1] <- pvu_phylumlist[i] #put the name of the phylum in the table
  res_prod[i,2] <- num_sigOTUs #put the number of sig OTUs in the table. 
  res_prod[i,3] <- "High-Nutrient"
  ## UNPRODUCTIVE
  df_oligo <- subset(pvu_oligo, Phylum == pvu_phylumlist[i]) #free phylum i
  num_sigOTUs <- length(unique(df_oligo$Species))
  res_oligo[i,1] <- pvu_phylumlist[i] #put the name of the phylum in the table.
  res_oligo[i,2] <- num_sigOTUs * -1 #put the number of sig OTUs in the table. 
  res_oligo[i,3] <- "Low-Nutrient"
  PVU_OTUresults <- rbind(res_prod,res_oligo)
  PVU_OTUresults$Comparison <- "High-Nutrient \n vs. \n Low-Nutrient"
}


#we need a df with the sums of significant OTUs in each phylum for TOP AND BOTTOM
### Top (negative) vs Bottom (positive)
#Select TOP
tvb_top <- subset(tvb, log2FoldChange < 0) 
length(unique(tvb_top$Species))  ### 51 OTUs are significantly overrepresented in FL across all habitats
tvb_top$Preference <- "Epilimnion"
#Select Bottom
tvb_bottom <- subset(tvb, log2FoldChange > 0)
length(unique(tvb_bottom$Species))  ### 169 OTUs are significantly overrepresented in FL across all habitats
tvb_bottom$Preference <- "Hypolimnion"
tvb2 <- rbind(tvb_top, tvb_bottom)
# Make a new data frame with each phylum and the number of sigs in each environment
tvb_phylumlist <- as.character(unique(tvb2$Phylum))
res_top <- data.frame(matrix(NA,length(tvb_phylumlist),3))
names(res_top) <- c("Phylum", "NumSigOTUs", "Preference")
res_bottom <- data.frame(matrix(NA,length(tvb_phylumlist),3))
names(res_bottom) <- c("Phylum", "NumSigOTUs", "Preference")
for(i in 1:length(tvb_phylumlist)){
  ## BOTTOM
  df_bottom <- subset(tvb_bottom, Phylum == tvb_phylumlist[i]) #free phylum i
  num_sigOTUs <- length(unique(df_bottom$Species))
  res_bottom[i,1] <- tvb_phylumlist[i] #put the name of the phylum in the table
  res_bottom[i,2] <- num_sigOTUs #put the number of sig OTUs in the table. 
  res_bottom[i,3] <- "Hypolimnion"
  ## TOP
  df_top <- subset(tvb_top, Phylum == tvb_phylumlist[i]) #free phylum i
  num_sigOTUs <- length(unique(df_top$Species))
  res_top[i,1] <- tvb_phylumlist[i] #put the name of the phylum in the table.
  res_top[i,2] <- num_sigOTUs * -1 #put the number of sig OTUs in the table. 
  res_top[i,3] <- "Epilimnion"
  TVB_OTUresults <- rbind(res_bottom,res_top)
  TVB_OTUresults$Comparison <- "Hypolimnion \n vs. \n Epilimnion"
}


### COMBINING ALL
allOTU_results <- rbind(TVB_OTUresults, PAFL_OTUresults, PVU_OTUresults)

allOTU_results$Preference <- factor(allOTU_results$Preference,
                                    levels = c("Particle-Associated","Free-Living","Low-Nutrient","High-Nutrient","Epilimnion", "Hypolimnion")) 

allOTU_results$Comparison <- factor(allOTU_results$Comparison,
                                    levels = c("Particle-Associated \nvs. \nFree-Living", "High-Nutrient \n vs. \n Low-Nutrient","Hypolimnion \n vs. \n Epilimnion")) 

allOTU_results$Phylum <- factor(allOTU_results$Phylum,levels = c("Bacteroidetes", "Cyanobacteria", "Verrucomicrobia", "Betaproteobacteria", "Actinobacteria", "Planctomycetes",  "Alphaproteobacteria", "Deltaproteobacteria",
                                                                 "Gammaproteobacteria", "Chloroflexi", "Lentisphaerae", "Chlorobi", "Armatimonadetes", "Firmicutes", "Acidobacteria", "Spirochaetae", "Candidate_division_OD1",
                                                                 "NPL-UPA2", "Deinococcus-Thermus", "Candidate_division_OP3", "TM6", "Chlamydiae", "Epsilonproteobacteria", "TA18", "Fibrobacteres", "Candidate_division_SR1", 
                                                                 "Candidate_division_BRC1", "BD1-5", "Gemmatimonadetes", "Fusobacteria", "Candidate_division_WS3", "Tenericutes", "Elusimicrobia", "WCHB1-60", "Deferribacteres", 
                                                                 "Candidate_division_TM7", "Candidate_division_OP8", "SPOTSOCT00m83", "Thermotogae", "Candidate_division_OP11", "Dictyoglomi", "unclassified"))

OTUtop20 <- c("Bacteroidetes", "Cyanobacteria", "Verrucomicrobia", "Betaproteobacteria", "Actinobacteria", "Planctomycetes",  "Alphaproteobacteria", "Deltaproteobacteria",
              "Gammaproteobacteria", "Chloroflexi", "Lentisphaerae", "Chlorobi", "Armatimonadetes", "Firmicutes", "Acidobacteria", "Spirochaetae","unclassified")

summed_20 <- allOTU_results[allOTU_results$Phylum %in% OTUtop20, ]

summed_20$Phylum = with(summed_20, factor(Phylum, levels = rev(levels(Phylum))))


#####  Plotting FIGURE 6  #####  Plotting FIGURE 6  #####  Plotting FIGURE 6  #####  Plotting FIGURE 6  #####  Plotting FIGURE 6
### This plot is flipped!
#tiff(filename="~/Final_PAFL_Trophicstate/Final_Figures/Fig.6_summed_otus_top17.tiff",  width= 30, height=25, units= "cm", pointsize= 8, res=200)
summed <- ggplot(summed_20, aes(y=NumSigOTUs, x=Phylum, fill=Preference)) + 
  geom_bar(stat="identity", position="identity") + coord_flip() + #ggtitle("Summed OTUs") +
  geom_bar(stat="identity", colour = "black", show_guide = FALSE, position="identity") +
  theme_gray() +# scale_y_continuous(breaks=seq(-200, 200, 25)) +
  facet_grid(. ~ Comparison, scales = "free_y", labeller = mf_labeller) + theme_bw() +
  ylab("Total Number of Significant OTUs") +  xlab("Phylum") + 
  scale_fill_manual(name = "", limits=c("Particle-Associated","Free-Living","Low-Nutrient","High-Nutrient","Epilimnion", "Hypolimnion"), 
                     values = c("firebrick1", "goldenrod1",  "turquoise3", "green4","palevioletred1","cornflowerblue"))+
  guides(fill = guide_legend(keywidth = 2, keyheight = 2)) + 
  theme(axis.title.x = element_text(face="bold", size=16),
        axis.text.x = element_text(colour = "black", size=14),
        axis.text.y = element_text(colour = "black", size=14),
        axis.title.y = element_text(face="bold", size=16),
        plot.title = element_text(face="bold", size = 20),
        legend.title = element_blank(),
        legend.text = element_text(size = 12),
        strip.text.y = element_blank(),
        strip.text.x = element_text(size = 14, face = "bold", colour = "black"),
        strip.background = element_blank(),
        legend.position="right"); summed
#dev.off()


####  Top 14!!!
phy14 <- c("Bacteroidetes", "Cyanobacteria", "Verrucomicrobia", "Betaproteobacteria", "Actinobacteria", "Planctomycetes",  "Alphaproteobacteria", "Deltaproteobacteria",
           "Gammaproteobacteria", "Chloroflexi", "Lentisphaerae", "Chlorobi", "Armatimonadetes", "Firmicutes")

otu14 <- allOTU_results[allOTU_results$Phylum %in% phy14, ]
otu14 <- subset(allOTU_results, Phylum %in% phy14)
otu14$Phylum = with(otu14, factor(Phylum, levels = rev(levels(Phylum))))

otu14$Phylum <- factor(otu14$Phylum, levels = c("Bacteroidetes", "Cyanobacteria", "Verrucomicrobia", "Betaproteobacteria", "Actinobacteria", "Planctomycetes",  "Alphaproteobacteria", 
                                                            "Deltaproteobacteria", "Gammaproteobacteria", "Chloroflexi", "Lentisphaerae", "Chlorobi", "Armatimonadetes", "Firmicutes"))

#otu14$Phylum = with(otu14, factor(Phylum, levels = rev(levels(Phylum)))) ### Reverse the order so its from the top down

###  PLOT ONLY THE PA AND FL!
otu14_PAFL <- subset(otu14, Comparison == "Particle-Associated \nvs. \nFree-Living")
otu14_PAFL <- droplevels(otu14_PAFL) # This function drops extra factors from the vector 


otu14_plot <- ggplot(otu14_PAFL, aes(y=NumSigOTUs, x=Phylum, fill=Preference)) + 
  geom_bar(stat="identity", position="identity") + #coord_flip() + #ggtitle("Summed OTUs") +
  #ggtitle("Number of Total Significant OTUs per Phylum in each Filter Fraction") + 
  geom_bar(stat="identity", colour = "black", show_guide = FALSE, position="identity") +
  ylab("Number of Significant OTUs") +  xlab("Phylum") + 
  scale_fill_manual(name = "", limits=c("Particle-Associated","Free-Living"), #,"Low-Nutrient","High-Nutrient","Epilimnion", "Hypolimnion"), 
                    values = c("firebrick1", "goldenrod1")) + #,  "turquoise3", "green4","palevioletred1","cornflowerblue"))+
  guides(fill = guide_legend(keywidth = 0.8, keyheight = 0.8)) + 
  theme(plot.title = element_text(face="bold", size=12),
        axis.title.x = element_text(face="bold", size=10),  #Set the x-axis title
        axis.title.y = element_text(face="bold", size=10, vjust=0.5),  #Set the y-axis title
        axis.text.x = element_text(colour = "black", size=8, angle = 30, hjust = 1, vjust = 1),  #Set the x-axis labels
        axis.text.y = element_text(colour = "black", size=8),  #Set the y-axis labels
        legend.title = element_text(size=8, face="bold"),  #Set the legend title 
        legend.text = element_text(size = 8),
        strip.text.x = element_text(size=8, face="bold"),  #Set the facet titles on x-axis 
        strip.text.y = element_text(size=8, face="bold"),  #Set the facet titles on x-axis 
        legend.position = c(0.8, 0.18),
        plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), "cm"),
        strip.background = element_blank()); otu14_plot
ggsave(otu14_plot, filename = "~/Final_PAFL_Trophicstate/otu14_summed_PA.pdf", height = 3, width = 4, dpi = 300)


##################################################################################### ABUNDANCE PLOTS ##################  Working up to FIGURE S3 and FIGURE S4
##################################################################################### ABUNDANCE PLOTS ##################  Working up to FIGURE S3 and FIGURE S4
##################################################################################### ABUNDANCE PLOTS ##################  Working up to FIGURE S3 and FIGURE S4
### Check lines 2261 for how sub_phy_melt_totals was created:
### Calculate the mean relative abundance based on ProdLevel + Quadrant for each PHYLUM 
sub_phy_melt_totals_nosherwin <- subset(sub_phy_melt_totals, Sample != "SHEE" & Sample != "SHEE3um" & Sample !="SHEH" & Sample != "SHEH3um")
prod_quad_phylum_stats <- ddply(sub_phy_melt_totals_nosherwin, c("ProdLevel","quadrant", "Phylum"), summarise, 
                                N = length(PercentAbund),
                                mean_abundance = mean(PercentAbund),
                                sd   = sd(PercentAbund),
                                se   = sd / sqrt(N))

abund_only_by_phylum <- ddply(prod_quad_phylum_stats, c("Phylum"), summarise, 
                              N = length(mean_abundance),
                              phylum_mean = mean(mean_abundance))

abund_only_by_phylum <-arrange(abund_only_by_phylum, desc(phylum_mean)) # abund_only_by_phylum[with(abund_only_by_phylum, order(-phylum_mean)), ]

phy_order <- as.character(abund_only_by_phylum$Phylum)  # THIS IS A VECTOR OF THE PHYLUM IN THE ORDER THAT WE WANT

### TOP 25 PHYLA!
top15phy <- prod_quad_phylum_stats[prod_quad_phylum_stats$Phylum %in% phy_order[1:15], ]

top15phy$Phylum <- factor(top15phy$Phylum,levels = c("Bacteroidetes", "Cyanobacteria", "Verrucomicrobia", "Betaproteobacteria", "Actinobacteria", "Planctomycetes", "Alphaproteobacteria", "unclassified",          
                                                     "Deltaproteobacteria", "Gammaproteobacteria", "Chloroflexi", "Lentisphaerae", "Armatimonadetes",  "Firmicutes",  "Chlorobi","Acidobacteria",         
                                                     "Spirochaetae", "Candidate_division_OD1",  "NPL-UPA2", "Deinococcus-Thermus", "Candidate_division_OP3", "TM6", "TA18", "Epsilonproteobacteria", "Chlamydiae"))


top15phy$quadrant <- factor(top15phy$quadrant,levels = c("Free Epilimnion", "Particle Epilimnion", "Free Hypolimnion", "Particle Hypolimnion"))

top15phy$Phylum = with(top15phy, factor(Phylum, levels = rev(levels(Phylum))))

abund_labeller <- function(var, value){
  value <- as.character(value)
  if (var=="quadrant") { 
    value[value=="Particle Epilimnion"] <- "Particle \n Epilimnion"
    value[value=="Free Epilimnion"]   <- "Free \n Epilimnion"
    value[value=="Particle Hypolimnion"]   <- "Particle \n Hypolimnion"
    value[value=="Free Hypolimnion"]   <- "Free \n Hypolimnion"
  }
  return(value)
}

#### ADDING NO SPACES FOR FREE LIVING and ADDING 2 SPACE FOR PARICLE ASSOCIATED AT THE END
top15phy$quadrant <- as.character(top15phy$quadrant)
top15phy$quadrant[top15phy$quadrant=="Particle Epilimnion"] <- "Epilimnion   "
top15phy$quadrant[top15phy$quadrant=="Free Epilimnion"]   <- "Epilimnion"
top15phy$quadrant[top15phy$quadrant=="Particle Hypolimnion"]   <- " Hypolimnion  "
top15phy$quadrant[top15phy$quadrant=="Free Hypolimnion"]   <- "Hypolimnion"

#### ADDING NO SPACES FOR PRODUCTIVE and ADDING 1 SPACE FOR UNPRODUCTIVE
top15phy$prod <- top15phy$ProdLevel
top15phy$prod <- as.character(top15phy$prod)
top15phy$prod[top15phy$prod =="Productive"]   <- ""
top15phy$prod[top15phy$prod =="Unproductive"]   <- " "

for(i in 1:length(top15phy$prod)){
  top15phy$prod_quad [i]<-paste(as.character(top15phy$prod[i]),as.character(top15phy$quadrant[i]))}

top15phy$prod_quad <- factor(top15phy$prod_quad,levels = c(" Epilimnion", " Hypolimnion", " Epilimnion   ", "  Hypolimnion  ", 
                                                           "  Epilimnion", "  Hypolimnion", "  Epilimnion   ",  "   Hypolimnion  "))

phy.colors.whew <- c(Acidobacteria = "grey26", Actinobacteria = "palevioletred2", Alphaproteobacteria = "steelblue4", Armatimonadetes = "red", Bacteroidetes = "darkorange",
                     "BD1-5" = "chartreuse", Betaproteobacteria = "royalblue", Caldiserica = "black","Candidate_division_BRC1" = "red",
                     "Candidate_division_JS1" = "aquamarine1",
                     "Candidate_division_OD1" = "#6DDE88", "Candidate_division_OP3" = "yellow1", "Candidate_division_OP8" = "goldenrod1", "Candidate_division_OP11" = "chocolate4",
                     "Candidate_division_SR1" = "tan3", "Candidate_division_TM7" = "skyblue1", "Candidate_division_WS3" = "magenta",
                     Chlamydiae="violet", Chlorobi="cyan2", Chloroflexi="darkgreen", Cyanobacteria = "chartreuse3", 
                     Deferribacteres = "slateblue3", "Deinococcus-Thermus" = "violetred", Dictyoglomi = "cornsilk4", Deltaproteobacteria = "deepskyblue", 
                     Elusimicrobia = "violetred4", Epsilonproteobacteria = "lightskyblue", Fibrobacteres = "hotpink", Firmicutes = "blue4", FGL7S = "palevioletred1",
                     Fusobacteria = "slateblue1", Gammaproteobacteria = "plum2", Gemmatimonadetes="black", GOUTA4 = "plum1", "Hyd24-12" = "sienna2", JTB23 = "seashell2",
                     Lentisphaerae = "yellow1", "NPL-UPA2"="royalblue", OC31 = "mediumpurple4", Planctomycetes = "mediumorchid3", Proteobacteria = "deepskyblue",
                     "SHA-109" = "lightsalmon3", SM2F11 = "lightskyblue2", SPOTSOCT00m83 = "orangered",
                     Spirochaetae = "gold3", Tenericutes="pink", Thermotogae = "chocolate1", TA06 = "lightslateblue",TA18 = "rosybrown3", TM6 = "olivedrab",
                     unclassified = "grey", Verrucomicrobia = "purple4", "WCHB1-60" = "palegreen")

top15plot <- ggplot(top15phy, aes(y=mean_abundance , x=Phylum, fill=Phylum)) + #coord_cartesian(xlim = c(0, 30)) + 
  guides(fill = guide_legend(reverse=TRUE)) + ggtitle("") + 
  geom_bar(stat="identity", position=position_dodge()) + #theme_classic() +
  geom_bar(stat="identity", position=position_dodge(), colour = "black", show_guide = FALSE) +
  facet_wrap(~ prod_quad, ncol = 8) + xlab("Top 15 Most Abundant Phyla") +
  scale_fill_manual(values = phy.colors.whew,name="Phylum") + 
  geom_errorbar(aes(ymin = mean_abundance -se, ymax = mean_abundance +se), width = 0.25, color = "black") + 
  ylab("Mean Percent Relative Abundance (%)") + coord_flip() + theme_bw() + 
  theme(axis.title.x = element_text(face="bold", size=16),
        axis.text.x = element_text(colour = "black", size=14),
        axis.text.y = element_text(colour = "black", size=14),
        axis.title.y = element_text(face="bold", size=16),
        plot.title = element_text(face="bold", size = 20),
        strip.text.x = element_text(size=12),
        strip.text.y = element_text(size=12),
        legend.title = element_text(size=12, face="bold"),
        legend.text = element_text(size = 12),
        strip.background = element_rect(fill="grey77", color = NA),
        legend.position="right"); 


library(gtable)
# get gtable object
z <- ggplot_gtable(ggplot_build(top15plot))

# add label for top strip
z <- gtable_add_rows(z, z$heights[[3]], 2)
z <- gtable_add_grob(z, list(rectGrob(gp = gpar(col = NA, linetype = 1, fill = gray(0.7))),
                             textGrob("Free-Living", gp = gpar(fontsize = 14, col = "black"))),
                     3, 4, 3, 7, name = paste(runif(2)))

z <- gtable_add_grob(z, list(rectGrob(gp = gpar(col = NA, linetype = 1, fill = gray(0.7))),
                             textGrob("Particle-Associated", gp = gpar(fontsize = 14, col = "black"))),
                     3, 9, 3, 13, name = paste(runif(2)))

z <- gtable_add_grob(z, list(rectGrob(gp = gpar(col = NA, linetype = 1, fill = gray(0.6))),
                             textGrob("High-Nutrient", gp = gpar(fontsize = 14, col = "black"))),
                     2, 4, 2, 13, name = paste(runif(2)))

z <- gtable_add_grob(z, list(rectGrob(gp = gpar(col = NA, linetype = 1, fill = gray(0.7))),
                             textGrob("Free-Living", gp = gpar(fontsize = 14, col = "black"))),
                     3, 15, 3, 19, name = paste(runif(2)))

z <- gtable_add_grob(z, list(rectGrob(gp = gpar(col = NA, linetype = 1, fill = gray(0.7))),
                             textGrob("Particle-Associated", gp = gpar(fontsize = 14, col = "black"))),
                     3, 22, 3, 25, name = paste(runif(2)))

z <- gtable_add_grob(z, list(rectGrob(gp = gpar(col = NA, linetype=1, fill = gray(0.6))),
                             textGrob("Low-Nutrient", gp = gpar(fontsize = 14, col = "black"))),
                     2, 15, 2, 25, name = paste(runif(2)))

# add margins
z <- gtable_add_cols(z, unit(1/8, "line"), 7)
z <- gtable_add_rows(z, unit(1/8, "line"), 3)

#####  Plotting FIGURE S3  #####  Plotting FIGURE S3  #####  Plotting FIGURE S3  #####  Plotting FIGURE S3  #####  Plotting FIGURE S3
#tiff(filename="~/Final_PAFL_Trophicstate/Final_Figures/Fig.S3_Abundance_top15_gtable.tiff", width= 40, height=25, units= "cm", pointsize= 10, res=200)
# draw it
grid.newpage()
grid.draw(z)
#dev.off()




###  MIDDLE 19 PHYLA 
mid20phy <- prod_quad_phylum_stats[prod_quad_phylum_stats$Phylum %in% phy_order[16:33], ]

mid20phy$Phylum <- factor(mid20phy$Phylum,levels = c("Acidobacteria", "Spirochaetae",  "Candidate_division_OD1", "NPL-UPA2", "Deinococcus-Thermus", "Candidate_division_OP3", 
                                                     "TM6", "TA18", "Epsilonproteobacteria", "Chlamydiae",  "Fibrobacteres", "Candidate_division_SR1",
                                                     "Candidate_division_BRC1", "BD1-5",  "Candidate_division_WS3", "Tenericutes", "Elusimicrobia","Gemmatimonadetes", "WCHB1-60", "Fusobacteria"))

#### ADDING NO SPACES FOR FREE LIVING and ADDING 2 SPACE FOR PARICLE ASSOCIATED AT THE END
mid20phy$quadrant <- as.character(mid20phy$quadrant)
mid20phy$quadrant[mid20phy$quadrant=="Particle Epilimnion"] <- "Epilimnion   "
mid20phy$quadrant[mid20phy$quadrant=="Free Epilimnion"]   <- "Epilimnion"
mid20phy$quadrant[mid20phy$quadrant=="Particle Hypolimnion"]   <- " Hypolimnion  "
mid20phy$quadrant[mid20phy$quadrant=="Free Hypolimnion"]   <- "Hypolimnion"

#### ADDING NO SPACES FOR PRODUCTIVE and ADDING 1 SPACE FOR UNPRODUCTIVE
mid20phy$prod <- mid20phy$ProdLevel
mid20phy$prod <- as.character(mid20phy$prod)
mid20phy$prod[mid20phy$prod =="Productive"]   <- ""
mid20phy$prod[mid20phy$prod =="Unproductive"]   <- " "

for(i in 1:length(mid20phy$prod)){
  mid20phy$prod_quad [i]<-paste(as.character(mid20phy$prod[i]),as.character(mid20phy$quadrant[i]))}

mid20phy$prod_quad <- factor(mid20phy$prod_quad,levels = c(" Epilimnion", " Hypolimnion", " Epilimnion   ", "  Hypolimnion  ", 
                                                           "  Epilimnion", "  Hypolimnion", "  Epilimnion   ",  "   Hypolimnion  "))


mid20phy$quadrant <- factor(mid20phy$quadrant,levels = c("Free Epilimnion", "Particle Epilimnion", "Free Hypolimnion", "Particle Hypolimnion"))

mid20phy$Phylum = with(mid20phy, factor(Phylum, levels = rev(levels(Phylum))))

mid20phy_plot <- ggplot(mid20phy, aes(y=mean_abundance , x=Phylum, fill=Phylum)) + #coord_cartesian(xlim = c(0, 30)) + 
  guides(fill = guide_legend(reverse=TRUE)) + ggtitle("") + 
  geom_bar(stat="identity", position=position_dodge()) + #theme_classic() +
  geom_bar(stat="identity", position=position_dodge(), colour = "black", show_guide = FALSE) +
  facet_wrap(~ prod_quad, ncol = 8) + xlab("Middle 18 Most Abundant Phyla") +
  scale_fill_manual(values = phy.colors.whew,name="Phylum") + 
  geom_errorbar(aes(ymin = mean_abundance -se, ymax = mean_abundance +se), width = 0.25, color = "black") + 
  ylab("Mean Percent Relative Abundance (%)") + coord_flip() + theme_bw() + 
  scale_y_continuous(breaks=seq(0, 1.6, 0.4), lim = c(0, 1.7)) + 
  theme(axis.title.x = element_text(face="bold", size=16),
        axis.text.x = element_text(colour = "black", size=14, vjust = 1, hjust = 1, angle = 60),
        axis.text.y = element_text(colour = "black", size=14),
        axis.title.y = element_text(face="bold", size=16),
        plot.title = element_text(face="bold", size = 20),
        strip.text.x = element_text(size=12),
        strip.text.y = element_text(size=12),
        legend.title = element_text(size=12, face="bold"),
        legend.text = element_text(size = 12),
        strip.background = element_rect(fill="grey77", color = NA),
        legend.position="right"); 



library(gtable)
# get gtable object
z <- ggplot_gtable(ggplot_build(mid20phy_plot))

# http://stackoverflow.com/questions/22818061/annotating-facet-title-as-strip-over-facet
# http://stackoverflow.com/questions/11353287/how-do-you-add-a-general-label-to-facets-in-ggplot2
# http://stackoverflow.com/questions/11442981/ggplot2-strip-text-labels-facet-wrap
# add label for Epilimnio strip
z <- gtable_add_rows(z, z$heights[[3]], 2)
z <- gtable_add_grob(z, list(rectGrob(gp = gpar(col = NA, linetype = 1, fill = gray(0.7))),
                             textGrob("Free-Living", gp = gpar(fontsize = 14, col = "black"))),
                     3, 4, 3, 7, name = paste(runif(2)))

z <- gtable_add_grob(z, list(rectGrob(gp = gpar(col = NA, linetype = 1, fill = gray(0.7))),
                             textGrob("Particle-Associated", gp = gpar(fontsize = 14, col = "black"))),
                     3, 9, 3, 13, name = paste(runif(2)))

z <- gtable_add_grob(z, list(rectGrob(gp = gpar(col = NA, linetype = 1, fill = gray(0.6))),
                             textGrob("High-Nutrient", gp = gpar(fontsize = 14, col = "black"))),
                     2, 4, 2, 13, name = paste(runif(2)))

z <- gtable_add_grob(z, list(rectGrob(gp = gpar(col = NA, linetype = 1, fill = gray(0.7))),
                             textGrob("Free-Living", gp = gpar(fontsize = 14, col = "black"))),
                     3, 15, 3, 19, name = paste(runif(2)))

z <- gtable_add_grob(z, list(rectGrob(gp = gpar(col = NA, linetype = 1, fill = gray(0.7))),
                             textGrob("Particle-Associated", gp = gpar(fontsize = 14, col = "black"))),
                     3, 22, 3, 25, name = paste(runif(2)))

z <- gtable_add_grob(z, list(rectGrob(gp = gpar(col = NA, linetype=1, fill = gray(0.6))),
                             textGrob("Low-Nutrient", gp = gpar(fontsize = 14, col = "black"))),
                     2, 15, 2, 25, name = paste(runif(2)))

# add margins
z <- gtable_add_cols(z, unit(1/8, "line"), 7)
z <- gtable_add_rows(z, unit(1/8, "line"), 3)

#####  Plotting FIGURE S4  #####  Plotting FIGURE S4  #####  Plotting FIGURE S4  #####  Plotting FIGURE S4  #####  Plotting FIGURE S4
#tiff(filename="~/Final_PAFL_Trophicstate/Final_Figures/Fig.S4_Abundance_mid20_gtable.tiff", width= 40, height=25, units= "cm", pointsize= 10, res=200)
# draw it
grid.newpage()
grid.draw(z)
#dev.off()
