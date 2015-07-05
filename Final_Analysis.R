#Load libraries
#library(phyloseq)  #must load the devtools version
#library("devtools")
#install_github("phyloseq", "joey711")
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

# Our phyloseq object!
merged_final 



####################################################  DEPTH PROFILE  ####################################################
####################################################  DEPTH PROFILE  ####################################################
####################################################  DEPTH PROFILE  ####################################################

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

### SUMMARIZING WITH DDPLY 
# STATS ON DO
ave_prof_DO <- ddply(profile, c("depth", "trophicstate"), summarise, 
                     N = length(DO),
                     mean = mean(DO),
                     sd   = sd(DO),
                     se   = sd / sqrt(N))
ave_prof_DO$Variable <- "DO"

# STATS ON TEMP
ave_prof_temp <- ddply(profile, c("depth", "trophicstate"), summarise, 
                       N = length(temp),
                       mean = mean(temp),
                       sd   = sd(temp),
                       se   = sd / sqrt(N))
ave_prof_temp$Variable <- "temp"

# STATS ON pH
ave_prof_pH <- ddply(profile, c("depth", "trophicstate"), summarise, 
                     N = length(pH),
                     mean = mean(pH),
                     sd   = sd(pH),
                     se   = sd / sqrt(N))
ave_prof_pH$Variable <- "pH"

# STATS ON SpC
ave_prof_SpC <- ddply(profile, c("depth", "trophicstate"), summarise, 
                      N = length(SpC),
                      mean = mean(SpC),
                      sd   = sd(SpC),
                      se   = sd / sqrt(N))
ave_prof_SpC$Variable <- "SpC"

mean_profile <- rbind(ave_prof_DO, ave_prof_temp, ave_prof_pH, ave_prof_SpC)

profile_labeller <- function(var, value){
  value <- as.character(value)
  if (var=="Variable") { 
    value[value=="temp"] <- "Temperature \n (Celsius)"
    value[value=="DO"]   <- "Dissolved Oxygen \n (mg/L)"
    value[value=="pH"]   <- "pH"
    value[value=="SpC"]   <- "Specific Conductivity \n (uS/cm)"
  }
  return(value)
}

mean_profile$Variable <-factor(mean_profile$Variable,levels=c("temp", "DO", "pH", "SpC"))
mean_profile$trophicstate <-factor(mean_profile$trophicstate,levels=c("Productive", "Unproductive", "Mixed"))
mean_profile13 <- subset(mean_profile, depth < 13.1)

## AVERAGE PLOT!
#jpeg(filename="~/Final_PAFL_Trophicstate/Figures/Fig.1_Average_PROD_profiles_13m_SE.jpeg", width= 40, height=30, units= "cm",pointsize= 18, res=500)
ggplot(mean_profile13, aes(x=mean, y = depth, color = trophicstate)) +   
  facet_grid(. ~ Variable, scales = "free", labeller = profile_labeller) +  
  geom_path(size=2, alpha = 0.8) + ylab("Depth (m)") + xlab("") + 
  theme_bw() +  geom_point(size=4, alpha = 0.8) + geom_hline(h=0) +
  scale_y_reverse(breaks=seq(0, 30, 2), lim = c(14,0), expand = c(0, 0)) + 
  geom_errorbarh(width=.1, aes(xmin = mean - se, xmax = mean + se)) +
  scale_color_manual(name = "", breaks = c("Productive", "Unproductive", "Mixed"), 
                     labels = c("Productive", "Unproductive", "Mixed"), 
                     values = c("deeppink","turquoise3", "blue")) +
  #scale_color_manual(name = "", breaks = c("Eutrophic", "Mesotrophic","Oligotrophic", "Mixed"), 
  #                   labels = c("Eutrophic", "Mesotrophic","Oligotrophic", "Mixed"), 
  #                   values = c("deeppink", "orange","turquoise3", "blue")) +
  theme(axis.title.x = element_text(face="bold", size=14),
        axis.text.x = element_text(colour = "black",size=14),
        axis.text.y = element_text(colour = "black", size=14),
        axis.title.y = element_text(face="bold", size=16),
        legend.title = element_text(size=12, face="bold"),
        legend.text = element_text(size = 12),
        strip.text.x = element_text(size = 16, face = "bold", colour = "black"),
        strip.background = element_blank(),
        legend.position = c(0.81, 0.08));
#dev.off()

# ALL LAKES 
library(tidyr)
library(dplyr)
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

## All LAKES
#jpeg(filename="~/Final_PAFL_Trophicstate/Figures/Fig.S1_AllLake_profiles.jpeg", width= 40, height=30, units= "cm",pointsize= 18, res=500)
ggplot(profile_all, aes(x=Value, y = depth, color = lakename)) +   
  facet_grid(. ~ Variable, scales = "free", labeller = profile_labeller) +  
  geom_path(size=2, alpha = 0.8) + ylab("Depth (m)") + xlab("") + 
  theme_bw() +  geom_point(size=4, alpha = 0.8) + geom_hline(h=0) +
  scale_y_reverse(breaks=seq(0, 30, 5), lim = c(30,0), expand = c(0, 0)) + 
  scale_color_manual(name = "Lake Name", 
                     labels = c("Baker", "Baseline", "Bassett", "Bristol", "Gull", "Lee", "Little Long", "Payne", "Sherman", "Sixteen", "Wintergreen"), 
                     values = c("blue", "red", "black", "forestgreen","violet","limegreen", "purple", "orange", "maroon1", "turquoise3", "thistle4")) +
  theme(axis.title.x = element_text(face="bold", size=14),
        axis.text.x = element_text(colour = "black",size=14),
        axis.text.y = element_text(colour = "black", size=14),
        axis.title.y = element_text(face="bold", size=16),
        legend.title = element_text(size=12, face="bold"),
        legend.text = element_text(size = 12),
        strip.text.x = element_text(size = 16, face = "bold", colour = "black"),
        strip.background = element_blank(),
        legend.position=c(0.81, 0.215));
#dev.off()


####################################################  ORDINATIONS  ####################################################
####################################################  ORDINATIONS  ####################################################
####################################################  ORDINATIONS  ####################################################
########  NMDS + PCoA Plots
### Get rid of the Wintergreen HYPOLIMNION Samples
nowin_merged <- subset_samples(merged_final, names != "WINH" & names != "WINH3um")

#########  PLOT ORDINATIONS FOR BOTH SCALED AND MANUAL
nowinOTU <- otu_table(nowin_merged)
#weighted
norm_bray <- vegdist(nowinOTU, method = "bray")  # calculates the Bray-Curtis Distances
bray_pcoa <- pcoa(norm_bray)
bray_pcoa2 <- bray_pcoa$vectors
bray_pcoa3 <- data.frame(bray_pcoa2[, 1:3])
bray_pcoa3$names <- row.names(bray_pcoa3)
bray_pcoa4 <- makeCategories_dups(bray_pcoa3)
bray_pcoa4$quadrant <- factor(bray_pcoa4$quadrant,levels = c("Free Epilimnion", "Free Mixed",  "Free Hypolimnion", "Particle Epilimnion", "Particle Mixed", "Particle Hypolimnion"))

pcoa_BC <- ggplot(bray_pcoa4, aes(Axis.1, Axis.2 * -1, color = quadrant, shape = trophicstate)) +
  xlab("PCoA1 [23.1%]") + ylab("PCoA2 [17.9%]") + #ggtitle("Bray-Curtis: Trophic State") +
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
        legend.position = "right"); pcoa_BC

##########  SORENSEN'S DISSIMILARITY
nowinOTU <- as.matrix(nowinOTU)
norm_soren <- vegdist(nowinOTU, method = "bray", binary = TRUE)
soren_pcoa <- pcoa(norm_soren)
soren_pcoa2 <- soren_pcoa$vectors
soren_pcoa3 <- data.frame(soren_pcoa2[, 1:3])
soren_pcoa3$names <- row.names(soren_pcoa3)
soren_pcoa4 <- makeCategories_dups(soren_pcoa3)
soren_pcoa4$quadrant <- factor(soren_pcoa4$quadrant,levels = c("Free Epilimnion", "Free Mixed",  "Free Hypolimnion", "Particle Epilimnion", "Particle Mixed", "Particle Hypolimnion"))

pcoa_soren <- ggplot(soren_pcoa4, aes(Axis.1, Axis.2 * -1, color = quadrant, shape = trophicstate)) +
  xlab("PCoA1 [23.1%]") + ylab("PCoA2 [17.9%]") + #ggtitle("Bray-Curtis: Trophic State") +
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


############## NMDS NMDS NMDS NMDS NDMS
nmds_bray <- metaMDS(normOTU, distance="bray")  #, autotransform = FALSE
nmds_bray <- data.frame(nmds_bray$points) #http://strata.uga.edu/software/pdf/mdsTutorial.pdf
nmds_bray$names<-row.names(nmds_bray) #new names column
nmds_bray <- makeCategories_dups(nmds_bray) #will add our categorical information:  lakenames, limnion, filter, quadrant and trophicstate
nmds_bray$quadrant <- factor(nmds_bray$quadrant,levels = c("Free Epilimnion", "Free Mixed",  "Free Hypolimnion", "Particle Epilimnion", "Particle Mixed", "Particle Hypolimnion"))

nmds_bc_quad <- ggplot(nmds_bray, aes(MDS1, MDS2, color = quadrant, shape = trophicstate)) +
  xlab("NMDS1") + ylab("NMDS2") + #ggtitle("Bray-Curtis: Scaled") +
  geom_point(size= 6, alpha=0.9) + theme_bw() + 
  annotate("text", label = " Stress = 0.17", x = (max(nmds_bray$MDS1) -0.25), y = (max(nmds_bray$MDS2) -0.02), size = 6, colour = "black") +
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
        legend.position = "right");  nmds_bc_quad


# UNWEIGHTED
nmds_soren <- metaMDS(normOTU, distance="bray", binary = TRUE)
nmds_soren <- data.frame(nmds_soren$points) #http://strata.uga.edu/software/pdf/mdsTutorial.pdf
nmds_soren$names<-row.names(nmds_soren) #new names column
nmds_soren <- makeCategories_dups(nmds_soren) #will add our categorical information:  lakenames, limnion, filter, quadrant and trophicstate
nmds_soren$quadrant <- factor(nmds_soren$quadrant,levels = c("Free Epilimnion", "Free Mixed",  "Free Hypolimnion", "Particle Epilimnion", "Particle Mixed", "Particle Hypolimnion"))

nmds_soren_quad <- ggplot(nmds_soren, aes(MDS1, MDS2, color = quadrant, shape = trophicstate)) +
  xlab("NMDS1") + ylab("NMDS2") + #ggtitle("Sorensen: Scaled") +
  geom_point(size= 6, alpha=0.9) + theme_bw() + 
  annotate("text", label = " Stress = 0.14", x = (max(nmds_soren$MDS1) -0.25), y = (max(nmds_soren$MDS2) -0.02), size = 6, colour = "black") +
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
        legend.position = "right");  nmds_soren_quad


multiplot(nmds_bc_quad, nmds_soren_quad, cols = 2)



####################################################  ALPHA DIVERSITY  ####################################################
####################################################  ALPHA DIVERSITY  ####################################################
####################################################  ALPHA DIVERSITY  ####################################################
# For alpha diversity we will RAREFY our phyloseq object.  Our raw-data phyloseq object was:
raw_merged  # Raw Merged Samples, we have 14,378 OTUs
raw_nowin <- subset_samples(raw_merged, names != "WINH" & names != "WINH3um")
raw_nowin <- prune_taxa(taxa_sums(raw_nowin) > 0, raw_nowin)

#Data Import for rarefied data 
raw_nowin  
min(sample_sums(raw_nowin))

# Following code from Michelle's Butterflygut website: http://rstudio-pubs-static.s3.amazonaws.com/25722_9fc9bdc48f0d4f13b05fa61baeda57a0.html#alpha-diversity
# Rarefy to 14937 reads with replacement 100 times to estimate species richness
# Since we are rarefying to 14937 reads we need to remove the two samples with less than 1000 reads
raredata14936 <- prune_samples(sample_sums(raw_nowin) > 14936, raw_nowin)

# Initialize matrices to store richness and evenness estimates
richness <- matrix(nrow = 41,ncol = 100)  # Store richness:  We have 41 samples 
row.names(richness) <- sample_names(raredata14936)
evenness <- matrix(nrow = 41,ncol = 100)  #Store evenness
row.names(evenness) <- sample_names(raredata14936)

# We want to be reproducible - so let's set the seed.
set.seed(3)

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
rich_stats$Test <- "ObsRich"

# Create a new matrix to hold the means and standard deviations of the evenness estimates
even_stats = matrix(nrow = nrow(evenness), ncol = 1)
even_stats[,1] = apply(evenness, 1, mean)
#even_stats[,2] = apply(evenness, 1, sd)
even_stats = data.frame(row.names(evenness), even_stats)
colnames(even_stats) = c("samples","mean")
even_stats$Test <- "InvSimpson"


###  Add SIMPSON'S EVENNESS!
#simps_even = as.data.frame(matrix(nrow = nrow(evenness), ncol = 2))
simps_even <- data.frame(matrix(nrow = nrow(evenness), ncol = 3))
colnames(simps_even) = c("samples","mean", "Test")
simps_even$samples <- even_stats$samples
simps_even$mean <-  even_stats$mean/rich_stats$mean
simps_even$Test <- "SimpsEven"

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



scaled_otu <- otu_table(merged_final)  
norm_scaled_bc <- vegdist(scaled_otu, method = "bray", binary = FALSE)   # calculates the Bray-Curtis Distances

