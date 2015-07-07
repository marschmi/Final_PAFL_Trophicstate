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
library(scales)

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
mean_profile13b <- subset(mean_profile13, Variable != "SpC")

## AVERAGE PLOT!
#jpeg(filename="~/Final_PAFL_Trophicstate/Figures/Fig.1_Average_PROD_profiles_13m_SE_NOspc.jpeg", width= 32, height=30, units= "cm",pointsize= 18, res=500)
ggplot(mean_profile13b, aes(x=mean, y = depth, color = trophicstate)) +   
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
        legend.text = element_text(size = 14),
        strip.text.x = element_text(size = 16, face = "bold", colour = "black"),
        strip.background = element_blank(),
        legend.position = c(0.1, 0.93));
        #legend.position = c(0.81, 0.08));
#dev.off()

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
### Clustering 
### Get rid of the Wintergreen HYPOLIMNION Samples
nowin_merged <- subset_samples(merged_final, names != "WINH" & names != "WINH3um")
nowin_merged <- prune_taxa(taxa_sums(nowin_merged) > 0, nowin_merged)
otu <- otu_table(nowin_merged)
otu_bray <- vegdist(otu, method = "bray")  # calculates the Bray-Curtis Distances
df_otu <- data.frame(otu)
otu_soren <- vegdist(df_otu, method = "bray", binary = TRUE)  # calculates the Bray-Curtis Distances


### Clustering based on 
#jpeg(filename="clustering_bray+soren.jpeg", width= 45, height=32, units= "cm", pointsize= 14, res=500)
par(mfrow = c(2,1))
plot(hclust(otu_bray), main = "Bray-Curtis Distance")
plot(hclust(otu_soren), main = "Sorensen Distance")
#dev.off()

# vector of colors
mypal = c("black", "red", "blue", "purple", "green3", "orange", "maroon1", "gold")
# cutting dendrogram in 5 clusters
hc <- hclust(otu_bray)
clus5 = cutree(hc, 5)
par(mfrow = c(1,1))
plot(as.phylo(hc), tip.color = mypal[clus5],  main = "Bray Curtis Dissimilarity") 

plot(as.phylo(hc), main = "", xlab = "Samples", sub = "", tip.color = mypal[clus5], ylab = "Bray Curtis Dissimilarity") 




########  NMDS + PCoA Plots
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
nowinOTU2 <- data.frame(otu_table(nowin_merged))
nowinOTU3 <- t(nowinOTU2)
norm_soren <- vegdist(nowinOTU3, method = "bray", binary = TRUE)

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
set.seed(3)
nmds_bray <- metaMDS(nowinOTU, distance="bray")  #, autotransform = FALSE
nmds_bray <- data.frame(nmds_bray$points) #http://strata.uga.edu/software/pdf/mdsTutorial.pdf
nmds_bray$names<-row.names(nmds_bray) #new names column
nmds_bray <- makeCategories_dups(nmds_bray) #will add our categorical information:  lakenames, limnion, filter, quadrant and trophicstate
nmds_bray$quadrant <- factor(nmds_bray$quadrant,levels = c("Free Epilimnion", "Free Mixed",  "Free Hypolimnion", "Particle Epilimnion", "Particle Mixed", "Particle Hypolimnion"))

nmds_bc_quad <- ggplot(nmds_bray, aes(MDS1*10, MDS2*10, color = quadrant, shape = trophicstate)) +
  xlab("NMDS1") + ylab("NMDS2") + ggtitle("Bray-Curtis Dissimilarity") +
  geom_point(size= 6, alpha=0.9) + theme_bw() +   geom_point(colour="white", size = 2) +
  annotate("text", label = " Stress = 0.17", x = (max(nmds_bray$MDS1*10) -0.2), y = (max(nmds_bray$MDS2*10) -0.02), size = 6, colour = "black") +
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
        plot.title = element_text(size = 16, face="bold"),
        ###LEGEND TOP RIGHT CORNER
        legend.position = "none");  nmds_bc_quad


# UNWEIGHTED
set.seed(3)
nmds_soren <- metaMDS(nowinOTU, distance="bray", binary = TRUE)
nmds_soren <- data.frame(nmds_soren$points) #http://strata.uga.edu/software/pdf/mdsTutorial.pdf
nmds_soren$names<-row.names(nmds_soren) #new names column
nmds_soren <- makeCategories_dups(nmds_soren) #will add our categorical information:  lakenames, limnion, filter, quadrant and trophicstate
nmds_soren$quadrant <- factor(nmds_soren$quadrant,levels = c("Free Epilimnion", "Free Mixed",  "Free Hypolimnion", "Particle Epilimnion", "Particle Mixed", "Particle Hypolimnion"))

nmds_soren_quad <- ggplot(nmds_soren, aes(MDS1, MDS2, color = quadrant, shape = trophicstate)) +
  xlab("NMDS1") + ylab("NMDS2") + ggtitle("Sorensen Dissimilarity") +
  geom_point(size= 6, alpha=0.9) + theme_bw() +   geom_point(colour="white", size = 2) +
  annotate("text", label = " Stress = 0.14", x = (max(nmds_soren$MDS1) -0.18), y = (max(nmds_soren$MDS2) -0.01), size = 6, colour = "black") +
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
        plot.title = element_text(size = 16, face="bold"),
        ###LEGEND TOP RIGHT CORNER
        legend.position = "right");  nmds_soren_quad


#multiplot(nmds_bc_quad, nmds_soren_quad, cols = 2)


#jpeg(filename="~/Final_PAFL_Trophicstate/Figures/Fig.2_NMDS_bc+soren.jpeg", width= 45, height=18, units= "cm", pointsize= 14, res=500)
grid.newpage()
pushViewport(viewport(layout=grid.layout(1,2,width=c(0.42,0.58))))
print(nmds_bc_quad, vp=viewport(layout.pos.row=1,layout.pos.col=1))
print(nmds_soren_quad, vp=viewport(layout.pos.row=1,layout.pos.col=2))
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


##  Calculate the BC dissimilarity and the Sorensen Dissimilarity
nosherwinOTU <- otu_table(nosherwin_merged)  # This is our OTU table that we will use for Adonis
#BCdist <- vegdist(nosherwinOTU, method = "bray", binary = FALSE)  # calculates the Bray-Curtis Distances
df_nosherwinOTU <- data.frame(nosherwinOTU)
otu_soren <- vegdist(df_nosherwinOTU, method = "bray", binary = TRUE)  ##SORENSEN DISTANCE --> Test's the presence/absence and makes 
BCdist <- otu_soren  ### This way we don't need to re-type the code for adonis, BUT - **BE CAREFUL** with this!  

##  Do you have the right object for BCdist (BC vs sorensen???)
##  Are you sure?!  

# #Run an ADONIS test!
adonis_PAFL_mult <- adonis(BCdist ~ limnion+filter +trophicstate+ DO+temp + pH, data=environ) # R2 = 0.36245
adonis_PAFL_quad <- adonis(BCdist~quadrant,data=environ) #  R2 = 
adonis_PAFL_filt <- adonis(BCdist~filter,data=environ)  #R2 = 
adonis_PAFL_limnion <- adonis(BCdist~limnion,data=environ) #R2 = 
adonis_PAFL_DO <- adonis(BCdist~DO,data=environ) # R2 = 
adonis_PAFL_trophicstate <- adonis(BCdist~trophicstate, data = environ) # R2 = 
adonis_PAFL_prod <- adonis(BCdist~ProdLevel, data = environ) # R2 = 
adonis_PAFL_temp <- adonis(BCdist~ temp,data=environ) # R2 = 
adonis_PAFL_pH <- adonis(BCdist~ pH,data=environ) # R2 = 


#adonis_PAFL_cross <- adonis(BCdist~limnion*DO,data=environ) #R2 = 
#adonis_PAFL_DO <- adonis(BCdist~DO,data=environ) # R2 = 
#adonis_PAFL_whew <- adonis(BCdist ~ limnion*filter *trophicstate* DO*temp * pH, data=environ) # R2 = 0.36245




# adonis_PAFL_troph <- adonis(BCdist~trophicstate,data=environ)

# Epilimnion!
epi <- subset_samples(nosherwin_merged, limnion == "Epilimnion")
epiOTU <- otu_table(epi)
epi_BC <- vegdist(epiOTU, method = "bray")

df_epiOTU <- data.frame(epi)
epi_soren <- vegdist(df_epiOTU, method = "bray", binary = TRUE)  ##SORENSEN DISTANCE --> Test's the presence/absence and makes 
BCdist <- epi_soren  ### This way we don't need to re-type the code for adonis, BUT - **BE CAREFUL** with this!  

##  Do you have the right object for BCdist (BC vs sorensen???)
##  Are you sure?! 

epi_env <- subset(environ, limnion == "Epilimnion")  #Load Environmental Data

epi_adon_mult<- adonis(epiOTU ~ filter+trophicstate+ DO + temp + pH, data=epi_env) #R2 = 
epi_adon_filt <- adonis(epiOTU ~ filter, data=epi_env) # R2 = 
epi_adon_troph <- adonis(epiOTU ~ trophicstate, data=epi_env) # R2 =  
epi_adon_prod <- adonis(epiOTU ~ ProdLevel, data=epi_env) # R2 =  
epi_adon_DO <- adonis(epiOTU ~ DO, data=epi_env) # R2 = 0.086
epi_adon_temp <- adonis(epiOTU ~ temp, data=epi_env) # R2 = 0.086
epi_adon_pH <- adonis(epiOTU ~ pH, data=epi_env) # R2 = 0.086


# hypolimnion!
hypo <- subset_samples(nosherwin_merged, limnion == "Hypolimnion")
hypoOTU <- otu_table(hypo)
hypo_BC <- vegdist(hypoOTU, method = "bray")
hypo_env <- subset(environ, limnion == "Hypolimnion")
hypo_adon_mult <- adonis(hypoOTU ~ filter+trophicstate+DO + temp + pH, data=hypo_env) # R2 = 
#hypo_adon_mult2 <- adonis(hypoOTU ~ trophicstate + filter, data=hypo_env) # 
hypo_adon_filt <- adonis(hypoOTU ~ filter, data=hypo_env) # R2 = 
hypo_adon_troph <- adonis(hypoOTU ~ trophicstate, data=hypo_env) # R2 = 
hypo_adon_prod <- adonis(hypoOTU ~ ProdLevel, data=hypo_env) # R2 = 
hypo_adon_DO <- adonis(hypoOTU ~ DO, data=hypo_env) # R2 = 
hypo_adon_temp <- adonis(hypoOTU ~ temp, data=hypo_env) # R2 = 
hypo_adon_pH <- adonis(hypoOTU ~ pH, data=hypo_env) # R2 = 


# Oligotrophic!!
oligo <- subset_samples(nosherwin_merged, trophicstate == "Oligotrophic")
oligoOTU <- otu_table(oligo)
oligo_BC <- vegdist(oligoOTU, method = "bray")
oligo_env <- subset(environ, trophicstate == "Oligotrophic")
oligo_adon_filt <- adonis(oligo_BC ~ filter, data=oligo_env) # R2 = 
oligo_adon_limnion <- adonis(oligo_BC ~ limnion, data=oligo_env) # R2 = 
oligo_adon_DO <- adonis(oligo_BC ~ DO, data=oligo_env) # R2 = 
oligo_adon_temp <- adonis(oligo_BC ~ temp, data=oligo_env) # R2 = 
oligo_adon_pH <- adonis(oligo_BC ~ pH, data=oligo_env) # R2 = 
oligo_adon_mult <- adonis(oligo_BC ~ limnion+filter+DO + temp + pH, data=oligo_env) # R2 = 


# MESO + EUTROPHIC!!
prod <- subset_samples(nosherwin_merged, ProdLevel == "Productive")
prodOTU <- otu_table(prod)
prod_BC <- vegdist(prodOTU, method = "bray")
prod_env <- subset(environ, ProdLevel == "Productive")
prod_adon_filt <- adonis(prod_BC ~ filter, data=prod_env) # R2 = 
prod_adon_limnion <- adonis(prod_BC ~ limnion, data=prod_env) # R2 = 
prod_adon_DO <- adonis(prod_BC ~ DO, data=prod_env) # R2 = 
prod_adon_temp <- adonis(prod_BC ~ temp, data=prod_env) # R2 = 
prod_adon_pH <- adonis(prod_BC ~ pH, data=prod_env) # R2 = 
prod_adon_mult <- adonis(prod_BC ~ limnion+filter+DO + temp + pH, data=prod_env) # R2 = 


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
sherm_OTU <- otu_table(sherman_merged)
#nowinOTU # This is our OTU table that we will use for Adonis 
sherm_BCdist <- vegdist(sherm_OTU, method = "bray", binary = FALSE)
# 
# #Run an ADONIS test!
adonis_sherm_mult <- adonis(sherm_BCdist ~ filter+limnion+DO+temp+pH, data=sherm_environ, nperm = 999) # R2 = 0.36245
adonis_sherm_filt <- adonis(sherm_BCdist ~ filter, data=sherm_environ, nperm = 999) # R2 = 0.36245
adonis_sherm_limnion <- adonis(sherm_BCdist ~ limnion, data=sherm_environ, nperm = 999) # R2 = 0.36245
adonis_sherm_temp <- adonis(sherm_BCdist ~ temp, data=sherm_environ, nperm = 999) # R2 = 0.36245
adonis_sherm_pH <- adonis(sherm_BCdist ~ pH, data=sherm_environ, nperm = 999) # R2 = 0.36245
adonis_sherm_DO <- adonis(sherm_BCdist ~ DO, data=sherm_environ, nperm = 999) # R2 = 0.36245




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



####   ALPHA DIVERSITY:  COMBINING EUTROPHIC AND MESOTROPHIC
####   ALPHA DIVERSITY:  COMBINING EUTROPHIC AND MESOTROPHIC
####   ALPHA DIVERSITY:  COMBINING EUTROPHIC AND MESOTROPHIC
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





####################################################  BETA DIVERSITY  ####################################################
####################################################  BETA DIVERSITY  ####################################################
####################################################  BETA DIVERSITY  ####################################################

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
ddply_beta <- ddply(nomix_beta2, c("troph_lim2", "trophicstate1", "filter1"), summarise, 
                    N = length(value),
                    mean = mean(value),
                    sd   = sd(value),
                    se   = sd / sqrt(N))

#jpeg(filename="~/Final_PAFL_Trophicstate/Figures/Fig.3c_beta_TROPH_SD.jpeg", width= 25, height=15, units= "cm", pointsize= 14, res=500)
beta_plot <- ggplot(ddply_beta, aes(x = troph_lim2, y = mean, color = troph_lim2)) + geom_point(size = 5) +
  geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=.2, position=position_dodge(.9)) +
  scale_color_manual(name = "", limits=c("Eutrophic Epilimnion Particle", "Eutrophic Epilimnion Free", "Eutrophic Hypolimnion Particle", "Eutrophic Hypolimnion Free",
                                        "Mesotrophic Epilimnion Particle", "Mesotrophic Epilimnion Free", "Mesotrophic Hypolimnion Particle", "Mesotrophic Hypolimnion Free",
                                        "Oligotrophic Epilimnion Particle", "Oligotrophic Epilimnion Free", "Oligotrophic Hypolimnion Particle", "Oligotrophic Hypolimnion Free"), 
                    values = c("deeppink", "deeppink", "deeppink", "deeppink","orange","orange","orange","orange", 
                               "turquoise3","turquoise3","turquoise3","turquoise3"))+
  scale_x_discrete(breaks=c("Eutrophic Epilimnion Particle", "Eutrophic Epilimnion Free", "Eutrophic Hypolimnion Particle", "Eutrophic Hypolimnion Free",
                            "Mesotrophic Epilimnion Particle", "Mesotrophic Epilimnion Free", "Mesotrophic Hypolimnion Particle", "Mesotrophic Hypolimnion Free",
                            "Oligotrophic Epilimnion Particle", "Oligotrophic Epilimnion Free", "Oligotrophic Hypolimnion Particle", "Oligotrophic Hypolimnion Free"),
                   labels=c("Epilimnion Particle", "Epilimnion Free", "Hypolimnion Particle", "Hypolimnion Free",
                            "Epilimnion Particle", "Epilimnion Free", "Hypolimnion Particle", "Hypolimnion Free",
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
        legend.position="none"); beta_plot
#dev.off()


################# BRAY CURTIS DISTANCE FOR PRODUCTIVE VS UNPRODUCTIVE 
################# BRAY CURTIS DISTANCE FOR PRODUCTIVE VS UNPRODUCTIVE 
################# BRAY CURTIS DISTANCE FOR PRODUCTIVE VS UNPRODUCTIVE 
prod_beta <- nomix_beta2
unique(prod_beta$troph_lim1)
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
                     values = c("deeppink", "deeppink", "deeppink", "deeppink",
                                "turquoise3","turquoise3","turquoise3","turquoise3"))+
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




####################################################  ALPHA + BETA COMBINED DIVERSITY  ####################################################
####################################################  ALPHA + BETA COMBINED DIVERSITY  ####################################################
####################################################  ALPHA + BETA COMBINED DIVERSITY  ####################################################
prod_even <- subset(prodalpha_stats, Test == "Simpson's Evenness")
prod_richobs <-subset(prodalpha_stats, Test == "Observed Richness") 

prod_invsimps <- ggplot(prod_even, aes(x = troph_lim, y = Meantroph_lim, color = troph_lim)) + geom_point(size = 5) +
  facet_grid(. ~ trophicstate, scales="free", space="free_x") + scale_y_continuous(breaks=seq(0, 0.10, 0.02), lim = c(0, 0.1)) +
  geom_errorbar(aes(ymin=Meantroph_lim-SDtroph_lim, ymax=Meantroph_lim+SDtroph_lim), width=.2, position=position_dodge(.9)) +
  xlab("Habitat") + ylab("Simpson's Evenness") + theme_bw() + 
  scale_color_manual(name = "", limits=c("Productive Epilimnion Particle", "Productive Epilimnion Free", "Productive Hypolimnion Particle", "Productive Hypolimnion Free",
                                         "Unproductive Epilimnion Particle", "Unproductive Epilimnion Free", "Unproductive Hypolimnion Particle", "Unproductive Hypolimnion Free",
                                         "Mixed Mixed Particle", "Mixed Mixed Free"), 
                     values = c("deeppink", "deeppink", "deeppink", "deeppink",
                                "turquoise3","turquoise3","turquoise3","turquoise3", "blue", "blue"))+
  scale_x_discrete(breaks=c("Productive Epilimnion Particle", "Productive Epilimnion Free", "Productive Hypolimnion Particle", "Productive Hypolimnion Free",
                            "Unproductive Epilimnion Particle", "Unproductive Epilimnion Free", "Unproductive Hypolimnion Particle", "Unproductive Hypolimnion Free",
                            "Mixed Mixed Particle", "Mixed Mixed Free"), 
                   labels=c("", "", "", "",
                            "", "", "", "",
                            "Particle", "Free")) +
  theme(axis.title.x = element_blank(), axis.ticks.x = element_blank(), axis.text.x = element_blank(),
        axis.text.y = element_text(colour = "black", size=14),
        axis.title.y = element_text(face="bold", size=16),
        plot.title = element_text(face="bold", size = 20),
        strip.text.x = element_text(size=12, face="bold"),
        strip.text.y = element_text(size=14, face="bold"),
        strip.background = element_blank(), 
        plot.margin = unit(c(0.1, 0.2, -0.44, 0.225), "cm"), #top, right, bottom, left
        strip.background = element_rect(colour="black"),
        legend.position="none");prod_invsimps

prod_obs <- ggplot(prod_richobs, aes(x = troph_lim, y = Meantroph_lim, color = troph_lim)) + geom_point(size = 5) +
  facet_grid(. ~trophicstate, scales="free", space="free_x") + 
  geom_errorbar(aes(ymin=Meantroph_lim-SDtroph_lim, ymax=Meantroph_lim+SDtroph_lim), width=.2, position=position_dodge(.9)) +
  xlab("  ") + ylab("Observed Richness") + theme_bw() + scale_y_continuous(breaks=seq(400, 1200, 200), lim = c(300,1225)) +
  scale_color_manual(name = "", limits=c("Productive Epilimnion Particle", "Productive Epilimnion Free", "Productive Hypolimnion Particle", "Productive Hypolimnion Free",
                                         "Unproductive Epilimnion Particle", "Unproductive Epilimnion Free", "Unproductive Hypolimnion Particle", "Unproductive Hypolimnion Free",
                                         "Mixed Mixed Particle", "Mixed Mixed Free"), 
                     values = c("deeppink", "deeppink", "deeppink", "deeppink",
                                "turquoise3","turquoise3","turquoise3","turquoise3", "blue", "blue"))+
  scale_x_discrete(breaks=c("Productive Epilimnion Particle", "Productive Epilimnion Free", "Productive Hypolimnion Particle", "Productive Hypolimnion Free",
                            "Unproductive Epilimnion Particle", "Unproductive Epilimnion Free", "Unproductive Hypolimnion Particle", "Unproductive Hypolimnion Free",
                            "Mixed Mixed Particle", "Mixed Mixed Free"), 
                   labels=c("", "", "", "",
                            "", "", "", "",
                            "Particle-\nAssociated", "Free-\nLiving")) +
  theme(axis.title.x = element_text(face="bold", size=16), axis.ticks.x = element_blank(), 
        axis.text.x = element_text(angle=60, colour = "black", vjust=1, hjust = 1, size=14),
        axis.text.y = element_text(colour = "black", size=14),
        axis.title.y = element_text(face="bold", size=16),
        plot.title = element_text(face="bold", size = 20),
        plot.margin = unit(c(-0.24, 0.2, -0.75, 0.1), "cm"), #top, right, bottom, left
        strip.background = element_blank(), strip.text = element_blank(),
        legend.position="none"); prod_obs


prodbeta_plot2 <- ggplot(ddply_prodbeta, aes(x = troph_lim1, y = mean, color = troph_lim1)) + geom_point(size = 5) +
  #geom_text(label = "C", x = "Productive Epilimnion Particle", y = max(ddply_prodbeta$mean)) +
  geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=.2, position=position_dodge(.9)) +
  scale_color_manual(name = "", limits=c("Productive Epilimnion Particle", "Productive Epilimnion Free", "Productive Hypolimnion Particle", "Productive Hypolimnion Free",
                                         "Unproductive Epilimnion Particle", "Unproductive Epilimnion Free", "Unproductive Hypolimnion Particle", "Unproductive Hypolimnion Free"), 
                     values = c("deeppink", "deeppink", "deeppink", "deeppink",
                                "turquoise3","turquoise3","turquoise3","turquoise3"))+
  scale_x_discrete(breaks=c("Productive Epilimnion Particle", "Productive Epilimnion Free", "Productive Hypolimnion Particle", "Productive Hypolimnion Free",
                            "Unproductive Epilimnion Particle", "Unproductive Epilimnion Free", "Unproductive Hypolimnion Particle", "Unproductive Hypolimnion Free"),
                   labels=c("Epilimnion \nParticle-Associated", "Epilimnion \nFree-Living", "Hypolimnion \nParticle-Associated", "Hypolimnion \nFree-Living",
                            "Epilimnion \nParticle-Associated", "Epilimnion \nFree-Living", "Hypolimnion \nParticle-Associated", "Hypolimnion \nFree-Living")) + 
  xlab("Habitat") + ylab("Bray Curtis Dissimilarity") + theme_bw() +  #scale_fill_brewer(palette="Paired") + 
  facet_grid(. ~ trophicstate1, scale = "free", space = "free") +
  theme(axis.title.x = element_text(face="bold", size=16),
        axis.text.x = element_text(angle=60, colour = "black", vjust=1, hjust = 1, size=14),
        axis.text.y = element_text(colour = "black", size=14),
        axis.title.y = element_text(face="bold", size=16),
        plot.title = element_text(face="bold", size = 20),
        axis.ticks.x = element_blank(),
        strip.background = element_blank(), strip.text = element_blank(),
        plot.margin = unit(c(-2.9, 3.95, 0.2, 0.5), "cm"),  #top, right, bottom, left   6.55
        legend.position="none"); prodbeta_plot2  

#jpeg(filename="~/Final_PAFL_Trophicstate/Figures/Fig.3.jpeg", width= 20, height=25, units= "cm", pointsize= 14, res=500)
grid.newpage()
pushViewport(viewport(layout=grid.layout(3,1,height=c(0.3,0.35, 0.35))))
print(prod_invsimps, vp=viewport(layout.pos.row=1,layout.pos.col=1))
print(prod_obs, vp=viewport(layout.pos.row=2,layout.pos.col=1))
print(prodbeta_plot2, vp=viewport(layout.pos.row=3,layout.pos.col=1))
#dev.off()








#########  Looking at Beta diversity within particles and free living 
#######  Creating geom_point plot for Beta Diverstiy.
# STATS ON BETA
betaPAFL <- subset(nomix_beta2, filter1 == filter2)

for(i in 1:length(betaPAFL$limnion1)){
  betaPAFL$troph_filt[i]<-paste(as.character(betaPAFL$trophicstate1[i]),
                            as.character(betaPAFL$filter1[i]))}

ddply_betaPAFL <- ddply(betaPAFL, c("troph_filt", "trophicstate1"), summarise, 
                        N = length(value),
                        mean = mean(value),
                        sd   = sd(value),
                        se   = sd / sqrt(N))

ddply_betaPAFL$troph_filt <-factor(ddply_betaPAFL$troph_filt,levels=c("Eutrophic Particle","Eutrophic Free","Mesotrophic Particle", "Mesotrophic Free", "Oligotrophic Particle", "Oligotrophic Free"))

#jpeg(filename="~/Final_PAFL_Trophicstate/Figures/Fig.3_beta_filter.jpeg", width= 25, height=15, units= "cm", pointsize= 14, res=500)
betaPAFL_plotSE <- ggplot(ddply_betaPAFL, aes(x = troph_filt, y = mean, color = troph_filt)) + geom_point(size = 5) +
  facet_grid(. ~ trophicstate1, scale = "free") + 
  scale_color_manual(name = "", limits=c("Eutrophic Particle","Eutrophic Free","Mesotrophic Particle", "Mesotrophic Free", "Oligotrophic Particle", "Oligotrophic Free"), 
                     values = c("deeppink", "deeppink","orange","orange","turquoise3","turquoise3"))+
  scale_x_discrete(breaks=c("Eutrophic Particle","Eutrophic Free","Mesotrophic Particle", "Mesotrophic Free", "Oligotrophic Particle", "Oligotrophic Free"),
                   labels=c("Particle-Associated", "Free-Living", "Particle-Associated", "Free-Living", "Particle-Associated", "Free-Living")) +
  geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=.2, position=position_dodge(.9)) +
  xlab("Filter Fraction") + ylab("Bray Curtis Dissimilarity: SD") + theme_bw() + #scale_fill_brewer(palette="Paired") + 
  theme(axis.title.x = element_text(face="bold", size=14),
        axis.text.x = element_text(angle=30, colour = "black", vjust=1, hjust = 1, size=14),
        axis.text.y = element_text(colour = "black", size=14),
        axis.title.y = element_text(face="bold", size=16),
        plot.title = element_text(face="bold", size = 20),
        strip.text.x = element_text(size=12, face = "bold", colour = "black"),
        strip.background = element_blank(),
        legend.position="none"); betaPAFL_plotSE
#dev.off()



####################################################  PHYLUM LEVEL LOG-2-FOLD RATIO ANALYSIS  ####################################################
####################################################  PHYLUM LEVEL LOG-2-FOLD RATIO ANALYSIS  ####################################################
####################################################  PHYLUM LEVEL LOG-2-FOLD RATIO ANALYSIS  ####################################################
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
    value[value=="Prod vs. Unprod"]   <- "Productive \n vs.\nUnproductive"
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
dfrat$Habitat[dfrat$Habitat == " Top Productive"] <- "Epilimnion Productive"
dfrat$Habitat[dfrat$Habitat == " Top Unproductive"] <- "Epilimnion Unproductive"
dfrat$Habitat[dfrat$Habitat == " Bottom Productive"] <- "Hypolimnion Productive"
dfrat$Habitat[dfrat$Habitat == " Bottom Unproductive"] <- "Hypolimnion Unproductive"
dfrat$Habitat[dfrat$Habitat == " PA Productive"] <- "Productive \nParticle-Associated"
dfrat$Habitat[dfrat$Habitat == " PA Unproductive"] <- "Unproductive \nParticle-Associated"
dfrat$Habitat[dfrat$Habitat == " FL Productive"] <- "Productive \nFree-Living"
dfrat$Habitat[dfrat$Habitat == " FL Unproductive"] <- "Unproductive \nFree-Living"
dfrat$Habitat[dfrat$Habitat == " PA Top"] <- "Particle-Associated \nEpilimnion"
dfrat$Habitat[dfrat$Habitat == " FL Top"] <- "Free-Living \nEpilimnion"
dfrat$Habitat[dfrat$Habitat == " FL Bottom"] <- "Free-Living \nHypolimnion"
dfrat$Habitat[dfrat$Habitat == " PA Bottom"] <- "Particle-Associated \nHypolimnion"

dfrat$Habitat <- factor(dfrat$Habitat,levels = c("Epilimnion Productive", "Hypolimnion Productive", "Mixed", "Epilimnion Unproductive", "Hypolimnion Unproductive",
                                                 "Particle-Associated \nHypolimnion", "Free-Living \nEpilimnion", "Free-Living \nHypolimnion", 
                                                 "Productive \nParticle-Associated", "Productive \nFree-Living", "Unproductive \nParticle-Associated", "Unproductive \nFree-Living"))

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
        plot.margin = unit(c(2, 2, 1.65, -1.25), "cm"), #top, right, bottom, left
        #panel.grid.minor=element_blank(), #panel.grid.major=element_blank(),
        legend.position="none"); #relabun_plot



#multiplot(relabun_plot, heat2, cols = 2)

#http://stackoverflow.com/questions/20817094/how-to-control-width-of-multiple-plots-in-ggplot2
#jpeg(filename="~/Final_PAFL_Trophicstate/Figures/Fig.4_heat+abund.jpeg", width= 35, height=35, units= "cm", pointsize= 8, res=250)
grid.newpage()
pushViewport(viewport(layout=grid.layout(1,2,width=c(0.8,0.2))))
print(heat, vp=viewport(layout.pos.row=1,layout.pos.col=1))
print(relabun_plot, vp=viewport(layout.pos.row=1,layout.pos.col=2))
#dev.off()


####################################################  OTU LEVEL LOG-2-FOLD RATIO ANALYSIS  ####################################################
####################################################  OTU LEVEL LOG-2-FOLD RATIO ANALYSIS  ####################################################
####################################################  OTU LEVEL LOG-2-FOLD RATIO ANALYSIS  ####################################################
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
#Surface productive
pafl1 <- subset(de_topProd, select = c(Phylum, Genus, Species, log2FoldChange, padj))
pafl1$Habitat <- "PA vs. FL: Epilimnion Productive"
#Surface unproductive
pafl2 <- subset(de_botProd, select = c(Phylum, Genus, Species, log2FoldChange, padj))
pafl2$Habitat <- "PA vs. FL: Epilimnion Unproductive"
#Bottom Productive
pafl3 <- subset(de_topUNProd, select = c(Phylum, Genus, Species, log2FoldChange, padj))
pafl3$Habitat <- "PA vs. FL: Hypolimnion Productive"
#Bottom Unproductive
pafl4 <- subset(de_botUNProd, select = c(Phylum, Genus, Species, log2FoldChange, padj))
pafl4$Habitat <- "PA vs. FL: Hypolimnion Unproductive"



###Top vs Bottom:  
otu_topbot1 <- subset(de_prodPA, select = c(Phylum, Genus, Species, log2FoldChange, padj))
otu_topbot1$Habitat <- "Top vs. Bottom: Productive \n Particle-Associated"
# Particle-Associated UNproductive
otu_topbot2 <- subset(de_unprodPA, select = c(Phylum, Genus, Species, log2FoldChange, padj))
otu_topbot2$Habitat <- "Top vs. Bottom: Unproductive \n Particle-Associated"
# Free-Living Productive 
otu_topbot3 <- subset(de_prodFL, select = c(Phylum, Genus, Species, log2FoldChange, padj))
otu_topbot3$Habitat <- "Top vs. Bottom: Productive \n Free-Living"
# Free-Living UNproductive
otu_topbot4 <- subset(de_unprodFL, select = c(Phylum, Genus, Species, log2FoldChange, padj))
otu_topbot4$Habitat <- "Top vs. Bottom: Unproductive \n Free-Living"



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
    value[value=="Prod vs. Unprod"]   <- "Productive \n vs. \nUnproductive"
  }
  return(value)
}

no_unclass_genus <- subset(otu_ratios, Genus != "unclassified")

#jpeg(filename="~/Final_PAFL_Trophicstate/Figures/genus_heat.jpeg", width= 30, height=60, units= "cm", pointsize= 8, res=250)
## GENUS LEVEL HEAT PLOT
ggplot(no_unclass_genus, aes(Habitat, Genus)) + geom_tile(aes(fill = log2FoldChange)) + 
  scale_fill_gradient2(name = "Odds-Ratio", mid = "gray", low = "darkorange", high = "blue4",  na.value = "white", guide = guide_colorbar(barwidth = 3, barheight = 18)) + #scale_y_reverse() + 
  theme_bw(base_size = 12) + scale_x_discrete(expand = c(0, 0)) + scale_y_discrete(expand = c(0, 0)) + ylab(NULL) + 
  #geom_text(aes(fill = splif2$Transformed, label = splif2$Transformed, size = 8)) +
  xlab("Habitat") + ylab("Genus") + 
  facet_grid(. ~ Comparison, scales = "free", space = "free", labeller=mf_labeller) + 
  theme(axis.text.x = element_text(colour="black", size=14, angle = 30, hjust = 1, vjust = 1), 
        axis.text.y = element_text(colour="black", vjust=0.5, size=14),
        axis.title.x = element_blank(), #text(face="bold", size=16),
        legend.title = element_text(face="bold", size=12),
        legend.text = element_text(size = 12),
        legend.position = c(0.1, 0.93), #"left",
        axis.title.y = element_text(face="bold", size=16),
        plot.margin = unit(c(0.5, 1, 0.5, 0.5), "cm"),
        strip.text.x = element_text(size=16, face = "bold", colour = "black"),
        strip.background = element_blank());  
#dev.off()




########   Attempting the Vanette Plot
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
  res_prod[i,3] <- "Productive"
  ## UNPRODUCTIVE
  df_oligo <- subset(pvu_oligo, Phylum == pvu_phylumlist[i]) #free phylum i
  num_sigOTUs <- length(unique(df_oligo$Species))
  res_oligo[i,1] <- pvu_phylumlist[i] #put the name of the phylum in the table.
  res_oligo[i,2] <- num_sigOTUs * -1 #put the number of sig OTUs in the table. 
  res_oligo[i,3] <- "Unproductive"
  PVU_OTUresults <- rbind(res_prod,res_oligo)
  PVU_OTUresults$Comparison <- "Productive \n vs. \n Unproductive"
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
                                    levels = c("Particle-Associated","Free-Living","Productive","Unproductive","Hypolimnion","Epilimnion")) 

allOTU_results$Comparison <- factor(allOTU_results$Comparison,
                                    levels = c("Particle-Associated \nvs. \nFree-Living", "Productive \n vs. \n Unproductive","Hypolimnion \n vs. \n Epilimnion")) 

allOTU_results$Phylum <- factor(allOTU_results$Phylum,levels = c("Bacteroidetes", "Cyanobacteria", "Verrucomicrobia", "Betaproteobacteria", "Actinobacteria", "Planctomycetes",  "Alphaproteobacteria", "Deltaproteobacteria",
                                                                 "Gammaproteobacteria", "Chloroflexi", "Lentisphaerae", "Chlorobi", "Armatimonadetes", "Firmicutes", "Acidobacteria", "Spirochaetae", "Candidate_division_OD1",
                                                                 "NPL-UPA2", "Deinococcus-Thermus", "Candidate_division_OP3", "TM6", "Chlamydiae", "Epsilonproteobacteria", "TA18", "Fibrobacteres", "Candidate_division_SR1", 
                                                                 "Candidate_division_BRC1", "BD1-5", "Gemmatimonadetes", "Fusobacteria", "Candidate_division_WS3", "Tenericutes", "Elusimicrobia", "WCHB1-60", "Deferribacteres", 
                                                                 "Candidate_division_TM7", "Candidate_division_OP8", "SPOTSOCT00m83", "Thermotogae", "Candidate_division_OP11", "Dictyoglomi", "unclassified"))

OTUtop20 <- c("Bacteroidetes", "Cyanobacteria", "Verrucomicrobia", "Betaproteobacteria", "Actinobacteria", "Planctomycetes",  "Alphaproteobacteria", "Deltaproteobacteria",
              "Gammaproteobacteria", "Chloroflexi", "Lentisphaerae", "Chlorobi", "Armatimonadetes", "Firmicutes", "Acidobacteria", "Spirochaetae","unclassified")

summed_20 <- allOTU_results[allOTU_results$Phylum %in% OTUtop20, ]

summed_20$Phylum = with(summed_20, factor(Phylum, levels = rev(levels(Phylum))))

### This plot is flipped!
#jpeg(filename="~/Final_PAFL_Trophicstate/Figures/summed_otus_top17.jpeg",  width= 25, height=35, units= "cm", pointsize= 8, res=250)
summed <- ggplot(summed_20, aes(y=NumSigOTUs, x=Phylum, fill=Preference)) + 
  geom_bar(stat="identity", position="identity") + coord_flip() + ggtitle("Summed OTUs") +
  geom_bar(stat="identity", colour = "black", show_guide = FALSE, position="identity") +
  theme_gray() +# scale_y_continuous(breaks=seq(-200, 200, 25)) +
  facet_grid(. ~ Comparison, scales = "free_y", labeller = mf_labeller) + theme_bw() +
  ylab("Number of Significant OTUs") +  xlab("Phylum") + 
  theme(axis.title.x = element_text(face="bold", size=16),
        axis.text.x = element_text(colour = "black", size=12),
        axis.text.y = element_text(colour = "black", size=12),
        axis.title.y = element_text(face="bold", size=16),
        plot.title = element_text(face="bold", size = 20),
        legend.title = element_blank(),
        legend.text = element_text(size = 12),
        strip.text.y = element_blank(),
        strip.text.x = element_text(size = 12, face = "bold", colour = "black"),
        strip.background = element_blank(),
        legend.position="none"); summed
#dev.off()


################# AVERAGE OTU PLOT ##############
################# AVERAGE OTU PLOT ##############
################# AVERAGE OTU PLOT ##############
###### Attempting to do the AVERAGE-OTU PLOT
#we need a df with the sums of unique & significant OTUs in each phylum for each HABITAT PA AND FL
# Make a new data frame with each phylum and the number of sigs in each environment
# PA VS FL:  FL (negative) PA (positive)
group_pvfPA <- group_by(pvf_PA, Phylum, Comparison, Habitat, Preference)
sum_group_pvfPA <- summarize(group_pvfPA,count=n())
group_pvfFL <- group_by(pvf_FL, Phylum, Comparison, Habitat, Preference)
sum_group_pvfFL <- summarize(group_pvfFL,count=n())
sum_group_pvfFL$count <- sum_group_pvfFL$count * -1
#library(plyr) #plyr way
#count_PA <-count(pvf_PA, c("Phylum", "Comparison", "Habitat", "Preference"))
#count_FL <-count(pvf_FL, c("Phylum", "Comparison", "Habitat", "Preference"))
#count_FL$freq <- count_FL$freq * -1
#count_PAFL <- rbind(count_PA, count_FL)

# PROD VS UNPROD  ###  Prod (POSITIVE) vs Unprod (NEGATIVE)
group_pvu_prod <- group_by(pvu_prod, Phylum, Comparison, Habitat, Preference)
sum_group_pvu_prod <- summarize(group_pvu_prod,count=n())
group_pvu_oligo <- group_by(pvu_oligo, Phylum, Comparison, Habitat, Preference)
sum_group_pvu_oligo <- summarize(group_pvu_oligo,count=n())
sum_group_pvu_oligo$count <- sum_group_pvu_oligo$count * -1
#plyr way
#count_prod <-count(pvu_prod, c("Phylum", "Comparison", "Habitat", "Preference"))
#count_oligo <-count(pvu_oligo, c("Phylum", "Comparison", "Habitat", "Preference"))
#count_oligo$freq <- count_oligo$freq * -1

### Top (negative) vs Bottom (positive)
group_tvb_bottom <- group_by(tvb_bottom, Phylum, Comparison, Habitat, Preference)
sum_group_tvb_bottom <- summarize(group_tvb_bottom,count=n())
group_tvb_top <- group_by(tvb_top, Phylum, Comparison, Habitat, Preference)
sum_group_tvb_top <- summarize(group_tvb_top,count=n())
sum_group_tvb_top$count <- sum_group_tvb_top$count * -1

#counts_ALL <- rbind(count_PA, count_FL, count_prod, count_oligo, count_top, count_bottom)
counts_ALL <- rbind(sum_group_pvfPA, sum_group_pvfFL, sum_group_pvu_prod, sum_group_pvu_oligo, sum_group_tvb_top, sum_group_tvb_bottom)
counts_ALL$Phylum <- factor(counts_ALL$Phylum,levels = c("Bacteroidetes", "Cyanobacteria", "Verrucomicrobia", "Betaproteobacteria", "Actinobacteria", "Planctomycetes",  "Alphaproteobacteria", "Deltaproteobacteria",
                                                         "Gammaproteobacteria", "Chloroflexi", "Lentisphaerae", "Chlorobi", "Armatimonadetes", "Firmicutes", "Acidobacteria", "Spirochaetae", "Candidate_division_OD1",
                                                         "NPL-UPA2", "Deinococcus-Thermus", "Candidate_division_OP3", "TM6", "Chlamydiae", "Epsilonproteobacteria", "TA18", "Fibrobacteres", "Candidate_division_SR1", 
                                                         "Candidate_division_BRC1", "BD1-5", "Gemmatimonadetes", "Fusobacteria", "Candidate_division_WS3", "Tenericutes", "Elusimicrobia", "WCHB1-60", "Deferribacteres", 
                                                         "Candidate_division_TM7", "Candidate_division_OP8", "SPOTSOCT00m83", "Thermotogae", "Candidate_division_OP11", "Dictyoglomi", "unclassified"))

counts_ALL$Phylum = with(counts_ALL, factor(Phylum, levels = rev(levels(Phylum))))



### SUMMARIZING WITH DDPLY 
#df_top16_positive <- df_top16
#df_top16_positive$freq <- abs(df_top16_positive$freq)
df_all <- ddply(counts_ALL, c("Phylum", "Preference"), summarise, 
                N = length(count),
                mean = mean(count),
                sd   = sd(count),
                se   = sd / sqrt(N))


# Add a Comparison Column 
df_all$Comparison <- ""  # Add an empty column named Comparison 
ddply_table <- data.table(df_all)
##  Add in the values
ddply_table[, Comparison := ifelse(Preference %in% c("Particle-Associated", "Free-Living"), "Particle-Associated \nvs. \nFree-Living",
                                   ifelse(Preference %in% c("Epilimnion", "Hypolimnion"), "Hypolimnion \n vs. \n Epilimnion",
                                          ifelse(Preference %in% c("Unproductive", "Productive"), "Productive \n vs. \n Unproductive", NA)))]
###  Order the values
ddply_table$Comparison <- factor(ddply_table$Comparison, levels = c("Particle-Associated \nvs. \nFree-Living", "Productive \n vs. \n Unproductive", "Hypolimnion \n vs. \n Epilimnion"))
ddply_table$Preference <- factor(ddply_table$Preference,levels = c("Particle-Associated", "Free-Living", "Productive", "Unproductive", "Hypolimnion", "Epilimnion"))

### PLOTTING THE AVERAGE! 
#jpeg(filename="~/Final_PAFL_Trophicstate/Figures/average_otus_all_SE.jpeg", width= 30, height=35, units= "cm", pointsize= 8, res=250)
ggplot(ddply_table, aes(y=mean, x=Phylum, fill=Preference)) + theme_bw() +
  geom_bar(stat="identity", position="identity") + coord_flip() + 
  geom_hline(yintercept = 0) + guides(fill = guide_legend(keywidth = 2, keyheight = 2)) +
  geom_bar(stat="identity", colour = "black", show_guide = FALSE, position="identity") +
  scale_y_continuous(breaks=seq(-20, 100, 20)) +#scale_y_continuous(breaks=seq(-200, 200, 25))+
  facet_grid(. ~ Comparison, scales = "free_y", labeller = mf_labeller) + 
  geom_errorbar(aes(ymin = mean - se, ymax = mean + se), width = 0.25) +
  ylab("Mean Number of Significant OTUs") +  xlab("Phylum") + 
  theme(axis.title.x = element_text(face="bold", size=16),
        axis.text.x = element_text(colour = "black", size=12),
        axis.text.y = element_text(colour = "black", size=14),
        axis.title.y = element_text(face="bold", size=16),
        plot.title = element_text(face="bold", size = 20),
        legend.title = element_blank(),
        legend.text = element_text(size = 12),
        strip.text.y = element_blank(),
        strip.text.x = element_text(size = 12, face = "bold", colour = "black"),
        strip.background = element_blank(),
        legend.position = "right")#c(0.88, 0.12))
#dev.off()



# Top 20 Phyla
OTUtop20 <- c("Bacteroidetes", "Cyanobacteria", "Verrucomicrobia", "Betaproteobacteria", "Actinobacteria", "Planctomycetes",  "Alphaproteobacteria", "Deltaproteobacteria",
              "Gammaproteobacteria", "Chloroflexi", "Lentisphaerae", "Chlorobi", "Armatimonadetes", "Firmicutes", "Acidobacteria", "Spirochaetae","unclassified")
top20_ddplyOTU <- subset(ddply_table, Phylum == OTUtop20)

OTU_20 <- ddply_table[ddply_table$Phylum %in% OTUtop20, ]

#jpeg(filename="~/Final_PAFL_Trophicstate/Figures/Fig.6_average_otus_top17_SE.jpeg", width= 25, height=20, units= "cm", pointsize= 8, res=250)
ave <- ggplot(OTU_20, aes(y=mean, x=Phylum, fill=Preference)) + theme_bw() +
  geom_bar(stat="identity", position="identity") + coord_flip() + 
  geom_hline(yintercept = 0) + guides(fill = guide_legend(keywidth = 2, keyheight = 2)) +
  geom_bar(stat="identity", colour = "black", show_guide = FALSE, position="identity") +
  scale_y_continuous(breaks=seq(-20, 100, 20)) +#scale_y_continuous(breaks=seq(-200, 200, 25))+
  facet_grid(. ~ Comparison, scales = "free_y", labeller = mf_labeller) + 
  geom_errorbar(aes(ymin = mean - se, ymax = mean + se), width = 0.25) +
  ylab("Mean Number of Significant OTUs") +  xlab("Phylum") + 
  theme(axis.title.x = element_text(face="bold", size=16),
        axis.text.x = element_text(colour = "black", size=12),
        axis.text.y = element_text(colour = "black", size=14),
        axis.title.y = element_text(face="bold", size=16),
        plot.title = element_text(face="bold", size = 20),
        legend.title = element_blank(),
        legend.text = element_text(size = 12),
        strip.text.y = element_blank(),
        strip.text.x = element_text(size = 12, face = "bold", colour = "black"),
        strip.background = element_blank(),
        legend.position = "right" ); ave # c(0.18, 0.12))
#dev.off()


average <- ave + ggtitle("Averaged OTUs") + xlab(" ") ; average


#jpeg(filename="~/Final_PAFL_Trophicstate/Figures/OTUs_summed+avg_top17.jpeg", width= 45, height=25, units= "cm", pointsize= 14, res=500)
grid.newpage()
pushViewport(viewport(layout=grid.layout(1,2,width=c(0.42,0.58))))
print(summed, vp=viewport(layout.pos.row=1,layout.pos.col=1))
print(average, vp=viewport(layout.pos.row=1,layout.pos.col=2))
#dev.off()


##################################################################################### ABUNDANCE PLOTS 
##################################################################################### ABUNDANCE PLOTS 
##################################################################################### ABUNDANCE PLOTS 
### Check lines 1211 for how sub_phy_melt_totals was created:
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

#jpeg(filename="Abundance_top15.jpeg", width= 28, height=25, units= "cm", pointsize= 10, res=250)
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
        legend.position="right"); top15plot
#dev.off()


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
                             textGrob("Productive", gp = gpar(fontsize = 14, col = "black"))),
                     2, 4, 2, 13, name = paste(runif(2)))

z <- gtable_add_grob(z, list(rectGrob(gp = gpar(col = NA, linetype = 1, fill = gray(0.7))),
                             textGrob("Free-Living", gp = gpar(fontsize = 14, col = "black"))),
                     3, 15, 3, 19, name = paste(runif(2)))

z <- gtable_add_grob(z, list(rectGrob(gp = gpar(col = NA, linetype = 1, fill = gray(0.7))),
                             textGrob("Particle-Associated", gp = gpar(fontsize = 14, col = "black"))),
                     3, 22, 3, 25, name = paste(runif(2)))

z <- gtable_add_grob(z, list(rectGrob(gp = gpar(col = NA, linetype=1, fill = gray(0.6))),
                             textGrob("Unproductive", gp = gpar(fontsize = 14, col = "black"))),
                     2, 15, 2, 25, name = paste(runif(2)))

# add margins
z <- gtable_add_cols(z, unit(1/8, "line"), 7)
z <- gtable_add_rows(z, unit(1/8, "line"), 3)

#jpeg(filename="~/Final_PAFL_Trophicstate/Figures/Fig.S2_Abundance_top15_gtable.jpeg", width= 40, height=25, units= "cm", pointsize= 10, res=250)
# draw it
grid.newpage()
grid.draw(z)
#dev.off()













###  MIDDLE 19 PHYLA 
mid20phy <- prod_quad_phylum_stats[prod_quad_phylum_stats$Phylum %in% phy_order[16:35], ]

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

#jpeg(filename="Abundance_middle19.jpeg", width= 28, height=25, units= "cm", pointsize= 10, res=250)
mid20phy_plot <- ggplot(mid20phy, aes(y=mean_abundance , x=Phylum, fill=Phylum)) + #coord_cartesian(xlim = c(0, 30)) + 
  guides(fill = guide_legend(reverse=TRUE)) + ggtitle("") + 
  geom_bar(stat="identity", position=position_dodge()) + #theme_classic() +
  geom_bar(stat="identity", position=position_dodge(), colour = "black", show_guide = FALSE) +
  facet_wrap(~ prod_quad, ncol = 8) + xlab("Middle 20 Most Abundant Phyla") +
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
        legend.position="right"); #mid20phy_plot
#dev.off()



library(gtable)
# get gtable object
z <- ggplot_gtable(ggplot_build(mid20phy_plot))

# http://stackoverflow.com/questions/22818061/annotating-facet-title-as-strip-over-facet
# http://stackoverflow.com/questions/11353287/how-do-you-add-a-general-label-to-facets-in-ggplot2
# http://stackoverflow.com/questions/11442981/ggplot2-strip-text-labels-facet-wrap
# add label for top strip
z <- gtable_add_rows(z, z$heights[[3]], 2)
z <- gtable_add_grob(z, list(rectGrob(gp = gpar(col = NA, linetype = 1, fill = gray(0.7))),
                             textGrob("Free-Living", gp = gpar(fontsize = 14, col = "black"))),
                     3, 4, 3, 7, name = paste(runif(2)))

z <- gtable_add_grob(z, list(rectGrob(gp = gpar(col = NA, linetype = 1, fill = gray(0.7))),
                             textGrob("Particle-Associated", gp = gpar(fontsize = 14, col = "black"))),
                     3, 9, 3, 13, name = paste(runif(2)))

z <- gtable_add_grob(z, list(rectGrob(gp = gpar(col = NA, linetype = 1, fill = gray(0.6))),
                             textGrob("Productive", gp = gpar(fontsize = 14, col = "black"))),
                     2, 4, 2, 13, name = paste(runif(2)))

z <- gtable_add_grob(z, list(rectGrob(gp = gpar(col = NA, linetype = 1, fill = gray(0.7))),
                             textGrob("Free-Living", gp = gpar(fontsize = 14, col = "black"))),
                     3, 15, 3, 19, name = paste(runif(2)))

z <- gtable_add_grob(z, list(rectGrob(gp = gpar(col = NA, linetype = 1, fill = gray(0.7))),
                             textGrob("Particle-Associated", gp = gpar(fontsize = 14, col = "black"))),
                     3, 22, 3, 25, name = paste(runif(2)))

z <- gtable_add_grob(z, list(rectGrob(gp = gpar(col = NA, linetype=1, fill = gray(0.6))),
                             textGrob("Unproductive", gp = gpar(fontsize = 14, col = "black"))),
                     2, 15, 2, 25, name = paste(runif(2)))

# add margins
z <- gtable_add_cols(z, unit(1/8, "line"), 7)
z <- gtable_add_rows(z, unit(1/8, "line"), 3)

#jpeg(filename="~/Final_PAFL_Trophicstate/Figures/Fig.S3_Abundance_mid20_gtable.jpeg", width= 40, height=25, units= "cm", pointsize= 10, res=250)
# draw it
grid.newpage()
grid.draw(z)
#dev.off()



####################################################################################  OTU SUM
####################################################################################  OTU SUM
####################################################################################  OTU SUM
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
PAFL_otu$SampleType <- c("Particle", "Particle", "Particle", "Free", "Free", "Free")
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
prod_otu$SampleType <- c("High Nutrient", "High Nutrient", "High Nutrient", "Low Nutrient", "Low Nutrient", "Low Nutrient")
prod_otu$Type <- c("Total", "High Nutrient Only", "Both", "Total", "Low Nutrient Only", "Both")
prod_otu$NumOTUs <- c(prod_total, prod_only, Prod_overlap, unprod_total, unprod_only, Prod_overlap)
prod_otu2 <- subset(prod_otu, Type != "Total")
#Fisher's test:  p = 0.0001


###  PUT IT ALL TOGETHER! 
PAFL_otu2$Comparison <- "Particle-Associated \n vs. Free-Living"
TB_otu2$Comparison <- "Epilimnion \n vs. Hypolimnion"
prod_otu2$Comparison <- "High  vs. Low \nNutrient"
otusums <- rbind(PAFL_otu2, TB_otu2, prod_otu2)

otusums$Comparison <- factor(otusums$Comparison,levels = c("Particle-Associated \n vs. Free-Living", "High  vs. Low \nNutrient", "Epilimnion \n vs. Hypolimnion"))
otusums$SampleType <- factor(otusums$SampleType,levels = c("Free", "Particle", "Low Nutrient", "High Nutrient", "Epilimnion", "Hypolimnion"))
otusums$Type <- factor(otusums$Type, levels = c("Both","Free-Living Only", "Particle-Associated Only", "Low Nutrient Only", "High Nutrient Only", "Epilimnion Only", "Hypolimnion Only"))
#otusums$Type <- factor(otusums$Type, levels = rev(levels(otusums$Type)))


label.df <- data.frame(Group = c("Productive \n vs. Unproductive", "Epilimnion \n vs. Hypolimnion"),
                       Value = c(8600, 8000))

#jpeg(filename="~/Final_PAFL_Trophicstate/Figures/Fig.5_DetectedOTUs_rarefied.jpeg", width= 23, height=18, units= "cm", pointsize= 10, res=250)
ggplot(otusums, aes(y=NumOTUs , x=SampleType, fill=Type, order=Type)) +
  facet_grid(. ~ Comparison,scales = "free") + #geom_text(x = 2, y = 8750, label = "***") +
  xlab("Sample Type ") + ylab("Number of Deteceted UniqueOTUs") + 
  geom_bar(stat="identity") +   geom_bar(stat="identity", colour="black", show_guide=FALSE) +  
  scale_y_continuous(expand = c(0,0), breaks=seq(0, 6000, 1000), lim = c(0, 6000)) + 
  scale_fill_manual(limits = c("Free-Living Only", "Particle-Associated Only", "Low Nutrient Only", "High Nutrient Only", "Epilimnion Only", "Hypolimnion Only", "Both"),
                    breaks = c("Free-Living Only", "Particle-Associated Only", "Low Nutrient Only", "High Nutrient Only", "Epilimnion Only", "Hypolimnion Only", "Both"),
                    values = c("orange", "red", "darkgreen", "limegreen", "deepskyblue", "blue4", "gray39")) +
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

