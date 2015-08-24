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
  scale_y_continuous(breaks=seq(0.01, 0.09, 0.02), lim = c(0.01,0.093)) +
  theme(axis.title.x = element_text(face="bold", size=14),
        axis.text.x = element_text(angle=30, colour = "black", vjust=1, hjust = 1, size=12),
        axis.text.y = element_text(colour = "black", size=12),
        axis.title.y = element_text(face="bold", size=14),
        plot.title = element_blank(),
        strip.text.x = element_blank(),
        strip.text.y = element_blank(),
        strip.background = element_blank(),
        legend.position="none",
        plot.margin = unit(c(-0.8, 0.1, 0.2, 0.23), "cm")); plot_even_sigs





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
  scale_y_continuous(breaks=seq(400, 1200, 200), lim = c(300,1300)) +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        axis.text.y = element_text(colour = "black", size=12),
        axis.title.y = element_text(face="bold", size=14),
        plot.title = element_blank(),
        strip.text.x = element_text(size=13, face="bold"),
        strip.text.y = element_blank(),
        strip.background = element_blank(),
        legend.position="none",
        plot.margin = unit(c(0.1, 0.1, 0.1, 0.13), "cm")); plot_richobs_sigs  #top, right, bottom, left)



#####  Plotting FIGURE 1  #####  Plotting FIGURE 1  #####  Plotting FIGURE 1  #####  Plotting FIGURE 1  #####  Plotting FIGURE 1
#tiff(filename="~/Final_PAFL_Trophicstate/Final_Figures/Fig.1_alpha_SIGS.tiff", width= 30, height=22, units= "cm", pointsize= 14, res=200)
pdf(file="~/Final_PAFL_Trophicstate/Final_Figures/Fig.1_alpha_SIGS_bigger2.pdf", width= 9, height=7)
grid.newpage()
pushViewport(viewport(layout=grid.layout(2,1,height=c(0.45,0.55))))
print(plot_richobs_sigs, vp=viewport(layout.pos.row=1,layout.pos.col=1))
print(plot_even_sigs, vp=viewport(layout.pos.row=2,layout.pos.col=1))
dev.off()

