## Phylogenetic conservation of freshwater lake habitat preferences varies between abundant bacterioplankton phyla.

###  **Authors:** Marian L. Schmidt, Jeffrey D. White, and Vincent J. Denef  



#### Accepted by _Environmental Microbiology_ on November 19th, 2015  

**********


####Information about this repository:  

##### **Original fastq files:**
The original Fastq files were submitted to the NCBI sequence read archive under BioProject PRJNA304344, SRA accession number SRP066777.


##### **Code Files:**
1. `SoMiLakes_PAFL_Analysis.Rmd` includes all of the data analysis!  
      - `SoMiLakes_PAFL_Analysis.html` is the rendered html of the above file.
2. `Functions_PAFL.R` has all of the functions written or sourced for the data analysis in the `SoMiLakes_PAFL_Analysis.Rmd` code file.  


##### **Directories/Directory Set up:**
1. **Nov6_FWDB_Silva:**  Includes the raw output files from the mothur pipline and the raw environmental data collected for the project.  Including: 
    -  **Mothur cons.taxonomy file:**  *SoMiLakes_copy_FWDB.files.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.pick.an.unique_list.0.03.cons.taxonomy*  
    -  **Mothur shared file:**  *SoMiLakes_copy_FWDB.files.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.pick.an.unique_list.shared*  
    -  **mothur batch file:**  *Nov4_FWDB_batch*  
2. **Nov6_FWDB_Silva/Nov_alpha_data:**  Includes the 100X sub-sampled (rarefied) data for within-sample (alpha) diversity calculations.  Also includes otu_melt, which is a large file containing the otu abundance counts for each sample in long-format.  
    -  **100X Sub-sampled Inverse Simpson:** *InvSimpson100_summed*  
    -  **100X Sub-sampled Observed Richness:** *ObservedRichness100_summed*  
    -  **Each OTU abundance per sample:** *otu_melt*
3. **FWDB_Figures:**  A home for all the figures created in the `SoMiLakes_PAFL_Analysis.Rmd` file.  
    - This is where all the figures will be outputted.



**Note:**  This project is under the general MIT License.
