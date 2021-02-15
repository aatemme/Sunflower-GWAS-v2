# Sunflower-GWAS-2.0

This is a R only rewrite and extention of the sunflower GWAS pipeline initiated by Rishi Masalia. (Masalia et al., Plos 2018) Please cite https://doi.org/10.1104/pp.20.00873 when using this code :)

NOTE: the GEMMA gwas algorithm used in this pipeline runs on Linux

This Sunflower-GWAS-2.0 pipeline adds:
- Streamlined folder structure of inputs/outputs/code/software
- GWAS calcuation using GEMMA
- Colocalization visualization using haplotype blocks
- Gene list per haplotype block
- Drawing of manhattan plots with haplotype blocks overlay

![picture](Overview.jpg)




**General idea for how using the haplotype blocks informs our GWAS procedure**
- The blocks are generated from the 1.5 million SNPs using plink. This divides up the genome in chunks that are co-inherited (according to some set thresholds). We now have condensed the genome from 1.5 million independent SNPs to ~20K independent chunks/regions/blocks. Now from the GEMMA output we get a list of SNPs and their p-values. Since not all SNPs are independent (as evidenced from the blocking procedure) the cutoff for significance is based on the 20K independent blocks we have.
- Another place where the block map we've constructed really shines is finding traits that colocate to the same region. If two or more traits have SNPs that lie in the same block we can say they colocate (possibly pleiotropic but most likely close linkage). So in determining if traits have significant SNPs in the same blocks (even if the SNPs are not the same) all we have to do is figure out if they have common blocks.
- A problem with the blocking procedure is that it's sensitive to missplaced SNPs. The algorithm grows a region by walking along the chromosome and stops if a certain fraction of SNPs don't fit the region it is growing (based on D'). What could happen is that two blocks that are adjacent in reality are broken up by a miss placed cluster of SNPs. If there are traits that hit to the first part of the "true" block and traits that hit to the second part of the "true" block we would be wrong in saying these traits are independent.
- To partially protect against this we can look at the SNPs that are significant for a trait at least once, line up those SNPs and re-do the blocking. If the SNPs from two or more blocks now fully go together in a single block we know that the intervening snps could have been missplaced and for the purposes of determining colocalization we should lump the blocks together.
In the triangle plots I'm showing the LD of the SNPs. The blocks the SNPs belong to in the "genome" row, and the new, significant only SNPs, blocks in the "significant" row. 

**Pipeline step by step**
* **Step zero: Preparation**
  * Make sure you have the required software and data files in the software folder.
  * Make sure you have the phenotype data in the correct format (following the example file) in the phenotype data folder.
  * Set the preferences in the scripts folder
  * Set the traits and environments you want to analyze in the top level “traits to run” and “environments to run” files.
* **Step one: GWAS**
  * Run script 1.
    * This scripts loops over each trait and environment combination. It takes the relevant phenotype data from the datafile and merges it with the .fam file in the software folder. Then it runs a call to the gemma program to do the actual association testing. Output files are generated each iteration of the loop and saved in “Tables/Assoc_files”
* **Step two: Manhattan plots**
   * Run script 2 (and/or 2b)
      * These scripts take the association testing output files generated in step 1 and draws simple Manhattan plots per trait. Output files are saved in “Plots/Manhattans”
* **Step three: Refine haplotype blocks based on significant SNPs**
  * Run script 3
    * This script loops over each association test file to identify the SNPs on the genome that are significant at least once. Then it reruns the haplotype blocking procedure to check if the haplotype blocks containing the significant SNPs can be combined. Subsequently it calculates D’ (D prime) between SNPs per chromosome (with a call to PLINK) and draws LD heatmaps per chromosome. The haplotype blocks as based on the full SNPset and the reblocked (only significant SNPs) set are overlaid on the linkage heatmaps. Output LD heatmaps are drawn per chromosome in “Plots/Colocalization”.
    * Additionally, this script (in the call to script 3b) generates a list of traits+environments and the haplotype blocks for which there are significant and/or suggestive SNPs. Output files are saved in “Tables/Blocks”
* **Step four: Draw colocalization of traits to similar genomic regions**
    * Run script 4 (and/or 4c)
      * This script takes the output from step three and visualizes colocalization of traits to genomic regions. Traits are ordered by hierarchical clustering and a dendrogram is drawn in step 4b. A plot per environment is drawn in “Plots/Colocalization”
      * Script 4c draws the same data visualization per chromosome in “Plots/Colocalization”
* **Step five: Determine genes in genome blocks**
  * Run script 5
    * This script finds the edges of the genomic regions that have significant SNPs and lists the genes contained therein. An output table is generated in “Tables/Genes”. A figure with the number of genes per region is drawn in “Plots/Colocalization”.
* **Step 6: Draw fancier Manhattan plots**
  * Run script 6
    * This script draws prettier Manhattan plots per trait+environment by overlaying the significant and suggestive blocks on the plot. Output figures are drawn in “Plots/Manhattan_regionhighlight”
* **Step 7: Output Relative Effect Size (RES) to table**
  * Run script 7
    * This script summarizes the RES per significant region per trait+environment and saves it as a table. Output table is saved in “Tables/blocks”
* **Step 8: Show significant blocks on haplotype map**
   * Run script 8
      * This script draws the haplotype block map of the genome as based on the SNP set. Blocks with significant SNPs are subtly noted by large black dots along the bottom. The output figure is saved in “Plots/Colocalization”

