# SNMF on LEA to get structure plot


## Convert VCF to Ped file
run interactive mode

            srun -p lotterhos -N 1 --pty /bin/bash

1. Create chromosome map

```bash
module load miniconda3
conda activate plink
conda activate bcftools

bcftools view -H populations.snps.vcf | cut -f 1 | uniq | awk '{print $0"\t"$0}' > populations.snps.vcf.chrom-map.txt
```

2. Make ped file using this chromosome map

```bash
module load vcftools

vcftools --vcf populations.snps.vcf --out BSB_r0.8_maf0.01 --plink --chrom-map populations.snps.vcf.chrom-map.txt
```

```bash
VCFtools - 0.1.17
(C) Adam Auton and Anthony Marcketta 2009

Parameters as interpreted:
	--vcf populations.snps.vcf
	--chrom-map populations.snps.vcf.chrom-map.txt
	--out BSB_r0.8_maf0.01
	--plink

After filtering, kept 117 out of 117 Individuals
Writing PLINK PED and MAP files ...
	Read 39 chromosome mapping file entries.
Done.
After filtering, kept 24896 out of a possible 24896 Sites
Run Time = 1.00 seconds
```
This should have created a .ped and .map file. Now download this to your local computer or wherever you will run R from

## Run SNMF on LEA (input file = .ped)

In R, install LEA and convert from a .ped file to a geno file to run SNMF (We will first need to convert to a lfmm file and then to a geno file)

http://membres-timc.imag.fr/Olivier.Francois/LEA/tutorial.html


```r
############
## SNMF ####
############
rm(list = ls())

#===== BSB =========
# 6 source populations
# 117 individuals
# 24,896 diploid loci
# ? env variables

#Set your working directory where you have your files
setwd("/Users/remygatins/GoogleDrive_gmail/Work/Projects/2021_Black\ Sea\ Bass/RADs/popmap_ref/BSB_all_r0.8_maf0.01_wsnp/")

###################
##### LEA #########
###################

#Install LEA
#if (!require("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")

#BiocManager::install("LEA")
library(LEA)


#Convert ped to lfmm and lfmm to geno
output = ped2lfmm("BSB_r0.8_maf0.01.ped")
#output = BSB_r0.8_maf0.01.lfmm (writen to the working directory)
# Now create a geno file:	"genotypes.geno".
output = lfmm2geno("BSB_r0.8_maf0.01.lfmm", "genotypes.geno")

library(vcfR)
vcf <- read.vcfR(file = "populations.snps.vcf")
#vcf <- read.vcfR(file = "filtered_ind_2.recode.vcf")
#vcf <- read.vcfR(file = "filtered.recode.vcf")
data <- vcfR2genlight(vcf)
data

pop_ID <-  c(rep("MA",20), rep("MD",20), rep("ME",18), rep("NC",13), rep("NJ",17), rep("RI",29))
ind_ID <- data$ind.names
library(stringr)
ind_ID <- str_sub(ind_ID,1,6) #keep only first 6 characters of the ind_ID

#Run SNMF
obj.snmf = snmf("genotypes.geno", K = 2, alpha = 100, project = "new") 
qmatrix = Q(obj.snmf, K = 2)
barplot(t(qmatrix), col = c("orange","violet","lightgreen"), border = NA, space = 0, xlab = "Individuals", ylab = "Admixture coefficients")

library("RColorBrewer")
cols <- brewer.pal(n = 6, name = "Paired")
barplot(t(qmatrix), col = cols, border = NA, space = 0, xlab = "Individuals", ylab = "Admixture coefficients")

# Choose the number of clusters
obj.snmf = snmf("genotypes.geno", K = 1:6, ploidy = 2, entropy = T, alpha = 100, project = "new")
plot(obj.snmf, col = "blue4", cex = 1.4, pch = 19)


#transform data to plot in ggplot
results <- as.data.frame(qmatrix)
results$pop <- pop_ID
results$indNames <- ind_ID

library(reshape2)
results_melt <- melt(results)
colnames(results_melt) <- c("Original_Pop","Sample","Ancestry_Pop","Ancestry")
#arrange population from north to south
results_melt$Original_Pop <- factor(results_melt$Original_Pop, levels = c("ME","MA","RI","NJ","MD","NC"))

library(ggplot2)
ggplot(results_melt, aes(x=Sample, y=Ancestry, fill=Ancestry_Pop)) +
  geom_bar(stat='identity') +
  scale_fill_manual(values = cols) +
  facet_grid(~Original_Pop, scales = "free", space = "free") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 3))+
  theme(axis.text.x = element_text(size = 4))+
  theme(axis.title.x = element_blank()) +
  theme(legend.position = "none")  

```




