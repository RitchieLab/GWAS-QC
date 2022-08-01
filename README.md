
<h1 align="center">
  <br>
  <a href="https://github.com/RitchieLab/GWAS-QC"><img src="https://user-images.githubusercontent.com/30478823/182218051-a9a58111-9b92-42e9-bdc1-c2ffad393fe2.png" alt="GWAS QC" width="200"></a>
  <br>
  GWAS Quality Control (QC)
  <br>
</h1>

<h4 align="center">A complementary tutorial to the <a href="https://pubmed.ncbi.nlm.nih.gov/21234875/" target="_blank">Truong, Woerner, Cherlin, et al. Paper (2022)</a>.</h4>

<p align="center">
  <a href="#basic-overview">Basic Overview</a> •
  <a href="#set-up">Set Up</a> •
  <a href="#workflow-steps">Workflow Steps</a> •
  <a href="#related-resources">Related Resources</a> •
  <a href="#license">License</a>
</p>




## Basic Overview
* We chose a publicly available dataset from the International Genome Sample Resource (IGSR) (www.internationalgenome.org). IGSR created and maintains the 1000 Genomes (1KG) Project to provide a public catalog of common human genetic variation and genotype data. The 1KG dataset has been kept up to date with current reference data sets, thus it is available for both GRCh37 and GRCh38. The latter is utilized here because the 2014 update increased the quantity of loci represented, resolved more than 1000 issues from the previous version of the assembly, and overall provides a better basis for alignment and subsequent analysis. Additionally, IGSR’s continued efforts will lead to the incorporation of various populations to the data which were not previously captured.

## Set Up
### Modules in BASH
* plink/1.9
* plink/2.0
* bcftools/1.9
* vcftools/0.1.12c
* tabix/0.2.6

### Notes on PLINK v1.9 and v2.0
* Not all commands are portable between PLINK version 1.9 and version 2.0. Since PLINK v2.0 is under heavy active development, the developers urge users to check certain results against an earlier, more widely-used version of PLINK. Some functions are available in v1.9 which are not in v2.0, and vice versa.
* Original version of PLINK: 1.07, https://zzz.bwh.harvard.edu/plink/plink2.shtml 
* Beta version: 1.90, https://www.cog-genomics.org/plink/1.9/
* Alpha version: 2.00, https://www.cog-genomics.org/plink/2.0/ 
* [Developer's comments](https://www.biostars.org/p/299855/) 1.9 and 2.0 serving as complementary resources
> The main difference is that plink 1.9 is essentially finished, while plink 2.0 is an alpha-stage program which will have significant unfinished components for a while to come. As a consequence, current development priorities for plink 2.0 are centered around things which are impossible to do with plink 1.9, such as handling multiallelic/phased variants and dosage data and reliably tracking REF/ALT alleles; while things that plink 1.9 already handles perfectly well, such as working with .ped/.map file pairs, have been deliberately deprioritized for now. So, you should stick to 1.9 as long as it's good enough for the jobs you need to perform. But once you need to do something outside 1.9's scope, you're likely to find that 2.0 already has the additional feature you need (or it'll be added quickly after you ask for it)


## Workflow Steps
* Need to add in Tess's updated workflow steps, https://ritchielab.org/blogs/tcherlin/2022/07/29/1-kg-gwas-tutorial-steps/ 
* Van deleted her pre-imputation data simulation steps since they're outdated now with the new dataset we used

## Related Resources
We recommend the following resources and tutorials developed for performing GWAS:
* Comphrehensive tutorial about GWAS and PRS by MareesAT: https://github.com/MareesAT/GWA_tutorial/
* GWAS Data Cleaning tutorial by the GENEVA Coordinating Center: https://www.bioconductor.org/packages/devel/bioc/vignettes/GWASTools/inst/doc/DataCleaning.pdf
* GWAS QC - theory and steps by the Pan African Bioinformatics Network for H3Africa: https://www.bioinf.wits.ac.za/courses/AIMS/QC_data.pdf 



## License
* To be decided. Suggestions welcome.







# STUFF BELOW THIS SECTION NEEDS TO BE EDITED & REORGANIZED INTO THE ABOVE SECTION HEADERS

## III. Imputation using TOPMed Imputation Server

* Imputation has become an essential component of GWAS quality control because it increases power, facilitates meta-analysis, and aids interpretation of signals. Genotype imputation is the statistical inference of unobserved genotype, which enables scientists to reconstruct the missing data in each genome and accurately evaluate the evidence for association at genetic markers that were not genotyped. Genotype imputation is achieved by comparing short stretches of an individual genome against stretches of previously characterized reference genomes.  It is usually performed on single nucleotide polymorphisms (SNPs), which are the most common type of genetic variation. 
* Several tools exist specifically for genotype imputation such as the Michigan and Trans-Omics for Precision Medicine (TOPMed) Imputation Servers where one uploads the phased or unphased GWAS genotypes in order to receive the imputed genomes in return. Each imputation server varies in terms of speed and accuracy. One of the most important considerations in imputation is the composition of the reference panel. For our study, we selected the TOPMed Imputation Reference panel  (version r2) because it is one of the most diverse reference panels available and contains information from 97,256 deeply sequenced human genomes containing 308,107085 genetic variants distributed across the 22 autosomes and the X chromosome. 
Theoretically, phased means that the two strands on each Chr are separated to identify which regions come from each parent whereas no phasing means that they are not separated. Essentially, for imputation phasing is the first step which is done in reference to the reference genome panel.

![image](https://user-images.githubusercontent.com/30478823/154597606-bc2f8b09-2741-493e-9c4a-dabdf238bd23.png)
![image](https://user-images.githubusercontent.com/30478823/154598027-ea78546d-e645-460a-bab7-b41fead62356.png)
![image](https://user-images.githubusercontent.com/30478823/154598503-691c4973-2858-4ff0-a9ac-95360317b405.png)
![image](https://user-images.githubusercontent.com/30478823/154600698-1d443de3-6691-4af9-a078-6bdb8f113e5a.png)
![image](https://user-images.githubusercontent.com/30478823/154738461-b951ab13-75b3-417f-bfc1-1e461dc4cf47.png)

## 19. Location of Imputed Data
```
# Download the completed imputation files using the wget commands provided by TOPMed to the location where you'll be working with it

# Location
~group/scratch/van/cphg-gwas-qc-imputed-data

# Unzip the files and enter the password that was emailed to you from TOPMed inside the quotes
module load p7zip
for file in *.zip; do 7z e $file -p"<password>"; done
```
![image](https://user-images.githubusercontent.com/30478823/154745163-97f3cb23-03db-487c-9638-63830eec92cc.png)


## IV. Post-Imputation QC

Much of the QC can be done in PLINK. For ease start by converting the output from the imputation from `vcf.gz` to bed/bim/fam file format.

```
#!/bin/bash

SOURCE1="postimp/chr"
SOURCE2=".dose.vcf.gz"
DEST1="postimp/chr"

for i in {2..22}
do
	SOURCE=${SOURCE1}$i${SOURCE2}
	DEST=${DEST1}$i
	plink --vcf $SOURCE --make-bed --out $DEST
done
```

Since this is a small dataset, it is possible to merge chromosome files into one bed/bim/fam file, so as to run fewer overall QC commands. This uses `mergelist.txt` which is just a file with all the chromosomes listed from 2-22 (withouut the filename extension), each on a separate line.

```
plink --bfile postimp/chr1 --merge-list postimp/mergelist.txt --make-bed --out postimp/merged
```

To actually run QC, we will run three PLINK commands on all the imputed data. The first includes only markers that have over 99% complete observations. The second command includes only individuals that have over 99% complete observations. With imputed data, these steps should not remove any variants or individuals. The third command only include variants where the minor allele frequency is over 5%, thus removing rare variants. The fourth command computes the Hardy-Weinberg equilibrium for alleles, which is important to check when assessing putatively influential variants/loci.

```
plink --bfile postimp/merged --geno 0.01 --make-bed --out postimp/merged1_geno
plink --bfile postimp/merged1_geno --mind 0.01 --make-bed--out postimp/merged2_mind
plink --bfile postimp/merged2_mind --maf 0.05 --make-bed --out postimp/merged3_maf
plink --bfile postimp/merged3_maf --hardy --out postimp/merged4_hwe
```

Sometimes it is valuable to remove related individuals to remove confounding variation such as relatedness. To calculate relatedness, run the following command to get pairs of individuals with kinship coefficients greater than 0.125.

```
plink --bfile postimp/merged3_maf --genome --min .125 --out postimp/merged5_related
```

This output can be read into R, where you can arbitrarily remove samples from the first column. Since this dataset is made up of related individuals, many were removed during this processing step.

```
related <- read.table("postimp/merged5_related.genome", header = TRUE)
write.table(unique(related[,1:2]), "postimp/related_IDs", append = FALSE, quote = FALSE, row.names = FALSE, col.names = FALSE)
```

PLINK can be run to actually remove these individuals.

```
plink --bfile postimp/merged3_maf --remove postimp/related_IDs --make-bed --out postimp/merged6_related
```

## Performing GWAS

## PLINK

With the phenotype/covariate file in the right format, here are the commands to perform the GWAS. The furst includes just sex as a covaraite, while the second command include sex and the firts 6 PCs as covariates.


```
plink --bfile postimp/merged6_related --logistic --pheno data/cov_pheno_rand.txt --pheno-name Pheno --covar data/cov_pheno_rand.txt --covar-name SEX --allow-no-sex --out postimp/merged9_rand_noPCA --hide-covar
```

```
plink --bfile postimp/merged6_related --logistic --pheno data/cov_pheno_rand.txt --pheno-name Pheno --covar data/cov_pheno_rand.txt --covar-name SEX-PC6 --allow-no-sex --out postimp/merged10_rand_wPCA --hide-covar
```

The GWAS results can be visualized with a QQ-plot and a manhattan plot. These can be generated in R.

```
library(ggplot2)
library(plyr)
library(gridExtra)
library(qqman)
library(RColorBrewer)
ppi <-300

d <- read.table("postimp/merged10_rand_wPCA.assoc.logistic", header=T) # takes a few minutes to read in
d$logP = -log10(d$P)
data_sorted <- d[order(d$CHR, d$BP), ]
data_sorted$order <- seq(1:length(data_sorted$CHR))

group.colors <- c("1"="black","2"="red","3"="blue","4"="green","5"="gray","6"="purple","7"="yellow","8"="pink","9"="black","10"="red","11"="blue","12"="green","13"="gray","14"="purple","15"="yellow","16"="pink","17"="black","18"="red","19"="blue","20"="green","21"="gray","22"="purple")

p <- ggplot(data_sorted, aes(order, logP, colour=factor(CHR),)) + geom_point(size=1) +scale_colour_manual(values=group.colors,name="Chromosome") +xlab("CHR") + ylab("-log10P") +theme(text=element_text(size=20))+
 theme(axis.ticks = element_blank(), axis.text.x = element_blank()) + ggtitle("CHPG: adjusted with PC1-6 and SEX")
ggsave(p, filename="postimp/figures/cphg_plink_PC1-6_SEX_manhattan.png",dpi=300, height=7, width=20, units="in", type="cairo-png")
png("postimp/figures/cphg_plink_PC1-6_SEX_qq.png", height=7, width=7, units="in", res=300, type="cairo-png")
qq(d$P, main="CHPG: adjusted with PC1-6 and SEX")
dev.off()
# Might take long (~10 minutes) to generate figures
```



## Regenie (run by Yuki)
Ok so here is the steps I took to create bgen sample file required to run Regenie:
*info.gz files from imputed data contains all the snp information, so I just want to extract 1 snp name out of it:
```
zcat ~/group/scratch/van/cphg-gwas-qc-data/imputed-data/chr22.info.gz | head -4
```
```
SNP	REF(0)	ALT(1)	ALT_Frq	MAF	AvgCall	Rsq	Genotyped	LooRsq	EmpR	EmpRsq	Dose0	Dose1
chr22:10581332:TATG:T	TATG	T	0.00013	0.00013	0.99988	0.42729	Imputed	-	-	-	-	-
chr22:10625878:T:C	T	C	0.00014	0.00014	0.99987	0.40740	Imputed	-	-	-	-	-
chr22:10670252:G:T	G	T	0.00010	0.00010	0.99990	0.51619	Imputed	-	-	-	-	-
```

So I am using the very first snp name (irrelevant which one you pick) to extract just 1 snp out of huge vcf file and create binary file with just 1 snp, so I can quickly check the *fam file.

```
plink2 --vcf chr22.dose.vcf.gz --id-delim _ --snp chr22:10581332:TATG:T --make-bed --out test22
```
This step create binary test22.bed/bim/fam PLINK files.
*fam file contains the sample information in this format:
```
HG00096 HG00096 0	0	0	-9
HG00097 HG00097 0	0	0	-9
HG00099 HG00099 0	0	0	-9
```
I need it to change this to this format to create bgen sample file:
```
ID_1 ID_2 missing sex
0 0 0 D
HG00096 HG00096 0 NA
HG00097 HG00097 0 NA
```
```
awk '{print $1,$2,"0","NA"}' test22.fam > cphg_sample.txt
```
And I manually added in 2 line header.

These specific steps are not all that important, but you need to get used to the idea of how to find the info you need from where.

```
<code>
```

## SAIGE

Probably won't show here...



