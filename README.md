
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


# Basic Overview
* We chose a publicly available dataset from the International Genome Sample Resource (IGSR) (www.internationalgenome.org). IGSR created and maintains the 1000 Genomes Project (1kGP) to provide a public catalog of common human genetic variation and genotype data. The 1kGP dataset has been kept up to date with current reference data sets, thus it is available for both GRCh37 and GRCh38. The latter is utilized here because the 2014 update increased the quantity of loci represented, resolved more than 1000 issues from the previous version of the assembly, and overall provides a better basis for alignment and subsequent analysis. Additionally, IGSR’s continued efforts will lead to the incorporation of various populations to the data which were not previously captured.

# Set Up
## Modules in BASH
<details> 
	<summary>👇 How to load modules </summary>
	<hr>
	
```
module load plink/1.9
module load plink/2.0
module load bcftools/1.9
module load vcftools/0.1.12c
module load tabix/0.2.6
module load liftOver/20180423
module load R
```
</details>


## Set up working directory
<details> 
	<summary>👇 How to set up your directory structure for the GWAS QC workflow </summary>
	<hr>

```
mkdir GWAS_QC
mkdir GWAS_QC/preImputation
mkdir GWAS_QC/preImutation/VCFfiles
mkdir GWAS_QC/Imputed
mkdir GWAS_QC/postImpuatation
```
</details>


## Notes on PLINK v1.9 and v2.0
* Not all commands are portable between PLINK version 1.9 and version 2.0. Since PLINK v2.0 is under heavy active development, the developers urge users to check certain results against an earlier, more widely-used version of PLINK. Some functions are available in v1.9 which are not in v2.0, and vice versa.
* Original version of PLINK: 1.07, https://zzz.bwh.harvard.edu/plink/plink2.shtml 
* Beta version: 1.90, https://www.cog-genomics.org/plink/1.9/
* Alpha version: 2.00, https://www.cog-genomics.org/plink/2.0/ 

<details> 
	<summary> 👇 PLINK developers' comments on 1.9 and 2.0 serving as complementary resources </summary>
	<hr>

> The main difference is that plink 1.9 is essentially finished, while plink 2.0 is an alpha-stage program which will have significant unfinished components for a while to come. As a consequence, current development priorities for plink 2.0 are centered around things which are impossible to do with plink 1.9, such as handling multiallelic/phased variants and dosage data and reliably tracking REF/ALT alleles; while things that plink 1.9 already handles perfectly well, such as working with .ped/.map file pairs, have been deliberately deprioritized for now. So, you should stick to 1.9 as long as it's good enough for the jobs you need to perform. But once you need to do something outside 1.9's scope, you're likely to find that 2.0 already has the additional feature you need (or it'll be added quickly after you ask for it). https://www.biostars.org/p/299855/
</details>

# Workflow Steps

## PART 1 -- Get Data

<details> 
	<summary>👇 Steps and code </summary>
	<hr>
	
* Enter your GWAS_QC directory 
```
cd GWAS_QC
```
* Download -- First, we need to download the publicly available dataset from the 1000 Genomes Project (1KGP). The data is Affy6.0 genotype data for 3,450 individuals with population-level data. 

* Files can be found here: http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/supporting/hd_genotype_chip/
* Two download options 
	- Download directly to your local computer by clicking the hyperlink for ALL.wgs.nhgri_coriell_affy_6.20140825.genotypes_has_ped.vcf.gz
	- Dowload using wget command: wget http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/supporting/hd_genotype_chip/ALL.wgs.nhgri_coriell_affy_6.20140825.genotypes_has_ped.vcf.gz

* At this point `ALL.wgs.nhgri_coriell_affy_6.20140825.genotypes_has_ped.vcf.gz` should be in your GWAS_QC directory.
```
gunzip ALL.wgs.nhgri_coriell_affy_6.20140825.genotypes_has_ped.vcf.gz

ls 

ALL.wgs.nhgri_coriell_affy_6.20140825.genotypes_has_ped.bed 
ALL.wgs.nhgri_coriell_affy_6.20140825.genotypes_has_ped.bim 
ALL.wgs.nhgri_coriell_affy_6.20140825.genotypes_has_ped.fam 
```

* The pedigree information can be downloaded here:
	- https://www.internationalgenome.org/faq/can-i-get-phenotype-gender-and-family-relationship-information-for-the-individuals/
</details>


## PART 2 -- Pre-Imputation
### Step 1 - Enter the preImpuatation directory

<details> 
	<summary>👇 Steps and code </summary>
	<hr>
	
```
cd preImputation
```
</details>


### Step 2 - Check heterogeneity and missingness 

<details> 
	<summary>👇 Steps and code </summary>
	<hr>
	
* First, run plink commands to calculate heterogenetiy and missingness for the data 
```
module load plink/1.9-20210416
plink --bfile ALL.wgs.nhgri_coriell_affy_6.20140825.genotypes_has_ped --het --out ALL.wgs.nhgri_coriell_affy_6.20140825.genotypes_has_ped_het
plink --bfile ALL.wgs.nhgri_coriell_affy_6.20140825.genotypes_has_ped --missing --out ALL.wgs.nhgri_coriell_affy_6.20140825.genotypes_has_ped_miss
```
* Then, plot in R

```

# Read in 1000 Genomes DATA
setwd("GWAS_QC")
het <- read.csv(file.path("GWAS_QC", "ALL.wgs.nhgri_coriell_affy_6.20140825.genotypes_has_ped_het.het"), sep="")
miss <- read.csv(file.path("GWAS_QC", "ALL.wgs.nhgri_coriell_affy_6.20140825.genotypes_has_ped_miss.smiss"), sep="")


# Organize the data
library(dplyr)
x <- miss %>% select(IID, F_MISS)
y <- het %>% mutate(HR = (`OBS_CT`-`O.HOM.`)/`OBS_CT`) %>% select(IID,HR) ## 1000 Genomes
to_plot <- inner_join(x,y, by = "IID")

# Calculate mean and sd of HR
HR_mean=mean(y$HR)
HR_mean
HR_SD=sd(y$HR)
HR_SD

# Calculate mean of Missingness
MISS_mean=mean(x$F_MISS)
MISS_mean

# Establish color scheme
rbPal <- colorRampPalette(c('#6baed6','#deebf7'))

to_plot$Col <- rbPal(200)[as.numeric(cut(to_plot$F_MISS,breaks = 200))]

png("1000_Genomes_Affy6_3450samples.png")
# Making the plot
plot(to_plot$F_MISS,to_plot$HR, 
     col = ifelse(to_plot$HR <= HR_mean-(3*HR_SD), "#c6dbef", 
                  ifelse(to_plot$HR >= HR_mean+(3*HR_SD), "#c6dbef", to_plot$Col2)),
     las = 1, 
     xlab = "Proportion of Missing Genotypes",
     ylab = "Heterozygosity rate", 
     pch = 19,
     cex = 0.5)
    # These commands make the threshold lines for HR and Missingness 
    # according to the H3ABioNet tutorial
    abline(h=HR_mean+(3*HR_SD),col=2,lty=3)
    abline(h=HR_mean-(3*HR_SD),col=2,lty=3)
    abline(v=MISS_mean,col=2,lty=3)
    abline(v=MISS_mean-(0.05*MISS_mean),col=1,lty=3)
    abline(v=MISS_mean+(0.05*MISS_mean),col=4,lty=3)
dev.off()

```

![1000 Genomes Affy6_3450samples](https://user-images.githubusercontent.com/66582523/183427858-1eaf5f79-4001-4dbf-8222-0c2d50fe7a74.png)

* Figure 2 uses the same code as above, on example internal data
</details>
	

### Step 3 - Update first column of file which has zeros

<details> 
	<summary>👇 Steps and code </summary>
	<hr>
	
```
cat ../ALL.wgs.nhgri_coriell_affy_6.20140825.genotypes_has_ped.fam | awk '{print $1,$2,$2,$2}' > ALL.wgs.nhgri_coriell_affy_6.20140825.genotypes_has_ped_toUpdate.txt

plink2 --bfile ALL.wgs.nhgri_coriell_affy_6.20140825.genotypes_has_ped --update-ids ALL.wgs.nhgri_coriell_affy_6.20140825.genotypes_has_ped_auto_toUpdate.txt --make-bed --out ALL.wgs.nhgri_coriell_affy_6.20140825.genotypes_has_ped_Updated
```
</details>
	
	
### Step 4 -- Add sex phenotype
	
<details> 
	<summary>👇 Steps and code </summary>
	<hr>
	
* First make sex file from 20130606_g1k.ped file

```
join -1 1 -2 1 <(cat ALL.wgs.nhgri_coriell_affy_6.20140825.genotypes_has_ped_Updated.fam |sort -k1,1) <(cat 20130606_g1k.ped |awk -F '\t' '{print $2,$5}' |sort -k1,1) |awk '{print $1,$2,$7}' > sex_file.txt
```
* Then update .fam file to include sex phenotype

```
plink2 --bfile ALL.wgs.nhgri_coriell_affy_6.20140825.genotypes_has_ped_Updated --update-sex sex_file.txt --make-bed --out ALL.wgs.nhgri_coriell_affy_6.20140825.genotypes_has_ped_Updated_withsex
```
* Next, check sex and remove sex inconsistencies

```
plink2 --bfile ALL.wgs.nhgri_coriell_affy_6.20140825.genotypes_has_ped_Updated_withsex --check-sex 
cat plink.sexcheck |  grep PROBLEM | sed 's/^ *//' > plink.sexcheck_list
plink2 --bfile ALL.wgs.nhgri_coriell_affy_6.20140825.genotypes_has_ped_Updated_withsex --remove plink.sexcheck_list --make-bed --out ALL.wgs.nhgri_coriell_affy_6.20140825.genotypes_has_ped_Updated_withsex_checked
```
</details>
	
	
### Step 5 -- Remove SNP variants that do not have SNP IDs
	
<details> 
	<summary>👇 Steps and code </summary>
	<hr>
	
```
echo . > noSNP
plink2 --bfile ALL.wgs.nhgri_coriell_affy_6.20140825.genotypes_has_ped_Updated_withsex_checked --exclude noSNP --make-bed --out ALL.wgs.nhgri_coriell_affy_6.20140825.genotypes_has_ped_Updated_withsex_checked_noDots
```
</details>
	

### Step 6 -- Run QC on data
	
<details> 
	<summary>👇 Steps and code </summary>
	<hr>
```
plink2 --bfile ALL.wgs.nhgri_coriell_affy_6.20140825.genotypes_has_ped_Updated_withsex_checked_noDots --geno 0.05 --mind 0.1 --make-bed --out ALL.wgs.nhgri_coriell_affy_6.20140825.genotypes_has_ped_Updated_withsex_checked_noDots_QC

head ALL.wgs.nhgri_coriell_affy_6.20140825.genotypes_has_ped_Updated_withsex_checked_noDots_QC.bim
```

</details>
	


