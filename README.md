
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


# NOTE: This entire repo is under active construction and code will be mismatched at this time.

# Basic Overview
* We chose a publicly available dataset from the International Genome Sample Resource (IGSR) (www.internationalgenome.org). IGSR created and maintains the 1000 Genomes Project (1kGP) to provide a public catalog of common human genetic variation and genotype data. The 1kGP dataset has been kept up to date with current reference data sets, thus it is available for both GRCh37 and GRCh38. The latter is utilized here because the 2014 update increased the quantity of loci represented, resolved more than 1000 issues from the previous version of the assembly, and overall provides a better basis for alignment and subsequent analysis. Additionally, IGSR’s continued efforts will lead to the incorporation of various populations to the data which were not previously captured.

# Set Up
## Modules in BASH
<details> 
	<summary>👇 How to load modules </summary>
	<hr>
	
```
module load plink/1.9-20210416
module load plink/2.0
module load bcftools/1.9
module load vcftools/0.1.12c # used by bcftools so may not need to be directly loaded
module load tabix/0.2.6 # this one might not work -- let me know if it's a problem use htslib instead (below)
module load htslib
module load perl
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
mkdir GWAS_QC/rawData
mkdir GWAS_QC/preImputation
mkdir GWAS_QC/preImputation/VCFfiles
mkdir GWAS_QC/Imputed
mkdir GWAS_QC/postImpuatation	
```	
</details>

```
.
├── GWAS_QC/ (current directory)
    ├── rawData/
    ├── preImputation/
    |    └── VCFfiles/
    ├── Imputed/
    └── postImputation/	
```


## Notes on PLINK v1.9 and v2.0
* HEREIN, ALL USAGE OF PLINK 1.9 INDICATED BY "plink" AND PLINK 2.0 BY "plink2"
* Not all commands are portable between PLINK version 1.9 and version 2.0. Since PLINK v2.0 is under heavy active development, the developers urge users to check certain results against an earlier, more widely-used version of PLINK. Some functions are available in v1.9 which are not in v2.0, and vice versa. Some of the same functions will produce different file formats and outputs as well.
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

### Step 1 - Retrieve 1000 Genomes Data
<details> 
	<summary>👇 Steps and code </summary>
	<hr>
	
* Enter your GWAS_QC directory 
```
cd GWAS_QC/rawData/
```
* Download -- First, we need to download the publicly available dataset from the 1000 Genomes Project (1KGP). The data is Affy6.0 genotype data for 3,450 individuals with population-level data. 

* Files can be found here: http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/supporting/hd_genotype_chip/
* Two download options 
	- Download directly to your local computer by clicking the hyperlink for `ALL.wgs.nhgri_coriell_affy_6.20140825.genotypes_has_ped.vcf.gz`
	- Dowload using wget command:

```
wget http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/supporting/hd_genotype_chip/ALL.wgs.nhgri_coriell_affy_6.20140825.genotypes_has_ped.vcf.gz
```

* At this point `ALL.wgs.nhgri_coriell_affy_6.20140825.genotypes_has_ped.vcf.gz` should be in your GWAS_QC directory.

```
gunzip ALL.wgs.nhgri_coriell_affy_6.20140825.genotypes_has_ped.vcf.gz	
```

* Then recode the `.vcf` file as `.bed`, `.bim`, `.bed`. 

```
plink2 --vcf ALL.wgs.nhgri_coriell_affy_6.20140825.genotypes_has_ped.vcf --make-bed --out ALL.wgs.nhgri_coriell_affy_6.20140825.genotypes_has_ped	
```

* When we check the directory, the following files should be unzipped:

```
ls 	
```
> ALL.wgs.nhgri_coriell_affy_6.20140825.genotypes_has_ped.bed 
> ALL.wgs.nhgri_coriell_affy_6.20140825.genotypes_has_ped.bim 
> ALL.wgs.nhgri_coriell_affy_6.20140825.genotypes_has_ped.fam 
> ALL.wgs.nhgri_coriell_affy_6.20140825.genotypes_has_ped.log
> ALL.wgs.nhgri_coriell_affy_6.20140825.genotypes_has_ped.vcf

</details>


### Step 2 - Download other relevant data
<details> 
	<summary>👇 Steps and code </summary>
	<hr>

* You will also want to dowload the pedigree information, which can be downloaded here:
	- https://www.internationalgenome.org/faq/can-i-get-phenotype-gender-and-family-relationship-information-for-the-individuals/
	-  20130606_g1k.ped
	
* You also need to dowload the TOPMed reference panel and code from the McCarthy Tools website
	- https://www.well.ox.ac.uk/~wrayner/tools/ 
	- Midway down the page there are instructions for Usage with the TOPMed reference panel
	- Download `CreateTOPMed.zip` from the McCarthy Tools website
	![image](https://user-images.githubusercontent.com/30478823/192685817-b8e3cc46-045d-4381-a55b-24e9e74bb64d.png)
	- Click on the link to navigate to the Bravo site where this can be downloaded from
	- Click on dbSNP in the top right panel and click "Download VCF" button to download `ALL.TOPMed_freeze5_hg38_dbSNP.vcf.gz` file
	![image](https://user-images.githubusercontent.com/30478823/192685864-a94ededb-84b0-4cd3-b921-ad627ff0f463.png)
	- Note: If you run the curl command that's given on the site, the filename will be different 
	- After downloading both files, move them to the GWAS_QC directory, and execute the following command to convert VCF to an HRC formatted reference legend according to the code provided by McCarthy Tools:
	```
	./CreateTOPMed.pl -i ALL.TOPMed_freeze5_hg38_dbSNP.vcf.gz
	```
	- By default this will create a file filtered for variants flagged as PASS only, if you wish to use all variants the -a flag overrides this. To override the default output file naming use -o filename.
	- If you get an error in the above step, try this variation. Both were run successfully on a local computer and server using perl/5.30.0. Depending on your setup, this may take a few hours to run:
	```
	perl CreateTOPMed.pl -i ALL.TOPMed_freeze5_hg38_dbSNP.vcf.gz	
	```
</details>

```
.
├── GWAS_QC/
    ├── rawData/ (current directory)
    |    ├── ALL.wgs.nhgri_coriell_affy_6.20140825.genotypes_has_ped.bed
    |    ├── ALL.wgs.nhgri_coriell_affy_6.20140825.genotypes_has_ped.bim
    |    ├── ALL.wgs.nhgri_coriell_affy_6.20140825.genotypes_has_ped.fam
    |    ├── ALL.wgs.nhgri_coriell_affy_6.20140825.genotypes_has_ped.log
    |    ├── ALL.wgs.nhgri_coriell_affy_6.20140825.genotypes_has_ped.vcf
    |    ├── CreateTOPMed.pl
    |    ├── CreateTOPMed.zip
    |    ├── LICENSE.txt
    |    └── PASS.Variants.TOPMed_freeze5_hg38_dbSNP.tab.gz
    ├──preImputation/
    |    └── VCFfiles/
    ├──Imputed/
    └── postImputation/	
```


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
	
* First, run plink commands to calculate heterogeneity and missingness for the data 

```
module load plink/2.0
plink --bfile ../rawData/ALL.wgs.nhgri_coriell_affy_6.20140825.genotypes_has_ped --het --out ALL.wgs.nhgri_coriell_affy_6.20140825.genotypes_has_ped_het	
```

* Output:
> PLINK v2.00a3LM 64-bit Intel (28 Oct 2020)     www.cog-genomics.org/plink/2.0/
> (C) 2005-2020 Shaun Purcell, Christopher Chang   GNU General Public License v3
> Logging to ALL.wgs.nhgri_coriell_affy_6.20140825.genotypes_has_ped_het.log.
> Options in effect:
>   --bfile ../rawData/ALL.wgs.nhgri_coriell_affy_6.20140825.genotypes_has_ped
>   --het
>   --out ALL.wgs.nhgri_coriell_affy_6.20140825.genotypes_has_ped_het
> 
> Start time: Wed Sep 28 13:22:52 2022
> 128235 MiB RAM detected; reserving 64117 MiB for main workspace.
> Using up to 56 threads (change this with --threads).
> 3450 samples (0 females, 0 males, 3450 ambiguous; 3450 founders) loaded from
> ../rawData/ALL.wgs.nhgri_coriell_affy_6.20140825.genotypes_has_ped.fam.
> 905788 variants loaded from
> ../rawData/ALL.wgs.nhgri_coriell_affy_6.20140825.genotypes_has_ped.bim.
> Note: No phenotype data present.
> Calculating allele frequencies... done.
> Excluding 37156 variants on non-autosomes from --het.
> --het: done.
> --het: Results written to
> ALL.wgs.nhgri_coriell_affy_6.20140825.genotypes_has_ped_het.het .

```	
plink --bfile ../rawData/ALL.wgs.nhgri_coriell_affy_6.20140825.genotypes_has_ped --missing --out ALL.wgs.nhgri_coriell_affy_6.20140825.genotypes_has_ped_miss	
```

* Output:
> PLINK v2.00a3LM 64-bit Intel (28 Oct 2020)     www.cog-genomics.org/plink/2.0/
> (C) 2005-2020 Shaun Purcell, Christopher Chang   GNU General Public License v3
> Logging to ALL.wgs.nhgri_coriell_affy_6.20140825.genotypes_has_ped_miss.log.
> Options in effect:
>   --bfile ../rawData/ALL.wgs.nhgri_coriell_affy_6.20140825.genotypes_has_ped
>   --missing
>   --out ALL.wgs.nhgri_coriell_affy_6.20140825.genotypes_has_ped_miss
> 
> 128235 MiB RAM detected; reserving 64117 MiB for main workspace.
> Using up to 56 threads (change this with --threads).
> 3450 samples (0 females, 0 males, 3450 ambiguous; 3450 founders) loaded from
> ../rawData/ALL.wgs.nhgri_coriell_affy_6.20140825.genotypes_has_ped.fam.
> 905788 variants loaded from
> ../rawData/ALL.wgs.nhgri_coriell_affy_6.20140825.genotypes_has_ped.bim.
> Note: No phenotype data present.
> Calculating sample missingness rates... done.
> Calculating allele frequencies... done.
> --missing: Sample missing data report written to
> ALL.wgs.nhgri_coriell_affy_6.20140825.genotypes_has_ped_miss.smiss .
> --missing: Variant missing data report written to
> ALL.wgs.nhgri_coriell_affy_6.20140825.genotypes_has_ped_miss.vmiss .

	
* Here's the current directory structure within GWAS_QC/preImputation/	
> ALL.wgs.nhgri_coriell_affy_6.20140825.genotypes_has_ped_het.het
> ALL.wgs.nhgri_coriell_affy_6.20140825.genotypes_has_ped_het.log
> ALL.wgs.nhgri_coriell_affy_6.20140825.genotypes_has_ped_miss.log
> ALL.wgs.nhgri_coriell_affy_6.20140825.genotypes_has_ped_miss.smiss
> ALL.wgs.nhgri_coriell_affy_6.20140825.genotypes_has_ped_miss.vmiss
```
.
├── GWAS_QC/
    ├── rawData/ 
    |    ├── ALL.wgs.nhgri_coriell_affy_6.20140825.genotypes_has_ped.bed
    |    ├── ALL.wgs.nhgri_coriell_affy_6.20140825.genotypes_has_ped.bim
    |    ├── ALL.wgs.nhgri_coriell_affy_6.20140825.genotypes_has_ped.fam
    |    ├── ALL.wgs.nhgri_coriell_affy_6.20140825.genotypes_has_ped.log
    |    ├── ALL.wgs.nhgri_coriell_affy_6.20140825.genotypes_has_ped.vcf
    |    ├── CreateTOPMed.pl
    |    ├── CreateTOPMed.zip
    |    ├── LICENSE.txt
    |    └── PASS.Variants.TOPMed_freeze5_hg38_dbSNP.tab.gz
    ├──preImputation/ (current directory)
    |    ├── ALL.wgs.nhgri_coriell_affy_6.20140825.genotypes_has_ped_het.het
    |    ├── ALL.wgs.nhgri_coriell_affy_6.20140825.genotypes_has_ped_het.log
    |    ├── ALL.wgs.nhgri_coriell_affy_6.20140825.genotypes_has_ped_miss.log
    |    ├── ALL.wgs.nhgri_coriell_affy_6.20140825.genotypes_has_ped_miss.smiss
    |    └── ALL.wgs.nhgri_coriell_affy_6.20140825.genotypes_has_ped_miss.vmiss
    ├──Imputed/
    └── postImputation/	
```

	
* Now, let's make some plots in R so first we need to enter into R in the terminal. You can run R code in Rstudio or your preferred IDE. 
	
```
R
```

* Read in the relevant data for the plots

```
# Read in 1000 Genomes DATA
setwd("~/GWAS_QC/preImputation/")
het <- read.csv(file.path("ALL.wgs.nhgri_coriell_affy_6.20140825.genotypes_has_ped_het.het"), sep="")     
```

> ```
> Output: 
> X.FID     IID O.HOM. E.HOM. OBS_CT         F
> 1     0 HG00096 634762 608326 867229 0.1021090
> 2     0 HG00097 635035 607916 866659 0.1048090
> 3     0 HG00098 629539 604988 861939 0.0955457
> 4     0 HG00099 635235 607948 866667 0.1054700
> 5     0 HG00100 633342 607067 865511 0.1016670
> 6     0 HG00101 633098 607641 866215 0.0984509
> ```      

```
miss <- read.csv(file.path("ALL.wgs.nhgri_coriell_affy_6.20140825.genotypes_has_ped_miss.smiss"), sep="")
```
> ```
>  X.FID     IID MISSING_CT OBS_CT     F_MISS
> 1     0 HG00096       1403 905531 0.00154937
> 2     0 HG00097       2016 905531 0.00222632
> 3     0 HG00098       6899 905531 0.00761873
> 4     0 HG00099       2042 905531 0.00225503
> 5     0 HG00100       3210 905531 0.00354488
> 6     0 HG00101       2428 905531 0.00268130
> ```

```
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
cat ../rawData/ALL.wgs.nhgri_coriell_affy_6.20140825.genotypes_has_ped.fam | awk '{print $1,$2,$2,$2}' > ALL.wgs.nhgri_coriell_affy_6.20140825.genotypes_has_ped_toUpdate.txt

plink2 --bfile ../rawData/ALL.wgs.nhgri_coriell_affy_6.20140825.genotypes_has_ped --update-ids ALL.wgs.nhgri_coriell_affy_6.20140825.genotypes_has_ped_toUpdate.txt --make-bed --out ALL.wgs.nhgri_coriell_affy_6.20140825.genotypes_has_ped_Updated
```
> PLINK v2.00a3LM 64-bit Intel (5 May 2021)      www.cog-genomics.org/plink/2.0/
> (C) 2005-2021 Shaun Purcell, Christopher Chang   GNU General Public License v3
> Logging to ALL.wgs.nhgri_coriell_affy_6.20140825.genotypes_has_ped_Updated.log.
> Options in effect:
>  --bfile ../rawData/ALL.wgs.nhgri_coriell_affy_6.20140825.genotypes_has_ped
>  --make-bed
>  --out ALL.wgs.nhgri_coriell_affy_6.20140825.genotypes_has_ped_Updated
>  --update-ids ALL.wgs.nhgri_coriell_affy_6.20140825.genotypes_has_ped_toUpdate.txt
> Start time: Mon Oct 10 16:25:54 2022
> 128235 MiB RAM detected; reserving 64117 MiB for main workspace.
> Using up to 56 threads (change this with --threads).
> 3450 samples (0 females, 0 males, 3450 ambiguous; 3450 founders) loaded from
> ../rawData/ALL.wgs.nhgri_coriell_affy_6.20140825.genotypes_has_ped.fam.
> 905788 variants loaded from
> ../rawData/ALL.wgs.nhgri_coriell_affy_6.20140825.genotypes_has_ped.bim.
> Note: No phenotype data present.
> --update-ids: 3450 samples updated.
> Writing ALL.wgs.nhgri_coriell_affy_6.20140825.genotypes_has_ped_Updated.fam ...
> done.
> Writing ALL.wgs.nhgri_coriell_affy_6.20140825.genotypes_has_ped_Updated.bim ...
> done.
> Writing ALL.wgs.nhgri_coriell_affy_6.20140825.genotypes_has_ped_Updated.bed ...
> done.

</details>

### Step 4 -- Add sex phenotype
	
<details> 
	<summary>👇 Steps and code </summary>
	<hr>
	
* First make sex file from 20130606_g1k.ped file

```
join -1 1 -2 1 <(cat ALL.wgs.nhgri_coriell_affy_6.20140825.genotypes_has_ped_Updated.fam |sort -k1,1) <(cat ../rawData/20130606_g1k.ped |awk -F '\t' '{print $2,$5}' |sort -k1,1) |awk '{print $1,$2,$7}' > sex_file.txt
```

* Then update .fam file to include sex phenotype
```
plink2 --bfile ALL.wgs.nhgri_coriell_affy_6.20140825.genotypes_has_ped_Updated --update-sex sex_file.txt --make-bed --out ALL.wgs.nhgri_coriell_affy_6.20140825.genotypes_has_ped_Updated_withsex
```
> PLINK v2.00a3LM 64-bit Intel (5 May 2021)      www.cog-genomics.org/plink/2.0/
> (C) 2005-2021 Shaun Purcell, Christopher Chang   GNU General Public License v3
> Logging to ALL.wgs.nhgri_coriell_affy_6.20140825.genotypes_has_ped_Updated_withsex.log.
> Options in effect:
>   --bfile ALL.wgs.nhgri_coriell_affy_6.20140825.genotypes_has_ped_Updated
>   --make-bed
>   --out ALL.wgs.nhgri_coriell_affy_6.20140825.genotypes_has_ped_Updated_withsex
>   --update-sex sex_file.txt
> Start time: Mon Oct 10 16:30:08 2022
> 128235 MiB RAM detected; reserving 64117 MiB for main workspace.
> Using up to 56 threads (change this with --threads).
> 3450 samples (0 females, 0 males, 3450 ambiguous; 3450 founders) loaded from
> ALL.wgs.nhgri_coriell_affy_6.20140825.genotypes_has_ped_Updated.fam.
> 905788 variants loaded from
> ALL.wgs.nhgri_coriell_affy_6.20140825.genotypes_has_ped_Updated.bim.
> Note: No phenotype data present.
> --update-sex: 3450 samples updated.
> Writing
> ALL.wgs.nhgri_coriell_affy_6.20140825.genotypes_has_ped_Updated_withsex.fam ...
> done.
> Writing
> ALL.wgs.nhgri_coriell_affy_6.20140825.genotypes_has_ped_Updated_withsex.bim ...
> done.
> Writing
> ALL.wgs.nhgri_coriell_affy_6.20140825.genotypes_has_ped_Updated_withsex.bed ...
> done.
> End time: Mon Oct 10 16:30:11 2022	
	
	
* Next, check sex and remove sex inconsistencies
```
plink --bfile ALL.wgs.nhgri_coriell_affy_6.20140825.genotypes_has_ped_Updated_withsex --check-sex 
```
> //' > plink.sexcheck_listPLINK v1.90p 64-bit (16 Apr 2021)              www.cog-genomics.org/plink/1.9/
> (C) 2005-2021 Shaun Purcell, Christopher Chang   GNU General Public License v3
> Logging to plink.log.
> Options in effect:
>   --bfile ALL.wgs.nhgri_coriell_affy_6.20140825.genotypes_has_ped_Updated_withsex
>   --check-sex
> 128235 MB RAM detected; reserving 64117 MB for main workspace.
> 905788 variants loaded from .bim file.
> 3450 people (1715 males, 1735 females) loaded from .fam.
> Using 1 thread (no multithreaded calculations invoked).
> Before main variant filters, 3450 founders and 0 nonfounders present.
> Calculating allele frequencies... done.
> Warning: 219657 het. haploid genotypes present (see plink.hh ); many commands
> treat these as missing.
> Warning: Nonmissing nonmale Y chromosome genotype(s) present; many commands
> treat these as missing.
> Total genotyping rate is 0.996711.
> 905788 variants and 3450 people pass filters and QC.
> Note: No phenotypes present.
> --check-sex: 36449 Xchr and 0 Ychr variant(s) scanned, 577 problems detected.
> Report written to plink.sexcheck .

	
```
cat plink.sexcheck |  grep PROBLEM | sed 's/^ *//' > plink.sexcheck_list
```

```
plink2 --bfile ALL.wgs.nhgri_coriell_affy_6.20140825.genotypes_has_ped_Updated_withsex --remove plink.sexcheck_list --make-bed --out ALL.wgs.nhgri_coriell_affy_6.20140825.genotypes_has_ped_Updated_withsex_checked
```
> PLINK v2.00a3LM 64-bit Intel (5 May 2021)      www.cog-genomics.org/plink/2.0/
> (C) 2005-2021 Shaun Purcell, Christopher Chang   GNU General Public License v3
> Logging to ALL.wgs.nhgri_coriell_affy_6.20140825.genotypes_has_ped_Updated_withsex_checked.log.
> Options in effect:
>   --bfile ALL.wgs.nhgri_coriell_affy_6.20140825.genotypes_has_ped_Updated_withsex
>   --make-bed
>   --out ALL.wgs.nhgri_coriell_affy_6.20140825.genotypes_has_ped_Updated_withsex_checked
>   --remove plink.sexcheck_list
> Start time: Mon Oct 10 16:35:55 2022
> 128235 MiB RAM detected; reserving 64117 MiB for main workspace.
> Using up to 56 threads (change this with --threads).
> 3450 samples (1735 females, 1715 males; 3450 founders) loaded from
> ALL.wgs.nhgri_coriell_affy_6.20140825.genotypes_has_ped_Updated_withsex.fam.
> 905788 variants loaded from
> ALL.wgs.nhgri_coriell_affy_6.20140825.genotypes_has_ped_Updated_withsex.bim.
> Note: No phenotype data present.
> --remove: 2873 samples remaining.
> 2873 samples (1160 females, 1713 males; 2873 founders) remaining after main
> filters.
> Writing
> ALL.wgs.nhgri_coriell_affy_6.20140825.genotypes_has_ped_Updated_withsex_checked.fam
> ... done.
> Writing
> ALL.wgs.nhgri_coriell_affy_6.20140825.genotypes_has_ped_Updated_withsex_checked.bim
> ... done.
> Writing
> ALL.wgs.nhgri_coriell_affy_6.20140825.genotypes_has_ped_Updated_withsex_checked.bed
> ... done.
> End time: Mon Oct 10 16:35:57 2022


</details>
	
### Step 5 -- Remove SNP variants that do not have SNP IDs
	
<details> 
	<summary>👇 Steps and code </summary>
	<hr>
	
```
echo . > noSNP
plink2 --bfile ALL.wgs.nhgri_coriell_affy_6.20140825.genotypes_has_ped_Updated_withsex_checked --exclude noSNP --make-bed --out ALL.wgs.nhgri_coriell_affy_6.20140825.genotypes_has_ped_Updated_withsex_checked_noDots
```
> PLINK v2.00a3LM 64-bit Intel (5 May 2021)      www.cog-genomics.org/plink/2.0/
> (C) 2005-2021 Shaun Purcell, Christopher Chang   GNU General Public License v3
> Logging to ALL.wgs.nhgri_coriell_affy_6.20140825.genotypes_has_ped_Updated_withsex_checked_noDots.log.
> Options in effect:
>   --bfile ALL.wgs.nhgri_coriell_affy_6.20140825.genotypes_has_ped_Updated_withsex_checked
>   --exclude noSNP
>   --make-bed
>   --out ALL.wgs.nhgri_coriell_affy_6.20140825.genotypes_has_ped_Updated_withsex_checked_noDots
> Start time: Mon Oct 10 16:38:29 2022
> 128235 MiB RAM detected; reserving 64117 MiB for main workspace.
> Using up to 56 threads (change this with --threads).
> 2873 samples (1160 females, 1713 males; 2873 founders) loaded from
> ALL.wgs.nhgri_coriell_affy_6.20140825.genotypes_has_ped_Updated_withsex_checked.fam.
> 905788 variants loaded from
> ALL.wgs.nhgri_coriell_affy_6.20140825.genotypes_has_ped_Updated_withsex_checked.bim.
> Note: No phenotype data present.
> --exclude: 887969 variants remaining.
> 887969 variants remaining after main filters.
> Writing
> ALL.wgs.nhgri_coriell_affy_6.20140825.genotypes_has_ped_Updated_withsex_checked_noDots.fam
> ... done.
> Writing
> ALL.wgs.nhgri_coriell_affy_6.20140825.genotypes_has_ped_Updated_withsex_checked_noDots.bim
> ... done.
> Writing
> ALL.wgs.nhgri_coriell_affy_6.20140825.genotypes_has_ped_Updated_withsex_checked_noDots.bed
> ... done.
> End time: Mon Oct 10 16:38:32 2022

</details>
	

### Step 6 -- Run QC on data
	
<details> 
	<summary>👇 Steps and code </summary>
	<hr>
	
```
plink2 --bfile ALL.wgs.nhgri_coriell_affy_6.20140825.genotypes_has_ped_Updated_withsex_checked_noDots --geno 0.05 --mind 0.1 --make-bed --out ALL.wgs.nhgri_coriell_affy_6.20140825.genotypes_has_ped_Updated_withsex_checked_noDots_QC
```
> PLINK v2.00a3LM 64-bit Intel (5 May 2021)      www.cog-genomics.org/plink/2.0/
> (C) 2005-2021 Shaun Purcell, Christopher Chang   GNU General Public License v3
> Logging to ALL.wgs.nhgri_coriell_affy_6.20140825.genotypes_has_ped_Updated_withsex_checked_noDots_QC.log.
> Options in effect:
>   --bfile ALL.wgs.nhgri_coriell_affy_6.20140825.genotypes_has_ped_Updated_withsex_checked_noDots
>   --geno 0.05
>   --make-bed
>   --mind 0.1
>   --out ALL.wgs.nhgri_coriell_affy_6.20140825.genotypes_has_ped_Updated_withsex_checked_noDots_QC
> Start time: Mon Oct 10 16:40:32 2022
> 128235 MiB RAM detected; reserving 64117 MiB for main workspace.
> Using up to 56 threads (change this with --threads).
> 2873 samples (1160 females, 1713 males; 2873 founders) loaded from
> ALL.wgs.nhgri_coriell_affy_6.20140825.genotypes_has_ped_Updated_withsex_checked_noDots.fam.
> 887969 variants loaded from
> ALL.wgs.nhgri_coriell_affy_6.20140825.genotypes_has_ped_Updated_withsex_checked_noDots.bim.
> Note: No phenotype data present.
> Calculating sample missingness rates... done.
> 0 samples removed due to missing genotype data (--mind).
> 2873 samples (1160 females, 1713 males; 2873 founders) remaining after main
> filters.
> Calculating allele frequencies... done.
> --geno: 7979 variants removed due to missing genotype data.
> 879990 variants remaining after main filters.
> Writing
> ALL.wgs.nhgri_coriell_affy_6.20140825.genotypes_has_ped_Updated_withsex_checked_noDots_QC.fam
> ... done.
> Writing
> ALL.wgs.nhgri_coriell_affy_6.20140825.genotypes_has_ped_Updated_withsex_checked_noDots_QC.bim
> ... done.
> Writing
> ALL.wgs.nhgri_coriell_affy_6.20140825.genotypes_has_ped_Updated_withsex_checked_noDots_QC.bed
> ... done.
> End time: Mon Oct 10 16:40:34 2022

```
head ALL.wgs.nhgri_coriell_affy_6.20140825.genotypes_has_ped_Updated_withsex_checked_noDots_QC.bim
```
> ```
> 1       rs10458597      0       564621  T       C
> 1       rs12565286      0       721290  C       G
> 1       rs12082473      0       740857  A       G
> 1       rs3094315       0       752566  A       G
> 1       rs2286139       0       761732  T       C
> 1       rs11240776      0       765269  G       A
> 1       rs2980319       0       777122  T       A
> 1       rs2980300       0       785989  C       T
> 1       rs2905036       0       792480  T       C
> 1       rs11240777      0       798959  A       G
> ```
</details>
	

### Step 7 -- LiftOver the data

<details> 
	<summary>👇 Steps and code </summary>
	<hr>
	
* We found out that the data was built using hg37 so we will to convert the data to hg38 using the liftOver module
	- First, make BED coordinate file

```
cat ALL.wgs.nhgri_coriell_affy_6.20140825.genotypes_has_ped_Updated_withsex_checked_noDots_QC.bim | awk '{print "chr"$1, $4, ($4+1), $4, $2}' > liftover_input.BED
```

* Check output
```
head liftover_input.BED

```
> ```
> chr1 564621 564622 564621 rs10458597
> chr1 721290 721291 721290 rs12565286
> chr1 740857 740858 740857 rs12082473
> chr1 752566 752567 752566 rs3094315
> chr1 761732 761733 761732 rs2286139
> chr1 765269 765270 765269 rs11240776
> chr1 777122 777123 777122 rs2980319
> chr1 785989 785990 785989 rs2980300
> chr1 792480 792481 792480 rs2905036
> chr1 798959 798960 798959 rs11240777 
> ```

```
sed -i 's/chr23/chrX/g' liftover_input.BED #I don't think this is a problem with this dataset
```

* Then, dowload the correct liftover file
	- http://hgdownload.cse.ucsc.edu/goldenpath/hg18/liftOver/
	- hg18ToHg38.over.chain
```
wget http://hgdownload.cse.ucsc.edu/goldenpath/hg18/liftOver/hg18ToHg38.over.chain
```	

* Finally, perform the actual liftOver
```
liftOver liftover_input.BED hg18ToHg38.over.chain liftover_newmap.txt liftover_exclude.txt

sed -i 's/chr//g' liftover_newmap.txt
awk '{print $5,$2}' liftover_newmap.txt > update_map.txt
cat liftover_exclude.txt | grep -v "#" | awk '{print $5}' > exclude_liftover.txt
```
	
</details>


### Step 8 -- Exclude data
	
<details> 
	<summary>👇 Steps and code </summary>
	<hr>
	
* Exclude any SNPs that do not liftOver and non-somatic chromosomes (X, Y)
```
plink2 --bfile ALL.wgs.nhgri_coriell_affy_6.20140825.genotypes_has_ped_Updated_withsex_checked_noDots_QC --exclude exclude_liftover.txt --update-map update_map.txt --not-chr X, Y --make-bed --out ALL.wgs.nhgri_coriell_affy_6.20140825.genotypes_has_ped_Updated_withsex_checked_noDots_QC_b38
```
> PLINK v2.00a3LM 64-bit Intel (5 May 2021)      www.cog-genomics.org/plink/2.0/
> (C) 2005-2021 Shaun Purcell, Christopher Chang   GNU General Public License v3
> Logging to ALL.wgs.nhgri_coriell_affy_6.20140825.genotypes_has_ped_Updated_withsex_checked_noDots_QC_b38.log.
> Options in effect:
>   --bfile ALL.wgs.nhgri_coriell_affy_6.20140825.genotypes_has_ped_Updated_withsex_checked_noDots_QC
>   --exclude exclude_liftover.txt
>   --make-bed
>   --not-chr X, Y
>   --out ALL.wgs.nhgri_coriell_affy_6.20140825.genotypes_has_ped_Updated_withsex_checked_noDots_QC_b38
>   --update-map update_map.txt
> Start time: Mon Oct 10 17:06:21 2022
> 128235 MiB RAM detected; reserving 64117 MiB for main workspace.
> Using up to 56 threads (change this with --threads).
> 2873 samples (1160 females, 1713 males; 2873 founders) loaded from
> ALL.wgs.nhgri_coriell_affy_6.20140825.genotypes_has_ped_Updated_withsex_checked_noDots_QC.fam.
> 843232 out of 879990 variants loaded from
> ALL.wgs.nhgri_coriell_affy_6.20140825.genotypes_has_ped_Updated_withsex_checked_noDots_QC.bim.
> Note: No phenotype data present.
> --update-map: 0 values updated.
> --exclude: 843232 variants remaining.
> 843232 variants remaining after main filters.
> Writing
> ALL.wgs.nhgri_coriell_affy_6.20140825.genotypes_has_ped_Updated_withsex_checked_noDots_QC_b38.fam
> ... done.
> Writing
> ALL.wgs.nhgri_coriell_affy_6.20140825.genotypes_has_ped_Updated_withsex_checked_noDots_QC_b38.bim
> ... done.
> Writing
> ALL.wgs.nhgri_coriell_affy_6.20140825.genotypes_has_ped_Updated_withsex_checked_noDots_QC_b38.bed
> ... done.
> End time: Mon Oct 10 17:06:23 2022


* Final pre-imputation variant count is 834,972 SNPs
```
wc -l <filename>
```
> * 905788 ALL.wgs.nhgri_coriell_affy_6.20140825.genotypes_has_ped_Updated.bim
> * 905788 ALL.wgs.nhgri_coriell_affy_6.20140825.genotypes_has_ped_Updated_withsex.bim
> * 905788 ALL.wgs.nhgri_coriell_affy_6.20140825.genotypes_has_ped_Updated_withsex_checked.bim
> * 887969 ALL.wgs.nhgri_coriell_affy_6.20140825.genotypes_has_ped_Updated_withsex_checked_noDots.bim
> * 879990 ALL.wgs.nhgri_coriell_affy_6.20140825.genotypes_has_ped_Updated_withsex_checked_noDots_QC.bim 
> * 834872 ALL.wgs.nhgri_coriell_affy_6.20140825.genotypes_has_ped_Updated_withsex_checked_noDots_QC_b38.bim 
>	* TO-DO Van got 843232 ALL.wgs.nhgri_coriell_affy_6.20140825.genotypes_has_ped_Updated_withsex_checked_noDots_QC_b38.bim

</details>
	

	
### Step 9 -- Principal Component Analysis (PCA)
	
<details> 
	<summary>👇 Steps and code </summary>
	<hr>

* First, run PCA on data using the plink_pca.sh script

```
module load plink_pca
plink_pca.sh -b ALL.wgs.nhgri_coriell_affy_6.20140825.genotypes_has_ped_Updated_withsex_checked_noDots_QC_b38 -g PCA_no1KG
```

> ```
> ...(some other printed output above this too)
> 128235 MB RAM detected; reserving 64117 MB for main workspace.
> 471514 variants loaded from .bim file.
> 2873 people (1713 males, 1160 females) loaded from .fam.
> Using 1 thread.
> Before main variant filters, 2873 founders and 0 nonfounders present.
> Calculating allele frequencies... done.
> Total genotyping rate is 0.997655.
> 471514 variants and 2873 people pass filters and QC.
> Note: No phenotypes present.
> Pruned 19116 variants from chromosome 1, leaving 18784.
> Pruned 19194 variants from chromosome 2, leaving 19285.
> Pruned 16293 variants from chromosome 3, leaving 16882.
> Pruned 15068 variants from chromosome 4, leaving 15207.
> Pruned 15643 variants from chromosome 5, leaving 15452.
> Pruned 16122 variants from chromosome 6, leaving 15251.
> Pruned 13104 variants from chromosome 7, leaving 13343.
> Pruned 12860 variants from chromosome 8, leaving 13062.
> Pruned 11610 variants from chromosome 9, leaving 11571.
> Pruned 13056 variants from chromosome 10, leaving 13071.
> Pruned 12506 variants from chromosome 11, leaving 11918.
> Pruned 11280 variants from chromosome 12, leaving 11950.
> Pruned 9480 variants from chromosome 13, leaving 9483.
> Pruned 7267 variants from chromosome 14, leaving 8155.
> Pruned 6011 variants from chromosome 15, leaving 7515.
> Pruned 6407 variants from chromosome 16, leaving 8114.
> Pruned 4564 variants from chromosome 17, leaving 6281.
> Pruned 6978 variants from chromosome 18, leaving 7804.
> Pruned 2692 variants from chromosome 19, leaving 3985.
> Pruned 5791 variants from chromosome 20, leaving 6533.
> Pruned 3195 variants from chromosome 21, leaving 3735.
> Pruned 2418 variants from chromosome 22, leaving 3478.
> Pruning complete.  230655 of 471514 variants removed.
> Marker lists written to
> plink_pca.ALL.wgs.nhgri_coriell_affy_6.20140825.genotypes_has_ped_Updated_withsex_checked_noDots_QC_b38.ld_list.prune.in
> and
> plink_pca.ALL.wgs.nhgri_coriell_affy_6.20140825.genotypes_has_ped_Updated_withsex_checked_noDots_QC_b38.ld_list.prune.out
> .
> PLINK v1.90b6.18 64-bit (16 Jun 2020)          www.cog-genomics.org/plink/1.9/
> (C) 2005-2020 Shaun Purcell, Christopher Chang   GNU General Public License v3
> Logging to plink_pca.ALL.wgs.nhgri_coriell_affy_6.20140825.genotypes_has_ped_Updated_withsex_checked_noDots_QC_b38.log.
> Options in effect:
>   --bfile preld.plink_pca.ALL.wgs.nhgri_coriell_affy_6.20140825.genotypes_has_ped_Updated_withsex_checked_noDots_QC_b38
>   --extract plink_pca.ALL.wgs.nhgri_coriell_affy_6.20140825.genotypes_has_ped_Updated_withsex_checked_noDots_QC_b38.ld_list.prune.in
>   --make-bed
>   --out plink_pca.ALL.wgs.nhgri_coriell_affy_6.20140825.genotypes_has_ped_Updated_withsex_checked_noDots_QC_b38
> 128235 MB RAM detected; reserving 64117 MB for main workspace.
> 471514 variants loaded from .bim file.
> 2873 people (1713 males, 1160 females) loaded from .fam.
> --extract: 240859 variants remaining.
> Using 1 thread (no multithreaded calculations invoked).
> Before main variant filters, 2873 founders and 0 nonfounders present.
> Calculating allele frequencies... done.
> Total genotyping rate is 0.99744.
> 240859 variants and 2873 people pass filters and QC.
> Note: No phenotypes present.
> --make-bed to
> plink_pca.ALL.wgs.nhgri_coriell_affy_6.20140825.genotypes_has_ped_Updated_withsex_checked_noDots_QC_b38.bed
> +
> plink_pca.ALL.wgs.nhgri_coriell_affy_6.20140825.genotypes_has_ped_Updated_withsex_checked_noDots_QC_b38.bim
> +
> plink_pca.ALL.wgs.nhgri_coriell_affy_6.20140825.genotypes_has_ped_Updated_withsex_checked_noDots_QC_b38.fam
> ... done.
> rm: cannot remove ‘preld.plink_pca.ALL.wgs.nhgri_coriell_affy_6.20140825.genotypes_has_ped_Updated_withsex_checked_noDots_QC_b38.nosex’: > No such file or directory
> rm: cannot remove ‘preld.plink_pca.ALL.wgs.nhgri_coriell_affy_6.20140825.genotypes_has_ped_Updated_withsex_checked_noDots_QC_b38.*.extracted.*’: No such file or directory
> rm: cannot remove ‘preld.plink_pca.ALL.wgs.nhgri_coriell_affy_6.20140825.genotypes_has_ped_Updated_withsex_checked_noDots_QC_b38.bed’: No such file or directory
> rm: cannot remove ‘preld.plink_pca.ALL.wgs.nhgri_coriell_affy_6.20140825.genotypes_has_ped_Updated_withsex_checked_noDots_QC_b38.bim’: No such file or directory
> rm: cannot remove ‘preld.plink_pca.ALL.wgs.nhgri_coriell_affy_6.20140825.genotypes_has_ped_Updated_withsex_checked_noDots_QC_b38.fam’: No such file or directory
> rm: cannot remove ‘preld.plink_pca.ALL.wgs.nhgri_coriell_affy_6.20140825.genotypes_has_ped_Updated_withsex_checked_noDots_QC_b38.log’: No such file or directory
> par_f=/tmp/tmp.2OQ8Ahb50A
> running smartpca...
> *> *> *Registered S3 methods overwritten by 'tibble':
>   method     from
>   format.tbl pillar
>   print.tbl  pillar
> Created Scree plot file plink_pca.ALL.wgs.nhgri_coriell_affy_6.20140825.genotypes_has_ped_Updated_withsex_checked_noDots_QC_b38_ScreePlot.png
> Created plot data file plink_pca.ALL.wgs.nhgri_coriell_affy_6.20140825.genotypes_has_ped_Updated_withsex_checked_noDots_QC_b38_variance.txt
> rm: cannot remove ‘preld.plink_pca.ALL.wgs.nhgri_coriell_affy_6.20140825.genotypes_has_ped_Updated_withsex_checked_noDots_QC_b38.*’: No such file or directory
> ```	
	
	
* The output should be 5 files that are made:
```
ls plink_pca.ALL.wgs.nhgri_coriell_affy_6.20140825.geno*
```

> * plink_pca.ALL.wgs.nhgri_coriell_affy_6.20140825.genotypes_has_ped_Updated_withsex_checked_noDots_QC_b38.eval
> * plink_pca.ALL.wgs.nhgri_coriell_affy_6.20140825.genotypes_has_ped_Updated_withsex_checked_noDots_QC_b38.evec
> * plink_pca.ALL.wgs.nhgri_coriell_affy_6.20140825.genotypes_has_ped_Updated_withsex_checked_noDots_QC_b38.excluded
> * plink_pca.ALL.wgs.nhgri_coriell_affy_6.20140825.genotypes_has_ped_Updated_withsex_checked_noDots_QC_b38_ScreePlot.png
> * plink_pca.ALL.wgs.nhgri_coriell_affy_6.20140825.genotypes_has_ped_Updated_withsex_checked_noDots_QC_b38_variance.txt 

* Then, create Scree and PCA plots in R
#### TODO JAKOB PLZ UPDATE ####
```
library(cowplot); library(tidyverse)

pca <- read_table("plink_pca.ALL.wgs.nhgri_coriell_affy_6.20140825.genotypes_has_ped_Updated_withsex_checked_noDots_QC_b38.evec", col_names = FALSE)
pve <- read_table("plink_pca.ALL.wgs.nhgri_coriell_affy_6.20140825.genotypes_has_ped_Updated_withsex_checked_noDots_QC_b38_variance.txt")
names(pve) <- c("PC","Variance")
ped <- read.csv("20130606_g1k.ped", sep = "\t")

scree_plot <- ggplot(pve, aes(PC, Variance)) + 
  geom_line() + 
  geom_point(color = "black", size = 3) +
  labs(x = "Principal Component", y = "Proportion variance explained") + 
  theme_light() +
  theme(axis.title.x = element_text(size = 15),
        axis.title.y = element_text(size = 14),
        axis.text = element_text(size = 14)) +
  scale_x_continuous(breaks=c(1,5,10,15,20))

pca_one_two <- pca[-1,1:3]
names(pca_one_two) <- c("Individual.ID", "PC1","PC2")
to_plot <- pca_one_two %>% left_join(ped %>% select(Individual.ID, Population), 
                                     by = "Individual.ID")
pc_plot <- to_plot %>% 
  ggplot(aes(PC1, PC2, color = Population)) +
  geom_point() + 
  theme_light() +
  theme(axis.text.y = element_text(angle = 90, size = 13, hjust = 0.5),
        axis.text.x = element_text(size = 13),
        axis.title = element_text(size = 15),
        legend.title = element_text(size = 15),
        legend.text = element_text(size = 12))

# Figure dimensions: 7x7
plot_grid(scree_plot,
          pc_plot, ncol = 1,
          labels = "AUTO", label_size = 17, 
          label_x = -0.005, label_y = 1.02)
```
#### This generates a plot similar to Figure 6 from the paper
![image](https://user-images.githubusercontent.com/66582523/183430931-e91733c5-8fde-4914-a931-945f6f297486.png)

</details>
	
	
### Step 10 -- Calculate frequency files and compare to TOPMed panel
	
<details> 
	<summary>👇 Steps and code </summary>
	<hr>
	
* First, calculate frequencies
	- !!! This needs to be done using Plink 1.9 !!! Check Plink version by typing "plink" and "plink2" on the command line.
```
plink --bfile ALL.wgs.nhgri_coriell_affy_6.20140825.genotypes_has_ped_Updated_withsex_checked_noDots_QC_b38 --freq --out ALL.wgs.nhgri_coriell_affy_6.20140825.genotypes_has_ped_Updated_withsex_checked_noDots_QC_b38_freq
```

* Examine the output
```
head ALL.wgs.nhgri_coriell_affy_6.20140825.genotypes_has_ped_Updated_withsex_checked_noDots_QC_b38_freq.frq
```

> ```
>   CHR          SNP   A1   A2          MAF  NCHROBS
>     1   rs10458597    T    C      0.00733     5730
>     1   rs12565286    C    G      0.03481     5716
>     1   rs12082473    A    G      0.07056     5726
>     1    rs3094315    G    A       0.3353     5732
>     1    rs2286139    C    T       0.4201     5732
>     1   rs11240776    G    A     0.008188     5740
>     1    rs2980319    A    T       0.2963     5740
>     1    rs2980300    T    C       0.4287     5678
>     1    rs2905036    C    T      0.03063     5746 
> ```

* Then, Compare variants to TOPMED panel
	- The perl script is created by the Wayner Tools group: https://www.well.ox.ac.uk/~wrayner/tools/
	- Make sure you've downloaded the following file: HRC-1000G-check-bim-v4.3.0.zip
	- Unzip and add the HRC-1000G-check-bim.pl script to your preImputation directory

```
perl HRC-1000G-check-bim.pl -r PASS.Variants.TOPMed_freeze5_hg38_dbSNP.tab -h -b ALL.wgs.nhgri_coriell_affy_6.20140825.genotypes_has_ped_Updated_withsex_checked_noDots_QC_b38.bim -f ALL.wgs.nhgri_coriell_affy_6.20140825.genotypes_has_ped_Updated_withsex_checked_noDots_QC_b38_freq.frq
```

* Tess's output
> ```
> Displayed
>  Matching to HRC
>  
>  Position Matches
>   ID matches HRC 0
>   ID Doesn't match HRC 832317
>   Total Position Matches 832317
>  ID Match
>   Position different from HRC 0
>  No Match to HRC 10845
>  Skipped (MT) 0
>  Total in bim file 843162
>  Total processed 843162 
> 
>   ID matches HRC 0
>   ID Doesn't match HRC 832317
>   Total Position Matches 832317
>  ID Match
>   Position different from HRC 0
>  No Match to HRC 10845
>  Skipped (MT) 0
>  Total in bim file 843162
>  Total processed 843162
> 
>  Indels 0
> 
>  SNPs not changed 1686
>  SNPs to change ref alt 795615
>  Strand ok 797298
>  Total Strand ok 797301
>  
>  Strand to change 10
>  Total checked 832317
>  Total checked Strand 797308
>  Total removed for allele Frequency diff > 0.2 163426
>  Palindromic SNPs with Freq > 0.4 19737
> 
>  Non Matching alleles 15272
>  ID and allele mismatching 15272; where HRC is . 0
>  Duplicates removed 0
> ```
 
* Run the Run-plink.sh script, which is generate by perl script 
* Pulls out chromosome info and makes VCF files
```
sed -i ‘s/plink/plink2/’ Run-plink.sh
 sed -i 's/--recode vcf/--recode vcf --output-chr chrM/g' Run-plink.sh
 chmod +x ./Run-plink.sh
 ./Run-plink.sh
> ```
</details>
	


<details> 
	<summary>👇 Steps and code </summary>
	<hr>
	
 * Flip files
 
 ```
export BCFTOOLS_PLUGINS=/appl/bcftools-1.9/libexec/bcftools/

for i in {1..22}; do bcftools +fixref ALL.wgs.nhgri_coriell_affy_6.20140825.genotypes_has_ped_Updated_withsex_checked_noDots_QC_b38-updated-chr$i'.vcf' -Ov -o ALL.wgs.nhgri_coriell_affy_6.20140825.genotypes_has_ped_Updated_withsex_checked_noDots_QC_b38-updated_flipped_chr$i'.vcf' -- -d -f ~/group/projects/PMBB/QC_Imputation/scripts/resources_broad_hg38_v0_Homo_sapiens_assembly38.fasta -m flip; done
```
</details>
	
### Step 11 - Sort and zip files to create VCF files for imputation
<details> 
	<summary>👇 Steps and code </summary>
	<hr>	
```
for i in {1..22}; do vcf-sort ALL.wgs.nhgri_coriell_affy_6.20140825.genotypes_has_ped_Updated_withsex_checked_noDots_QC_b38-updated_flipped_chr$i.vcf | bgzip -c > VCFfiles/ALL.wgs.nhgri_coriell_affy_6.20140825_ImputationInput_TOPMED_chr$i.vcf.gz; done
```
</details>
	
	
### Step 12 -- If not already on local computer, copy VCF files to local computer in order to upload to TOPMed impuatation server

<details> 
	<summary>👇 Steps and code </summary>
	<hr>
	
```
scp -r login@serveraddress:/~/GWAS_QC/VCFfiles/ ~/Desktop/
```
</details>
	

### Step 13 -- Last Pre-Impuation step - Calculate relateds (need later, for GWAS)

<details> 
	<summary>👇 Steps and code </summary>
	<hr>
	
```
module load drop_relateds.sh
drop_relateds.sh -b ALL.wgs.nhgri_coriell_affy_6.20140825.genotypes_has_ped_Updated_withsex_checked_noDots_QC_b38 -i ALL.wgs.nhgri_coriell_affy_6.20140825.genotypes_has_ped_pruned10_genome_updated.genome -p remove_related
```
</details>
	

## PART 3 -- Genotype Imputation
### Step 14 -- Imputation using TOPMed Imputation Server
* Imputation has become an essential component of GWAS quality control because it increases power, facilitates meta-analysis, and aids interpretation of signals. Genotype imputation is the statistical inference of unobserved genotype, which enables scientists to reconstruct the missing data in each genome and accurately evaluate the evidence for association at genetic markers that were not genotyped. Genotype imputation is achieved by comparing short stretches of an individual genome against stretches of previously characterized reference genomes.  It is usually performed on single nucleotide polymorphisms (SNPs), which are the most common type of genetic variation. 
* Several tools exist specifically for genotype imputation such as the Michigan and Trans-Omics for Precision Medicine (TOPMed) Imputation Servers where one uploads the phased or unphased GWAS genotypes in order to receive the imputed genomes in return. Each imputation server varies in terms of speed and accuracy. One of the most important considerations in imputation is the composition of the reference panel. For our study, we selected the TOPMed Imputation Reference panel  (version r2) because it is one of the most diverse reference panels available and contains information from 97,256 deeply sequenced human genomes containing 308,107085 genetic variants distributed across the 22 autosomes and the X chromosome. 
* Theoretically, phased means that the two strands on each Chr are separated to identify which regions come from each parent whereas no phasing means that they are not separated. Essentially, for imputation phasing is the first step which is done in reference to the reference genome panel.

<details> 
	<summary>👇 Steps and code </summary>
	<hr>
	
### TODO - NEED TO UPDATE IMAGES & CODE

![image](https://user-images.githubusercontent.com/30478823/154597606-bc2f8b09-2741-493e-9c4a-dabdf238bd23.png)
![image](https://user-images.githubusercontent.com/30478823/154598027-ea78546d-e645-460a-bab7-b41fead62356.png)
![image](https://user-images.githubusercontent.com/30478823/154598503-691c4973-2858-4ff0-a9ac-95360317b405.png)
![image](https://user-images.githubusercontent.com/30478823/154600698-1d443de3-6691-4af9-a078-6bdb8f113e5a.png)
![image](https://user-images.githubusercontent.com/30478823/154738461-b951ab13-75b3-417f-bfc1-1e461dc4cf47.png)
</details>
	
	
### Step 15 -- Download Imputed Data to Working Directory
	
<details> 
	<summary>👇 Steps and code </summary>
	<hr>
	
```
# Download the completed imputation files using the wget commands provided by TOPMed to the location where you'll be working with it

# Location
~group/scratch/van/cphg-gwas-qc-imputed-data

# Unzip the files and enter the password that was emailed to you from TOPMed inside the quotes
module load p7zip
for file in *.zip; do 7z e $file -p"<password>"; done
```
![image](https://user-images.githubusercontent.com/30478823/154745163-97f3cb23-03db-487c-9638-63830eec92cc.png)
</details>

## PART 4 -- Post-Imputation QC
### Step 16 -- 
### TODO - NEED TO UPDATE IMAGES & CODE
	
<details> 
	<summary>👇 Steps and code </summary>
	<hr>
	
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
</details>
	
	
## PART 5 Performing GWAS
### Step 17 -- GWAS with PLINK or Regenie
### TODO - NEED TO UPDATE IMAGES & CODE
   
### PLINK

<details> 
	<summary>👇 Steps and code </summary>
	<hr>
	
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
</details>
	
	
### Regenie

<details> 
	<summary>👇 Steps and code </summary>
	<hr>
	
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
</details> 
	
	
# Related Resources
We recommend the following resources and tutorials developed for performing GWAS:
* Comphrehensive tutorial about GWAS and PRS by MareesAT: https://github.com/MareesAT/GWA_tutorial/
* GWAS Data Cleaning tutorial by the GENEVA Coordinating Center: https://www.bioconductor.org/packages/devel/bioc/vignettes/GWASTools/inst/doc/DataCleaning.pdf
* GWAS QC - theory and steps by the Pan African Bioinformatics Network for H3Africa: https://www.bioinf.wits.ac.za/courses/AIMS/QC_data.pdf 
* The International Sample Genome Resource (IGSR) GitHub: https://github.com/igsr 


# License
* To be decided. Suggestions welcome.
