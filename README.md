
<h1 align="center">
  <br>
  <a href="https://github.com/RitchieLab/GWAS-QC"><img src="https://user-images.githubusercontent.com/30478823/182218051-a9a58111-9b92-42e9-bdc1-c2ffad393fe2.png" alt="GWAS QC" width="200"></a>
  <br>
  GWAS Quality Control (QC)
  <br>
</h1>

<h4 align="center">A complementary tutorial to the <a href="https://currentprotocols.onlinelibrary.wiley.com/doi/10.1002/cpz1.603" target="_blank">Truong, Woerner, Cherlin, et al. Paper (2022)</a>.</h4>

<p align="center">
  <a href="#basic-overview">Basic Overview</a> •
  <a href="#set-up">Set Up</a> •
  <a href="#workflow-steps">Workflow Steps</a> •
  <a href="#related-resources">Related Resources</a>
</p>

# To Be Updated
* Upload plink_pca code needs to uploaded 
* ✅ Upload heterogeneity and missingness code 
* Upload scree plot / pca plot code added as an R script that you can run easily
* Upload drop_relateds.sh but we need to adapt the dependencies / paths
* Fix any hardcoded paths
* Add file with simulated phenotypes
* Repeat to clarify specific Plink versions
* Make a note of ideal imputed graph vs the tutorial imputed R-sq graph
* Locate Jakob's .genome file: *~/group/personal/jakob/gwas/cphg/affy/ALL.wgs.nhgri_coriell_affy_6.20140825.genotypes_has_ped_pruned10_genome.genome* 
* ✅ Upload all the steps for .genome file at the end of the tutorial


# Basic Overview
* We chose a publicly available dataset from the International Genome Sample Resource (IGSR) (www.internationalgenome.org). IGSR created and maintains the 1000 Genomes Project (1kGP) to provide a public catalog of common human genetic variation and genotype data. The 1kGP dataset has been kept up to date with current reference data sets, thus it is available for both GRCh37 and GRCh38. The latter is utilized here because the 2014 update increased the quantity of loci represented, resolved more than 1000 issues from the previous version of the assembly, and overall provides a better basis for alignment and subsequent analysis. Additionally, IGSR’s continued efforts will lead to the incorporation of various populations to the data which were not previously captured.

# Set Up
## Set up working directory
<details> 
	<summary>👇 How to set up your directory structure for the GWAS QC workflow </summary>
	<hr>	

```
# use -p option to create multiple directories at once
mkdir -p GWAS_QC/rawData GWAS_QC/preImputation GWAS_QC/preImputation/VCFfiles GWAS_QC/Imputed GWAS_QC/postImpuatation	
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

## Modules in BASH
<details> 
	<summary>👇 How to load modules </summary>
	<hr>
	
```
module load plink/1.9-20210416
module load bcftools/1.9
module load vcftools/0.1.12c # used by bcftools so may not need to be directly loaded
module load tabix/0.2.6 # this one might not work -- let me know if it's a problem use htslib instead (below)
module load htslib
module load perl
module load liftOver/20180423
module load R
```
</details>

## Notes on PLINK v1.9 and v2.0
* HEREIN, ALL USAGE OF PLINK 1.9 IS INDICATED BY "plink" AND PLINK 2.0 BY "plink2"
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
	
	
* Download -- First, we need to download the publicly available dataset from the 1000 Genomes Project (1KGP). The data is Affy6.0 genotype data for 3,450 individuals with population-level data. 

* Files can be found here: http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/supporting/hd_genotype_chip/
* Two download options 
	- Download directly to your local computer by clicking the hyperlink for `ALL.wgs.nhgri_coriell_affy_6.20140825.genotypes_has_ped.vcf.gz` and move it to your working directory
	- OR Dowload using wget command, for example:

* Enter your GWAS_QC directory 
```
cd GWAS_QC/rawData/
```
	
```
wget http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/supporting/hd_genotype_chip/ALL.wgs.nhgri_coriell_affy_6.20140825.genotypes_has_ped.vcf.gz
```

* At this point `ALL.wgs.nhgri_coriell_affy_6.20140825.genotypes_has_ped.vcf.gz` should be in your GWAS_QC/rawData/ directory.

```
gunzip ALL.wgs.nhgri_coriell_affy_6.20140825.genotypes_has_ped.vcf.gz	
```

* Then recode the `.vcf` file as `.bed`, `.bim`, `.bed`. 

```
plink --vcf ALL.wgs.nhgri_coriell_affy_6.20140825.genotypes_has_ped.vcf --make-bed --out ALL.wgs.nhgri_coriell_affy_6.20140825.genotypes_has_ped	
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

```
wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/working/20130606_sample_info/20130606_g1k.ped
```	
	
* You also need to dowload the TOPMed reference panel and code from the McCarthy Tools website
	- https://www.well.ox.ac.uk/~wrayner/tools/ 
	- Midway down the page there are instructions for Usage with the TOPMed reference panel
	- Download `CreateTOPMed.zip` from the McCarthy Tools website
	![image](https://user-images.githubusercontent.com/30478823/192685817-b8e3cc46-045d-4381-a55b-24e9e74bb64d.png)
	- Click on the link to navigate to the Bravo site where `ALL.TOPMed_freeze5_hg38_dbSNP.vcf.gz` can be downloaded from
	- Click on dbSNP in the top right panel and click "Download VCF" button to download `ALL.TOPMed_freeze5_hg38_dbSNP.vcf.gz` file
	![image](https://user-images.githubusercontent.com/30478823/192685864-a94ededb-84b0-4cd3-b921-ad627ff0f463.png)
	- Note: If you run the curl command that's given on the site, the filename will be different 
	- After downloading both files, move them to the GWAS_QC/rawData/ directory, and unzip the  `CreateTOPMed.zip` file
	- Then execute the following command to convert VCF to an HRC formatted reference legend according to the code provided by McCarthy Tools:
	```
	./CreateTOPMed.pl -i ALL.TOPMed_freeze5_hg38_dbSNP.vcf.gz
	```
	- If the shabang "./" doesn't work, you may need to run the command with "perl" instead
	- By default this will create a file filtered for variants flagged as PASS only, if you wish to use all variants the -a flag overrides this. To override the default output file naming use -o filename.
	- If you get an error in the above step, try this variation. Both were run successfully on a local computer and server using perl/5.30.0. Depending on your setup, this may take a few hours to run:
	```
	perl CreateTOPMed.pl -i ALL.TOPMed_freeze5_hg38_dbSNP.vcf.gz	
	```
	- Since it takes awhile to run, you may want to submit a job to the queue if you're working on a shared cluster
</details>

* This is how the file system should look at the end of Step 2
```
.
├── GWAS_QC/
    ├── rawData/ (current directory)
    |    ├── 20130606_g1k.ped
    |    ├── ALL.TOPMed_freeze5_hg38_dbSNP.vcf.gz
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
### Step 1 - Enter the preImputation directory

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

To calculate heterogeneity, you can use the --het command in Plink. The command would look like this: *plink --bfile [file_name] --het*

This will calculate the heterozygosity for each individual in the dataset and write the results to a file called "[file_name].het".
	
```
module load plink/1.9
plink --bfile ../rawData/ALL.wgs.nhgri_coriell_affy_6.20140825.genotypes_has_ped --het --out ALL.wgs.nhgri_coriell_affy_6.20140825.genotypes_has_ped_het	
```

* Output:
> ```
> PLINK v1.90p 64-bit (16 Apr 2021)              www.cog-genomics.org/plink/1.9/
> (C) 2005-2021 Shaun Purcell, Christopher Chang   GNU General Public License v3
> Logging to ALL.wgs.nhgri_coriell_affy_6.20140825.genotypes_has_ped_het.log.
> Options in effect:
>   --bfile ../rawData/ALL.wgs.nhgri_coriell_affy_6.20140825.genotypes_has_ped
>   --het
>   --out ALL.wgs.nhgri_coriell_affy_6.20140825.genotypes_has_ped_het
> 
> 128235 MB RAM detected; reserving 64117 MB for main workspace.
> 905788 variants loaded from .bim file.
> 3450 people (0 males, 0 females, 3450 ambiguous) loaded from .fam.
> Ambiguous sex IDs written to
> ALL.wgs.nhgri_coriell_affy_6.20140825.genotypes_has_ped_het.nosex .
> Using 1 thread (no multithreaded calculations invoked).
> Before main variant filters, 3450 founders and 0 nonfounders present.
> Calculating allele frequencies... done.
> Warning: Nonmissing nonmale Y chromosome genotype(s) present; many commands
> treat these as missing.
> Total genotyping rate is 0.996712.
> 905788 variants and 3450 people pass filters and QC.
> Note: No phenotypes present.
> --het: 868601 variants scanned, report written to
> ALL.wgs.nhgri_coriell_affy_6.20140825.genotypes_has_ped_het.het .
> ```

To calculate missingness, you can use the --missing command in Plink. The command would look like this: *plink --bfile [file_name] --missing* 
	
This will calculate the proportion of missing genotypes for each individual in the dataset and write the results to a file called "[file_name].lmiss"

You can also use the --missing2 command which will report the frequency of missing genotypes at each variant, the command is: *plink --bfile [file_name] --missing2*

This will give a file with extension .lmiss2 which will have frequency of missing genotypes for each variant	


```	
plink --bfile ../rawData/ALL.wgs.nhgri_coriell_affy_6.20140825.genotypes_has_ped --missing --out ALL.wgs.nhgri_coriell_affy_6.20140825.genotypes_has_ped_miss	
```

* Output:
> ```
> PLINK v1.90p 64-bit (16 Apr 2021)              www.cog-genomics.org/plink/1.9/
> (C) 2005-2021 Shaun Purcell, Christopher Chang   GNU General Public License v3
> Logging to ALL.wgs.nhgri_coriell_affy_6.20140825.genotypes_has_ped_miss.log.
> Options in effect:
>   --bfile ../rawData/ALL.wgs.nhgri_coriell_affy_6.20140825.genotypes_has_ped
>   --missing
>   --out ALL.wgs.nhgri_coriell_affy_6.20140825.genotypes_has_ped_miss
> 
> 128235 MB RAM detected; reserving 64117 MB for main workspace.
> 905788 variants loaded from .bim file.
> 3450 people (0 males, 0 females, 3450 ambiguous) loaded from .fam.
> Ambiguous sex IDs written to
> ALL.wgs.nhgri_coriell_affy_6.20140825.genotypes_has_ped_miss.nosex .
> Using 1 thread (no multithreaded calculations invoked).
> Before main variant filters, 3450 founders and 0 nonfounders present.
> Calculating allele frequencies... done.
> Warning: Nonmissing nonmale Y chromosome genotype(s) present; many commands
> treat these as missing.
> Total genotyping rate is 0.996712.
> --missing: Sample missing data report written to
> ALL.wgs.nhgri_coriell_affy_6.20140825.genotypes_has_ped_miss.imiss, and
> variant-based missing data report written to
> ALL.wgs.nhgri_coriell_affy_6.20140825.genotypes_has_ped_miss.lmiss.
> ```
	
* Here's the current directory structure within GWAS_QC/preImputation/	
> ALL.wgs.nhgri_coriell_affy_6.20140825.genotypes_has_ped_het.het
> ALL.wgs.nhgri_coriell_affy_6.20140825.genotypes_has_ped_het.log
> ALL.wgs.nhgri_coriell_affy_6.20140825.genotypes_has_ped_miss.log
> ALL.wgs.nhgri_coriell_affy_6.20140825.genotypes_has_ped_miss.imiss
> ALL.wgs.nhgri_coriell_affy_6.20140825.genotypes_has_ped_miss.lmiss
```
.
├── GWAS_QC/
    ├── rawData/
    |    ├── 20130606_g1k.ped
    |    ├── ALL.TOPMed_freeze5_hg38_dbSNP.vcf.gz
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
    |    ├── ALL.wgs.nhgri_coriell_affy_6.20140825.genotypes_has_ped_het.nosex
    |    ├── ALL.wgs.nhgri_coriell_affy_6.20140825.genotypes_has_ped_miss.imiss
    |    ├── ALL.wgs.nhgri_coriell_affy_6.20140825.genotypes_has_ped_miss.lmiss
    |    ├── ALL.wgs.nhgri_coriell_affy_6.20140825.genotypes_has_ped_miss.log
    |    ├── ALL.wgs.nhgri_coriell_affy_6.20140825.genotypes_has_ped_miss.nosex
    |    └── sex_file.txt
    ├──Imputed/
    └── postImputation/	
```

	

* Plot heterogeneity vs. missingness
* Needed: 
	- Code_Heterogeneity_Missingness.R
	- ALL.wgs.nhgri_coriell_affy_6.20140825.genotypes_has_ped_miss.imiss 
	- ALL.wgs.nhgri_coriell_affy_6.20140825.genotypes_has_ped_het.het
* Input file: ALL.wgs.nhgri_coriell_affy_6.20140825.genotypes_has_ped
* Output file: ALL.wgs.nhgri_coriell_affy_6.20140825.genotypes_has_ped_HM_plot.pdf

```
Rscript Code_Heterogeneity_Missingness.R ALL.wgs.nhgri_coriell_affy_6.20140825.genotypes_has_ped ALL.wgs.nhgri_coriell_affy_6.20140825.genotypes_has_ped_HM_plot.pdf
```
	
<img width="686" alt="Screen Shot 2022-11-03 at 1 45 56 PM" src="https://user-images.githubusercontent.com/66582523/199796370-91260dee-bedb-4bd1-b65d-f6244c08df1f.png">

* Figure 2 uses the same code as above, on example internal data
</details>
	

### Step 3 - Update first column of file which has zeros

<details> 
	<summary>👇 Steps and code </summary>
	<hr>
	
```
cat ../rawData/ALL.wgs.nhgri_coriell_affy_6.20140825.genotypes_has_ped.fam | awk '{print $1,$2,$2,$2}' > ALL.wgs.nhgri_coriell_affy_6.20140825.genotypes_has_ped_toUpdate.txt

plink --bfile ../rawData/ALL.wgs.nhgri_coriell_affy_6.20140825.genotypes_has_ped --update-ids ALL.wgs.nhgri_coriell_affy_6.20140825.genotypes_has_ped_toUpdate.txt --make-bed --out ALL.wgs.nhgri_coriell_affy_6.20140825.genotypes_has_ped_Updated
```
> ```
PLINK v1.90p 64-bit (16 Apr 2021)
Options in effect:
  --bfile ../rawData/ALL.wgs.nhgri_coriell_affy_6.20140825.genotypes_has_ped
  --make-bed
  --out ALL.wgs.nhgri_coriell_affy_6.20140825.genotypes_has_ped_Updated
  --update-ids ALL.wgs.nhgri_coriell_affy_6.20140825.genotypes_has_ped_toUpdate.txt
Hostname: superman
Working directory: /project/ritchie07/personal/ikuzwa/PROJECTS/GWAS_QC/preImputation2.0
Start time: Mon Jan 23 13:55:46 2023
Random number seed: 1674500146
128235 MB RAM detected; reserving 64117 MB for main workspace.
905788 variants loaded from .bim file.
3450 people (0 males, 0 females, 3450 ambiguous) loaded from .fam.
Ambiguous sex IDs written to
ALL.wgs.nhgri_coriell_affy_6.20140825.genotypes_has_ped_Updated.nosex .
--update-ids: 3450 people updated.
Using 1 thread (no multithreaded calculations invoked).
Before main variant filters, 3450 founders and 0 nonfounders present.
Calculating allele frequencies... done.
Warning: Nonmissing nonmale Y chromosome genotype(s) present; many commands
treat these as missing.
Total genotyping rate is 0.996712.
905788 variants and 3450 people pass filters and QC.
Note: No phenotypes present.
--make-bed to
ALL.wgs.nhgri_coriell_affy_6.20140825.genotypes_has_ped_Updated.bed +
ALL.wgs.nhgri_coriell_affy_6.20140825.genotypes_has_ped_Updated.bim +
ALL.wgs.nhgri_coriell_affy_6.20140825.genotypes_has_ped_Updated.fam ... done.
> ```
	
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
> ```
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
> ```	
	
* Next, check sex and remove sex inconsistencies
```
plink --bfile ALL.wgs.nhgri_coriell_affy_6.20140825.genotypes_has_ped_Updated_withsex --check-sex 
```

> ```
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
> ```
	
```
cat plink.sexcheck |  grep PROBLEM | sed 's/^ *//' > plink.sexcheck_list
```

```
plink2 --bfile ALL.wgs.nhgri_coriell_affy_6.20140825.genotypes_has_ped_Updated_withsex --remove plink.sexcheck_list --make-bed --out ALL.wgs.nhgri_coriell_affy_6.20140825.genotypes_has_ped_Updated_withsex_checked
```
> ```
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
> ```

</details>
	
### Step 5 -- Remove SNP variants that do not have SNP IDs
	
<details> 
	<summary>👇 Steps and code </summary>
	<hr>
	
```
echo . > noSNP
plink2 --bfile ALL.wgs.nhgri_coriell_affy_6.20140825.genotypes_has_ped_Updated_withsex_checked --exclude noSNP --make-bed --out ALL.wgs.nhgri_coriell_affy_6.20140825.genotypes_has_ped_Updated_withsex_checked_noDots
```
> ```
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
> ```
	
</details>
	

### Step 6 -- Run QC on data
	
<details> 
	<summary>👇 Steps and code </summary>
	<hr>
	
```
plink2 --bfile ALL.wgs.nhgri_coriell_affy_6.20140825.genotypes_has_ped_Updated_withsex_checked_noDots --geno 0.05 --mind 0.1 --make-bed --out ALL.wgs.nhgri_coriell_affy_6.20140825.genotypes_has_ped_Updated_withsex_checked_noDots_QC
```
> ```
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
> ```

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
cat ALL.wgs.nhgri_coriell_affy_6.20140825.genotypes_has_ped_Updated_withsex_checked_noDots_QC.bim | awk '{print "chr"$1, $4, ($4+1), $4, $2}' > liftover_input.bed
```

* Check output
```
head liftover_input.bed
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
sed -i 's/chr23/chrX/g' liftover_input.bed #I don't think this is a problem with this dataset
```

* Then, download the correct liftover file
	- http://hgdownload.cse.ucsc.edu/goldenpath/hg19/liftOver/
	- dowload `hg19ToHg38.over.chain.gz` to your local computer
	- If doing analysis on server, transfer file to working directory using scp
	- unzip the data
	
* Transfer file to working directory on server
```
scp /home/directory/path/hg19ToHg38.over.chain.gz /server/directory/path/
```	

* Unzip file
```
gunzip 	hg19ToHg38.over.chain.gz
```

* Finally, perform the actual liftOver
```
liftOver liftover_input.bed hg19ToHg38.over.chain liftover_newmap.txt liftover_exclude.txt
```
* You will end up with two output files `liftover_newmap.txt` and `liftover_exclude.txt`
* `liftover_newmap.txt` will contain the coordinates of each SNP in assembly `GRCh38` after the liftover in column 2
```
head liftover_newmap.txt
```
```
chr1    629241  629242  564621  rs10458597
chr1    785910  785911  721290  rs12565286
chr1    805477  805478  740857  rs12082473
chr1    817186  817187  752566  rs3094315
chr1    826352  826353  761732  rs2286139
chr1    829889  829890  765269  rs11240776
chr1    841742  841743  777122  rs2980319
chr1    850609  850610  785989  rs2980300
chr1    857100  857101  792480  rs2905036
chr1    863579  863580  798959  rs11240777
```

* To double check that this worked correctly take the first SNP `rs10458597` and look it up in dbSNP https://www.ncbi.nlm.nih.gov/snp/
	
	
<img width="808" alt="Screen Shot 2023-01-31 at 3 13 02 PM" src="https://user-images.githubusercontent.com/66582523/215872229-be837816-fb34-4b3b-a1b3-7559b9474ce2.png">

* You should see that SNP `rs10458597` is in position `629241` in assembly `GRCh38` but position `564621` in assembly `GRCh37`
* `liftover_exclude.txt` will contain the SNPs that were not able to be aligned during the liftOver process
```
head exclude_liftover.txt
```
```
rs613882
rs615185
rs17457350
rs17423513
rs10492935
rs17423547
rs10915404
rs557477
rs560335
rs565941
```	
	
* Some formatting is required for both files before moving on to Step 8	
```
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
	
> ```
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
> ```

* Final pre-imputation variant count is 834,872 SNPs
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
	
	
* The command should have made an output of 5 files:
```
ls plink_pca.ALL.wgs.nhgri_coriell_affy_6.20140825.geno*
```

> * plink_pca.ALL.wgs.nhgri_coriell_affy_6.20140825.genotypes_has_ped_Updated_withsex_checked_noDots_QC_b38.eval
> * plink_pca.ALL.wgs.nhgri_coriell_affy_6.20140825.genotypes_has_ped_Updated_withsex_checked_noDots_QC_b38.evec
> * plink_pca.ALL.wgs.nhgri_coriell_affy_6.20140825.genotypes_has_ped_Updated_withsex_checked_noDots_QC_b38.excluded
> * plink_pca.ALL.wgs.nhgri_coriell_affy_6.20140825.genotypes_has_ped_Updated_withsex_checked_noDots_QC_b38_ScreePlot.png
> * plink_pca.ALL.wgs.nhgri_coriell_affy_6.20140825.genotypes_has_ped_Updated_withsex_checked_noDots_QC_b38_variance.txt 

* Then, create Scree and PCA plots in R
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
names(pca_one_two)  <- c("Individual.ID", "PC1","PC2")

pc_plot <- pca_one_two %>% 
  left_join(ped %>% 
              select(Individual.ID, Population), 
            by = "Individual.ID") %>% 
  # Collapse to Superpopulations
  mutate(Population = fct_collapse(Population, 
                                   AFR = c("YRI","LWK","GWD","MSL","ESN"),
                                   AMR = c("ASW","ACB","MXL","PUR","CLM","PEL"),
                                   EAS = c("CHB","JPT","CHS","CDX","KHV"), #,"CHD"), # Denver Chinese not in our data set
                                   EUR = c("CEU","TSI","GBR","FIN","IBS"),
                                   SAS = c("GIH","PJL","BEB","STU","ITU")),
         Population = factor(Population, levels = c("AFR","AMR","EAS","EUR","SAS"))) %>% 
  ggplot(aes(PC1, PC2, color = Population)) +
  geom_point(alpha = 0.8) + 
  theme_light() +
  theme(axis.text.y = element_text(angle = 90, size = 13, hjust = 0.5),
        axis.text.x = element_text(size = 13),
        axis.title = element_text(size = 15),
        legend.title = element_text(size = 15),
        legend.text = element_text(size = 12)) +
  scale_color_manual(values = palette.colors(palette = "set1"))

# Figure dimensions: 7x7
plot_grid(scree_plot,
          pc_plot, ncol = 1,
          labels = "AUTO", label_size = 17, 
          label_x = -0.005, label_y = 1.02)
```
#### This generates a plot similar to Figure 6 from the paper
![Figure_6](https://user-images.githubusercontent.com/86791191/195680420-c6a84bb3-d02f-4b53-90bd-4cfc5b22ddb5.png)

</details>
	
	
### Step 10 -- Calculate frequency files and compare to TOPMed panel
	
<details> 
	<summary>👇 Steps and code </summary>
	<hr>
	
* First, calculate frequencies
	- Note: This needs to be done using Plink 1.9. Check Plink version by typing "plink" and "plink2" on the command line.
```
plink --bfile ALL.wgs.nhgri_coriell_affy_6.20140825.genotypes_has_ped_Updated_withsex_checked_noDots_QC_b38 --freq --out ALL.wgs.nhgri_coriell_affy_6.20140825.genotypes_has_ped_Updated_withsex_checked_noDots_QC_b38_freq
```

> ```
> PLINK v1.90p 64-bit (16 Apr 2021)              www.cog-genomics.org/plink/1.9/
> (C) 2005-2021 Shaun Purcell, Christopher Chang   GNU General Public License v3
> Logging to ALL.wgs.nhgri_coriell_affy_6.20140825.genotypes_has_ped_Updated_withsex_checked_noDots_QC_b38_freq.log.
> Options in effect:
>   --bfile ALL.wgs.nhgri_coriell_affy_6.20140825.genotypes_has_ped_Updated_withsex_checked_noDots_QC_b38
>   --freq
>   --out ALL.wgs.nhgri_coriell_affy_6.20140825.genotypes_has_ped_Updated_withsex_checked_noDots_QC_b38_freq
> 
> 128235 MB RAM detected; reserving 64117 MB for main workspace.
> 843232 variants loaded from .bim file.
> 2873 people (1713 males, 1160 females) loaded from .fam.
> Using 1 thread (no multithreaded calculations invoked).
> Before main variant filters, 2873 founders and 0 nonfounders present.
> Calculating allele frequencies... done.
> Total genotyping rate is 0.997609.
> --freq: Allele frequencies (founders only) written to
> ALL.wgs.nhgri_coriell_affy_6.20140825.genotypes_has_ped_Updated_withsex_checked_noDots_QC_b38_freq.frq
> .

	
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

* Then, compare variants to TOPMED panel
	- The perl script is created by the Wayner Tools group: https://www.well.ox.ac.uk/~wrayner/tools/
	- Make sure you've downloaded the following file: `HRC-1000G-check-bim-v4.3.0.zip`
	- Unzip and add the `HRC-1000G-check-bim.pl` script to your rawData/ directory
	- Remember that we already have `PASS.Variants.TOPMed_freeze5_hg38_dbSNP.tab` in rawData/ as well

```
perl ../rawData/HRC-1000G-check-bim.pl -r ../rawData/PASS.Variants.TOPMed_freeze5_hg38_dbSNP.tab -h -b ALL.wgs.nhgri_coriell_affy_6.20140825.genotypes_has_ped_Updated_withsex_checked_noDots_QC_b38.bim -f ALL.wgs.nhgri_coriell_affy_6.20140825.genotypes_has_ped_Updated_withsex_checked_noDots_QC_b38_freq.frq
```

* Tess's output
> ```
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
 	
	
* Van's output
> ```
> Matching to HRC
> 
> Position Matches
>  ID matches HRC 0
>  ID Doesn't match HRC 139530
>  Total Position Matches 139530
> ID Match
>  Position different from HRC 0
> No Match to HRC 703702
> Skipped (MT) 0
> Total in bim file 843232
> Total processed 843232
> 
> Indels 0
> 
> SNPs not changed 15188
> SNPs to change ref alt 37739
> Strand ok 38265
> Total Strand ok 52927
> 
> Strand to change 28960
> Total checked 139530
> Total checked Strand 67225
> Total removed for allele Frequency diff > 0.2 45531
> Palindromic SNPs with Freq > 0.4 291
> 
> Non Matching alleles 72014
> ID and allele mismatching 72014; where HRC is . 0
> Duplicates removed 0
> ```
 
* Run the Run-plink.sh script, which is generate by perl script 
	- Pulls out chromosome info and makes VCF files
```
sed -i 's/--recode vcf/--recode vcf --output-chr chrM/g' Run-plink.sh
chmod +x ./Run-plink.sh
./Run-plink.sh
```

* Example output:
	
> ```
> Start time: Tue Oct 11 17:18:24 2022
> 128235 MiB RAM detected; reserving 64117 MiB for main workspace.
> Using up to 56 threads (change this with --threads).
> 2873 samples (1160 females, 1713 males; 2873 founders) loaded from
> /project/ritchie02/projects/cphg-gwas-qc/GWAS_QC/preImputation/ALL.wgs.nhgri_coriell_affy_6.20140825.genotypes_has_ped_Updated_withsex_checked_noDots_QC_b38-updated.fam.
> 236 out of 21694 variants loaded from
> /project/ritchie02/projects/cphg-gwas-qc/GWAS_QC/preImputation/ALL.wgs.nhgri_coriell_affy_6.20140825.genotypes_has_ped_Updated_withsex_checked_noDots_QC_b38-updated.bim.
> Note: No phenotype data present.
> Writing
> /project/ritchie02/projects/cphg-gwas-qc/GWAS_QC/preImputation/ALL.wgs.nhgri_coriell_affy_6.20140825.genotypes_has_ped_Updated_withsex_checked_noDots_QC_b38-updated-chr22.fam
> ... done.
> Writing
> /project/ritchie02/projects/cphg-gwas-qc/GWAS_QC/preImputation/ALL.wgs.nhgri_coriell_affy_6.20140825.genotypes_has_ped_Updated_withsex_checked_noDots_QC_b38-updated-chr22.bim
> ... done.
> Writing
> /project/ritchie02/projects/cphg-gwas-qc/GWAS_QC/preImputation/ALL.wgs.nhgri_coriell_affy_6.20140825.genotypes_has_ped_Updated_withsex_checked_noDots_QC_b38-updated-chr22.bed
> ... done.
> End time: Tue Oct 11 17:18:24 2022
> PLINK v2.00a3LM 64-bit Intel (5 May 2021)      www.cog-genomics.org/plink/2.0/
> (C) 2005-2021 Shaun Purcell, Christopher Chang   GNU General Public License v3
> Error: Duplicate --output-chr flag.
> ```

</details>
	

	
### Step 11 - Sort and zip files to create VCF files for imputation
<details> 
	<summary>👇 Steps and code </summary>
	<hr>

* Flip files
 ```
export BCFTOOLS_PLUGINS=/appl/bcftools-1.9/libexec/bcftools/

for i in {1..22}; do bcftools +fixref ALL.wgs.nhgri_coriell_affy_6.20140825.genotypes_has_ped_Updated_withsex_checked_noDots_QC_b38-updated-chr$i'.vcf' -Ov -o ALL.wgs.nhgri_coriell_affy_6.20140825.genotypes_has_ped_Updated_withsex_checked_noDots_QC_b38-updated_flipped_chr$i'.vcf' -- -d -f ~/group/projects/PMBB/QC_Imputation/scripts/resources_broad_hg38_v0_Homo_sapiens_assembly38.fasta -m flip; done
```

* Example Output:
> ```
> ...repeated multiple times...
> # SC, guessed strand convention
> SC      TOP-compatible  0
> SC      BOT-compatible  0
> # ST, substitution types
> ST      A>C     2       0.8%
> ST      A>G     28      11.9%
> ST      A>T     1       0.4%
> ST      C>A     3       1.3%
> ST      C>G     3       1.3%
> ST      C>T     81      34.3%
> ST      G>A     73      30.9%
> ST      G>C     3       1.3%
> ST      G>T     5       2.1%
> ST      T>A     0       0.0%
> ST      T>C     35      14.8%
> ST      T>G     2       0.8%
> # NS, Number of sites:
> NS      total           236
> NS      ref match       236     100.0%
> NS      ref mismatch    0       0.0%
> NS      flipped         0       0.0%
> NS      swapped         0       0.0%
> NS      flip+swap       0       0.0%
> NS      unresolved      7       3.0%
> NS      fixed pos       0       0.0%
> NS      skipped         0
> NS      non-ACGT        0
> NS      non-SNP         0
> NS      non-biallelic   0
> ```
	
	
* Sort VCF and zip 
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
	

### Step 13 -- Last Pre-Imputation step to calculate relateds which will be needed later for GWAS

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
* Several tools exist specifically for genotype imputation such as the Michigan and Trans-Omics for Precision Medicine (TOPMed) Imputation Servers where one uploads the phased or unphased GWAS genotypes in order to receive the imputed genomes in return. Each imputation server varies in terms of speed and accuracy. One of the most important considerations in imputation is the composition of the reference panel. For our study, we selected the TOPMed Imputation Reference panel  (version r2) because it is one of the most diverse reference panels available and contains information from 97,256 deeply sequenced human genomes containing 308,107085 genetic variants distributed across the 22 autosomes and the X chromosome. 

<details> 
	<summary>👇 Steps for form entry </summary>
	<hr>

![image](https://user-images.githubusercontent.com/30478823/195203144-bf28dd0a-78e7-41bb-a6f6-19729ff5f48d.png)


</details>
		
<details> 
	<summary>👇 TOPMed Server Output </summary>
	<hr>	

* Click on the "Jobs" tab at the top to view your current jobs
* Select the relevant job run to view your results and download the completed data
* The data will expire within a certain timeframe
* There will be sections for: Input Validation, Quality Control, Quality Control (Report), Pre-phasing and Imputation, and Data Compression and Encryption

![image](https://user-images.githubusercontent.com/30478823/197097052-5936e080-732d-4678-b6a4-f15058c2a029.png)

![image](https://user-images.githubusercontent.com/30478823/197097105-cc631c01-7931-4b66-8c07-9fd4a2267b0b.png)

</details>
	
	
### Step 15 -- Download Imputed Data to Working Directory
	
<details> 
	<summary>👇 Steps and code </summary>
	<hr>
	
* Download the completed imputation files using the wget commands provided by TOPMed to the location where you'll be working with it
	
* Change to the appropriate location
```
cd GWAS_QC/Imputed/
```
	
* Copy imputed files to your working directory from the TOPMed Imputation Server. The output will be a curl command that you can copy & paste to run in your active directory. Below is an example:
```
curl -sL <URL> | bash	
```
	
* Unzip the files and enter the password that was emailed to you from TOPMed inside the quotes module load p7zip
```
for file in *.zip; do 7z e $file -p"<password>"; done
```

</details>

	
## PART 4 -- Post-Imputation QC
### Step 16 -- Merge Chromosomes
	
	
<details> 
	<summary>👇 Steps and code </summary>
	<hr>

	
* After unzipping we need to merge the files to do downstream analsysis 
	
* First make mergelist.txt file
	
```	
seq 2 22 | sed 's/^/chr/' > mergelist.txt
```

	
* NEED PLINK1.9 for merge
```
module load plink/1.9
```
* Merge chromosome files into 1 bim/bed/bam file
**THIS TAKES A LONG TIME**
```
plink --bfile chr1 --merge-list mergelist.txt --make-bed --out merged	
```

</details>
	
### Step 17 -- Run pre-QC Commands
	
	
<details> 
	<summary>👇 Steps and code </summary>
	<hr>	

NEED PLINK 2 for these QC Steps	
#Run QC Commands

	
*  First update the .fam files to have the correct ID orientation	
```
head merged.fam	
```
```	
0 HG00096_HG00096 0 0 0 -9
0 HG00098_HG00098 0 0 0 -9
0 HG00101_HG00101 0 0 0 -9
0 HG00102_HG00102 0 0 0 -9
0 HG00105_HG00105 0 0 0 -9
0 HG00106_HG00106 0 0 0 -9
0 HG00107_HG00107 0 0 0 -9
0 HG00108_HG00108 0 0 0 -9
0 HG00109_HG00109 0 0 0 -9
0 HG00111_HG00111 0 0 0 -9
```

* This code will reorient that .fam files to be usable for plink tools downstream analysis
* Col 1 = FID
* Col 2 = IID
* Col 3 = XXX
* Col 4 = XXX
* Col 5 = XXX
* Col 6 = XXX
```
cat merged.fam | awk '{print $1,$2}' | sed s'/_/ /g' | awk '{print $1,$2"_"$2,$2,$3}' > merged_toUpdate.txt
plink2 --bfile merged --update-ids merged_toUpdate.txt --make-bed --out merged_Updated
head merged_Updated.fam
```
```
HG00096	HG00096	0	0	0	-9
HG00098	HG00098	0	0	0	-9
HG00101	HG00101	0	0	0	-9
HG00102	HG00102	0	0	0	-9
HG00105	HG00105	0	0	0	-9
HG00106	HG00106	0	0	0	-9
HG00107	HG00107	0	0	0	-9
HG00108	HG00108	0	0	0	-9
HG00109	HG00109	0	0	0	-9
HG00111	HG00111	0	0	0	-9
```
```
ls 
```
* Should have the following files:
```
merged_Updated.bed
merged_Updated.bim
merged_Updated.fam
merged_Updated.log
```
* If you don't, open the merged_Updated.log file and try to discern the error. It's possible you ran the code in the wrong directory.
</details>	

### Step 18 -- Run QC Commands
	
	
<details> 
	<summary>👇 Steps and code </summary>
	<hr>	

* Want to QC on 
	- geno 0.01 (geno filters out all variants with missing call rates exceeding the provided value to be removed)
```	
plink2 --bfile merged_Updated --geno 0.01 --make-bed --out merged_Updated_1_geno
ls merged_Updated_1_geno.*
```
```
merged_Updated_1_geno.bed
merged_Updated_1_geno.bim
merged_Updated_1_geno.fam
merged_Updated_1_geno.log
```
* Want to QC on 
	- mind 0.01 (filters out al lthe samples with missing call rates exceeding the provided value to be removed)
	- maf 0.05 (filters out all variants with minor allele frequency below the provided threshold)
	- hardy (writes a list of genotype counts and Hardy-Weinberg equilibrium exact test statistics to plink.hwe) 
```
plink2 --bfile merged_Updated_1_geno --mind 0.01 --maf 0.05 --hwe 0.05 --make-bed --out merged_Updated_2_QC
ls merged_Updated2_QC.*
```
```
merged_Updated_2_QC.bed
merged_Updated_2_QC.bim
merged_Updated_2_QC.fam
merged_Updated_2_QC.log
```
	
</details>	
	

### Step 19 --Remove Relateds
	
<details> 
	<summary>👇 Steps and code </summary>
	<hr>
	
```	
# Get sex inconsistent samples (from raw)
plink --bfile ~/group/personal/tess/GWAS_Tutorial/ALL.wgs.nhgri_coriell_affy_6.20140825.genotypes_has_ped --check-sex --out prefix_checksex

awk '{ if ($4 == 0) { print $2} }' prefix_checksex.sexcheck > nosex.txt

# Remove sex inconsistent samples (from raw)
plink --bfile ~/group/personal/tess/GWAS_Tutorial/ALL.wgs.nhgri_coriell_affy_6.20140825.genotypes_has_ped --remove nosex.txt --make-bed --out prefix_nosex

# Call rate (95% snp call rate / 90% sample call rate)
plink --bfile prefix_nosex --geno 0.05 --mind 0.1 --make-bed --out prefix_call

# MAF
plink --bfile prefix_call --maf 0.05 --make-bed --out prefix_maf

# Pruning down to 99,999
plink --bfile prefix_maf --indep-pairwise 50 5 0.12831 --out prune_100k

# Keep those samples (Around 97k are autosomal)
plink --bfile prefix_maf --extract prune_100k.prune.in --make-bed --out prefix_pruned

# IBD
plink --bfile prefix_pruned --genome --out prefix_pruned_100k_genome


cat  ~/group/personal/jakob/gwas/cphg/affy/ALL.wgs.nhgri_coriell_affy_6.20140825.genotypes_has_ped_pruned10_genome.genome | awk '{print $2,$2,$4,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14}' > ALL.wgs.nhgri_coriell_affy_6.20140825.genotypes_has_ped_pruned10_genome_updated.genome
```
	
* Manually change the header so it is FID1 IID1 FID2 IID2
	- Steps
	- i to initiate insert mode
	- Delete headers 1-4
	- Replace with FID1, IID1, FID2, IID2 (Space separated)
	- esc to exit from the insert mode
	- ZZ to close and save the file
```
vi ALL.wgs.nhgri_coriell_affy_6.20140825.genotypes_has_ped_pruned10_genome_updated.genome
```

* Actually remove the relateds
```
drop_relateds.sh -b merged_Updated_2_QC -i ALL.wgs.nhgri_coriell_affy_6.20140825.genotypes_has_ped_pruned10_genome_updated.genome -p merged_Updated_2_QC_remove_related
```
</details>
	
	
# GWAS-Related Resources
We recommend the following resources and tutorials developed for performing GWAS. Due to time constraints, we were not able to fully vet every resource listed below.
* Genome-wide association studies review paper: https://www.nature.com/articles/s43586-021-00056-9
	- Table 1 outlines open access tools that can be applied at each stage of GWAS
* Methods and Tools in Genome-wide Association Studies: https://link.springer.com/protocol/10.1007/978-1-4939-8618-7_5
* Comphrehensive tutorial about GWAS and PRS by MareesAT: https://github.com/MareesAT/GWA_tutorial/
* GWAS Data Cleaning tutorial by the GENEVA Coordinating Center: https://www.bioconductor.org/packages/devel/bioc/vignettes/GWASTools/inst/doc/DataCleaning.pdf
* GWAS QC - theory and steps by the Pan African Bioinformatics Network for H3Africa: https://www.bioinf.wits.ac.za/courses/AIMS/QC_data.pdf 
* The International Sample Genome Resource (IGSR) GitHub: https://github.com/igsr 

### PLINK
* PLINK Documentation for Association Analysis: https://www.cog-genomics.org/plink/1.9/assoc
* Using PLINK for Genome-Wide Association Study: https://lsl.sinica.edu.tw/Activities/class/files/20210506821.pdf
* Plink Tutorial from BIOS25328 Cancer Genomics Class: https://bios25328.hakyimlab.org/post/2021/04/09/plink-tutorial/
* GWAS and PLINK Practical: https://ibg.colorado.edu/cdrom2019/colodro_grasby/GWAS_QC_part2/GWAS_QC_part2_practical.pdf
### SAIGE
* SAIGE GWAS Walkthrough: https://documentation.dnanexus.com/science/scientific-guides/saige-gwas-walkthrough
* SAIGE Github: https://github.com/weizhouUMICH/SAIGE
GWAS in large-scale biobanks and cohorts: https://www.colorado.edu/ibg/sites/default/files/attached-files/boulderworkshop-saige-part2.pdf

### Regenie
* Regenie Documentation: https://rgcgithub.github.io/regenie/options/
* Glow tutorial which implements a distribute version of the Regenie method: https://glow.readthedocs.io/en/latest/tutorial.html
* DNAnexus GWAS on the Research Analysis Platform using regenie: https://www.youtube.com/watch?v=762PVlyZJ-U

