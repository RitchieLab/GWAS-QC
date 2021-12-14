# GWAS-QC
## Workflow is implemented on 1000Genomes on Build 38 (abbreviated as 1KG here)


## Data locations
* Raw files: /project/ritchie00/datasets/1KG_Phase3/plink_files/plink_raw_files/b38
* Project output diretory: ~group/projects/cphg-gwas-qc
  * Illumina GSA Manifest file (cleaned & converted to bed format): GSA-24v3-0_A2_cleaned.bed

## Browsing raw data
* We go into the directory containing our raw genotype files and inspect to get a sense of them
* For 1KG, we have .bed, .bim, and .fam files
  * *.fam files contain information about samples, one sample per line
  * *.bim files contain information about markers, one marker per line
  * *.bed contain binary genotype information (don’t view this file directly)
* If you're unfamiliar with genomic file types, we will end up needing all of these files because they each contain some form of information we need
* It looks like we have files for chromosome 1-22

## Imputation Preparation (Pre-Imputation QC)
* Data is already lifted-over to build38??
* Split data if its size maxes out the TOPMed Server
* Calculating freq
* Check SNPs 
* Create VCF files aligned to build38 reference alleles
* Sort VCF and zip files
* Run VCF check

```
## Modules we'll need
module load plink
module load bcftools/1.9
export BCFTOOLS_PLUGINS=/appl/bcftools-1.9/libexec/bcftools/
module load vcftools/0.1.12c
module load tabix/0.2.6
```

```
## Preparing files for TOPMed imputation

## Is data lifted over?

## Extract SNPs from GSA Manifest file
in_path=/project/ritchie00/datasets/1KG_Phase3/plink_files/
plink_raw_files/b38/
out_path=~/group/projects/cphg-gwas-qc/

 for i in {1..22}; do \
 plink --make-bed \
 --bfile ${in_path}ALL.chr${i}_GRCh38.genotypes.20170504.genotypes \
 --extract range GSA-24v3-0_A2_cleaned.bed \
 --out ${out_path}prepare/chr${i}_GSA-filtered; done
![image](https://user-images.githubusercontent.com/30478823/146036594-4590d0eb-0753-4bd7-9404-32de174f3c89.png)

## Calculating freq
plink --bfile <input_files> --freq --out <output_filename>

## Checking snps against TOPMed
(http://www.well.ox.ac.uk/~wrayner/tools/)
## Checks Strand, alleles, position, Ref/Alt assignments and frequency differences
## Updates: Strand, position, ref/alt assignment
## Removes: A/T & G/C SNPs if MAF > 0.4, SNPs with differing alleles, SNPs with > 0.2 allele frequency difference, SNPs not in reference panel
## Using HRC-1000G-check-bim-v4.3.0
<use Yuki's perl script?>

<make sure the chromosomes are in the right chr# format>
plink --recode vcf --output-chr

## Creating VCF file aligned with build38 reference alleles (downloaded from: https://console.cloud.google.com/storage/browser/genomics-public-data/resources/broad/hg38/)

#Sorting VCF and zipping files using VCFtools and tabix (make sure the module are loaded first)

## Run VCF check (downloaded from https://github.com/zhanxw/checkVCF)

```
* 


## Imputation using TOPMed Imputation Server
XYZ

```
<code>
```

## Post-Imputation QC
XYZ
```
<code>
```

## GWAS
XYZ

```
<code>
```


