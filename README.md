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

```
## Preparing files for TOPMed imputation

## Is data lifted over?

## Extract SNPs from GSA Manifest file
plink --make-bed \
 --bfile <input_files> \
 --extract range <GSA bed file> \
 --out <output_filenames>

## Calculating freq
plink --bfile <input_files> --freq --out <output_filename>

## Checking snps against TOPMed

```
* 


## Imputation using TOPMed Imputation Server


## Post-Imputation QC

## GWAS


