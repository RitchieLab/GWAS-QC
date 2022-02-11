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


## Starting out
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
in_path=/project/ritchie00/datasets/1KG_Phase3/plink_files/plink_raw_files/b38/
out_path=~/group/projects/cphg-gwas-qc/

 for i in {1..22}; do \
 plink --make-bed \
 --bfile ${in_path}ALL.chr${i}_GRCh38.genotypes.20170504.genotypes \
 --extract range GSA-24v3-0_A2_cleaned.bed \
 --exclude 1KG_GSA-filtered_merged-merge.missnp
 --out ${out_path}prepare/chr${i}_GSA-filtered; done
```
![image](https://user-images.githubusercontent.com/30478823/146036594-4590d0eb-0753-4bd7-9404-32de174f3c89.png)

## Merge per-chromosome files into whole-genome input file
```
## create a text file which is our merge list containing a line for each chromosome we are working with
cat > mergelist.txt
# then type in all the chr file prefixes separated by line minus the first chromosome, chr1
```
![image](https://user-images.githubusercontent.com/30478823/146052632-032cc0b4-f0b3-4b10-a369-ae58ffce5a4c.png)


```
## cd into prepare/ then merge the chr into 1 file
plink --bfile chr1_GSA-filtered --merge-list ../mergelist.txt --make-bed --out 1KG_GSA-filtered_merged
 
## if it finds error: # variants with 3+ alleles present, go back and re-generate chromosomes from raw input with the aberrant snps using the --exclude command with the .missnp file
for i in {4..4}; do \
 plink --make-bed \
 --bfile ${in_path}ALL.chr${i}_GRCh38.genotypes.20170504.genotypes \
 --extract range GSA-24v3-0_A2_cleaned.bed \
 --exclude prepare/1KG_GSA-filtered_merged-merge.missnp \
 --out ${out_path}prepare/chr${i}_GSA-filtered; done
 
 ## now repeat the merge again
 plink --bfile chr1_GSA-filtered --merge-list ../mergelist.txt --make-bed --out 1KG_GSA-filtered_merged
```
![image](https://user-images.githubusercontent.com/30478823/146056961-3a72cfa8-aaa4-4dc2-86fe-63f5c5947565.png)

Based on the above image, we see that there are 2504 ambigious individuals with 0 males and 0 females. 

## create sex-file (FID, IID, sex (coded as 1 or 2) and pheno-file (FID, IID, pheno)
20130606_g1k.ped = File with Family ID, Individual ID, Gender etc to use to create our sex-file
igsr_samples.tsv = File with ancestryto use to create our pheno-file

```

```

## update-sex and update-pheno in the same command on the
```
plink --bfile 1KG_GSA-filtered_merged --update-sex <sex-file> --update-pheno <pheno-file> --make-bed --out 1KG_GSA-filtered_merged-with-sex-pheno
```
This is what it looks like after we only update sex.
![image](https://user-images.githubusercontent.com/30478823/153644627-ef3dd315-c113-45ee-85ab-b2d075210ba8.png)


## Run Pre-imputation QC before you start your freq command
Yuki's usual PMBB pre-imputation QC criteria are: 95% snp call rate / 90% sample call rate / sexcheck-failed samples removed.
```
```

## Prepare files to upload to TOPMed Imputation Server
```
## Calculating freq
plink --bfile 1KG_GSA-filtered_merged --freq --out 1KG_GSA-filtered_merged_freq
```
![image](https://user-images.githubusercontent.com/30478823/146056904-ab43216b-7d42-4469-acfa-580d963856fa.png)

![image](https://user-images.githubusercontent.com/30478823/146056829-047a7afc-44c0-4435-8afd-d0e8975f9863.png)

```
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


