# GWAS-QC
## Workflow is implemented on 1000Genomes on Build 38 (abbreviated as 1KG here)

## GWAS QC Workflow on Example Dataset
We chose a publicly available dataset from the International Genome Sample Resource (IGSR) (www.internationalgenome.org). IGSR created and maintains the 1000 Genomes (1KG) Project to provide a public catalog of common human genetic variation and genotype data. The 1KG dataset has been kept up to date with current reference data sets, thus it is available for both GRCh37 and GRCh38. The latter is utilized here because the 2014 update increased the quantity of loci represented, resolved more than 1000 issues from the previous version of the assembly, and overall provides a better basis for alignment and subsequent analysis. Additionally, IGSR’s continued efforts will lead to the incorporation of various populations to the data which were not previously captured.

## Data locations
* Raw files: /project/ritchie00/datasets/1KG_Phase3/plink_files/plink_raw_files/b38
* Project output diretory: ~group/projects/cphg-gwas-qc
  * Illumina GSA Manifest file (cleaned & converted to bed format): ~metadata/GSA-24v3-0_A2_cleaned.bed

## Browsing raw data
* We go into the directory containing our raw genotype files and inspect to get a sense of them
* For 1KG, we have .bed, .bim, and .fam files
  * *.fam files contain information about samples, one sample per line
  * *.bim files contain information about markers, one marker per line
  * *.bed contain binary genotype information (don’t view this file directly)
* If you're unfamiliar with genomic file types, we will end up needing all of these files because they each contain some form of information we need
* It looks like we have files for chromosome 1-22

## Pre-Imputation QC
During the Pre-Imputation steps, we will filter the 1000 Genomes raw files to retain only the SNPs which match the Illumina Infinium Global Screening Array v3.0 Manifest File for our build. Although the 1KG dataset is complete, we will sample it down in order to perform the GWAS QC steps. Pre-imputation QC is done mostly in PLINK. If we are running GWAS on unrelated samples using a generalized regression, then we would use PLINK. Instead, if we are running analyses on all samples using a mixed linear model or bayesian approach, we would use SAIGE or Regenie, respectively.

For manuscript text: “For example, we took 1000 genomes data and we got ### variants. Write how many variants at each step for Shefali to check my work”
* Sex check
* Missingness
* Flipped

First, download the Global Screening Array (GSA) Manifest File from [Illumina’s website] (https://support.illumina.com/downloads/infinium-global-screening-array-v3-0-product-files.html)
* Note: SNP IDs are not always portable so extract the Chr position without having to input SNP ID which is a command like “extract range”
* 1000 Genomes is the whole genome sequence data from which we want to extra the SNPs that are in the GSA manifest file. We are doing this so the filtered data matches a common genotyping chip.

## Download the GSA Manifest File from Illumina’s site and copy it to the project location
> scp GSA-24v3-0_A2.csv heyvan@superman.pmacs.upenn.edu:/home/heyvan/group/projects/cphg-gwas-qc

## Look at the first few lines to see what fields it has
```
head -10 GSA-24v3-0_A2.csv
> output
Illumina, Inc.
[Heading]
Descriptor File Name,GSA-24v3-0_A2.bpm
Assay Format,Infinium HTS
Date Manufactured,8/7/2019
Loci Count ,654027
[Assay]
IlmnID,Name,IlmnStrand,SNP,AddressA_ID,AlleleA_ProbeSeq,AddressB_ID,AlleleB_ProbeSe
q,GenomeBuild,Chr,MapInfo,Ploidy,Species,Source,SourceVersion,SourceStrand,SourceSe
q,TopGenomicSeq,BeadSetID,Exp_Clusters,Intensity_Only,RefStrand
1:103380393-0_B_R_2346041316,1:103380393,BOT,[T/C],0009663149,AATAAACTTTTATGCAAAACT
TGTAAGATAACTCTTCTTTCCTTCTTCTT,,,38,1,102914837,diploid,Homo sapiens,1000genomes,0,T
OP,GCTTCCCCTTTCTCTCCTCTTTCTCCTTTGGGACCCTAAACAATGTTAAAAAAAAAAAAA[A/G]AAGAAGAAGGAAAGA
AGAGTTATCTTACAAGTTTTGCATAAAAGTTTATTAACCTTGGCA,GCTTCCCCTTTCTCTCCTCTTTCTCCTTTGGGACCCT
AAACAATGTTAAAAAAAAAAAAA[A/G]AAGAAGAAGGAAAGAAGAGTTATCTTACAAGTTTTGCATAAAAGTTTATTAACCT
TGGCA,1895,3,0,-
1:109439680-0_T_F_2348625138,1:109439680,TOP,[A/G],0026641362,GTTGGGCAATGCTTATTTCTA
TTTGCATGATTATGCCAAAGCATTAGAAT,,,38,1,108897058,diploid,Homo sapiens,1000genomes,0,T
OP,TTTACAGCCAGTTGGGCAATGCTTATTTCTATTTGCATGATTATGCCAAAGCATTAGAAT[A/G]TCACCATCATGATTT
AACCCTTGCAAGGTAATTAATTTAAGCTTTTAAATATTCTTCCTT,TTTACAGCCAGTTGGGCAATGCTTATTTCTATTTGCA
TGATTATGCCAAAGCATTAGAAT[A/G]TCACCATCATGATTTAACCCTTGCAAGGTAATTAATTTAAGCTTTTAAATATTCT
TCCTT,1964,3,0,+
```

## Parse the csv to grab the fields we need (chr, start position & stop position from mapinfo)
```
awk -F ',' '{print "chr"$10, $11, $11+1, NR}' GSA-24v3-0_A2.csv > GSA-24v3-0_A2.bed
```

## Check the file lengths
```
wc -l <filename>
```

## There’s some weird metadata in that .csv file so here we’ll grab just the lines with 4 columns
```
awk '{if (NF == 4) {print $0}}' GSA-24v3-0_A2_cleaned.bed
```

## Starting out
* Data is already Build 38 so it doesn't need ot be lifted over from Build 37 to 38 using liftOver
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

## Extract SNPs from GSA Manifest file
```
in_path=/project/ritchie00/datasets/1KG_Phase3/plink_files/plink_raw_files/b38/
out_path=~/group/projects/cphg-gwas-qc/

 for i in {1..22}; do \
 plink --make-bed \
 --bfile ${in_path}ALL.chr${i}_GRCh38.genotypes.20170504.genotypes \
 --extract range ~metadata/GSA-24v3-0_A2_cleaned.bed \
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
plink --bfile ~unmerged/chr1_GSA-filtered --merge-list ../mergelist.txt --make-bed --out 1KG_GSA-filtered_merged
 
## if it finds error: # variants with 3+ alleles present, go back and re-generate chromosomes from raw input with the aberrant snps using the --exclude command with the .missnp file
for i in {4..4}; do \
 plink --make-bed \
 --bfile ${in_path}ALL.chr${i}_GRCh38.genotypes.20170504.genotypes \
 --extract range ~metadata/GSA-24v3-0_A2_cleaned.bed \
 --exclude prepare/1KG_GSA-filtered_merged-merge.missnp \
 --out ${out_path}prepare/unmerged/chr${i}_GSA-filtered; done
 
 ## now repeat the merge again
 plink --bfile ~unmerged/chr1_GSA-filtered --merge-list ../mergelist.txt --make-bed --out 1KG_GSA-filtered_merged
```
![image](https://user-images.githubusercontent.com/30478823/146056961-3a72cfa8-aaa4-4dc2-86fe-63f5c5947565.png)

Based on the above image, we see that there are 2504 ambigious individuals with 0 males and 0 females. 

Checking the .fam file from our merge, we see that column 5 has 0's when it should have 1's or 2's to represent female/male. 
![image](https://user-images.githubusercontent.com/30478823/153645848-a03b98f3-d997-4857-bcb3-350d99bc2a2c.png)


## Create sex-file (FID, IID, sex (coded as 1 or 2) and pheno-file (FID, IID, pheno)
20130606_g1k.ped = File with Family ID, Individual ID, Gender etc to use to create our sex-file
![image](https://user-images.githubusercontent.com/30478823/153914145-636c8e0c-1fe2-4f1d-b702-a8ef0ec7ee99.png)

igsr_samples.tsv = File with ancestryto use to create our pheno-file
![image](https://user-images.githubusercontent.com/30478823/153914187-f5420695-a71c-4b62-a331-e8231eccb095.png)

```
# my query to grab the family ID, individual ID, and sex info (not good)
awk -F '\t' '{print $1, $2, $5}'  20130606_g1k.ped > sex-file.txt
```
![image](https://user-images.githubusercontent.com/30478823/153914094-00a1cd19-59ec-4a30-943f-6c3c92052496.png)


```
# Shefali's query to grab the items (this one worked best)
join -1 1 -2 1 <(cat 1KG_GSA-filtered_merged.fam |sort -k1,1) <(cat ../20130606_g1k.ped |awk -F '\t' '{print $2,$5}' |sort -k1,1) |awk '{print $1,$2,$7}' >sex_file.txt
```
![image](https://user-images.githubusercontent.com/30478823/153914059-e3e680ba-9326-4b74-8fc9-7cb6430d3c5d.png)


```
# do the same for phenotype to have pheno_file on hand for GWAS
#  Extract the population info from igsr_samples.tsv file
# then replace AFR=2 and others (EUR, AMR, EAS, SAS)=1
```

## Update-sex to update our 1KG files with the appropriate sex info
```
plink --bfile 1KG_GSA-filtered_merged --update-sex sex_file.txt --make-bed --out 1KG_GSA-filtered_merged_withsex
```
This is what it looks like after we update sex.
![image](https://user-images.githubusercontent.com/30478823/153913980-85bb7a89-8107-4daa-b9be-32e7c86d8663.png)


## Run pre-imputation QC before you start your freq command
Yuki's usual PMBB pre-imputation QC criteria are: 95% snp call rate / 90% sample call rate / sexcheck-failed samples removed.
```
plink --bfile 1KG_GSA-filtered_merged_withsex --geno 0.05 --mind 0.1 --make-bed --out 1KG_GSA-filtered_merged_withsex_QC
```
![image](https://user-images.githubusercontent.com/30478823/153922085-c9e51727-315c-4a29-bc9d-fc682067b3ae.png)

## Perform heterozygosity check in R to generate plots
```
```

## Prepare files to upload to TOPMed Imputation Server
```
## Calculating freq
plink --bfile 1KG_GSA-filtered_merged_withsex_QC --freq --out 1KG_GSA-filtered_merged_withsex_QC_freq
```
![image](https://user-images.githubusercontent.com/30478823/146056904-ab43216b-7d42-4469-acfa-580d963856fa.png)

```
## Checking snps against TOPMed
(http://www.well.ox.ac.uk/~wrayner/tools/)
## Checks Strand, alleles, position, Ref/Alt assignments and frequency differences
## Updates: Strand, position, ref/alt assignment
## Removes: A/T & G/C SNPs if MAF > 0.4, SNPs with differing alleles, SNPs with > 0.2 allele frequency difference, SNPs not in reference panel
## (TO-DO) Using HRC-1000G-check-bim-v4.3.0
perl /home/yub7/group/projects/PMBB/QC_Imputation/scripts/HRC-1000G-check-bim.pl \
-r /home/yub7/group/projects/PMBB/QC_Imputation/scripts/PASS.Variants.TOPMed_freeze5_hg38_dbSNP.tab -h \
-b 1KG_GSA-filtered_merged_withsex_QC.bim \
-f 1KG_GSA-filtered_merged_withsex_QC_freq.frq
```
![image](https://user-images.githubusercontent.com/30478823/154357815-632c6c56-eac7-4ce2-ae26-efaea8ac2de5.png)
![image](https://user-images.githubusercontent.com/30478823/154360757-f09aabb9-2298-45eb-b6b0-ebe7e5722b78.png)

## Coding Chr from "1","23" into "Chr1","ChrX" format
```
sed -i 's/--recode vcf/--recode vcf --output-chr chrM/g' Run-plink.sh
chmod +x ./Run-plink.sh
./Run-plink.sh
```

## TODO Creating VCF file aligned with build38 reference alleles (downloaded from: https://console.cloud.google.com/storage/browser/genomics-public-data/resources/broad/hg38/)
```
module load bcftools/1.9
export BCFTOOLS_PLUGINS=/appl/bcftools-1.9/libexec/bcftools/

for i in {1..22}; do \
bcftools +fixref 1KG_GSA-filtered_merged_withsex_QC-updated-chr$i'.vcf' \
-Ov -o 1KG_GSA-filtered_merged_withsex_QC-updated_flipped_chr$i'.vcf' \
-- -d -f ~/group/projects/PMBB/QC_Imputation/scripts/resources_broad_hg38_v0_Homo_sapiens_assembly38.fasta -m flip; done
```

## TODO Sorting VCF and zipping files using VCFtools and tabix (make sure the module are loaded first)
```
module load vcftools/0.1.12c
module load tabix/0.2.6
for i in {1..22}; do \
vcf-sort 1KG_GSA-filtered_merged_withsex_QC-updated_flipped_chr$i.vcf | \
bgzip -c > 1KG_ImputationInput_TOPMED_chr$i.vcf.gz; done
```

## TODO Run VCF check (downloaded from https://github.com/zhanxw/checkVCF)
```
for i in {1..23}; do python ~/group/projects/PMBB/QC_Imputation/scripts/checkVCF.py \
-r ~/group/projects/PMBB/QC_Imputation/scripts/resources_broad_hg38_v0_Homo_sapiens_assembly38.fasta \
 -o test 1KG_ImputationInput_TOPMED_chr$i.vcf.gz; mv test.check.log test.check_chr$i.log >> test_check_set1.log; done
```


## Compute Principal Component Analyses on Pre-Imputed Data
```
```

## Imputation using TOPMed Imputation Server
Imputation has become an essential component of GWAS quality control because it increases power, facilitates meta-analysis, and aids interpretation of signals. Genotype imputation is the statistical inference of unobserved genotype, which enables scientists to reconstruct the missing data in each genome and accurately evaluate the evidence for association at genetic markers that were not genotyped. Genotype imputation is achieved by comparing short stretches of an individual genome against stretches of previously characterized reference genomes.  It is usually performed on single nucleotide polymorphisms (SNPs), which are the most common type of genetic variation. 
Several tools exist specifically for genotype imputation such as the Michigan and Trans-Omics for Precision Medicine (TOPMed) Imputation Servers where one uploads the phased or unphased GWAS genotypes in order to receive the imputed genomes in return. Each imputation server varies in terms of speed and accuracy. One of the most important considerations in imputation is the composition of the reference panel. For our study, we selected the TOPMed Imputation Reference panel  (version r2) because it is one of the most diverse reference panels available and contains information from 97,256 deeply sequenced human genomes containing 308,107085 genetic variants distributed across the 22 autosomes and the X chromosome. 
Theoretically, phased means that the two strands on each Chr are separated to identify which regions come from each parent whereas no phasing means that they are not separated. Essentially, for imputation phasing is the first step which is done in reference to the reference genome panel.

## Post-Imputation QC
XYZ
```
<code>
```

## Performing GWAS

## PLINK
```
<code>
```

## Regenie
```
<code>
```

## SAIGE
```
<code>
```


## Additional material
We recommend the following resources and tutorials developed for performing GWAS:
* Comphrehensive tutorial about GWAS and PRS by MareesAT: https://github.com/MareesAT/GWA_tutorial/
* GWAS Data Cleaning tutorial by the GENEVA Coordinating Center: https://www.bioconductor.org/packages/devel/bioc/vignettes/GWASTools/inst/doc/DataCleaning.pdf
