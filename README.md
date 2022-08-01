
<h1 align="center">
  <br>
  <a href="https://github.com/RitchieLab/GWAS-QC"><img src="https://user-images.githubusercontent.com/30478823/182218051-a9a58111-9b92-42e9-bdc1-c2ffad393fe2.png" alt="GWAS QC" width="200"></a>
  <br>
  GWAS Quality Control (QC)
  <br>
</h1>

<h4 align="center">A complementary tutorial to the <a href="https://pubmed.ncbi.nlm.nih.gov/21234875/" target="_blank">Truong, Woerner, Cherlin, et al 2022 Paper</a>.</h4>

<p align="center">
  <a href="#basic-overview">Basic Overview</a> •
  <a href="#set-up">Set Up</a> •
  <a href="#download">Download</a> •
  <a href="#related-resources">Related Resources</a> •
  <a href="#license">License</a>
</p>




## Basic Overview


## Set Up

## Download

## Related Resources

## License



# GWAS-QC
## Workflow is implemented on 1000Genomes on Build 38 (abbreviated as 1KG here)

## GWAS QC Workflow on Example Dataset
We chose a publicly available dataset from the International Genome Sample Resource (IGSR) (www.internationalgenome.org). IGSR created and maintains the 1000 Genomes (1KG) Project to provide a public catalog of common human genetic variation and genotype data. The 1KG dataset has been kept up to date with current reference data sets, thus it is available for both GRCh37 and GRCh38. The latter is utilized here because the 2014 update increased the quantity of loci represented, resolved more than 1000 issues from the previous version of the assembly, and overall provides a better basis for alignment and subsequent analysis. Additionally, IGSR’s continued efforts will lead to the incorporation of various populations to the data which were not previously captured.

## Setting up your environment
### Modules in BASH
* plink/1.9
* plink/2.0
* bcftools/1.9
* vcftools/0.1.12c
* tabix/0.2.6
* 

### Notes on PLINK v1.9 and v2.0
* Not all commands are portable between PLINK version 1.9 and version 2.0. Since PLINK v2.0 is under heavy active development, the developers urge users to check certain results against an earlier, more widely-used version of PLINK. Some functions are available in v1.9 which are not in v2.0, and vice versa.
* Original version of PLINK: 1.07, https://zzz.bwh.harvard.edu/plink/plink2.shtml 
* Beta version: 1.90, https://www.cog-genomics.org/plink/1.9/
* Alpha version: 2.00, https://www.cog-genomics.org/plink/2.0/ 
* [Developer's comments](https://www.biostars.org/p/299855/) 1.9 and 2.0 serving as complementary resources
> The main difference is that plink 1.9 is essentially finished, while plink 2.0 is an alpha-stage program which will have significant unfinished components for a while to come. As a consequence, current development priorities for plink 2.0 are centered around things which are impossible to do with plink 1.9, such as handling multiallelic/phased variants and dosage data and reliably tracking REF/ALT alleles; while things that plink 1.9 already handles perfectly well, such as working with .ped/.map file pairs, have been deliberately deprioritized for now. So, you should stick to 1.9 as long as it's good enough for the jobs you need to perform. But once you need to do something outside 1.9's scope, you're likely to find that 2.0 already has the additional feature you need (or it'll be added quickly after you ask for it)

## Data locations
* New raw files: ~/group/personal/tess/GWAS_Tutorial
* Old raw files: /project/ritchie00/datasets/1KG_Phase3/plink_files/plink_raw_files/b38
* Metadata: ~/group/projects/cphg-gwas-qc/metadata
  * 1000 Genome's pedigree file: ~/metadata/20130606_g1k.ped
  * Illumina GSA Manifest file (cleaned & converted to bed format): ~metadata/GSA-24v3-0_A2_cleaned.bed
* Project output diretory: TBD
* QC'd Pre-Imputation Data: ~group/projects/cphg-gwas-qc/prepare
* Imputated Data: ~group/scratch/van/cphg-gwas-qc/imputed-data

## Browsing raw data
* We go into the directory containing our raw genotype files and inspect to get a sense of them
* For 1KG, we have .bed, .bim, and .fam files
  * *.fam files contain information about samples, one sample per line
  * *.bim files contain information about markers, one marker per line
  * *.bed contain binary genotype information (don’t view this file directly)
* If you're unfamiliar with genomic file types, we will end up needing all of these files because they each contain some form of information we need
* It looks like we have files for chromosome 1-22

# I. Subsampling our complete data for this tutorial
During the Pre-Imputation steps, we will filter the 1000 Genomes raw files to retain only the SNPs which match the Illumina Infinium Global Screening Array v3.0 Manifest File for our build. Although the 1KG dataset is complete, we will sample it down in order to perform the GWAS QC steps. Pre-imputation QC is done mostly in PLINK. If we are running GWAS on unrelated samples using a generalized regression, then we would use PLINK. Instead, if we are running analyses on all samples using a mixed linear model or bayesian approach, we would use SAIGE or Regenie, respectively.

For manuscript text: “For example, we took 1000 genomes data and we got ### variants. Write how many variants at each step for Shefali to check my work”
* Sex check
* Missingness
* Flipped

First, download the Global Screening Array (GSA) Manifest File from [Illumina’s website] (https://support.illumina.com/downloads/infinium-global-screening-array-v3-0-product-files.html)
* Note: SNP IDs are not always portable so extract the Chr position without having to input SNP ID which is a command like “extract range”
* 1000 Genomes is the whole genome sequence data from which we want to extra the SNPs that are in the GSA manifest file. We are doing this so the filtered data matches a common genotyping chip.

## 1. Download the GSA Manifest File from Illumina’s site and copy it to the project location
> scp GSA-24v3-0_A2.csv heyvan@superman.pmacs.upenn.edu:/home/heyvan/group/projects/cphg-gwas-qc/metadata

## 2. Look at the first few lines to see what fields it has
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

## 3. Parse the csv to grab the fields we need (chr, start position & stop position from mapinfo)
```
awk -F ',' '{print "chr"$10, $11, $11+1, NR}' GSA-24v3-0_A2.csv > GSA-24v3-0_A2.bed
```

## 4. Check the file lengths
```
wc -l <filename>
```

## 5. There’s some weird metadata in that .csv file so here we’ll grab just the lines with 4 columns
```
awk '{if (NF == 4) {print $0}}' GSA-24v3-0_A2_cleaned.bed
```

## 6. Starting out
* Data is already Build 38 so it doesn't need to be lifted over from Build 37 to 38 using liftOver
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

## 7. Extract SNPs from GSA Manifest file
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

## 8. Merge per-chromosome files into whole-genome input file
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


## 9. Create sex-file (FID, IID, sex (coded as 1 or 2) and pheno-file (FID, IID, pheno)
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

## 10. Update-sex to update our 1KG files with the appropriate sex info
```
plink --bfile 1KG_GSA-filtered_merged --update-sex sex_file.txt --make-bed --out 1KG_GSA-filtered_merged_withsex
```
This is what it looks like after we update sex.
![image](https://user-images.githubusercontent.com/30478823/153913980-85bb7a89-8107-4daa-b9be-32e7c86d8663.png)

# II. Pre-Imputation QC

## 11. Run pre-imputation QC before you start your freq command
Yuki's usual PMBB pre-imputation QC criteria are: 95% snp call rate / 90% sample call rate / sexcheck-failed samples removed.
```
plink --bfile 1KG_GSA-filtered_merged_withsex --geno 0.05 --mind 0.1 --make-bed --out 1KG_GSA-filtered_merged_withsex_QC
```
![image](https://user-images.githubusercontent.com/30478823/153922085-c9e51727-315c-4a29-bc9d-fc682067b3ae.png)

## 12. Perform heterozygosity check in R to generate plots
```
```

## 13. Prepare files to upload to TOPMed Imputation Server
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

## 14. Coding Chr from "1","23" into "Chr1","ChrX" format
This will recode the chromosomes into the correct syntax that Plink v2.00 can use it. Plink v1.9 and earlier does not write vcf in the right format for TopMED. 

```
sed -i 's/--recode vcf/--recode vcf --output-chr chrM/g' Run-plink.sh
chmod +x ./Run-plink.sh
./Run-plink.sh
```

## 15. Creating VCF file aligned with build38 reference alleles (downloaded from: https://console.cloud.google.com/storage/browser/genomics-public-data/resources/broad/hg38/)
```
module load bcftools/1.9
export BCFTOOLS_PLUGINS=/appl/bcftools-1.9/libexec/bcftools/

for i in {1..22}; do \
bcftools +fixref 1KG_GSA-filtered_merged_withsex_QC-updated-chr$i'.vcf' \
-Ov -o 1KG_GSA-filtered_merged_withsex_QC-updated_flipped_chr$i'.vcf' \
-- -d -f ~/group/projects/PMBB/QC_Imputation/scripts/resources_broad_hg38_v0_Homo_sapiens_assembly38.fasta -m flip; done
```
![image](https://user-images.githubusercontent.com/30478823/154592044-bd1e7f13-ffff-4b95-8efb-e21e02c46f02.png)

## 16. Sorting VCF and zipping files using VCFtools and tabix (make sure the module are loaded first)
```
module load vcftools/0.1.12c
module load tabix/0.2.6
for i in {1..22}; do \
vcf-sort 1KG_GSA-filtered_merged_withsex_QC-updated_flipped_chr$i.vcf | \
bgzip -c > 1KG_ImputationInput_TOPMED_chr$i.vcf.gz; done
```
![image](https://user-images.githubusercontent.com/30478823/154592893-04983705-37db-4aa1-b761-7dbe0da27a7b.png)
![image](https://user-images.githubusercontent.com/30478823/154592906-7c88fb93-79b8-451a-b4e7-acf925d5bb08.png)

## 17. Run VCF check (downloaded from https://github.com/zhanxw/checkVCF)
```
for i in {1..23}; do python ~/group/projects/PMBB/QC_Imputation/scripts/checkVCF.py \
-r ~/group/projects/PMBB/QC_Imputation/scripts/resources_broad_hg38_v0_Homo_sapiens_assembly38.fasta \
 -o test 1KG_ImputationInput_TOPMED_chr$i.vcf.gz; mv test.check.log test.check_chr$i.log >> test_check_set1.log; done
```
![image](https://user-images.githubusercontent.com/30478823/154593256-cc5a4372-421f-402c-ba02-907fbe854424.png)


## 18. Compute Principal Component Analyses on Pre-Imputed Data
```
```

## III. Imputation using TOPMed Imputation Server
Imputation has become an essential component of GWAS quality control because it increases power, facilitates meta-analysis, and aids interpretation of signals. Genotype imputation is the statistical inference of unobserved genotype, which enables scientists to reconstruct the missing data in each genome and accurately evaluate the evidence for association at genetic markers that were not genotyped. Genotype imputation is achieved by comparing short stretches of an individual genome against stretches of previously characterized reference genomes.  It is usually performed on single nucleotide polymorphisms (SNPs), which are the most common type of genetic variation. 

Several tools exist specifically for genotype imputation such as the Michigan and Trans-Omics for Precision Medicine (TOPMed) Imputation Servers where one uploads the phased or unphased GWAS genotypes in order to receive the imputed genomes in return. Each imputation server varies in terms of speed and accuracy. One of the most important considerations in imputation is the composition of the reference panel. For our study, we selected the TOPMed Imputation Reference panel  (version r2) because it is one of the most diverse reference panels available and contains information from 97,256 deeply sequenced human genomes containing 308,107085 genetic variants distributed across the 22 autosomes and the X chromosome. 
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


## Additional material
We recommend the following resources and tutorials developed for performing GWAS:
* Comphrehensive tutorial about GWAS and PRS by MareesAT: https://github.com/MareesAT/GWA_tutorial/
* GWAS Data Cleaning tutorial by the GENEVA Coordinating Center: https://www.bioconductor.org/packages/devel/bioc/vignettes/GWASTools/inst/doc/DataCleaning.pdf
