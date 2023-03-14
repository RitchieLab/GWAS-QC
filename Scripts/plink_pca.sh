#!/bin/bash

# bed_pca.sh ( or should we call it plink_pca.sh)

# Uses PLINK 1.9 and eigenstrat 6.0 to calculate principal components from PLINK bed/bim/fam files.

scriptname="plink_pca.sh"

usage(){
  echo ""
  echo "	Usage:  $scriptname: -b binary_file -i ibd_file [-gsfT] [-m MAF] [-p prefix] [-s marker_selection] [-l pruning_args] [-e merge_args] -n [num_evecs] -r [ldregress] -N [removal_iterations] -P [eigensoft_opts] -T [num_threads]"
  echo "		-b:  (Required) Prefix for bed, bim, fam file (plink)"
  echo "		-R:  (Required) Directory for 1000 genomes plink format reference files"
  echo "		     Build 37 files are at /project/ritchie00/datasets/1KG_Phase3/plink_files/plink_raw_files/biallelic_autosomes"
  echo "                     Build 38 files are at /project/ritchie00/datasets/1KG_Phase3/plink_files/plink_raw_files/b38/biallelic_autosomes"
  echo "		-m:  Minimum MAF (default=0.05)"
  echo " 		-s:  PLINK Marker selection arguments (default=\"--autosome --hwe 0.000001\")" 
  echo "		-l:  PLINK LD pruning arguments (default=\"--indep-pairwise 50 5 0.5\") -- to turn off pruning pass an empty string (-l \"\")"
  echo "		-e:  Plink Merge arguments (default=\"\")"
  echo "		-g:  Turns off projection onto 1000 Genome samples" 
  echo " 		-S:  Calculate projection using super-populations (EUR, etc) (default=false)"
  echo "		-f:  Do not use Fast Eigenvalue Approx. (Non fast may not work) (default=use)"
  echo "		-T:  Do not calculate Tracy-Widom Stats when not using Fast EigenValue Approx (default=calculate)"
  echo "		-n:  Number of Eigenvectors (default=20)?"
  echo "		-p:  PREFIX for files (default=PLINK_PCA)"
  echo "		-r:  LD regression parameters (default=0)"
  echo "		-N:  Number of Outlier Removal Iterations (default=5)"
  echo "		-P:  Other Eigensoft Options (comma-separated) (default is none)"
  echo "		-t:  Number of threads to use (default=1)"
  echo 
  exit 0
}

module unload plink
module load plink/1.90Beta6.18

# set defaults
maf=0.05
sel_args="--autosome --hwe 0.000001"
ld_args="--indep-pairwise 50 5 0.5"
merge_args=""
project_1kg=true
project_superop=false
fast_pca=true
twstats=true
num_evec=20
ldregress=0
numoutlier=5
pca_opts=""
nthreads=1

plink_exec=plink

# 1kg sample file  ALL.chr3.snp.biallelic.bed
#onekg_dir=/project/ritchie00/datasets/1KG_Phase3/plink_files/plink_raw_files/biallelic_autosomes/

POP_FILE=/project/ritchie00/datasets/1KG_Phase3/info_files/integrated_call_samples_v3.20130502.ALL.panel

while getopts ":b:R:m:s:l:p:e:n:p:r:N:P:t:gSfT" opt; do
  case ${opt} in
    b ) binary_prefix=$OPTARG
     ;;
    R ) onekg_dir=$OPTARG
     ;;
    m ) maf=$OPTARG
     ;; 
    s ) sel_args=$OPTARG
     ;;
    l ) ld_args=$OPTARG
     ;;
    p ) prefix=$OPTARG
     ;;
    e ) merge_args=$OPTARG
     ;;
    g ) project_1kg=false
     ;;
    S ) project_superpop=true
     ;;
    f ) fast_pca=false
     ;;
    T ) twstats=false
     ;;
    n ) num_evec=$OPTARG
     ;;
    p ) prefix=$OPTARG
     ;;
    r ) ldregress=$OPTARG
     ;;
    N ) numoutlier=$OPTARG
     ;;
    P ) pca_opts=$OPTARG
     ;;
    t ) nthreads=$OPTARG
     ;;
    \? ) usage
     ;;
    : )
      echo "Invalid option: $OPTARG requires an argument" 1>&2
     ;;
  esac
done


if [ -z "$binary_prefix" ]
then
	usage
fi

if [ -z "$onekg_dir" ]
then
  if test "$project_1kg" = true; then
	usage
  fi
fi

module load gsl/2.5
module load eigensoft/6.0.1

if test -z "$prefix"; then
  prefix=plink_pca.$binary_prefix
fi

preld=preld.$prefix
finalprefix=$prefix

onekg_dir="${onekg_dir}/"

# run downsample on the bed/bim/fam file
eval $plink_exec --bfile $binary_prefix --maf $maf "$sel_args" --out $preld -allow-no-sex --threads $nthreads --make-bed || touch $preld.bed $preld.bim $preld.fam

if test -s ${preld}.bed; then
	pre_runld=$preld
	
	if test "$project_1kg" = true; then
		SNP_LIST=$(mktemp)
		cut -f2 ${preld}.bim > $SNP_LIST
		
		for c in $(sed 's/  */\t/g' $preld.bim | cut -f1 | sort -u); do
			binary_1kg=${onekg_dir}ALL.chr${c}.snp.biallelic
			binary_out=${preld}.ALL.chr${c}.extracted

			ncolons=$(perl -n -e '$colons=split(":",(split (" ", $_))[1]);print $colons-1;last;' $binary_prefix.bim)
			onekg_bimfile=${onekg_dir}ALL.chr${c}.snp.biallelic.bim
			delete_bim=false

			if [[ $ncolons = 1 ]]; then
				onekg_colons=$(perl -n -e '$colons=split(":",(split (" ", $_))[1]);print $colons-1;last;' $binary_1kg.bim)
				if [[ $onekg_colons > 1 ]] ; then
					newbim_file=tmp.ALL.chr${c}.snp.biallelic.bim
					perl -n -p -e '@id=split " ";@pcs=split(":",$id[1]);$id[1]="$pcs[0]:$pcs[1]";$_=join("\t",@id);$_.="\n";' $onekg_bimfile > $newbim_file
					onekg_bimfile=$newbim_file
					delete_bim=true
				fi
			fi

			# extract the markers in preld
			$plink_exec --bfile $binary_1kg --bim $onekg_bimfile --extract $SNP_LIST --out $binary_out --make-bed --allow-no-sex
			if [ -f  $binary_out.bed ]; then
				if [ -z $MERGE_FILE ]; then
					MERGE_FILE=$(mktemp)
				fi
				echo -e "${binary_out}.bed\t${binary_out}.bim\t${binary_out}.fam" >> $MERGE_FILE
			fi
			if test "$delete_bim" = true ; then
				echo "deleting $onekg_bimfile"
				rm tmp.ALL.chr*.snp.biallelic.bim
			fi
		done
		
		if [ -n $MERGE_FILE ]; then
			# now, extract the markers overlapping the 1kg data
			GEN_SNPS=$(mktemp)
			for f in $(cut -f2 $MERGE_FILE); do
				cut -f2 $f >> $GEN_SNPS
			done
		
			# merge the files together
			premerge=premerge_${prefix}
			$plink_exec --bfile $preld --extract $GEN_SNPS --out $premerge --make-bed --allow-no-sex 
		
			pre_runld=pre_runld_${prefix}
			$plink_exec --bfile $premerge --merge-list $MERGE_FILE --out $pre_runld --allow-no-sex --make-bed
			if [ -e ${pre_runld}-merge.missnp ] ; then
				# remove problem snps from each file
				NEW_MERGE_FILE=$(mktemp)
				for b in $(cut -f1 $MERGE_FILE | sed 's/\.bed//') ; do
					outname=tmp.$b
					$plink_exec --bfile $b --exclude ${pre_runld}-merge.missnp --make-bed --out $outname
					echo -e "${outname}.bed\t${outname}.bim\t${outname}.fam" >> $NEW_MERGE_FILE
					# rm $b.*
				done
				# remove from user input file
				$plink_exec --bfile $premerge --exclude ${pre_runld}-merge.missnp --make-bed --out tmp.$premerge
				$plink_exec --bfile tmp.$premerge --merge-list $NEW_MERGE_FILE --out $pre_runld --allow-no-sex --make-bed
				rm tmp.$premerge.* tmp.preld.* $premerge.*
				mv ${pre_runld}-merge.missnp $prefix-merge.missnp
			fi
		fi
	fi
	#run_ld "$OUTDIR" "$PREFIX" $PRELD_NAME "$ld_args"
	# LD prune (if needed)
	if test -n "$ld_args"; then
		eval $plink_exec --bfile $pre_runld "$ld_args" --threads $nthreads -allow-no-sex --out $prefix.ld_list
		$plink_exec --bfile $pre_runld --extract $prefix.ld_list.prune.in --make-bed --out $finalprefix
		rm $prefix.ld_list.*  $preld.bim $preld.fam $preld.bed $preld.log $preld.nosex $preld.*.extracted.* $pre_runld.*
	else
		for ext in bed bim fam; do 
			mv $pre_runld.$ext $finalprefix.$ext
		done
		rm $pre_runld.log $pre_runld.nosex
	fi
		
 else
	# use the original bed/bim/fam files created above
	for ext in bed bim fam; do
		mv $preld.$ext $finalprefix.$ext
	done
fi


# run pca on the bed/bim/fam created (with $prefix as name)
#FAM_OVERALL=$(mktemp)
#MAX_N=$(cat $finalprefix.fam | wc -l)
#sed 's/[ \t][ \t]*/\t/g' f_$i.fam | cut -f1-2 | sort -t'\0' > $FAM_OVERALL
#N_F=1

input_fn=$finalprefix
# Allow some sample-level dropping to happen here (i.e. geno) if specified
SAMPLE_DROPPED=$(mktemp)

if  [ "$merge_args" ]; then
	input_fn=$finalprefix.in
	eval $plink_exec --bfile $finalprefix "$merge_args" --out $input_fn --make-bed -allow-no-sex
	join -v1 -t'\0' $finalprefix.fam <(sed 's/[ \t][ \t]*/\t/g' $input_fn.fam | cut -f1-2 | sort -t'\0') > $SAMPLE_DROPPED
	rm $finalprefix.bed $finalprefix.fam $finalprefix.bim
fi

PAR_F=$(mktemp)

# if projecting, add populations
if test "$project_1kg" = true; then
	FAM_LINENO=$(mktemp)
	nl -ba -nln -w1 $input_fn.fam | sed -e 's/  */\t/g' > $FAM_LINENO
	POP_COL=2
	if test "$project_superpop" = true; then
		POP_COL=3
	fi
	
#	join -t$'\t' -1 3 -2 1 -a 1 <(sort -t$'\t' -k3,3 $FAM_LINENO | cut -f1-6) <(cut -f1,$POP_COL $POP_FILE | sort -t$'\t' -k1,1) \
#		| sort -k2,2n -t$'\t' | cut -f1,3- \
#		| awk '{if (NF<6) print $0 "\tUNK"; else print $0;}' > $input_fn.fam.new

	join -t$'\t' -1 3 -2 1 -a 1 -o 1.1,1.2,1.3,1.4,1.5,1.6,2.2 <(sort -t$'\t' -k3,3 $FAM_LINENO | cut -f1-6) <(cut -f1,$POP_COL $POP_FILE | sort -t$'\t' -k1,1) \
	| sort -k1,1n -t$'\t' | cut -f2,3- \
	| awk '{if (NF<6) print $0 "UNK"; else print $0;}' > $input_fn.fam.new
	
	mv $input_fn.fam.new $input_fn.fam
	
	POPLIST=$(mktemp)
	cut -f $POP_COL $POP_FILE | tail -n+2 | sort -u > $POPLIST
	echo "poplist: $POPLIST" >> $PAR_F
	
else
	# remove the -9 in the last column of the fam file
	sed -i 's/-9$/UNK/' $input_fn.fam
fi


echo "genotypename: $input_fn.bed" >> $PAR_F
echo "snpname: $input_fn.bim" >> $PAR_F
echo "indivname: $input_fn.fam" >> $PAR_F
echo "evecoutname: $prefix.evec" >> $PAR_F
echo "numoutevec: $num_evec" >> $PAR_F
echo "ldregress: $ldregress" >> $PAR_F
echo "familynames: NO" >> $PAR_F

if test "$fast_pca" = false; then
	echo "evaloutname: $prefix.eval" >> $PAR_F
	echo "numthreads: $nthreads" >> $PAR_F
	echo "numoutlieriter: $numoutlier" >> $PAR_F
	echo "outlieroutname: $prefix.outlier" >> $PAR_F
	echo "numoutlierevec: $num_evec" >> $PAR_F
else
	echo "fastmode: YES" >> $PAR_F
fi

echo "$pca_opts" | sed 's/"//g' | tr ',' '\n' >> $PAR_F
echo par_f=$PAR_F

ulimit -c unlimited
echo "running smartpca..."
smartpca -p $PAR_F > $prefix.eval

rm $input_fn.bed $input_fn.bim $input_fn.fam $input_fn.log
if [ -e $input_fn.nosex ] ; then
	rm $input_fn.nosex
fi

if [ -e $prefix.outlier ]; then
  cat $SAMPLE_DROPPED <(awk '{print $3 "\t" $3}' $prefix.outlier) > $prefix.excluded
else
  cat $SAMPLE_DROPPED > $prefix.excluded
fi

if test "$fast_pca" = false; then
	if test "$twstats" = true; then
		# do the twstats calculation here
          	echo "running twstats..." 
		twstats -t /usr/share/twtable -i $prefix.eval -o $prefix.twstats
	fi
fi

if [ -s "$prefix.evec" ]; then
  pca_scree.R $prefix.evec $prefix 
fi

rm $preld.*

