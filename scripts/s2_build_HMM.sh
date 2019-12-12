#!/bin/bash -l
#SBATCH -N 1 -n 1
#SBATCH --mem=25G
#SBATCH --time=0-03:00:00     # 3 hr
#SBATCH --mail-user=lotus.lofgren@ucr.edu
#SBATCH --mail-type=ALL

###This script builds an HMM model for each class.

cd ~/bigdata/GPCR/HMM_model_seqs/expanded_GPCR_set/expanded_fastas

#load modules
module load hmmer
module load clustalo
module load parallel

#concatonate the expanded and known input files
for i in *.expanded.fasta
do
	b=$(basename "$i" .expanded.fasta)
	cat $i ./seqs/"$b"_whole.fasta > $b.HMMinput.fasta
done

#align the concatenated files with clustal omega
for i in *.HMMinput.fasta
do
	b=$(basename "$i" .HMMinput.fasta)
	clustalo -i $i -o $b.fasaln
done

#build the profile HMM in the form:
#hmmbuild [-options] hmmout.file align.file
for i in *.fasaln
do
	b=$(basename "$i" .fasaln)
	hmmbuild $b.hmm $i
done


##search the Suillus et all. sequences against each model 
#first clean to sequences - sometimes JGI puts a weird lowercase "x" at the end of the lines, before the stop codon.
perl -i -p -e 'unless( /^>/ ) {  s/x// }' *_wo_stops.fasta

#concatenate the HMM files 
cat *.hmm > gpcr_models.hmm

#compress and index the concatenated flatfile with hmmpress
hmmpress gpcr_models.hmm

#search each species against each HMM model 
parallel -j 4 hmmscan --cpu 4 --domtblout {.}.domtblout gpcr_models.hmm {} \> {.}.hmmscan ::: $(ls ./seqs/*_wo_stops.fasta)




#compress and index the flatfile with hmmpress
#for i in *.hmm
#do
#	hmmpress $i
#done

#search each species against each HMM model 
#for i in ./seqs/*_wo_stops.fasta
#do
#   for j in *.hmm
#   do
#	b=$(basename "$i" _wo_stops.fasta)
#	hmmscan --tblout $b.hmm.output.jgi.600 $j $i
#   done
#done