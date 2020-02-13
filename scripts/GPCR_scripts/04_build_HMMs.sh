#!/bin/bash -l
#SBATCH -N 1 -n 1
#SBATCH --mem=25G
#SBATCH --time=0-05:00:00     # 5 hr
#SBATCH --mail-user=lotus.lofgren@ucr.edu
#SBATCH --mail-type=ALL

##This script builds an HMM model from each of the trimmed, padded 7tm databases for each GPCR class 

#load modules
module load hmmer
module load clustalo
module load parallel

#concatonate the expanded and known input files
for i in *.expanded_trimmed.fasta
do
	b=$(basename "$i" .expanded_trimmed.fasta)
	cat $i ./seqs/"$b".fasta > $b.HMMinput.fasta
done


#align the concatenated files with clustal omega
for i in *.HMMinput.fasta
do
	b=$(basename "$i" .HMMinput.fasta)
	clustalo -i $i -o $b.fasaln
done

#rm *.HMMinput.fasta

#build the profile HMM in the form:
#hmmbuild [-options] hmmout.file align.file
for i in *.fasaln
do
	b=$(basename "$i" .fasaln)
	hmmbuild $b.hmm $i
done

# rm *.fasaln
