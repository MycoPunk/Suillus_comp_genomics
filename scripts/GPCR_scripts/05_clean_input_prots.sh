#!/bin/bash -l
#SBATCH -N 1 -n 16
#SBATCH --mem=25G
#SBATCH --time=0-01:00:00     #1 hr
#SBATCH --mail-user=lotus.lofgren@ucr.edu
#SBATCH --mail-type=ALL
#SBATCH -p batch

#this script cleans the input (whole genome protein) files and preps them for the HMM matching 

cd seqs/proteins

#format input seq headers

for i in *.fasta
do
	b=$(basename "$i" .fasta)
	perl -p -e 's/>jgi\|(\S+)\|(\d+)\|/>$1|$1_$2 /' $i > $b.clean.fa
done

rm *.fasta

#remove stop codons
for i in *.fa
do
	b=$(basename "$i" .fa)
	sed 's/*//g' $i > $b.fasta
done


#remove weird x chs. 
perl -i -p -e 'unless( /^>/ ) {  s/x// }' ~/bigdata/GPCRs/seqs/proteins/*.fasta
