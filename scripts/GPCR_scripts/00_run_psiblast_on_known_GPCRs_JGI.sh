#!/bin/bash -l
#SBATCH -N 1 -n 16
#SBATCH --mem=25G
#SBATCH --time=2-01:00:00     #2 day
#SBATCH --mail-user=lotus.lofgren@ucr.edu
#SBATCH --mail-type=ALL
#SBATCH -p batch

#GPCR code to pad the database
#module load db-ncbi/20190805
module load ncbi-blast/2.9.0+ 


#to blast protein seqs to jgi database, for each family
#tell the system where to find the database
#NCBI_DB=/srv/projects/db/ncbi/preformatted/20190805/
JGI_DB=/rhome/llofgren/bigdata/GPCRs/seqs/db/1kfg_2019.aa.short.fasta


for i in ./seqs/*whole.fasta
do
	a=$(basename "$i" .whole.fasta)
	psiblast -db $JGI_DB -num_threads 16 -num_iterations 5 -num_alignments 800 -num_descriptions 800 -query $i -out $a.whole.psiblastout
done
