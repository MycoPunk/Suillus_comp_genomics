#!/bin/bash -l
#SBATCH -N 1 -n 16 -p short
#SBATCH --mem=25G
#SBATCH --time=0-2:00:00     # 2 hr
#SBATCH --mail-user=lotuslofgren@gmail.com
#SBATCH --mail-type=ALL	
module load ncbi-blast/2.9.0+ 
module load ncbi_tools

#get the fasta names for all sequences that have 7 tm domains in each padded database	
for i in *_phobius.tab
do
	b=$(basename "$i" _phobius.tab)
	awk '($2 == "6" || $2 == "7" || $2 == "8") {print $1}' $i > $b.fasta.names
	#awk '$2 == "7" {print $1}' $i > $b.fasta.names
done

#get set fasta files for the names above
	#grab from <>expanded.wo.stops.fasta


#set blast db
JGI_DB=/rhome/llofgren/bigdata/GPCRs/seqs/db/*1kfg_2019.aa.short.fasta


#use the blast command line applications (blastdbcmd) to retrieve the fasta files from the whichever database you used to blast against
for i in *.fasta.names
do
	b=$(basename "$i" .fasta.names)
	blastdbcmd -entry_batch $i -db $JGI_DB > $b.expanded_trimmed.fasta
done

rm *.fasta.names
