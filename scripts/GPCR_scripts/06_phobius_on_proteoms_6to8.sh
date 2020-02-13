#!/bin/bash -l
#SBATCH -N 1 -n 20
#SBATCH --mem=25G
#SBATCH --time=0-2:00:00     # 2 days
#SBATCH --mail-user=lotuslofgren@gmail.com
#SBATCH --mail-type=ALL

#load modules
module load ncbi-blast/2.9.0+ 
module load ncbi_tools

cd seqs/proteins

#CPU=$SLURM_CPUS_ON_NODE
##if [ -z $CPU ]; then
## CPU=1
#fi
##module load parallel
##module load phobius
##run_phobius () {
## name=$1
## base=$(basename "$name" .fasta) 
##if [ ! -s ${base}_phobius.tab ]; then
##        phobius -short $name > ${base}_phobius.tab
## fi
##}
##export -f run_phobius

##time parallel -j $CPU run_phobius ::: $(ls *.fasta)


#get the fasta names for all sequences that have 7 tm domains in each padded database	
for i in *_phobius.tab
do
	b=$(basename "$i" _phobius.tab)
	awk '($2 == "6" || $2 == "7" || $2 == "8") {print $1}' $i > $b.fasta.names
	#awk '$2 == "7" {print $1}' $i > $b.fasta.names
done

#set blast db
JGI_DB=/rhome/llofgren/bigdata/GPCRs/seqs/db/*_updated.fasta

#use the blast command line applications (blastdbcmd) to retrieve the fasta files from the whichever database you used to blast against
for i in *.fasta.names
do
	b=$(basename "$i" .fasta.names)
	blastdbcmd -entry_batch $i -db $JGI_DB > $b.trimmed.fasta
done

rm *.fasta.names
rm *_phobius.tab
