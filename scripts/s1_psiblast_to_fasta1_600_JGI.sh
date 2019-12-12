#!/bin/bash -l
#SBATCH -N 1 -n 16
#SBATCH --mem=25G
#SBATCH --time=0-02:00:00     #2 hr
#SBATCH --mail-user=lotus.lofgren@ucr.edu
#SBATCH --mail-type=ALL

#this script turn psi_blast_output into fasta files to build HMM models 

cd /bigdata/GPCR/HMM_model_seqs/expanded_GPCR_set/expanded_fastas

#load modules
module load ncbi-blast/2.9.0+ 
module load ncbi_tools

#limit sequence ID length to 49 chrs
#perl -p -e 's/^>(\S{1,49})\S*.+/>$1/' 1kfg_2019.aa.fasta > 1kfg_2019.aa.short.fasta

#index the headers
#makeblastdb -in 1kfg_2019.aa.short.fasta -parse_seqids -dbtype prot

#set database
JGI_DB=/rhome/llofgren/bigdata/GPCR/HMM_model_seqs/expanded_GPCR_set/expanded_fastas/1kfg_2019.aa.short.fasta

#extract just the final iteration
for i in GPCR_*
do
	b=$(basename "$i" .psiblastout)
	tac $i | grep 'Results' -m 1 -B 999999 | tac > $b.final_iteration.txt
done

#grep all lines that start with ">"  (non-redundant seq hits), remove ">" and everything after the first white space
for i in *final_iteration.txt
do
	b=$(basename "$i" .final_iteration.txt)
	grep '^>' $i | perl -p -e 's/>(\S+)\s+.+/$1/' > $b.ID_lines.txt
done


#use the blast command line applications (blastdbcmd) to retrieve the fasta files from the whichever database you used to blast against
for i in *ID_lines.txt
do
	b=$(basename "$i" .ID_lines.txt)
	blastdbcmd -entry_batch $i -db $JGI_DB > $b.expanded.fasta
done
