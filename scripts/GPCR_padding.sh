#!/bin/bash -l
#SBATCH -N 1 -n 16
#SBATCH --mem=25G
#SBATCH --time=1-00:00:00     # 1 day
#SBATCH --mail-user=lotus.lofgren@ucr.edu
#SBATCH --mail-type=ALL

#GPCR code to pad the database
module load db-ncbi/20190805
module load ncbi-blast/2.9.0+ 
#psiblast -help

##blast protein seqs to ncbi database, for each family
#tell the system where to find the database
NCBI_DB=/srv/projects/db/ncbi/preformatted/20190805/

psiblast -db $NCBI_DB/nr -num_threads 16 -query GPCR_class2_whole.fasta -out class2_GPCR_output.txt

##after it runs, retrieve the fasta sequences like so:
module load ncbi_tools
#blastdbcmd -h

#first, extract just the final iteration
tac class2_GPCR_output.txt | grep 'Results' -m 1 -B 99999 | tac > class2_third_iteration.txt

#grep all lines that start with ">"  (non-redundant seq hits), remove ">" and everything after the first white space
grep '^>' class2_third_iteration.txt | perl -p -e 's/>(\S+)\s+.+/$1/' > ID_lines_class2.txt

#then use the blast command line applications (blastdbcmd) to retrieve the fasta files from the whichever database you used to blast against
blastdbcmd -entry_batch ID_lines_class2.txt -db $NCBI_DB/nr > class2_expanded.fasta