#!/bin/bash -l
#SBATCH -N 1 -n 1
#SBATCH --mem=25G
#SBATCH --time=0-01:00:00     # 1 hr
#SBATCH --mail-user=lotus.lofgren@ucr.edu
#SBATCH --mail-type=ALL

#load modules
module load R 
module load hmmer
module load clustalo

#initialize the output csv with headers
echo "gene_name	protein_is	call_is	false_pos	false_neg	correct_call" > GPCR_itteration_totals.csv

#set the number of itterations
for x in {1..3} ; do
Rscript s1_leave_one_out.R
wait

#align the reference -1 with clustal omega

clustalo -i minus_one.fasta -o minus_one.fasaln

#build the profile HMM in the form:
#hmmbuild [-options] hmmout.file align.file
hmmbuild minus_one_profile.hmm minus_one.fasaln

#compress and index the flatfile with hmmpress
hmmpress minus_one_profile.hmm

#search the left out sequence against the profile HMM in the form:
#hmmscan --tblout output_file P-fam_database your_protein_sequences
hmmscan --tblout hmm_output_file minus_one_profile.hmm one_out.fasta

Rscript s2_analyze_output.R

#clean up intermediate files
rm slurm*
rm minus_one_profile.hmm*
rm *.fasaln
rm one_out.fasta
rm minus_one.fasta
rm hmm_output_file
done
