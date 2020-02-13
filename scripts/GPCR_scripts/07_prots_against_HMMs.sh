#!/bin/bash -l
#SBATCH -N 1 -n 1
#SBATCH --mem=25G
#SBATCH --time=0-03:00:00     # 3 hr
#SBATCH --mail-user=lotus.lofgren@ucr.edu
#SBATCH --mail-type=ALL


##This script searches the proteomes of each species of interest against each the HMM model for each GPCR class.

#load modules
module load hmmer
module load clustalo
module load parallel

#concatenate the HMM files 
#cat *7TM_only.hmm > gpcr_models_7TM_only.hmm
cat *whole.hmm > gpcr_models_whole.hmm


#compress and index the concatenated flatfile with hmmpress
#hmmpress gpcr_models_7TM_only.hmm
hmmpress gpcr_models_whole.hmm

#search each species against each HMM model 
#parallel -j 4 hmmscan --cpu 4 --domtblout {.}.domtblout gpcr_models_7TM_only.hmm {} \> {.}.hmmscan ::: $(ls ./seqs/proteins/*.trimmed.fasta)
parallel -j 4 hmmscan --cpu 4 --domtblout {.}.domtblout gpcr_models_whole.hmm {} \> {.}.hmmscan ::: $(ls ./seqs/proteins/*.trimmed.fasta)
