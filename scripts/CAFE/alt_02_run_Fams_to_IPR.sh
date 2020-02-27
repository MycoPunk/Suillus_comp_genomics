#!/bin/bash -l
#SBATCH --mem=25G
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH --time=0-02:00:00     # 2 hrs
#SBATCH --mail-user=lotus.lofgren@ucr.edu
#SBATCH --mail-type=ALL

for i in test.contractions*
do
b=$(basename "$i" .txt)
while IFS= read -r line 
do
grep $line IPR_files.tab | awk '{ print $2 }' > IPR_anno_output.$b.txt
done < $i
done 

for h in test.expansions*
do
b=$(basename "$h" .txt)
while IFS= read -r line 
do
grep $line IPR_files.tab | awk '{ print $2 }' > IPR_anno_output.$b.txt
done < $h
done 