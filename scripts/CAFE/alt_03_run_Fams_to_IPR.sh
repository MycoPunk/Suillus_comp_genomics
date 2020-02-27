#!/bin/bash -l
#SBATCH --mem=25G
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH --time=0-02:00:00     # 2 hrs
#SBATCH --mail-user=lotus.lofgren@ucr.edu
#SBATCH --mail-type=ALL

##print the IPR annotations for each family 

for i in *contractions.clean.txt
do
b=$(basename "$i" .txt)
while IFS= read -r line 
do
grep $line IPR_files.tab | awk '{ print $2 }' > $b.IPR_anno_output.txt
done < $i
done 


for i in *expansions.clean.txt
do
b=$(basename "$i" .txt)
while IFS= read -r line 
do
grep $line IPR_files.tab | awk '{ print $2 }' > $b.IPR_anno_output.txt
done < $i
done 


