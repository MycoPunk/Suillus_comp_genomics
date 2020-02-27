#!/bin/bash -l
#SBATCH --mem=25G
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH --time=0-02:00:00     # 2 hrs
#SBATCH --mail-user=lotus.lofgren@ucr.edu
#SBATCH --mail-type=ALL

#split the files
awk '{for(i=2;i<=NF;i++){printf "%s ", $i >> $1"_expansions.txt"};printf "\n" >> $1"_expansions.txt"; close($1"_expansions.txt")}' genes_in_expanded_fams.tab
awk '{for(i=2;i<=NF;i++){printf "%s ", $i >> $1"_contractions.txt"};printf "\n" >> $1"_contractions.txt"; close($1"_contractions.txt")}' genes_in_contracted_fams.tab

#clean up the files - grep for  white space, trans to new line and remove NAs'
for i in *.txt
do
b=$(basename "$i" .txt)
sed 's/\tNA//g' $i > $b.clean.txt
sed -i 's/\t/\n/g' $b.clean.txt 
sed -i '/^$/d' $b.clean.txt
done


##print the GO terms for each family 
#expansions
for i in *contractions.clean.txt
do
b=$(basename "$i" .txt)
while IFS= read -r line 
do
grep $line GO_files.tab | awk '{ print $2 }' > $b.GO_terms_output.txt
done < $i
done 


for i in *expansions.clean.txt
do
b=$(basename "$i" .txt)
while IFS= read -r line 
do
grep $line GO_files.tab | awk '{ print $2 }' > $b.GO_terms_output.txt
done < $i
done 
