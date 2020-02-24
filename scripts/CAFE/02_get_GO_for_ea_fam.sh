#!/bin/bash -l
#SBATCH --mem=25G
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH --time=0-02:00:00     # 2 hrs
#SBATCH --mail-user=lotus.lofgren@ucr.edu
#SBATCH --mail-type=ALL

#split the files
awk '{filename = sprintf("expansions_%d.txt", NR); print >filename; close(filename)}' genes_in_expanded_fams.tab
awk '{filename = sprintf("contractions_%d.txt", NR); print >filename; close(filename)}' genes_in_contracted_fams.tab
 
#clean up the files - grep for  white space, trans to new line and remove NAs'
for i in *.txt
do
b=$(basename "$i" .txt)
sed 's/\tNA//g' $i > test.$b.txt
sed -i 's/\t/\n/g' test.$b.txt 
sed -i '/^$/d' test.$b.txt
done

##print the GO terms for each family 
#expansions
for i in test.contractions*
do
b=$(basename "$i" .txt)
while IFS= read -r line 
do
grep $line GO_files.tab | awk '{ print $2 }' > GO_terms_output.$b.txt
done < $i
done 

#contractions
for h in test.expansions*
do
b=$(basename "$h" .txt)
while IFS= read -r line 
do
grep $line GO_files.tab | awk '{ print $2 }' > GO_terms_output.$b.txt
done < $h
done 
