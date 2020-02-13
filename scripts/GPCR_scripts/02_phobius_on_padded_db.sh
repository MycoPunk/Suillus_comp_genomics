#!/bin/bash -l
#SBATCH -N 1 -n 16 -p short
#SBATCH --mem=25G
#SBATCH --time=2-0:00:00     # 2 days
#SBATCH --mail-user=lotuslofgren@gmail.com
#SBATCH --mail-type=ALL

CPU=$SLURM_CPUS_ON_NODE
if [ -z $CPU ]; then
 CPU=1
fi
module load parallel
module load phobius
run_phobius () {
 name=$1
 base=$(basename "$name" .expanded.wo.stops.fasta) 
#if [ ! -s ${base}_phobius.tab ]; then
        phobius -short $name > ${base}_phobius.tab
 #fi
}
export -f run_phobius

time parallel -j $CPU run_phobius ::: $(ls *.expanded.wo.stops.fasta)
