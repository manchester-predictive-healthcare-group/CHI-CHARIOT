#!/bin/bash --login
#SBATCH -p himem
#SBATCH -t 7-0
#SBATCH -a 1-40
#SBATCH -n 7

module purge
module load apps/gcc/R/4.5.0

### Map task ID (1-40) to parameter combinations

genders=(1 2)
splits=(devel valid)
chains=(1 2 3 4 5 6 7 8 9 10)

### Convert SLURM_ARRAY_TASK_ID (1-based) to 0-based index
idx=$((SLURM_ARRAY_TASK_ID - 1))

### Extract each parameter from the index
n_split=${#splits[@]}       # 2
n_chain=${#chains[@]}       # 10

split_idx=$((idx % n_split))
idx=$((idx / n_split))

chain_idx=$((idx % n_chain))
idx=$((idx / n_chain))

gender_idx=$((idx))

### Assign parameter values
gender=${genders[$gender_idx]}
split=${splits[$split_idx]}
chain=${chains[$chain_idx]}

### Define output file
outfile="p4_impute_cohort_prototype3_split_${split}_${gender}_${chain}.out"

echo "Running: gender=${gender} split=${split} chain=${chain}"
echo "Output:  ${outfile}"

### Run script
Rscript p4_impute_cohort_prototype3_split.R ${chain} ${gender} ${split} > ${outfile}
