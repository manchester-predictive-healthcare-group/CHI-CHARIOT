#!/bin/bash --login
#SBATCH -p himem
#SBATCH -t 1-06:00:00
#SBATCH -a 1-2000
#SBATCH -n 1
#SBATCH --mem=250G

module purge
module load apps/gcc/R/4.5.0

### Map task ID to parameter combinations

genders=(1 2)
iters=({1..10000})

### Convert SLURM_ARRAY_TASK_ID (1-based) to 0-based index
idx=$((SLURM_ARRAY_TASK_ID - 1))

### Extract each parameter from the index
n_gender=${#genders[@]}

gender_idx=$((idx % n_gender))
idx=$((idx / n_gender))

iter_idx=$((idx))

### Assign parameter values
gender=${genders[$gender_idx]}
iter=${iters[$iter_idx]}

### Define output file
outfile="p2_optimism_loop_${gender}_${iter}.out"

echo "Running: gender=${gender} iter=${iter}"
echo "Output:  ${outfile}"

### Run script
Rscript p2_optimism_loop.R ${gender} ${iter} > ${outfile}
