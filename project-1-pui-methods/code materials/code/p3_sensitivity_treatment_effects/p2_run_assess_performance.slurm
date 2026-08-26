#!/bin/bash --login
#SBATCH -p multicore
#SBATCH -t 12:00:00
#SBATCH -a 1-16
#SBATCH -n 4

module purge
module load apps/gcc/R/4.4.2

### Map task ID (1-16) to parameter combinations

genders=(1 2)
scenarios=(1 2 3 4 5 6 7 8)

### Convert SLURM_ARRAY_TASK_ID (1-based) to 0-based index
idx=$((SLURM_ARRAY_TASK_ID - 1))

### Extract each parameter from the index
n_scenario=${#scenarios[@]}       # 8

scenario_idx=$((idx % n_scenario))
idx=$((idx / n_scenario))

gender_idx=$((idx))

### Assign parameter values
gender=${genders[$gender_idx]}
scenario=${scenarios[$scenario_idx]}

### Define output file
outfile="p2.1_assess_performance_${gender}_${scenario}.out"

echo "Running: gender=${gender} scenario=${scenario}"
echo "Output:  ${outfile}"

Rscript p2.1_assess_performance.R \
  ${gender} ${scenario} > ${outfile}
