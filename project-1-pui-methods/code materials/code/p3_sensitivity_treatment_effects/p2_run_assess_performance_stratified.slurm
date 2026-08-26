#!/bin/bash --login
#SBATCH -p multicore
#SBATCH -t 12:00:00
#SBATCH -a 1-72
#SBATCH -n 4

module purge
module load apps/gcc/R/4.4.2

### Map task ID (1-64) to parameter combinations

genders=(1 2)
meds=(statins antihypertensives)
group_types=(0 1)
scenarios=(1 2 3 4 5 6 7 8 9)

### Convert SLURM_ARRAY_TASK_ID (1-based) to 0-based index
idx=$((SLURM_ARRAY_TASK_ID - 1))

### Extract each parameter from the index
n_med=${#meds[@]}          # 2
n_group_types=${#group_types[@]}          # 2
n_scenario=${#scenarios[@]}       # 8

med_idx=$((idx % n_med))
idx=$((idx / n_med))

group_type_idx=$((idx % n_group_types))
idx=$((idx / n_group_types))

scenario_idx=$((idx % n_scenario))
idx=$((idx / n_scenario))

gender_idx=$((idx))

### Assign parameter values
gender=${genders[$gender_idx]}
med=${meds[$med_idx]}
group_type=${group_types[group_type_idx]}
scenario=${scenarios[$scenario_idx]}

### Define output file
outfile="p2.2_assess_performance_stratified_treatment_${gender}_${med}${group_type}_${scenario}.out"

echo "Running: gender=${gender} med=${med} group_type=${group_type} scenario=${scenario}"
echo "Output:  ${outfile}"

Rscript p2.2_assess_performance_stratified_treatment.R \
  ${gender} ${med} ${group_type} ${scenario} > ${outfile}
