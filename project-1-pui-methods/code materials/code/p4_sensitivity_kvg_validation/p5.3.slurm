#!/bin/bash --login
#SBATCH -p multicore
#SBATCH -t 24:00:00
#SBATCH -a 1-32
#SBATCH -n 7

module purge
module load apps/gcc/R/4.4.2

### Map task ID (1-32) to parameter combinations
### Order: gender > med > t_fup > model_flex
### 2 genders x 2 meds x 4 t_fups x 2 flexes = 32

genders=(1 2)
meds=(statins antihypertensives)
t_fups=(1 3 5 10)
flexes=(flex linear)

### Convert SLURM_ARRAY_TASK_ID (1-based) to 0-based index
idx=$((SLURM_ARRAY_TASK_ID - 1))

### Extract each parameter from the index
n_flex=${#flexes[@]}       # 2
n_t=${#t_fups[@]}          # 4
n_med=${#meds[@]}          # 2

flex_idx=$((idx % n_flex))
idx=$((idx / n_flex))

t_idx=$((idx % n_t))
idx=$((idx / n_t))

med_idx=$((idx % n_med))
idx=$((idx / n_med))

gender_idx=$((idx))

### Assign parameter values
gender=${genders[$gender_idx]}
med=${meds[$med_idx]}
t_fup=${t_fups[$t_idx]}
flex=${flexes[$flex_idx]}

### Define output file
outfile="p5.1_estimate_intervention_effects_create_cohorts_${med}${gender}_t${t_fup}_${flex}.out"

echo "Running: gender=${gender} med=${med} t_fup=${t_fup} flex=${flex}"
echo "Output:  ${outfile}"

Rscript p5.3_estimate_intervention_effects_and_diagnostics.R \
  ${gender} ${med} ${t_fup} ${flex} > p5.3_estimate_intervention_effects_and_diagnostics_${med}${gender}_t${t_fup}_${flex}.out