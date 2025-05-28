#!/bin/bash

missed_jobs=(
156644 156657
)

for job_id in "${missed_jobs[@]}"; do
  echo "Submitting job for index $job_id"
  qsub -v PBS_ARRAY_INDEX=$job_id /srv/scratch/z5394590/pglmm/set1/pbs/run_pglmm_sim_missed.pbs
done
