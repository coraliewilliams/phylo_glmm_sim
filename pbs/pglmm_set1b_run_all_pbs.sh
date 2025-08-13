#!/usr/bin/env bash
# submit_all_pbs.sh
# Generate and qsub 20 PBS array jobs to cover 1–24,000 in chunks of 24,000

# (1) adjust this to your actual working directory:
WORKDIR="/srv/scratch/z5394590/pglmm/set1/"

# (2) create a place to store all the little PBS scripts
PBS_DIR="${WORKDIR}/pbs_scripts"
mkdir -p "${PBS_DIR}"

# (3) for each of the 3 blocks:
for (( chunk=0; chunk<3; chunk++ )); do
  start=$(( chunk*10000 + 1 ))
  end=$(( start + 10000 - 1 ))
  jobname="pglmm_set1b_${start}_${end}"
  script="${PBS_DIR}/${jobname}.pbs"

  cat > "${script}" <<EOF
#!/bin/bash
#PBS -N ${jobname}
#PBS -l select=1:ncpus=4:mem=6gb
#PBS -l walltime=00:30:00
#PBS -J ${start}-${end}

cd ${WORKDIR}

module purge
module load r/4.3.1
module load udunits/2.2.28
module load geos/3.9.1
module load udunits/2.2.28
module load gdal/3.5.3
module load sqlite/3.39.4
module load openssl/1.1.1s
module load proj/8.2.1
module load gsl/2.7.1

Rscript pglmm_sim_set1b.R
EOF

  # submit it
  qsub "${script}"
done