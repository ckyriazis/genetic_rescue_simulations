#! /bin/bash
#$ -wd /u/home/c/ckyriazi/project-klohmuel/genetic_rescue_simulations/scripts/rescue
#$ -l h_rt=140:00:00,h_data=48G,highp
#$ -o /u/home/c/ckyriazi/project-klohmuel/genetic_rescue_simulations/output/rescue
#$ -e /u/home/c/ckyriazi/project-klohmuel/genetic_rescue_simulations/output/rescue
#$ -N rescue_10000Ns_50Ts_25Nb_5nM_1nR_h0.0_het_5319
#$ -m bea
#$ -t 1:25


# MAKE SURE TO CHANGE NAME OF JOB

source /u/local/Modules/default/init/modules.sh
module load gcc/4.9.3

SLIMDIR=/u/home/c/ckyriazi/project-klohmuel/software/slim_build

${SLIMDIR}/slim /u/home/c/ckyriazi/project-klohmuel/genetic_rescue_simulations/scripts/rescue/slim_rescue_10000Ns_50Ts_25Nb_5nM_1nR_h0.0_het_5319.slim
