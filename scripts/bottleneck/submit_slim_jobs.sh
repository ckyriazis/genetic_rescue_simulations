#! /bin/bash
#$ -wd /u/home/c/ckyriazi/project-klohmuel/genetic_rescue_simulations/scripts/bottleneck
#$ -l h_rt=200:00:00,h_data=48G,highp
#$ -o /u/home/c/ckyriazi/project-klohmuel/genetic_rescue_simulations/output/bottleneck
#$ -e /u/home/c/ckyriazi/project-klohmuel/genetic_rescue_simulations/output/bottleneck
#$ -N bottleneck_10000Na_25Nb_h0.0_3619
#$ -m bea
#$ -t 1:25

# MAKE SURE TO CHANGE NAME OF JOB

source /u/local/Modules/default/init/modules.sh
module load gcc/4.9.3

SLIMDIR=/u/home/c/ckyriazi/project-klohmuel/software/slim_build

${SLIMDIR}/slim /u/home/c/ckyriazi/project-klohmuel/genetic_rescue_simulations/scripts/bottleneck/run_scripts/slim_bottleneck_10000Na_25Nb_h0.0_3619.slim 
