#! /bin/bash
#$ -wd /u/home/c/ckyriazi/project-klohmuel/genetic_rescue_simulations/scripts/bottleneck
#$ -l h_rt=312:00:00,h_data=96G,highp
#$ -o /u/home/c/ckyriazi/project-klohmuel/genetic_rescue_simulations/output/bottleneck
#$ -e /u/home/c/ckyriazi/project-klohmuel/genetic_rescue_simulations/output/bottleneck
#$ -N bottleneck_15000Na_100Nb_100nF_h0.0_4319
#$ -m bea
#$ -t 1:10

# MAKE SURE TO CHANGE NAME OF JOB

source /u/local/Modules/default/init/modules.sh
module load gcc/4.9.3

SLIMDIR=/u/home/c/ckyriazi/project-klohmuel/software/slim_build

${SLIMDIR}/slim /u/home/c/ckyriazi/project-klohmuel/genetic_rescue_simulations/scripts/bottleneck/slim_bottleneck_15000Na_100Nb_100nF_h0.0_4319.slim
