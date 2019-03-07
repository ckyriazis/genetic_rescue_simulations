#! /bin/bash
#$ -wd /u/home/c/ckyriazi/project-klohmuel/genetic_rescue_simulations/scripts/bottleneck
#$ -l h_rt=32:00:00,h_data=32G,highp
#$ -o /u/home/c/ckyriazi/project-klohmuel/genetic_rescue_simulations/output
#$ -e /u/home/c/ckyriazi/project-klohmuel/genetic_rescue_simulations/output
#$ -N bottleneck_1000Na_50Nb_h0.0_3619
#$ -m bea
#$ -t 1:10

# MAKE SURE TO CHANGE NAME OF JOB

source /u/local/Modules/default/init/modules.sh
module load gcc/4.9.3

SLIMDIR=/u/home/c/ckyriazi/project-klohmuel/software/slim_build

${SLIMDIR}/slim /u/home/c/ckyriazi/project-klohmuel/genetic_rescue_simulations/scripts/bottleneck/slim_bottleneck_1000Na_50Nb_h0.0_3619.slim

