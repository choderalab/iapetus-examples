#!/usr/bin/env bash
# Set walltime limit
#BSUB -W 48:00
#
# Set output file
#BSUB -o  iapetus.%J.%I.log
#
# Specify node group
#BSUB -m "ls-gpu lt-gpu lp-gpu lg-gpu"
#BSUB -q gpuqueue
#
# nodes: number of nodes and GPU request
#BSUB -n 1 -gpu "num=1:j_exclusive=yes:mode=shared" -R "rusage[mem=12]"
#
# job name (default = name of script file)
#BSUB -J "iapetus[1-10]"

jobnames="unused arg comp7 comp7_nowat comp8 comp8_nowat glu his_neut his_posit imi mero"
arr=($jobnames)
jobname=${arr[$LSB_JOBINDEX]}
echo $jobname
iapetus --gromacs $jobname --ligseq 423 --output ${jobname}.nc --verbose
