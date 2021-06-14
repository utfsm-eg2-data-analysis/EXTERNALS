#!/bin/sh                                                   
echo 'starting'

sbatch --export=opt="Fe_liq_RD" -J "FeD2" job.slurm
sbatch --export=opt="Pb_liq_RD" -J "PbD2" job.slurm
sbatch --export=opt="Ca_liq_RD" -J "CaD2" job.slurm
sbatch --export=opt="Fe_sol_RD" -J "Fe" job.slurm
sbatch --export=opt="Pb_sol_RD" -J "Pb" job.slurm
sbatch --export=opt="Ca_sol_RD" -J "Ca" job.slurm

echo 'done'

