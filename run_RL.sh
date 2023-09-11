#!/usr/bin/env bash
#SBATCH -n 10
#SBATCH --gres=gpu:1
#SBATCH --mem-per-cpu=2048
#SBATCH --time=4-00:00:00
#SBATCH --mail-type=END
#SBATCH --nodelist gnode039

#eval "$(conda shell.bash hook)"
#conda activate peptides

python3 RL_run.py

# while true; do continue; done;
