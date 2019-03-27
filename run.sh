#!/bin/bash

sbatch -v slurm-exec-PtDep-noWeight -J toyFlow-PtDep-noWeight -o errors-PtDep-noWeight.txt
sbatch -v slurm-exec-PtDep-Weight -J toyFlow-PtDep-Weight -o errors-PtDep-Weight.txt
sbatch -v slurm-exec-noPtDep-noWeight -J toyFlow-noPtDep-noWeight -o errors-noPtDep-noWeight.txt
sbatch -v slurm-exec-noPtDep-Weight -J toyFlow-noPtDep-Weight -o errors-noPtDep-Weight.txt
