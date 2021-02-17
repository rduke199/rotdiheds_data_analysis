#!/bin/bash
#SBATCH -t 14-00:00:00   			      #Time for the job to run
#SBATCH --job-name=RunRotDiheds           #Name of the job
#SBATCH -N 1 	        			       	#Number of nodes required
#SBATCH -n 8				                #Number of cores needed for the job
#SBATCH -p SAN16M64_L 
#SBATCH --account=col_cmri235_uksr  #Name of account to run under

#SBATCH --mail-type ALL			      	#Send email on start/end
#SBATCH --mail-user rdu230@uky.edu	#Where to send email
#SBATCH --error=SLURM_RunRotDiheds_%j.err		#Name of error file
#SBATCH --output=SLURM_RunRotDiheds_%j.out 	#Name of output file

#module load ccs/gaussian/g16-A.03/g16-sandybridge
ulimit -u 32768
export OPENBLAS_NUM_THREADS=1
echo "Job $SLURM_JOB_ID running on SLURM NODELIST: $SLURM_NODELIST "

. ~/.bashrc
source $(conda info --base)/etc/profile.d/conda.sh
conda activate myenv

python rotDiheds2.py

