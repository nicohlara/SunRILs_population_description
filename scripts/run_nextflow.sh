#!/bin/bash
#SBATCH --account=guedira_seq_map
#SBATCH --time=24:00:00  # walltime limit (HH:MM:SS)
#SBATCH --partition=atlas
#SBATCH --nodes=1
#SBATCH --ntasks=48  
#SBATCH --job-name="SunRILs_nextflow"
#SBATCH --mail-user=nalara@ncsu.edu  # email address
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL

module load miniconda3
module load nextflow
conda activate /home/nicolas.lara/.conda/envs/imputation

cd /project/guedira_seq_map/nico/SunRILs_population_description

export _JAVA_OPTIONS="-Xmx350G"

nextflow scripts/00_vcf_preprocess_combine.nf
