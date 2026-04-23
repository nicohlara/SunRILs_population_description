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
conda activate /home/nicolas.lara/.conda/envs/imputation

ref_gen_dir=/90daydata/guedira_seq_map/nico/Hilliard_ref
#mkdir ${ref_gen_dir}

cd ${ref_gen_dir}

#cp /project/guedira_seq_map/mwillman/Genome_Assembly_MW/Taes_Hilliard_1.2.fasta .

#samtools faidx Taes_Hilliard_1.2.fasta
bwa index Taes_Hilliard_1.2.fasta
