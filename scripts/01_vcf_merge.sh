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
conda activate /home/nicolas.lara/.conda/envs/imp_2

cd /project/guedira_seq_map/nico/SunRILs_population_description

export _JAVA_OPTIONS="-Xmx350G"

cd data/processed_vcf
# find -name "*_imp.vcf.gz" > list.txt
# cat list.txt | while read file; do
#   bcftools index "${file/input_file/sorted_input_file}"
#   bcftools sort $file -Oz -o "${file/input_file/sorted_input_file}"
# done
# bcftools merge *_imp.vcf.gz --force-samples -Oz -o SunRILs_imp.vcf.g



beagle gt=SunRILs_imp.vcf.gz out=SunRILs_mergeimp map=../../linkage_map/monotonic/consensus_GBS_monotonic.map nthreads=40 window=350
