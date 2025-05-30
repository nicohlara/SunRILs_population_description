#!/bin/bash
#SBATCH --job-name="beagle_impute"                   
#SBATCH --qos=normal
#SBATCH -p atlas                                 
#SBATCH -A guedira_seq_map                       
#SBATCH -N 1                                      
#SBATCH -n 48                                     
#SBATCH -t 7-00:00:00                             
#SBATCH --mail-user=nalara@ncsu.edu                                      
#SBATCH --mail-type=END                           
#SBATCH --mail-type=FAIL                          
#SBATCH -o "stdout.%x.%j.%N"                      
#SBATCH -e "stderr.%x.%j.%N"  

##impute using beagle5
##adapted from Z. Winn's atlas pipeline
##Nico Lara, 2025-3-14

module load beagle


cd /90daydata/guedira_seq_map/nico/SunFilt/output/

echo "Begin imputation: $(date +"%Y-%m-%d %H:%M:%S")"
#bgzip SunRILs_prod_filt2.vcf

beagle gt=SunRILs_prod_filt2.vcf.gz
            out=SunRILs_prod_filt_imp
	    map=/90daydata/guedira_seq_map/nico/SunRILs/HPC-GBS-Pipeline/SynOp_RIL906_v1.0_GBS_monotonic.map \
            nthreads=40 \
            window=350
