#!/bin/bash
#SBATCH --job-name="bcf_filtering"
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


##UNAPLIED_FILTERS:
##line_het=0.05


##BCFTOOLS:
#biallelic SNP
#unaligned chromosomes
#MAF=0.0012
#depth=1
#var_miss=0.2

##R/gaston:
#var_het=0.4
#line_miss=0.7
#rename lines
#aggregate duplicated lines

module load miniconda3
conda activate imp_2

DIR=/90daydata/guedira_seq_map/nico/SunRILs/raw_VCF

cd $DIR
VCF_IN=SunRILs_production.vcf.gz
VCF_OUT=SunRILs_prod_filt.vcf.gz

##filter by read depth, minimum MAF, and maximum missingness
##filter by maf
##smallest biparental has 144 RILs, total pop has 3.8k.
##Took 144/3800*0.03125 (smallest rare allele likelihood in this pop) = threshold
bcftools view -i 'FORMAT/DP>1 && MAF>0.0012 && F_MISSING<0.1' ${VCF_IN} -Oz -o  DP_MAF_MISS_filt.vcf.gz

##filter to biallelic SNP                                                                                                           
bcftools view -m2 -M2 -v snps DP_MAF_MISS_filt.vcf.gz -Oz -o  biallelic.vcf.gz

##remove unaligned reads
bcftools view -t "^UNKNOWN" biallelic.vcf.gz -Oz -o  ${VCF_OUT}


vcf_files=("${VCF_IN}" "DP_MAF_MISS_filt.vcf.gz" "biallelic.vcf.gz" "${VCF_OUT}")  # Add your VCF file names here

# Create an output table with headers
echo -e "VCF File\tTotal Markers" > comparison_table.txt

# Loop through each VCF file and get the total number of markers
for vcf in "${vcf_files[@]}"
	do
	# Get the total number of records from bcftools stats
	total_markers=$(bcftools stats "$vcf" | grep "number of records" | awk '{print $6}')
			    
	# Append the file name and total markers count to the comparison table
	echo -e "${vcf}\t${total_markers}" >> comparison_table.txt
done
