#!/bin/bash

base_dir=/mnt/c/Users/nalara/Documents/GitHub/SunRILs_population/data
RAW_VCF=${base_dir}/processed_vcf_20250815/SunRILs_raw.vcf.gz
pop_table=${base_dir}/cross_info.csv
output_dir=${base_dir}/biparental_vcf

##basic universal filtering
#bcftools view -i 'FORMAT/DP > 3 && FORMAT/DP < 100 && 'FORMAT/GQ' >= 20' ${RAW_VCF} -Oz -o DP_filter.vcf.gz
#bcftools view -m2 -M2 -v snps DP_filter.vcf.gz -Oz -o biallelic.vcf.gz
#bcftools index -c biallelic.vcf.gz
#bcftools view -t "^UNKNOWN" biallelic.vcf.gz -Oz -o "filtered.vcf.gz"

## Index vcf
#bcftools index filtered.vcf.gz

## Extract sample names from the VCF file header
#bcftools query -l filtered.vcf.gz > all_samples.txt
grep '^UX[0-9]\{4\}-' all_samples.txt | sed 's/-.*//' | sort -u > cross_ids_present.txt
awk -F',' 'NR==FNR {keep[$1]; next} $1 in keep' cross_ids_present.txt ${pop_table} > subset_pop_table.txt


#mkdir ${output_dir}

##set internal field separator to comma
IFS=","
# For each population, subset out the population samples and parents
# Filter after subsetting
while read Cross_ID Parent_1 Parent_2; do
    grep "^${Cross_ID}-" all_samples.txt > temp_list.txt
    echo "${Parent_1}" >> temp_list.txt
    echo "${Parent_2}" >> temp_list.txt

    grep -Fxf temp_list.txt all_samples.txt > "${Cross_ID}_list.txt"

    bcftools view -S ${Cross_ID}_list.txt filtered.vcf.gz -Oz -o SunRILs_${Cross_ID}_subset.vcf.gz
    bcftools view -e "MAF < 0.3 || F_MISSING > 0.5" SunRILs_${Cross_ID}_subset.vcf.gz -Oz -o ${output_dir}/SunRILs_${Cross_ID}_filt.vcf.gz
done < subset_pop_table.txt

ls -l *.vcf.gz || echo "No matching VCFs created."

#rm DP_filter.vcf.gz
#rm biallelic.vcf.gz*
#rm filtered.vcf.gz*
#rm all_samples.txt
#rm cross_ids_present.txt
#rm subset_pop_table.txt
#rm temp_list.txt

