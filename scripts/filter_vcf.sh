



module load miniconda3
conda activate imp_2

VCF_IN=a.vcf
VCF_OUT=b.vcf

var_miss=0.2
line_miss=0.7
var_het=0.4
line_het=0.05
MAF=0.0012
depth=1




##remove unknown chromosome
bcftools -t "^Unknown" $VCF_IN > vcf1.vcf




##filter by maf
##smallest biparental has 144 RILs, total pop has 3.8k. 
##Took 144/3800*0.03125 (smallest rare allele likelihood in this pop) = threshold
#bcftools view -i 'MAF>0.0012' vcf.vcf > vcf.vcf


bcftools view \
  --max-missing $var_miss \  # Filter variants by missing data
  --exclude 'F_MISSING > $line_miss' \  # Filter individuals (samples) by missing data
  --exclude 'F_HET > $var_het' \  # Filter variants by heterozygous calls
  --exclude 'F_HET < $line_het' \  # Filter individuals by heterozygous calls
  -t "^Unknown" \ #filter out unknown chromosome
  --min-af $MAF \  # Filter variants by minor allele frequency
  --minDP $depth \  # Filter based on sequencing depth
  -Oz -o ${VCF_OUT} ${VCF_IN}  # Output the filtered VCF in compressed format
