##############################
## Beacon pipeline 
##############################

## creates a beacon friendly vcf 

###############################
## Filter Samples by SNP 
###############################
 
module load tabix/1.2.1
module load bcftools/1.20
module load vt/0.5772

echo "[FILTER]"

bcftools filter -i 'TYPE="snp"' /results/dengen_2211_merged.v2.tags.vcf.gz > /results/dengen.v2.tags.SNP_filtered.vcf.gz

echo "[FILTER]"


###############################
## Decompose Samples
###############################

echo "[DECOMPOSE]"

vt decompose /results/dengen_2211_merged.v2.tags.SNP_filtered.vcf.gz -o /results/dengen_2211_merged.v2.tags.SNP_filtered.decomposed.vcf.gz

echo "[DECOMPOSE]"

################################
## Rename Chromossome 
################################

echo "[RENAME]"

bcftools annotate --rename-chrs chr_name_conv.txt /results/dengen_2211_merged.v2.tags.SNP_filtered.decomposed.vcf.gz -Oz -o /results/dengen_2211_merged.v2.tags.SNP_filtered.decomposed.renamed.vcf.gz

tabix -p vcf /results/dengen_2211_merged.v2.tags.SNP_filtered.decomposed.renamed.vcf.gz

echo "[RENAME]"

