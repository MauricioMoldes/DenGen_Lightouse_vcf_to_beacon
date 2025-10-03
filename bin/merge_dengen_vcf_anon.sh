
#################################
## Merge DenGen VCFs into project VCF 
#################################

module load bcftools/1.20

bcftools merge -l /tmp/dengen_normalized.paths --missing-to-ref -m none -Oz -o dengen_2211_merged.v2.vcf.gz

