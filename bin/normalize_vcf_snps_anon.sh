###########################################
## Normalize vcf 
###########################################

module load bcftools/1.20

# Normalize each DenGen sample 

#!/bin/bash
input="dengen_2211_list.paths"
while IFS= read -r line
do
  echo "$line"
  filename=$(basename -- "$line")
  extension="${filename##*.}"
  filename="${filename%.*.*}"
   
	bcftools norm --fasta-ref /ngc/shared/resources/h.sapiens/hg38/genomes/GCA_000001405.15_GRCh38_no_alt_analysis_set/20210411/GCA_000001405.15_GRCh38_no_alt_analysis_set_maskedGRC_exclusions.fasta --multiallelics - --threads 8 -Oz -o /results/normalization/"$filename".normalized.vcf.gz "$line"


done < "$input"




