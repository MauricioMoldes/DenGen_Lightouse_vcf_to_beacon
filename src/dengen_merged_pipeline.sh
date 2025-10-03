
################
## DenGen Beacon Project VCF pipeline 
################

# Executed as an array for the possibility of a project like vcf, or multiple samples 

declare -a arr=("dengen_2211_merged.v2.tags")

# now loop through the above array
for sample in "${arr[@]}"
do

echo "$sample"

###########################
## Filter Sample by SNP
###########################

echo [File checker : Filtering]


if [ ! -f "$sample".SNP_filtered.vcf.gz ]; then
    
	echo "[Filtering Start]"

	bcftools filter -i 'TYPE="snp"' "$sample".vcf.gz > "$sample".SNP_filtered.vcf.gz

	echo "[Filtering End"]

fi

###########################
## Decompose using VT
###########################

echo [File checker : Decompose]

	if [ ! -f "$sample".SNP_filtered.decomposed.vcf.gz ]; then


	echo "[Decompose Start]"

	vt decompose "$sample".SNP_filtered.vcf.gz -o "$sample".SNP_filtered.decomposed.vcf.gz

	echo "[Decompose End]"

	fi 

###########################
## Rename Chromossome
###########################

echo [File checker : Rename ]

	if [ ! -f "$sample".SNP_filtered.decomposed.renamed.vcf.gz ]; then


	echo "[Rename Start]"


	bcftools annotate --rename-chrs chr_name_conv.txt "$sample".SNP_filtered.decomposed.vcf.gz -Oz -o "$sample".SNP_filtered.decomposed.renamed.vcf.gz

        tabix -p vcf "$sample".SNP_filtered.decomposed.renamed.vcf.gz


	echo "[Rename End]"

	fi


###########################
## Add to ri-tools container
###########################

cp "$sample".SNP_filtered.decomposed.renamed.vcf.gz /dengen_beacon/beacon2-pi-api/ri-tools/files/vcf/files_to_read


###########################
## Run ri-tools container
###########################

cd /dengen_beacon/beacon2-pi-api/

docker compose restart beacon-ri-tools 
docker exec -it ri-tools python genomicVariations_vcf.py


#####################
## Remove VCF from beacon ri tools mountpoint
#####################


cd /dengen_beacon/beacon2-pi-api/ri-tools/files/vcf/files_to_read
rm "$sample".SNP_filtered.decomposed.renamed.vcf.gz

####################
## Clean up on each sample intermediary files
####################

cd /dengen_beacon/beacon2-pi-api/ri-tools/files/vcf/data

rm *.SNP_filtered.* 
rm *.head
rm *.variants


done
