
###################
## Inserts DenGen Cohort into Mongo DB 
###################

docker cp analyses.json mongoprod:tmp/analyses.json
docker cp biosamples.json mongoprod:tmp/biosamples.json
docker cp cohorts.json mongoprod:tmp/cohorts.json
docker cp datasets.json mongoprod:tmp/datasets.json
docker cp individuals.json mongoprod:tmp/individuals.json
docker cp runs.json mongoprod:tmp/runs.json

docker exec mongoprod mongoimport --jsonArray --uri "mongodb://root:example@127.0.0.1:27017/beacon?authSource=admin" --file /tmp/datasets.json --collection datasets
docker exec mongoprod mongoimport --jsonArray --uri "mongodb://root:example@127.0.0.1:27017/beacon?authSource=admin" --file /tmp/analyses.json --collection analyses
docker exec mongoprod mongoimport --jsonArray --uri "mongodb://root:example@127.0.0.1:27017/beacon?authSource=admin" --file /tmp/biosamples.json --collection biosamples
docker exec mongoprod mongoimport --jsonArray --uri "mongodb://root:example@127.0.0.1:27017/beacon?authSource=admin" --file /tmp/cohorts.json --collection cohorts
docker exec mongoprod mongoimport --jsonArray --uri "mongodb://root:example@127.0.0.1:27017/beacon?authSource=admin" --file /tmp/individuals.json --collection individuals
docker exec mongoprod mongoimport --jsonArray --uri "mongodb://root:example@127.0.0.1:27017/beacon?authSource=admin" --file /tmp/runs.json --collection runs
docker exec mongoprod mongoimport --jsonArray --uri "mongodb://root:example@127.0.0.1:27017/beacon?authSource=admin" --file /tmp/genomicVariations.json --collection genomicVariations


#########################
## Re-index
##########################

docker exec beaconprod python /beacon/connections/mongo/reindex.py

############################
## Extract filtering terms 
############################

docker exec beaconprod python beacon/connections/mongo/extract_filtering_terms.py

##############################
## Get descendents 
##############################

docker exec beaconprod python beacon/connections/mongo/get_descendants.py
