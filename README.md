# Metagenomics
Repository for work on the diet partitioning study at Portal

## File Structure

### Data

#### CollectionData

Files:
* _vial_id.csv_: contains data on the vial ID code, the sample ID code (SKME_### for plants, PIT tag for rodents), and whether the sample was plant material or fecal 
* _plant_voucher_collection.csv_: data on plant vouchers collected at the Portal Project to be identified and barcoded
* _fecal_sample_collection.csv_: data for each fecal sample collected. Contains species, sex, PIT tag, sample type, and whether it was from an experiment (trap_and_bait) or not
* _Portal_plant_species.csv_: current plant list from Portal at the time. Likely needs to be updated

#### SequencedData

Insects: we ran a handful of samples to see if any insects popped out in the fecal samples, but they did not

Plants: contains both data on the DNA sequences for both plant vouchers and fecal samples
