# Accurate Determination of Tumor Sample-Patient Sample Pairs for T-Cell Immunotherapy Quality Assurance
***Created by Alexandre R. Sathler for the Earle A. Chiles Research Institute***


## Project Background
This project and its code draws from genome fingerprinting and fingerprint comparison methods outlined in 'Ultrafast Comparison of Personal Genomes via Precomputed Genome Fingerprints,' published in *Frontiers in Genetics* ([Glusman, 2017](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5623000/)). The paper posits that fingerprints can be created from any genome by comparing the genome to the average genome sequence of a population. When making this comparison, Single Nucleotide Variants (SNVs) at every locus become apparent, and can be stored in orderly fashion in the Variant Call Format (`.vcf`). Because every genome will have unique SNVs when compared to a population genome, these SNVs can be used to create a unique genome fingerprint. When compared, fingerprints can be used to identify individuals and determine the similarity between two genomes much more quickly than full genome analysis.

Glusman also published Perl scripts that create fingerprints from VCF files and compare any two fingerprints. These scripts can be found [here](https://github.com/gglusman/genome-fingerprints/), and will be called regularly during this project.

## Project Abstract
The Earle A. Chiles Cancer Research Institute (EACRI) works closely with Providence Saint-Joseph Health and Services (PSJH) to provide innovative immunotherapy cancer treatments to the hospital network's patients. One such immunotherapy treatment involves harvesting a cancer patient's tumor tissue for analysis and sequencing. The epitomes from the tumor are used to train and grow transformed T-cells from the patient that are able to recognize and attack the patient's tumor cells. This tumor is both patient-specific and tumor-specific, however, and therefore requires that the tumor samples are properly paired to the patient to ensure effective treatment.

As such, this project aims to implement genome fingerprinting and comparison into the tissue processing pipeline as a quality measure that conclusively determines whether a tumor sample comes from a patient or not. This task will be accomplished in the following steps:
1. Create and compare fingerprints from the hundreds of tumor tissue samples and normal tissue samples already processed by PHS
2. Analyze the pairwise fingerprint comparison data, data and statistical analysis tools, and a tumor-normal pair key to determine what correlation between two compared fingerprints ensures tumor-normal pairing.
3. Implement fingerprinting and fingerprint comparison into the PHS tissue processing pipeline.


## Key Files:
File Name | Description
--------- | -----------
[create_fingerprints_from_VCFs](https://git.providence.org/BWG-2019/Sathler-BWG-Project-2019/blob/master/Scripts/python/create_fingerprints_from_VCFs.py) | Searches a designated folder for `.vcf.gz` files and creates fingerprints from them using the `createDMFs.pl` script. The corresponding fingerprints that are created are moved into a second designated folder.
[compare_fingerprints_in_folder](https://git.providence.org/BWG-2019/Sathler-BWG-Project-2019/blob/master/Scripts/python/compare_fingerprints_in_folder.py) | Completes all pairwise comparisons for `.outn.gz` files in a folder (fingerprint files) using the `compareDMFs.pl` script. The comparison `.tsv` files are saved to the `/Comparisons/misc_tsvs` folder, which is created if necessary.
[merge_CSVs](https://git.providence.org/BWG-2019/Sathler-BWG-Project-2019/blob/master/Scripts/python/merge_CSVs.py) | Merges all of the `.tsv` comparison files in the `/misc_tsvs` folder into a large `all_comparisons.csv` file. The script also creates the intermediate `/row_removed_tsvs` and `/tsvs_with_name_rows` folders.
[filter_VCFs](https://git.providence.org/BWG-2019/Sathler-BWG-Project-2019/blob/master/Scripts/python/filter_VCFs.py) | Finds all of the `.vcf.gz` files in a folder and sorts it's contents by line (each line represents a base-pair locus). Lines are written to an output file in a designated folder if they demonstrate sufficient average read depth, average variant allele frequency, and number of positive reads. Each of these parameters can be specified when the script is called.

## References:
1. Glusman G, Mauldin DE, Hood LE, Robinson M. Ultrafast Comparison of Personal Genomes via Precomputed Genome Fingerprints. Front Genet. 2017;8:136. Published 2017 Sep 26. [doi:10.3389/fgene.2017.00136](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5623000/)