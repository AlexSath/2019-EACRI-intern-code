# Accurate Determination of Tumor Sample-Patient Sample Pairs for T-Cell Immunotherapy Quality Assurance
***Created by Alexandre R. Sathler for the Earle A. Chiles Research Institute***

## Project Abstract
The Earle A. Chiles Cancer Research Institute (EACRI) works closely with Providence Saint-Joseph Health and Services (PSJH) to provide innovative immunotherapy cancer treatments to the hospital network's patients. One such immunotherapy treatment involves harvesting a cancer patient's tumor tissue for analysis and sequencing. The epitomes from the tumor are used to train and grow transformed T-cells from the patient that are able to recognize and attack the patient's tumor cells. This tumor is both patient-specific and tumor-specific, however, and therefore requires that the tumor samples are properly paired to the patient to ensure effective treatment.

As such, this project aims to implement genome fingerprinting and comparison into the tissue processing pipeline as a quality measure that conclusively determines whether a tumor sample comes from a patient or not. This task will be accomplished in the following steps:
1. Create and compare fingerprints from the hundreds of tumor tissue samples and normal tissue samples already processed by PHS
2. Analyze the pairwise fingerprint comparison data, data and statistical analysis tools, and a tumor-normal pair key to determine what correlation between two compared fingerprints ensures tumor-normal pairing.
3. Implement fingerprinting and fingerprint comparison into the PHS tissue processing pipeline.

## Project Background
This project and its code draws from genome fingerprinting and fingerprint comparison methods outlined in 'Ultrafast Comparison of Personal Genomes via Precomputed Genome Fingerprints,' published in *Frontiers in Genetics* ([Glusman, 2017](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5623000/)). The paper posits that fingerprints can be created from any genome by comparing the genome to the average genome sequence of a population. When making this comparison, Single Nucleotide Variants (SNVs) at every locus become apparent, and can be stored in orderly fashion in the Variant Call Format (`.vcf`). Because every genome will have unique SNVs when compared to a population genome, these SNVs can be used to create a unique genome fingerprint. When compared, fingerprints can be used to identify individuals and determine the similarity between two genomes much more quickly than full genome analysis.

Glusman also published Perl scripts that create fingerprints from VCF files and compare any two fingerprints. These scripts can be found [here](https://github.com/gglusman/genome-fingerprints/), and will be called regularly during this project.


# General Notes:
### Analyses Included in this Repository:
Some analyses and graphs have been included in this repository for understanding of some output files. Note that comparison files such as `.ttnnp.csv` (see below) cannot be included because they include sensitive personal health information. Instead PDFs of Graphs and TSVs of Data Summaries will be included to show what is possible with this process/package.
- v050: A 758 sample dataset of WES data from EACRI/PSJH cancer patients. One filter_fingerprint_compare process takes several days.
- 50v050: A 108 sample subset (54 tumor-normal-pairs) of the v050 dataset - better for quick analyses. One filter_fingerprint_compare process takes a few hours.
- 'Dump': A small 26 sample dataset (13 tumor-normal-pairs). Some samples are shared with the v050 dataset. One filter_fingerprint_compare process takes an hour at most.

### Naming Scheme:
This project has a fairly unique naming scheme to help identify the results of different filters. This section aims to help explain this naming scheme.

Lets take the following as an example:
> 50v050c1fGatkFrbys_50_10_2_1

It should be interpreted like so:
> *[Identifier]* _ *[min-allele-depth]* _ *[min-variant-allele-frequency]* _ *[min-number-of-SNV-callers-called]* _ *[max-ref-population-allele-freq]*

In this case, the identifier is '50v050c1fGatkFrbys', meaning that the 50v050 dataset was used, the variable c in `computeDMFs.pl` was set to 1, and the VCF file was filtered only taking callers GATK and Freebayes into account. The remaining numbers indicate the other filters. For a SNV to pass the filter, it had to have:
- a minimum allele depth of 50
- a minimum allele frequency of 10%
- a minimum number of 2 callers calling the SNV
- a maximum reference population allele frequency of 1%. So fewer than 1% of the reference population (the script looks for EXAC_all and 1000g) must have the SNV for it to pass the filter.

### File Types:
- `"All Comparisons CSV"`:
- `.cp.csv`:
- `.ttnnp.csv`:
- `Comp_Type_SNV_Findx-[max/min]_[dec. filter].tsv`:
- `Smpl_Type_SNV_Findx-[max/min]_[dec. filter].tsv`:


# Script Descriptions and Uses:
### Files From Glusman et al.'s Paper:
File Name | Description
--------- | -----------
[computeDMFs.pl](https://github.com/AlexSath/2019-EACRI-intern-code/blob/master/Scripts/filtering_processing/fingerprintScripts/computeDMF.pl) | Exactly the same as found in Glusman's repository, except default C values are set to 1 instead of twenty. This prevents close-together SNVs from being filtered out of the Fingerprinting process as this was found to be beneficial for Whole Exome Sequencing Data.
[compareDMFs.pl](https://github.com/AlexSath/2019-EACRI-intern-code/blob/master/Scripts/filtering_processing/fingerprintScripts/compareDMF.pl) | This file has not been changed from Glusman's repository.

### Core Process Files:
File Name | Description
--------- | -----------
[filter_fingerprint_compare.py](https://github.com/AlexSath/2019-EACRI-intern-code/blob/master/Scripts/filter_fingerprint_compare.py) | The primary file for this process. When provided a root directory that contains gzipped vcfs (it is assumed that all VCFs are 7 levels deep within the root directory) and an output folder, a new folder is created and named using a given identifier and the VCF filters used. Then, using a Tumor-Normal-Pair key provided by the user, all VCFs are fingerprinted and compared and all resulting comparison data is graphed.
[filter_comparisons_and_graph.py](https://github.com/AlexSath/2019-EACRI-intern-code/blob/master/Scripts/filter_comparisons_and_graph.py) | This file should be used when filtering and comparing a set of VCFs has already been done, but a subset of all pairwise comparisons must be graphed. This folder will take a comparison csv as input and filter the comparisons by samples that are found in a provided key file. The subset of comparisons will then be graphed in a new sub-folder. Mostly obsolete in this version of the process.
[find_shared_SNVs_FingIndex.py](https://github.com/AlexSath/2019-EACRI-intern-code/blob/master/Scripts/find_shared_SNVs_FingIndex.py) | This script, when provided a `.ttnnp.csv` file as input, will find all comparisons with correlations above or below a desired value for all provided comparison types (e.g. tumor-normal-pair), and output one .tsv with individual shared SNVs and Fingerprint Indexes for all comparisons, and two data summary .tsvs for better understanding of the fingerprinting process with your data.
[processing_graphing.py](https://github.com/AlexSath/2019-EACRI-intern-code/blob/master/Scripts/comparison_analysis/processing_graphing.py) | This file is automatically called by `filter_fingerprint_compare.py` but can also be called manually. When provided an all_comparisons.csv file and tumor_normal_key file as input, it will output a tumor-normal-pair true/false `.cp.tnp.csv` and a comparison_type `.ttnnp.csv` as output and create violin plots from these csv types.

### Files for Filtering and Processing:
File Name | Description
--------- | -----------
[create_fingerprints_from_VCFs.py](https://github.com/AlexSath/2019-EACRI-intern-code/blob/master/Scripts/filtering_processing/python/create_fingerprints_from_VCFs.py) | Searches a designated folder for `.vcf.gz` files and creates fingerprints from them using the `createDMFs.pl` script. The corresponding fingerprints that are created are moved into a second designated folder.
[compare_fingerprints_in_folder.py](https://github.com/AlexSath/2019-EACRI-intern-code/blob/master/Scripts/filtering_processing/python/compare_fingerprints_in_folder.py) | Completes all pairwise comparisons for `.outn.gz` files in a folder (fingerprint files) using the `compareDMFs.pl` script. The comparison `.tsv` files are saved to the `/Comparisons/misc_tsvs` folder, which is created if necessary.
[merge_CSVs.py](https://github.com/AlexSath/2019-EACRI-intern-code/blob/master/Scripts/filtering_processing/python/merge_CSVs.py) | Merges all of the `.tsv` comparison files in the `/misc_tsvs` folder into a large `all_comparisons.csv` file. The script also creates the intermediate `/row_removed_tsvs` and `/tsvs_with_name_rows` folders.
[filter_VCFs.py](https://github.com/AlexSath/2019-EACRI-intern-code/blob/master/Scripts/filtering_processing/python/filter_VCF.py) | Finds all of the `.vcf.gz` files in a folder and sorts it's contents by line (each line represents a base-pair locus). Lines are written to an output file in a designated folder if they demonstrate sufficient average read depth, average variant allele frequency, and number of positive reads. Each of these parameters can be specified when the script is called.

### Files for Data Processing and Graphing:

## References:
1. Glusman G, Mauldin DE, Hood LE, Robinson M. Ultrafast Comparison of Personal Genomes via Precomputed Genome Fingerprints. Front Genet. 2017;8:136. Published 2017 Sep 26. [doi:10.3389/fgene.2017.00136](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5623000/)