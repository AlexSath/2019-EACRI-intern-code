#!/usr/bin/python
import os
import sys
import itertools
import subprocess
#import filter_VCF as fVCF
#import create_fingerprints_from_VCFs as cffV
#import compare_fingerprints_in_folder as cffr
#import merge_CSVs as mCSV

current_file_path = os.path.abspath(__file__)
current_dir_path = os.path.dirname(current_file_path)
alg_range = True
if '--help' in sys.argv:
    with open(os.path.join(current_dir_path, "helptext", "filter_fingerprint_compare.txt"), 'r') as fin:
        for line in fin:
            print(line)
    sys.exit()
else:
    vcf_folder_path = sys.argv[1]
    if os.path.isdir(vcf_folder_path) is not True:
        sys.exit(f"{vcf_folder_path} is not a valid directory.")
    output_folder = sys.argv[2]
    if os.path.isdir(output_folder) is not True:
        sys.exit(f"{output_folder} is not a valid directory")
    categorizer = sys.argv[3]
    min_depth = sys.argv[4]
    min_freq = sys.argv[5]
    min_reads = sys.argv[6]
    max_pop_freq = sys.argv[7]
    tnp_path = sys.argv[8]
    alg_num1 = sys.argv[9]
    alg_num2 = ''
    try:
        alg_num2 = sys.argv[10]
    except:
        if alg_num2 is not '':
            sys.exit(f"{alg_num2} is type {type(alg_num2)}. It must be an Int or NoneType")
        else:
            alg_range = False

GCP_VCF_paths = []
os_file_paths = []

def main():
    proj_dir = os.path.join(output_folder, f"{categorizer}_{min_depth}_{min_freq}_{min_reads}_{max_pop_freq}")
    fAnalyses_path = os.path.join(proj_dir, "fAnalyses")
    fingerprints_path = os.path.join(proj_dir, "Fingerprints")
    Comparisons_path = os.path.join(proj_dir, "Comparisons")
    paths_to_create = [proj_dir, fAnalyses_path, fingerprints_path, Comparisons_path]
    for folder in paths_to_create:
        os.makedirs(folder, exist_ok=True)
    if alg_range is True:
        subprocess.call(['python3', os.path.join(current_dir_path, "filtering_processing", "python", 'filter_VCF.py'), \
                         vcf_folder_path, fAnalyses_path, min_depth, min_freq, \
                         min_reads, max_pop_freq, alg_num1, alg_num2])
    else:
        subprocess.call(['python3', os.path.join(current_dir_path, "filtering_processing", "python", "filter_VCF.py"), \
                         vcf_folder_path, fAnalyses_path, min_depth, min_freq, \
                         min_reads, max_pop_freq, alg_num1])
    subprocess.call(['python3', os.path.join(current_dir_path, "filtering_processing", "python", "create_fingerprints_from_VCFs.py"), \
                     fAnalyses_path, fingerprints_path, 'true'])
    subprocess.call(['python3', os.path.join(current_dir_path, "filtering_processing", "python", "compare_fingerprints_in_folder.py"), \
                     fingerprints_path, Comparisons_path])
    subprocess.call(['python3', os.path.join(current_dir_path, "filtering_processing", "python", "merge_CSVs.py"), Comparisons_path])

    all_comp_path = os.path.join(Comparisons_path, 'all_fingerprint_comparisons.csv')
    updated_csv_path = os.path.join(Comparisons_path, f"{categorizer}_{min_depth}_{min_freq}_{min_reads}_{max_pop_freq}.csv")
    subprocess.call(['mv', f"{all_comp_path}", f"{updated_csv_path}"])

    subprocess.call(['python3', os.path.join(current_dir_path, 'comparison_analysis', 'processing_graphing.py'), updated_csv_path, tnp_path])

if __name__ == "__main__":
    main()
