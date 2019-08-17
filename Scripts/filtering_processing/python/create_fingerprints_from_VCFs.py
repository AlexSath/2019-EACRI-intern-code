#!/usr/bin/python
#Created by Alex R. Sathler for the EACRI
#Perl Script for computing genome fingerprints found at /gglusman/genome-fingerprints/ on GitHub
import os
import sys
import subprocess
import itertools
import fingerprint_processing_tools as fpt

current_file_path = os.path.abspath(__file__)
current_dir_path = os.path.dirname(current_file_path)
vcf_directory = sys.argv[1] #First argument is where the .vcf files can be found (nested folders are ok)
output_directory = sys.argv[2] # Second argument is where the .vcf files will be placed. Subfolders will be created to organize some files
filtered_or_not = 'False'
filtered_or_not = sys.argv[3].lower().strip() #Third argument tells the function whether or not the .vcfs are filtered or not
vcf_file_paths = []
vcf_file_names = []

def main():
    print('Retrieving VCF file paths and file names...')
    vcf_names_paths = retrieve_vcf_files(vcf_directory, ".vcf.gz")

    if len(vcf_names_paths[0]) != len(vcf_names_paths[1]):
        print(f"Fatal Error: vcf_file_paths and vcf_file_names are not of equal length")
        print(f"Number of File Paths: {len(vcf_names_paths[1])}")
        sys.exit(f"Number of File Names: {len(vcf_names_paths[0])}")

    if os.path.isdir(os.path.join(output_directory, 'outn')):
        sys.exit(f"Folder outn already exists. Exiting under the assumption that fingerprints have already been made...")

    #Iterate through both names and paths while pairing names and paths together:
    total_fingerprints = len(vcf_names_paths[0])
    current_fingerprint = 1
    for (n, p) in zip(vcf_names_paths[0], vcf_names_paths[1]):
        print(f'Creating fingerprint for {n}... [{current_fingerprint}/{total_fingerprints}]')
        #Find and call the computeDMF script, with a file name output and filepath input
        subprocess.call(["perl", fpt.find_file_path(current_dir_path, "computeDMF.pl", 2), os.path.join(output_directory, "..", n), p])
        current_fingerprint += 1

    #Find and call the bash script that organizes the .out, .outn, and .out.close files into their respective folders
    subprocess.call(["bash", fpt.find_file_path(current_dir_path, "organizeFingerprints.sh", 2), output_directory])


'''
a function that will retrieve all of the file names or file directories with
.vcf.gz extension within the path directory and all of its sub-directories

vcf_path [string]: the path to the directory that will be searched (and
therefore sub-directories as well)
is_path [integer]: if 1, the function will output file paths to the files array.
Otherwise, the function will output file names to the files array
'''
def retrieve_vcf_files(vcf_path, string1):
    names_paths = [[], []]
    # r=root, d=directories, f=files
    for r, d, f in os.walk(vcf_path):
        for file in f:
            file_path = os.path.join(r, file)
            file_path_arr = file_path[1:].split('/')
            if '_filtered_' in file and string1 in file:
                if (file_path_arr[-2] == "exome") or (file_path_arr[-2] == "Analysis") or (".vcf.gz" in file_path_arr[-2]):
                    names_paths[0].append(file[:-7])
                    names_paths[1].append(os.path.join(r, file))
            if string1 in file and '_filtered_' not in file:
                if file_path_arr[-7] not in names_paths[0]:
                    names_paths[0].append(file_path_arr[-7])
                else:
                    names_paths[0].append(f"{file_path_arr[-7]}~2")
                names_paths[1].append(os.path.join(r, file))

    return names_paths

if __name__ == "__main__":
    main()
