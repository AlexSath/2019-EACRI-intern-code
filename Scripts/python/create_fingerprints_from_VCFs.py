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
filtered_or_not = sys.argv[3] #Third argument tells the function whether or not the .vcfs are filtered or not
vcf_file_paths = []
vcf_file_names = []

def main():
    print('Retrieving VCF file paths...')
    vcf_file_paths = retrieve_vcf_files(vcf_directory, 0) #An array of absolute .vcf file paths
    print('Retrieving VCF file names...')
    #If the .vcfs are filtered, then use the second case of the funciton, if not, use the third
    if (filtered_or_not == 'True'):
        vcf_file_names = retrieve_vcf_files(vcf_directory, 1) #An array of .vcf file names without file extensions
    else:
        vcf_file_names = retrieve_vcf_files(vcf_directory, 2)

    #Iterate through both names and paths while pairing names and paths together:
    total_fingerprints = len(vcf_file_names)
    current_fingerprint = 1
    for (n, p) in zip(vcf_file_names, vcf_file_paths):
        print(f'Creating fingerprint for {n}... [{current_fingerprint}/{total_fingerprints}]')
        #Find and call the computeDMF script, with a file name output and filepath input
        subprocess.call(["perl", fpt.find_file_path(current_dir_path, "computeDMF.pl", 2), n, p])
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
def retrieve_vcf_files(vcf_path, type):
    files = []
    # r=root, d=directories, f=files
    if (type == 0): #if finding paths, go here
        for r, d, f in os.walk(vcf_path):
            for file in f:
                if '.vcf.gz' in file: #if it is a vcf.gz file
                    file_path = os.path.join(r, file) #the absolute path of the file
                    file_no_slash = file_path.split('/') #the file path split into an array by '/'
                    num_path_terms = len(file_no_slash) #the length of the file path array
                    if (file_no_slash[num_path_terms - 2] == "exome"): #if the file's parent directory is called exome, then append.
                        files.append(file_path)
    elif (type == 1): #if finding filtered file names, go here;
        for r, d, f in os.walk(vcf_path):
            for file in f:
                if '_filtered_' in file and '.vcf.gz' in file:
                    file_path = os.path.join(r, file) #the absolute path of the file
                    file_no_slash = file_path.split('/') #the file path split into an array by '/'
                    num_path_terms = len(file_no_slash) #the length of the file path array
                    if (file_no_slash[num_path_terms - 2] == "exome"): #if the file's parent directory is called exome, then append.
                        files.append(file[:-7])
    else:
        for r, d, f in os.walk(vcf_path):
            for dir in d:
                #If the parent folder is the Comparisons folder, then use the current folder as the .vcf name
                folder_parent = os.path.abspath(os.path.join(r, dir, ".."))
                if (folder_parent == os.path.abspath(vcf_path)):
                    files.append(dir)
    return files

if __name__ == "__main__":
    main()
