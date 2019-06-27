import os
import sys
import subprocess
import itertools
import fingerprint_processing_tools

current_file_path = os.path.abspath(__file__)
current_dir_path = os.path.dirname(current_file_path)
vcf_directory = ""

def main():
    vcf_directory = fingerprint_processing_tools.ask_for_path('containing desired .vcf.gz files')
    #vcf_directory = sys.argv[1]
    output_directory = fingerprint_processing_tools.ask_for_path('where fingerprints will output')
    #output_directory = sys.argv[2]

    print('Retrieving VCF file paths...')
    vcf_file_paths = retrieve_vcf_files(vcf_directory, 1)
    print('Retrieving VCF file names...')
    vcf_file_names = retrieve_vcf_files(vcf_directory, 0)

    for (n, p) in zip(vcf_file_names, vcf_file_paths):
        print('Creating fingerprint for {}...'.format(n))
        subprocess.call(["perl", fingerprint_processing_tools.find_file_path(current_dir_path, "computeDMF.pl", 2), n, p])

    print(fingerprint_processing_tools.find_file_path(current_dir_path, "computeDMF.pl", 2))
    subprocess.call(["bash", fingerprint_processing_tools.find_file_path(current_dir_path, "organizeFingerprints.sh", 2), output_directory])


'''
a function that will retrieve all of the file names or file directories with
.vcf.gz extension within the path directory and all of its sub-directories

vcf_path [string]: the path to the directory that will be searched (and
therefore sub-directories as well)
is_path [integer]: if 1, the function will output file paths to the files array.
Otherwise, the function will output file names to the files array
'''
def retrieve_vcf_files(vcf_path, is_path):
    files = []
    # r=root, d=directories, f=files
    for r, d, f in os.walk(vcf_path):
        if (is_path == 1): #if finding paths, go here
            for file in f:
                if '.vcf.gz' in file: #if it is a vcf.gz file
                    file_path = os.path.join(r, file) #the absolute path of the file
                    file_no_slash = file_path.split('/') #the file path split into an array by '/'
                    num_path_terms = len(file_no_slash) #the length of the file path array
                    if (file_no_slash[num_path_terms - 2] == "exome"): #if the file's parent directory is called exome, then append.
                        files.append(file_path)
        else: #if finding file names, go here;
            for dir in d:
                folder_parent = os.path.abspath(os.path.join(r, dir, ".."))
                if (folder_parent == os.path.abspath(vcf_path)):
                    files.append(dir)
    return files

main()
