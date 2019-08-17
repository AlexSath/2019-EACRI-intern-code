#!/usr/bin/python
import os
import sys
import itertools
import subprocess
import fingerprint_processing_tools as fpt

current_file_path = os.path.abspath(__file__)
current_dir_path = os.path.dirname(current_file_path)
vcf_folder_path = sys.argv[1]
dir_to_populate = sys.argv[2]

GCP_VCF_paths = []
os_file_paths = []
def main():
    process = subprocess.Popen(["gsutil", "ls", "gs://phs-exome/Analyses/**/ann.vcf.gz"], stdout=subprocess.PIPE, encoding='utf8')
    stdout, stderr = process.communicate()

    temp_paths = stdout.split('\n')
    num_temp_paths = len(temp_paths)
    print(num_temp_paths)

    current_path = 0
    while (current_path != num_temp_paths):
        if "v0.5.0" in temp_paths[current_path]:
            print("Appending " + temp_paths[current_path] + "...")
            GCP_VCF_paths.append(temp_paths[current_path])
            current_path += 1
        else:
            print("v0.5.0 not found in: " + temp_paths[current_path] + "...")
        current_path += 1

    num_paths = len(GCP_VCF_paths)

    current_path = 0
    while (current_path != num_paths):
        print(str(current_path) + " " + GCP_VCF_paths[current_path][15:-11])
        path_to_create = os.path.join(dir_to_populate, GCP_VCF_paths[current_path][15:-11])
        fpt.create_dir_if_absent(path_to_create)

        os_file_paths.append(os.path.join(dir_to_populate, GCP_VCF_paths[current_path][15:]))
        current_path += 1

    current_path = 0
    while (current_path != num_paths):
        p = subprocess.Popen(["gsutil", "cp", GCP_VCF_paths[current_path], os_file_paths[current_path]])
        p.wait()
        current_path += 1

    return False

if __name__ == "__main__":
    main()
