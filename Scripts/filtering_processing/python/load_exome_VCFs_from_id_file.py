#!/usr/bin/python
import os
import sys
import itertools
import subprocess
import fingerprint_processing_tools as fpt

current_file_path = os.path.abspath(__file__)
current_dir_path = os.path.dirname(current_file_path)
gs_bucket_string = sys.argv[1]
# "gs://phs-exome/Analyses/**/ann.vcf.gz"
dir_to_populate = sys.argv[2]
root_dir = os.path.dirname(dir_to_populate)
ids_to_find_file = sys.argv[3]

GCP_VCF_paths = []
os_file_paths = []
id_list = []

with open(ids_to_find_file) as fin:
    for line in fin:
        id_list.append(line.replace('\n',''))

def main():
    process = subprocess.Popen(["gsutil", "ls", gs_bucket_string], stdout=subprocess.PIPE, encoding='utf8')
    stdout, stderr = process.communicate()

    temp_paths = stdout.split('\n')[:-1]
    num_temp_paths = len(temp_paths)
    print(num_temp_paths)

    with open(os.path.join(root_dir, 'copied_exome_ids.txt'), 'w+') as fout:
        for path in temp_paths:
            for id in id_list:
                mod_id = f"{id[0:2]}-{id[2:5]}-{id[5:10]}"
                if (id in path) or (mod_id in path):
                    if ('v0.4.0' in path) or ('v0.5.0' in path):
                        if ('QC' not in path) and ('POS' not in path) and ('NEC' not in path):
                            GCP_VCF_paths.append(path)
                            fout.write(f"{id}\t{path}\n")
                            break

    current_path = 0
    for path in GCP_VCF_paths:
        print(str(current_path) + " " + path)
        temp_path = path[24:]
        temp_arr = temp_path.split('/')

        smpl_name = temp_arr[0].replace('_DNA','').split('_')[-1].replace('-','')
        path_to_create = os.path.join(dir_to_populate, smpl_name, '/'.join(temp_arr[1:-1]))
        print(path_to_create)
        os.makedirs(path_to_create, exist_ok = True)

        os_file_paths.append(os.path.join(dir_to_populate, smpl_name, '/'.join(temp_arr[1:-1])))
        current_path += 1

    if len(GCP_VCF_paths) != len(os_file_paths):
        sys.exit(f"Number of GCP Paths {len(GCP_VCF_paths)} does not equal Number of OS Paths P{len(os_file_paths)}")

    for gcp_path, os_path in zip(GCP_VCF_paths, os_file_paths):
        p = subprocess.Popen(["gsutil", "cp", gcp_path, os_path])
        p.wait()
        current_path += 1

    return False

if __name__ == "__main__":
    main()
