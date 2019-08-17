#!/usr/bin/python
import os
import sys
import itertools
import subprocess

current_file_path = os.path.abspath(__file__)
current_dir_path = os.path.dirname(current_file_path)
str_to_search = sys.argv[1]
output_folder = sys.argv[2]
root_dir = os.path.dirname(output_folder)
ids_to_find_file = sys.argv[3]
# "gs://phs-exome/Analyses/**/ann.vcf.gz"
# "gs://phs-exome/PPMP-TST170/**/variants.annotated.vcf"

id_list = []
with open(ids_to_find_file) as fin:
    for line in fin:
        id_list.append(line.replace('\n',''))

GCP_VCF_paths = []
os_file_paths = []
def main():
    process = subprocess.Popen(["gsutil", "ls", str_to_search], stdout=subprocess.PIPE, encoding='utf8')
    stdout, stderr = process.communicate()

    temp_paths = stdout.split('\n')[:-1]
    num_temp_paths = len(temp_paths)
    print(num_temp_paths)

    with open(os.path.join(root_dir,'copied_tst170_ids.txt'), 'w+') as fout:
        for path in temp_paths:
            path_smpl = path.split('/')[5]
            for id in id_list:
                mod_id = f"{id[0:2]}-{id[2:5]}-{id[5:10]}"
                if (id in path_smpl) or (mod_id in path_smpl):
                    if ('QC' not in path) and ('POS' not in path) and ('NEC' not in path):
                        GCP_VCF_paths.append(path)
                        fout.write(f"{id}\t{path}\n")
                        break

    current_path = 0
    for path in GCP_VCF_paths:
        print(str(current_path) + " " + path)
        temp_path = GCP_VCF_paths[current_path][27:]
        temp_arr = temp_path.split('/')

        new_arr = [temp_arr[1].replace('_DNA','').split('_')[-1].replace('-','')]
        new_arr.append(temp_arr[0])
        new_arr.append('hg19')
        curr = 2
        while curr <= len(temp_arr) - 1:
            new_arr.append(temp_arr[curr])
            curr += 1

        path_to_create = os.path.join(output_folder, '/'.join(new_arr))
        os_file_paths.append(path_to_create)
        os.makedirs('/'.join(path_to_create.split('/')[:-1]), exist_ok = True)
        current_path += 1

    if len(GCP_VCF_paths) != len(os_file_paths):
        sys.exit(f"Number of GCP Paths {len(GCP_VCF_paths)} does not equal Number of OS Paths P{len(os_file_paths)}")

    for gcp_path, os_path in zip(GCP_VCF_paths, os_file_paths):
        p = subprocess.Popen(["gsutil", "cp", gcp_path, os_path])
        p.wait()
        p2 = subprocess.Popen(['gzip', os_path])
        p2.wait()

    return False

if __name__ == "__main__":
    main()
