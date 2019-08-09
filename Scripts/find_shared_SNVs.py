#!/usr/bin/python

import io
import os
import sys
import gzip
import ntpath
import itertools
import subprocess
import numpy as np

ttnnp_file = sys.argv[1]
if os.path.isfile(ttnnp_file) is not True:
    sys.exit(f"The first argument must be a valid file path. You entered '{ttnnp_file}'")
max_or_min = sys.argv[2]
if (max_or_min != 'max') and (max_or_min != 'min'):
    sys.exit(f"The second argument must be 'max' or 'min'. You entered {max_or_min}")
filter = sys.argv[3]
key_file = sys.argv[4]
rawAnalyses_path = sys.argv[5]
comp_types = sys.argv[6].split(',')
try:
    float(filter)
except:
    sys.exit(f"The fourth argument failed to become a float. It must be type float.")

def main():
    tn = tumor_normal_arrs(key_file)

    print(comp_types)
    data_anal_dir = os.path.join(os.path.dirname(ttnnp_file), 'comp_analysis', f'{max_or_min}_{filter}')
    os.makedirs(data_anal_dir, exist_ok=True)
    imp_pairs_path = os.path.join(data_anal_dir, f"{ntpath.basename(ttnnp_file)[:-4]}.{max_or_min}_{filter}.csv")
    if os.path.isfile(imp_pairs_path) is not True:
        for comp_type in comp_types:
            with open(ttnnp_file, 'r') as input:
                with open(imp_pairs_path, 'a+') as fout:
                    for line in input:
                        line_arr = line.replace('\n','').split(',')
                        first_sample = line_arr[0].split('_filtered')[0]
                        second_sample = line_arr[1].split('_filtered')[0]
                        first_sample_type = 'unknown'
                        second_sample_type = 'unknown'
                        line_comp = line_arr[-1]
                        l_200 = line_arr[-2]

                        if line_arr[-1] == comp_type:
                            if first_sample in tn[0]:
                                first_sample_type = 'tumor'
                            elif first_sample in tn[1]:
                                first_sample_type = 'normal'
                            if second_sample in tn[0]:
                                second_sample_type = 'tumor'
                            elif second_sample in tn[1]:
                                second_sample_type = 'normal'

                            if float(l_200) != 1.0:
                                if max_or_min == 'max':
                                    if l_200 < filter:
                                        fout.write(f"{first_sample},{second_sample},{first_sample_type},{second_sample_type},{l_200},{comp_type}\n")
                                elif max_or_min == 'min':
                                    if l_200 > filter:
                                        fout.write(f"{first_sample},{second_sample},{first_sample_type},{second_sample_type},{l_200},{comp_type}\n")
    else:
        print(f"Skipping finding important samples because 'imp_pairs' CSV already exists...")

    imp_comps = nf_f_paths(imp_pairs_path, rawAnalyses_path)
    subprocess.call(['wc', '-l', f"{imp_pairs_path}"])
    print(f"Number of pairs in imp_pairs array: {len(imp_comps)}")

    _snvs_folder = os.path.join(data_anal_dir, f'_snvs_{max_or_min}_{filter}')
    os.makedirs(_snvs_folder, exist_ok=True)
    starting_num_snvs = []
    filtered_num_snvs = []
    num_shared_in_comp = []
    awk_path = os.path.join(_snvs_folder, "awk.sh")
    with open(awk_path, 'w+') as fout:
        fout.write(
        """#!/bin/bash\n
        awk -F'\t' '{if ($3 != "") print $0}'""")
    curr_comp = 1
    for comp in imp_comps:
        print(f'\n-----------[{curr_comp}/{len(imp_comps)}]')
        first_f_snvs = os.path.join(_snvs_folder, f"{comp[0][0][0]}_filtered.snvs")
        second_f_snvs = os.path.join(_snvs_folder, f"{comp[1][0][0]}_filtered.snvs")
        if os.path.isfile(first_f_snvs) is not True:
            print(f"Processing {first_f_snvs[1:].split('/')[-1]}...")
            with gzip.open(comp[0][2], 'rt') as fin:
                with open(first_f_snvs, 'w+') as fout:
                    for line in fin:
                        if '#' not in line:
                            line_arr = line.replace('\n', '').split('\t')
                            fout.write(f"{line_arr[0]}-{line_arr[1]}-{line_arr[3]}-{line_arr[4]}\n")
        else:
            print(f"{first_f_snvs[1:].split('/')[-1]} already exists. Skipping...")

        print(f"Counting SNVs in {first_f_snvs[1:].split('/')[-1]}...")
        p = subprocess.Popen(['wc', '-l', first_f_snvs], stdout = subprocess.PIPE)
        num = p.stdout.read()
        filtered_num_snvs.append([int(num.split(b'/')[0].strip()), comp[0][0][1]])

        print(f"Counting SNVs in {'/'.join(comp[0][1][1:].split('/')[-7:])}...")
        p1 = subprocess.Popen(['gunzip', '-c', comp[0][1]], stdout = subprocess.PIPE)
        p2 = subprocess.Popen(['wc', '-l'], stdin = p1.stdout, stdout = subprocess.PIPE)
        num = p2.stdout.read()
        starting_num_snvs.append([int(num.split(b'/')[0].strip()), comp[0][0][1]])

        if os.path.isfile(second_f_snvs) is not True:
            print(f"Processing {second_f_snvs[1:].split('/')[-1]}...")
            with gzip.open(comp[1][2], 'rt') as fin:
                with open(second_f_snvs, 'w+') as fout:
                    for line in fin:
                        if '#' not in line:
                            line_arr = line.replace('\n', '').split('\t')
                            fout.write(f"{line_arr[0]}-{line_arr[1]}-{line_arr[3]}-{line_arr[4]}\n")
        else:
            print(f"{second_f_snvs[1:].split('/')[-1]} already exists. Skipping...")

        print(f"Counting SNVs in {second_f_snvs[1:].split('/')[-1]}...")
        p = subprocess.Popen(['wc', '-l', second_f_snvs], stdout = subprocess.PIPE)
        num = p.stdout.read()
        filtered_num_snvs.append([int(num.split(b'/')[0].strip()), comp[1][0][1]])

        print(f"Counting SNVs in {'/'.join(comp[1][1][1:].split('/')[-7:])}...")
        p1 = subprocess.Popen(['gunzip', '-c', comp[1][1]], stdout = subprocess.PIPE)
        p2 = subprocess.Popen(['wc', '-l'], stdin = p1.stdout, stdout = subprocess.PIPE)
        num = p2.stdout.read()
        starting_num_snvs.append([int(num.split(b'/')[0].strip()), comp[1][0][1]])

        print(f"\nComparing {comp[0][0][0]} with {comp[1][0][0]}")
        p1 = subprocess.Popen(['comm', first_f_snvs, second_f_snvs], stdout = subprocess.PIPE)
        p2 = subprocess.Popen(['bash', awk_path], stdin = p1.stdout, stdout = subprocess.PIPE)
        p3 = subprocess.Popen(['wc', '-l'], stdin = p2.stdout, stdout = subprocess.PIPE)
        num_shared = p3.stdout.read()
        num_shared = int(num_shared.strip())

        if (comp[0][0][1] == "tumor" and comp[1][0][1] == "normal") \
            or (comp[0][0][1] == "normal" and comp[1][0][1] == "tumor"):
            if f"{comp[0][0][1]},{comp[1][0][1]}" or f"{comp[1][0][1]},{comp[0][0][1]}" in open(key_file).read():
                print(f"num shared SNVs: {num_shared} | comp type: {comp[0][0][1]}-{comp[1][0][1]}-pair")
                num_shared_in_comp.append([num_shared, f"{comp[0][0][1]}-{comp[1][0][1]}-pair"])
            else:
                print(f"num shared SNVs: {num_shared} | comp type: {comp[0][0][1]}-{comp[1][0][1]}-unpaired")
                num_shared_in_comp.append([num_shared, f"{comp[0][0][1]}-{comp[1][0][1]}-unpaired"])
        else:
            print(f"num shared SNVs: {num_shared} | comp type: {comp[0][0][1]}-{comp[1][0][1]}")
            num_shared_in_comp.append([num_shared, f"{comp[0][0][1]}-{comp[1][0][1]}"])

        curr_comp += 1

    os.remove(awk_path)
    print(num_shared_in_comp)

    tumor_pre_filter = get_sample_type_nums(starting_num_snvs, 'tumor')
    normal_pre_filter = get_sample_type_nums(starting_num_snvs, 'normal')
    tumor_post_filter = get_sample_type_nums(filtered_num_snvs, 'tumor')
    normal_post_filter = get_sample_type_nums(filtered_num_snvs, 'normal')

    tumor_pre_sum = get_data_sum(tumor_pre_filter)
    normal_pre_sum = get_data_sum(normal_pre_filter)
    tumor_post_sum = get_data_sum(tumor_post_filter)
    normal_post_sum = get_data_sum(normal_post_filter)

    comp_SNV_array = get_comp_comm_SNVs(num_shared_in_comp)

    with open(os.path.join(data_anal_dir, f"tnp_SNV_anal_{max_or_min}_{filter}.tsv"), 'a+') as fout:
        fout.write(f"##Source: {ttnnp_file}\n")
        fout.write(f"##Filter: {max_or_min} correlation of {filter} for {comp_types}\n")
        fout.write(f"#sample_type:\tnumsmpl\tmean\tmed.\tstd.\tmin\tmax\n")
        fout.write(f"Tumor (pre):\t{len(tumor_pre_filter)}\t{tumor_pre_sum[0]}\t{tumor_pre_sum[1]}\t{tumor_pre_sum[2]}\t{tumor_pre_sum[3]}\t{tumor_pre_sum[4]}\n")
        fout.write(f"Tumor (post):\t{len(tumor_post_filter)}\t{tumor_post_sum[0]}\t{tumor_post_sum[1]}\t{tumor_post_sum[2]}\t{tumor_post_sum[3]}\t{tumor_post_sum[4]}\n")
        fout.write(f"Normal (pre):\t{len(normal_pre_filter)}\t{normal_pre_sum[0]}\t{normal_pre_sum[1]}\t{normal_pre_sum[2]}\t{normal_pre_sum[3]}\t{normal_pre_sum[4]}\n")
        fout.write(f"Normal (post):\t{len(normal_post_filter)}\t{normal_post_sum[0]}\t{normal_post_sum[1]}\t{normal_post_sum[2]}\t{normal_post_sum[3]}\t{normal_post_sum[4]}\n")

    with open(os.path.join(data_anal_dir, f"comp_type_anal_{max_or_min}_{filter}.tsv"), 'a+') as fout:
        fout.write(f"##Source: {ttnnp_file}\n")
        fout.write(f"##Filter: {max_or_min} correlation of {filter} for {comp_types}\n")
        fout.write(f"## _ = 'shared SNVs'\n")
        fout.write(f"#type\tnumcomp\tmean _\tmed. _\tstd. _\tmin _\tmax _\n")
        for comp in comp_SNV_array:
            fout.write(f"{comp[0]}\t{comp[1]}\t{comp[2][0]}\t{comp[2][1]}\t{comp[2][2]}\t{comp[2][3]}\t{comp[2][4]}\n")

def tumor_normal_arrs(filepath):
    tumor = []
    normal = []
    with open(filepath, 'r') as fin:
        curr_line = 1
        for line in fin:
            line_arr = line.replace('\n', '').split(',')
            if curr_line == 1:
                if 'tumor_id' == line_arr[0]:
                    tumor_pos = 0
                    normal_pos = 1
                else:
                    normal_pos = 0
                    tumor_pos = 1
            else:
                tumor.append(line_arr[tumor_pos])
                normal.append(line_arr[normal_pos])
            curr_line += 1

    return [tumor, normal]


def nf_f_paths(name_file, rawAnalyses):
    all_samples_array = []
    with open(name_file, 'r') as fin:
        current_name = -2
        current_line = 1
        for line in fin:
            if current_line != 'header':
                line_arr = line.replace('\n', '').split(',')
                all_samples_array.append([[[line_arr[0], line_arr[2]], \
                                         get_vcf_path(line_arr[0], rawAnalyses), \
                                         get_vcf_path(line_arr[0], '')],
                                         [[line_arr[1], line_arr[3]], \
                                         get_vcf_path(line_arr[1], rawAnalyses), \
                                         get_vcf_path(line_arr[1], '')]])
            current_name += 1
            current_line += 1
    return all_samples_array


def get_vcf_path(name, dir):
    path = ''
    if dir == '':
        dir = os.path.abspath(os.path.join(os.path.dirname(ttnnp_file), '..', '..', 'fAnalyses'))
        if os.path.isdir(dir) is not True:
            sys.exit(f"Was not able to find fAnalyses using the ttnnp filepath")
        for root, dirs, files in os.walk(dir):
            for f in files:
                if '.vcf' and 'filtered' in f:
                    if name in root:
                        if name in f:
                            path = os.path.join(root, f)
                            break
    else:
        for root, dirs, files in os.walk(dir):
            for f in files:
                if '.vcf.gz' in f:
                    if name in root:
                        path = os.path.join(root, f)
                        break
    if (path == '') or (os.path.isfile(path) is False):
         sys.exit(f"Was not able to find {name} in {dir}")
    return path

def get_sample_type_nums(array, sample_type):
    num_arr = []
    curr_sample = 0
    for num_type in array:
        if num_type[1] == sample_type:
            num_arr.append(array[curr_sample][0])
        curr_sample += 1
    return num_arr

def get_arr_mean(array):
    sum = 0
    for num in array:
        sum += num
    try:
        mean = sum / len(array)
    except:
        mean = None
    return mean

def get_data_sum(arr):
    return [int(round(get_arr_mean(arr))), int(round(np.median(arr))), \
            int(round(np.std(arr))), int(round(np.min(arr))), \
            int(round(np.max(arr)))]

def get_comp_comm_SNVs(arr):
    tt = [0]
    nn = [0]
    tnp = [0]
    tnup = [0]
    for num_type in arr:
        if num_type[1] == "tumor-tumor":
            if tt == [0]:
                tt = []
            tt.append(num_type[0])
        if num_type[1] == "normal-normal":
            if nn == [0]:
                nn = []
            nn.append(num_type[0])
        if num_type[1] == "tumor-normal-pair" or num_type[1] == "normal-tumor-pair":
            if tnp == [0]:
                tnp = []
            tnp.append(num_type[0])
        if num_type[1] == "tumor-normal-unpaired" or num_type[1] == "normal-tumor-unpaired":
            if tnup == [0]:
                tnup = []
            tnup.append(num_type[0])
    return [['tt', len(tt), get_data_sum(tt)], \
            ['nn', len(nn), get_data_sum(nn)], \
            ['tnp', len(tnp), get_data_sum(tnp)], \
            ['tnup', len(tnup), get_data_sum(tnup)]]

if __name__ == "__main__":
    main()
