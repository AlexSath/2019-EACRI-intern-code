#!/usr/bin/python

import io
import os
import sys
import time
import gzip
import ntpath
import itertools
import subprocess
import numpy as np

file_path = os.path.abspath(__file__)
current_dir_path = os.path.dirname(file_path)

if '--help' in sys.argv:
    with open(os.path.join(current_dir_path, "helptext", "find_shared_SNVs_FingIndex.txt"), 'r') as fin:
        for line in fin:
            print(line)
    sys.exit()
else:
    ttnnp_file = sys.argv[1]
    fAnalyses_path = os.path.abspath(os.path.join(os.path.dirname(ttnnp_file), '..', '..', 'fAnalyses'))
    if os.path.isdir(fAnalyses_path) is not True:
        sys.exit(f"Was not able to find fAnalyses using the ttnnp filepath")
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
    imp_pairs_file = os.path.join(data_anal_dir, f"{ntpath.basename(ttnnp_file)[:-10]}.{max_or_min}{filter}.csv")

    if os.path.isfile(imp_pairs_file) is not True:
        with open(ttnnp_file, 'r') as input:
            with open(imp_pairs_file, 'a+') as fout:
                fout.write(f"first_smpl,second_smpl,first_smpl_type,sec_smpl_type,corr_l_200,comp_type\n")
                for line in input:
                    line_arr = line.replace('\n','').split(',')
                    first_sample = line_arr[0].split('_filtered')[0]
                    second_sample = line_arr[1].split('_filtered')[0]
                    first_sample_type = 'unknown'
                    second_sample_type = 'unknown'
                    line_comp = line_arr[-1]
                    l_200 = line_arr[-2]

                    if line_comp in comp_types:
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
                                    fout.write(f"{first_sample},{second_sample},{first_sample_type},{second_sample_type},{l_200},{line_comp}\n")
                            elif max_or_min == 'min':
                                if l_200 > filter:
                                    fout.write(f"{first_sample},{second_sample},{first_sample_type},{second_sample_type},{l_200},{line_comp}\n")
    else:
        print(f"Skipping finding important samples because 'imp_pairs' CSV already exists...")

    _snvs_folder = os.path.join(data_anal_dir, f'_snvs_{max_or_min}_{filter}')
    os.makedirs(_snvs_folder, exist_ok=True)

    _findexes_folder = os.path.join(data_anal_dir, f'_findexes_{max_or_min}_{filter}')
    os.makedirs(_findexes_folder, exist_ok=True)

    awk_path = os.path.join(_snvs_folder, "awk.sh")
    with open(awk_path, 'w+') as fout:
        fout.write("""#!/bin/bash\n awk -F'\t' '{if ($3 != "") print $0}'""")

    curr_line = 1
    curr_comp = 1
    total_comps = get_line_count(imp_pairs_file) - 1
    #               Both sample types contain:
    #               Each with RAW and FILTERED sub-arrays
    #               SNV data          Findex data
    tumor_smpl_data = [[[], []],      [[], []]]
    normal_smpl_data = [[[], []],     [[], []]]
    #           each has comm_snvs and comm_findexes
    #              tt       nn       tnp      tnup
    comp_type_data = [[[],[]], [[],[]], [[],[]], [[],[]]]
    with open(imp_pairs_file) as fin:
        for line in fin:
            if curr_line == 1:
                curr_line += 1
                continue
            print(f"\n-----------[{curr_comp}/{total_comps}]")
            line_arr = line.replace("\n",'').split(',')
            first_pre_path = get_vcf_path(line_arr[0], rawAnalyses_path)
            first_post_path = get_vcf_path(line_arr[0], fAnalyses_path)
            sec_pre_path = get_vcf_path(line_arr[1], rawAnalyses_path)
            sec_post_path = get_vcf_path(line_arr[1], fAnalyses_path)

            first_f_snvs = os.path.join(_snvs_folder, f"{line_arr[0]}_filtered.snvs")
            sec_f_snvs = os.path.join(_snvs_folder, f"{line_arr[1]}_filtered.snvs")

            print(f"Counting SNVs in {'/'.join(first_pre_path[1:].split('/')[-7:])}...")
            num_first_Rsnvs = get_gzipped_line_count(first_pre_path)

            print(f"Finding SNVs for {ntpath.basename(first_post_path)}...")
            find_filtered_snvs(first_f_snvs, first_post_path)

            print(f"Counting SNVs in {ntpath.basename(first_f_snvs)}")
            num_first_Fsnvs = get_line_count(first_f_snvs)

            print(f"Counting SNVs in {'/'.join(sec_pre_path[1:].split('/')[-7:])}")
            num_sec_Rsnvs = get_gzipped_line_count(sec_pre_path)

            print(f"Finding SNVs for {ntpath.basename(sec_post_path)}")
            find_filtered_snvs(sec_f_snvs, sec_post_path)

            print(f"Counting SNVs in {ntpath.basename(sec_f_snvs)}")
            num_sec_Fsnvs = get_line_count(sec_f_snvs)

            first_r_findexes = os.path.join(_findexes_folder, f"{line_arr[0]}_raw.findexes")
            first_f_findexes = os.path.join(_findexes_folder, f"{line_arr[0]}_filtered.findexes")
            sec_r_findexes = os.path.join(_findexes_folder, f"{line_arr[1]}_raw.findexes")
            sec_f_findexes = os.path.join(_findexes_folder, f"{line_arr[1]}_filtered.findexes")

            print(f"\nFinding findexes for {'/'.join(first_pre_path[1:].split('/')[-7:])}")
            find_findexes(first_r_findexes, first_pre_path)

            print(f"Counting findexes for {ntpath.basename(first_r_findexes)}")
            num_first_Rfindexes = get_line_count(first_r_findexes)

            print(f"Finding findexes for {ntpath.basename(first_post_path)}")
            find_findexes(first_f_findexes, first_post_path)

            print(f"Counting findexes for {ntpath.basename(first_f_findexes)}")
            num_first_Ffindexes = get_line_count(first_f_findexes)

            print(f"Finding findexes for {'/'.join(sec_pre_path[1:].split('/')[-7:])}")
            find_findexes(sec_r_findexes, sec_pre_path)

            print(f"Counting findexes for {ntpath.basename(sec_r_findexes)}")
            num_sec_Rfindexes = get_line_count(sec_r_findexes)

            print(f"Finding findexes for {ntpath.basename(sec_post_path)}")
            find_findexes(sec_f_findexes, sec_post_path)

            print(f"Counting findexes for {ntpath.basename(sec_f_findexes)}")
            num_sec_Ffindexes = get_line_count(sec_f_findexes)

            comm_path = f"{imp_pairs_file[:-4]}.comm.tsv"
            with open(comm_path, 'a+') as fout:
                if os.stat(comm_path).st_size == 0:
                    fout.write(f"first_smpl        \tsec_smpl        \tfirst\tsnvs\tfindxs\tsecond\tsnvs\tfindxs\tl_200\tcomparison_type        \tcmnsnv\t%first\t%sec\tcmnfdx\t%first\t%sec\n")
                print(f"\nComparing {line_arr[0]} SNVs with {line_arr[1]} SNVs...")
                p1 = subprocess.Popen(['comm', first_f_snvs, sec_f_snvs], stdout = subprocess.PIPE)
                p2 = subprocess.Popen(['bash', awk_path], stdin = p1.stdout, stdout = subprocess.PIPE)
                p3 = subprocess.Popen(['wc', '-l'], stdin = p2.stdout, stdout = subprocess.PIPE)
                num_shared_snvs = int(p3.stdout.read().strip())

                print(f"Comparing {line_arr[0]} FingIndexes with {line_arr[1]} FingIndexes...")
                p1a = subprocess.Popen(['comm', first_f_findexes, sec_f_findexes], stdout = subprocess.PIPE)
                p2b = subprocess.Popen(['bash', awk_path], stdin = p1a.stdout, stdout = subprocess.PIPE)
                p3c = subprocess.Popen(['wc', '-l'], stdin = p2b.stdout, stdout = subprocess.PIPE)
                num_shared_findexes = int(p3c.stdout.read().strip())

                temp = '\t'.join(line_arr[0:3])
                temp2 = '\t'.join(line_arr[-2:])
                fout.write(f"{temp}\t{num_first_Fsnvs}\t" + \
                           f"{num_first_Ffindexes}\t{line_arr[3]}\t{num_sec_Fsnvs}\t" + \
                           f"{num_sec_Ffindexes}\t{temp2}\t" + \
                           f"{num_shared_snvs}\t{get_percent(num_shared_snvs,num_first_Fsnvs)}\t" + \
                           f"{get_percent(num_shared_snvs,num_sec_Fsnvs)}\t" + \
                           f"{num_shared_findexes}\t{get_percent(num_shared_findexes,num_first_Ffindexes)}\t" + \
                           f"{get_percent(num_shared_findexes,num_sec_Ffindexes)}\n")

            #Appending data to relevant arrays...
            if line_arr[2] == 'tumor':
                tumor_smpl_data[0][0].append(num_first_Rsnvs)
                tumor_smpl_data[0][1].append(num_first_Fsnvs)
                tumor_smpl_data[1][0].append(num_first_Rfindexes)
                tumor_smpl_data[1][1].append(num_first_Ffindexes)
            if line_arr[2] == 'normal':
                normal_smpl_data[0][0].append(num_first_Rsnvs)
                normal_smpl_data[0][1].append(num_first_Fsnvs)
                normal_smpl_data[1][0].append(num_first_Rfindexes)
                normal_smpl_data[1][1].append(num_first_Ffindexes)
            if line_arr[3] == 'tumor':
                tumor_smpl_data[0][0].append(num_sec_Rsnvs)
                tumor_smpl_data[0][1].append(num_sec_Fsnvs)
                tumor_smpl_data[1][0].append(num_sec_Rfindexes)
                tumor_smpl_data[1][1].append(num_sec_Ffindexes)
            if line_arr[3] == 'normal':
                normal_smpl_data[0][0].append(num_sec_Rsnvs)
                normal_smpl_data[0][1].append(num_sec_Fsnvs)
                normal_smpl_data[1][0].append(num_sec_Rfindexes)
                normal_smpl_data[1][1].append(num_sec_Ffindexes)

            if line_arr[5] == 'normal-normal':
                comp_type_data[0][0].append(num_shared_snvs)
                comp_type_data[0][1].append(num_shared_findexes)
            elif line_arr[5] == 'tumor-tumor':
                comp_type_data[1][0].append(num_shared_snvs)
                comp_type_data[1][1].append(num_shared_findexes)
            elif line_arr[5] == 'tumor-normal-pair':
                comp_type_data[2][0].append(num_shared_snvs)
                comp_type_data[2][1].append(num_shared_findexes)
            elif line_arr[5] == 'tumor-normal-unpaired':
                comp_type_data[3][0].append(num_shared_snvs)
                comp_type_data[3][1].append(num_shared_findexes)
            else:
                print(f"Comp_Type {line_arr[5]} did not match any known comp_type_data")

            curr_comp += 1

        tumor_Rsnv_summ = get_data_sum(tumor_smpl_data[0][0])
        tumor_Fsnv_summ = get_data_sum(tumor_smpl_data[0][1])
        normal_Rsnv_summ = get_data_sum(normal_smpl_data[0][0])
        normal_Fsnv_summ = get_data_sum(normal_smpl_data[0][1])

        tumor_Rfindx_summ = get_data_sum(tumor_smpl_data[1][0])
        tumor_Ffindx_summ = get_data_sum(tumor_smpl_data[1][1])
        normal_Rfindx_summ = get_data_sum(normal_smpl_data[1][0])
        normal_Ffindx_summ = get_data_sum(normal_smpl_data[1][1])

        tt_sC_summ = get_data_sum(comp_type_data[0][0])
        nn_sC_summ = get_data_sum(comp_type_data[1][0])
        tnp_sC_summ = get_data_sum(comp_type_data[2][0])
        tnup_sC_summ = get_data_sum(comp_type_data[3][0])

        tt_fC_summ = get_data_sum(comp_type_data[0][1])
        nn_fC_summ = get_data_sum(comp_type_data[1][1])
        tnp_fC_summ = get_data_sum(comp_type_data[2][1])
        tnup_fC_summ = get_data_sum(comp_type_data[3][1])

        with open(os.path.join(data_anal_dir, f"Smpl_Type_SNV_Findx-{max_or_min}_{filter}.tsv"), 'w+') as fout:
            fout.write(f"##Source: {ttnnp_file}\n")
            fout.write(f"##Filter: {max_or_min} correlation of {filter} for {comp_types}\n")
            fout.write(f"##First Number: SNV data, Second Number: Fingerprint Index data\n")
            fout.write(f"## _ = 'SNVs/Findxs'\n")
            fout.write(f"#sample_type:\tnumsmpl\tmean _\tmed. _\tstd. _\tmin. _\tmax. _\n")
            fout.write(f"Tumor (pre):\t{tumor_Rsnv_summ[0]}\t{tumor_Rsnv_summ[1]}\t{tumor_Rsnv_summ[2]}\t{tumor_Rsnv_summ[3]}\t{tumor_Rsnv_summ[4]}\t{tumor_Rsnv_summ[5]}\n")
            fout.write(f"            \t{tumor_Rfindx_summ[0]}\t{tumor_Rfindx_summ[1]}\t{tumor_Rfindx_summ[2]}\t{tumor_Rfindx_summ[3]}\t{tumor_Rfindx_summ[4]}\t{tumor_Rfindx_summ[5]}\n")
            fout.write(f"Tumor (post):\t{tumor_Fsnv_summ[0]}\t{tumor_Fsnv_summ[1]}\t{tumor_Fsnv_summ[2]}\t{tumor_Fsnv_summ[3]}\t{tumor_Fsnv_summ[4]}\t{tumor_Fsnv_summ[5]}\n")
            fout.write(f"             \t{tumor_Ffindx_summ[0]}\t{tumor_Ffindx_summ[1]}\t{tumor_Ffindx_summ[2]}\t{tumor_Ffindx_summ[3]}\t{tumor_Ffindx_summ[4]}\t{tumor_Ffindx_summ[5]}\n")
            fout.write(f"Normal (pre):\t{normal_Rsnv_summ[0]}\t{normal_Rsnv_summ[1]}\t{normal_Rsnv_summ[2]}\t{normal_Rsnv_summ[3]}\t{normal_Rsnv_summ[4]}\t{normal_Rsnv_summ[5]}\n")
            fout.write(f"             \t{normal_Rfindx_summ[0]}\t{normal_Rfindx_summ[1]}\t{normal_Rfindx_summ[2]}\t{normal_Rfindx_summ[3]}\t{normal_Rfindx_summ[4]}\t{normal_Rfindx_summ[5]}\n")
            fout.write(f"Normal (post):\t{normal_Fsnv_summ[0]}\t{normal_Fsnv_summ[1]}\t{normal_Fsnv_summ[2]}\t{normal_Fsnv_summ[3]}\t{normal_Fsnv_summ[4]}\t{normal_Fsnv_summ[5]}\n")
            fout.write(f"              \t{normal_Ffindx_summ[0]}\t{normal_Ffindx_summ[1]}\t{normal_Ffindx_summ[2]}\t{normal_Ffindx_summ[3]}\t{normal_Ffindx_summ[4]}\t{normal_Ffindx_summ[5]}\n")

        with open(os.path.join(data_anal_dir, f"Comp_Type_SNV_Findx-{max_or_min}_{filter}.tsv"), 'w+') as fout:
            fout.write(f"##Source: {ttnnp_file}\n")
            fout.write(f"##Filter: {max_or_min} correlation of {filter} for {comp_types}\n")
            fout.write(f"##First Number: SNV data, Second Number: Fingerprint Index data\n")
            fout.write(f"## _ = 'shared SNVs/Findxs'\n")
            fout.write(f"#type\tnumcomp\tmean _\tmed. _\tstd. _\tmin. _\tmax. _\n")
            fout.write(f"tt:\t{tt_sC_summ[0]}\t{tt_sC_summ[1]}\t{tt_sC_summ[2]}\t{tt_sC_summ[3]}\t{tt_sC_summ[4]}\t{tt_sC_summ[5]}\n")
            fout.write(f"   \t{tt_fC_summ[0]}\t{tt_fC_summ[1]}\t{tt_fC_summ[2]}\t{tt_fC_summ[3]}\t{tt_fC_summ[4]}\t{tt_fC_summ[5]}\n")
            fout.write(f"nn:\t{nn_sC_summ[0]}\t{nn_sC_summ[1]}\t{nn_sC_summ[2]}\t{nn_sC_summ[3]}\t{nn_sC_summ[4]}\t{nn_sC_summ[5]}\n")
            fout.write(f"   \t{nn_fC_summ[0]}\t{nn_fC_summ[1]}\t{nn_fC_summ[2]}\t{nn_fC_summ[3]}\t{nn_fC_summ[4]}\t{nn_fC_summ[5]}\n")
            fout.write(f"tnp:\t{tnp_sC_summ[0]}\t{tnp_sC_summ[1]}\t{tnp_sC_summ[2]}\t{tnp_sC_summ[3]}\t{tnp_sC_summ[4]}\t{tnp_sC_summ[5]}\n")
            fout.write(f"    \t{tnp_fC_summ[0]}\t{tnp_fC_summ[1]}\t{tnp_fC_summ[2]}\t{tnp_fC_summ[3]}\t{tnp_fC_summ[4]}\t{tnp_fC_summ[5]}\n")
            fout.write(f"tnup:\t{tnup_sC_summ[0]}\t{tnup_sC_summ[1]}\t{tnup_sC_summ[2]}\t{tnup_sC_summ[3]}\t{tnup_sC_summ[4]}\t{tnup_sC_summ[5]}\n")
            fout.write(f"     \t{tnup_fC_summ[0]}\t{tnup_fC_summ[1]}\t{tnup_fC_summ[2]}\t{tnup_fC_summ[3]}\t{tnup_fC_summ[4]}\t{tnup_fC_summ[5]}\n")

    os.remove(awk_path)


def get_percent(num, denom):
    if denom == 0:
        percent = '!DIV0'
    else:
        percent = round((num/denom) * 100, 2)
    return percent


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


def get_vcf_path(name, dir):
    path = ''
    file_found = False
    for root, dirs, files in os.walk(dir):
        if file_found is not True:
            for f in files:
                if '.vcf.gz' in f:
                    if name in root:
                        path = os.path.join(root, f)
                        file_found = True
                        break
        else:
            break
    if (path == '') or (os.path.isfile(path) is False):
         sys.exit(f"Was not able to find {name} in {dir}")
    return path


def get_gzipped_line_count(input_gzipped_vcf):
    p1 = subprocess.Popen(['gunzip', '-c', input_gzipped_vcf], stdout = subprocess.PIPE)
    p2 = subprocess.Popen(['wc', '-l'], stdin = p1.stdout, stdout = subprocess.PIPE)
    num_snvs = int(p2.stdout.read().strip()) - 294
    return num_snvs


def get_line_count(file_path):
    fin = open(file_path, 'r')
    pr = subprocess.Popen(['wc', '-l'], stdin = fin, stdout = subprocess.PIPE)
    raw_out = pr.stdout.read()
    print(raw_out)
    num_lines = int(raw_out)
    print(num_lines)
    return num_lines


def find_filtered_snvs(output_file, input_path):
    if os.path.isfile(output_file) is not True:
        with gzip.open(input_path, 'rt') as fin:
            with open(output_file, 'w+') as fout:
                for line in fin:
                    if '#' not in line:
                        line_arr = line.replace('\n', '').split('\t')
                        output = f"{line_arr[0]}-{line_arr[1]}-{line_arr[3]}-{line_arr[4]}\n"
                        fout.write(output)
    else:
        print(f"{ntpath.basename(output_file)} already exists. Skipping...")


def find_findexes(output_file, input_path):
    if os.path.isfile(output_file) is not True:
        os.makedirs(os.path.join(os.path.dirname(output_file), "temp"), exist_ok=True)
        tmp_output = os.path.join(os.path.dirname(output_file), "temp", f"{ntpath.basename(output_file)}.temp")
        with gzip.open(input_path, 'rt') as fin:
            with open(tmp_output, 'w+') as fout:
                prev = []
                for line in fin:
                    if '#' not in line:
                        line_arr = line.replace('\n','').split('\t')
                        curr = [line_arr[1], f"{line_arr[3]}{line_arr[4]}"]
                        if prev != []:
                            output = f"{round((int(curr[0]) - int(prev[0])) % 200)}-{prev[1] + curr[1]}\n"
                            fout.write(output)
                        prev = curr
        p1 = subprocess.Popen(["sort", tmp_output], stdout = subprocess.PIPE)
        fout = open(output_file, 'w+')
        p2 = subprocess.Popen(['uniq', '-u'], stdin = p1.stdout, stdout = fout)
        p2.wait()
    else:
        print(f"{ntpath.basename(output_file)} already exists. Skipping...")


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
    if arr == []:
        arr = [0]
    return [len(arr), int(round(get_arr_mean(arr))), \
            int(round(np.median(arr))), int(round(np.std(arr))), \
            int(round(np.min(arr))), int(round(np.max(arr)))]


if __name__ == "__main__":
    main()
