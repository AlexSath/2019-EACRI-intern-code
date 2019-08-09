#!/usr/bin/python3

import os
import sys
import itertools
import ntpath

script_path = os.path.abspath(__file__)
input_path = sys.argv[1]
input_name = ntpath.basename(input_path)
tnp_key_path = sys.argv[2]
output_folder = sys.argv[3]
output_path = os.path.join(output_folder, input_name[:-4] + ".ttnnp.csv")
sep = '_filtered'

def main():
    input_header = get_header_arr(input_path)
    tnp_header = get_header_arr(tnp_key_path)
    tumor_arr = get_col_data('tumor_id', tnp_header, tnp_key_path)
    normal_arr = get_col_data('normal_id', tnp_header, tnp_key_path)
    print(f"{tumor_arr}\n{normal_arr}")

    fin = open(input_path, 'r')
    fout = open(output_path, 'w+')

    f_sample_num = get_col_num('First_Sample', input_header)
    s_sample_num = get_col_num('Second_Sample', input_header)
    pair_col_num = get_col_num('is_pair', input_header)

    curr_line = 1
    for line in fin:
        line = line.replace("\n", "")
        line_arr = line.replace("\n", "").split(",")
        if curr_line == 1:
            fout.write(f"{line},comp_type\n")
        else:
            if (line_arr[f_sample_num].split(sep)[0] in tumor_arr) \
               and (line_arr[s_sample_num].split(sep)[0] in tumor_arr):
                fout.write(f"{line},tumor-tumor\n")
            elif (line_arr[f_sample_num].split(sep)[0] in normal_arr) \
                 and (line_arr[s_sample_num].split(sep)[0] in normal_arr):
                fout.write(f"{line},normal-normal\n")
            elif ((line_arr[f_sample_num].split(sep)[0] in tumor_arr) \
                 and (line_arr[s_sample_num].split(sep)[0] in normal_arr)) \
                 or ((line_arr[f_sample_num].split(sep)[0] in normal_arr) \
                 and (line_arr[s_sample_num].split(sep)[0] in tumor_arr)):
                if (line_arr[pair_col_num] == 'true'):
                    fout.write(f"{line},tumor-normal-pair\n")
                else:
                    fout.write(f"{line},tumor-normal-unpaired\n")
        curr_line += 1

    fin.close()
    fout.close()


def get_header_arr(file_path):
    fin = open(file_path, 'r')
    header_arr = []
    curr_line = 1
    for line in fin:
        if curr_line == 1:
            header_arr = line.replace("\n", "").split(",")
            break
    return header_arr


def get_col_num(colname, header_arr):
    col_num = 0
    found_col = False
    for col in header_arr:
        if col == colname:
            found_col = True
            break
        col_num += 1
    if not found_col:
        print(f"Column {colname} not found in {header_arr}, please try again...")
        exit()
    return col_num


def get_col_data(colname, header_arr, filepath):
    fin = open(filepath, 'r')
    col_data = []

    col_num = get_col_num(colname, header_arr)

    curr_line = 1
    for line in fin:
        line = line.replace("\n", "").split(",")
        if curr_line != 1:
            col_data.append(line[col_num])
        curr_line += 1
    return col_data


if __name__ == "__main__":
    main()
