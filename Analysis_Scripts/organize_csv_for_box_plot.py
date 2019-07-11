#!/usr/bin/python3
# ---! RENDERED OBSOLETE BY THE 'MELT' FUNCTION !--- #

import os
import sys
import ntpath
import itertools

input_path = sys.argv[1]
output_folder = os.path.dirname(input_path)
output_path = os.path.join(output_folder, ntpath.basename(input_path)[:-4] + ".boxed.csv")


def main():
    fin = open(input_path, 'r')
    hdr_arr = []
    key_col_arr = []
    key_col_index = []

    current_line = 1
    for line in fin:
        if (current_line == 1):
            hdr_arr = line.replace("\n", "").split(",")
            break
    fin.close()

    num_cols = len(hdr_arr) - 1
    pssd_bin = False
    pssd_L200 = False
    bin_col_num = 0
    L200_col_num = 0

    curr_col = 0
    while not pssd_L200:
        if (hdr_arr[curr_col] == "binary"):
            bin_col_num = curr_col
            pssd_bin = True
        if pssd_bin:
            key_col_arr.append(hdr_arr[curr_col])
            key_col_index.append(curr_col)
        if (hdr_arr[curr_col] == "L_200"):
            L200_col_num = curr_col
            pssd_L200 = True
        curr_col += 1

    print(key_col_arr)
    print(key_col_index)
    fin = open(input_path, 'r')
    fout = open(output_path, 'w+')
    offset = bin_col_num - 1

    current_line = 1
    for line in fin:
        if (current_line == 1):
            fout.write("l_num,correlation\n")
        else:
            line = line.replace("\n", "").split(",")
            curr_col2 = 0
            for column, index in zip(key_col_arr, key_col_index):
                fout.write(f"{column},{line[index]}\n")
        current_line += 1
    fin.close()
    fout.close()


if __name__ == "__main__":
    main()
