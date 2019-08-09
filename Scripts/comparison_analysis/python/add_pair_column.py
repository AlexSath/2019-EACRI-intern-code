#!/usr/bin/python3

import os
import sys
import ntpath
import itertools
import re

comp_csv = sys.argv[1]
tnp_key = sys.argv[2]
swp = '_filtered'

def main():
    first_sample_location = 0
    second_sample_location = 0
    tumor_array = []
    normal_array = []

    current_line = 1
    with open(tnp_key, 'r') as fin:
        for line in fin:
            line_arr = line.replace('\n', '').split(",")
            if (current_line == 1):
                if "tumor_id" not in line and "normal_id" not in line:
                    print("Was not able to find columns 'tumor_id' and/or 'normal_id' in the first line of the file")
                else:
                    current_column = 0
                    for column in line_arr:
                        if (column == "tumor_id"):
                            tumor_col_id = current_column
                        elif (column == "normal_id"):
                            normal_col_id = current_column
                        current_column += 1
            else:
                tumor_array.append(line_arr[tumor_col_id])
                normal_array.append(line_arr[normal_col_id])
            current_line += 1

    print(tumor_array)
    print(normal_array)

    with open(comp_csv, 'r') as fin:
        with open(f"{os.path.dirname(comp_csv)}/_f{ntpath.basename(tnp_key[:-4])}/{ntpath.basename(comp_csv)[:-4]}.pc.csv", 'w+') as pc:
            with open(f"{os.path.dirname(comp_csv)}/_f{ntpath.basename(tnp_key[:-4])}/{ntpath.basename(comp_csv)[:-4]}.ttnnp.csv", 'w+') as ttnnp:
                current_line = 1
                for line in fin:
                    line = line.replace('\n', '')
                    if current_line == 1:
                        pc.write(f"{line},is_pair\n")
                        ttnnp.write(f"{line},comp_type\n")
                    else:
                        pair_not_found = True
                        line = line.replace("\n", "")
                        line_arr = line.split(",")
                        f = line_arr[0].split(swp)[0]
                        s = line_arr[1].split(swp)[0]
                        if ((s in tumor_array) or (s in normal_array)) and ((f in tumor_array) or (f in normal_array)):
                            for t, n in zip(tumor_array, normal_array):
                                if ((f == t) and (s == n)) or ((f == n) and (s == t)):
                                    pc.write(f"{line},true\n")
                                    ttnnp.write(f"{line},tumor-normal-pair\n")
                                    pair_not_found = False
                                    break
                            if pair_not_found:
                                pc.write(f"{line},false\n")
                                if ((f in tumor_array) and (s in normal_array)) \
                                     or ((f in normal_array) and (s in tumor_array)):
                                    ttnnp.write(f"{line},tumor-normal-unpaired\n")
                                elif (f in tumor_array) and (s in tumor_array):
                                    ttnnp.write(f"{line},tumor-tumor\n")
                                elif (f in normal_array) and (s in normal_array):
                                    ttnnp.write(f"{line},normal-normal\n")
                    current_line += 1

if __name__ == "__main__":
    main()
