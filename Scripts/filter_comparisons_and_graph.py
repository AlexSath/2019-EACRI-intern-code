#!/usr/bin/python

import os
import sys
import ntpath
import itertools
import subprocess

script_path = os.path.abspath(__file__)
script_dir = os.path.dirname(script_path)

if '--help' in sys.argv:
    with open(os.path.join(script_dir, "helptext", "filter_comparisons_and_graph.txt"), 'r') as fin:
        for line in fin:
            print(line)
    sys.exit()
else:
    tnp_key = sys.argv[1]
    tnp_name = ntpath.basename(tnp_key)[:-4]
    dir_to_search = sys.argv[2]
    swp = '_filtered'

def main():
    tumor_arr = []
    normal_arr = []

    with open(tnp_key, 'r') as fin:
        for line in fin:
            if '#' not in line:
                tumor_arr.append(line.split(',')[0])
                normal_arr.append(line.split(',')[1])

    for root, dir, file in os.walk(dir_to_search):
        for f in file:
            if '.csv' in f \
               and 'tnp' not in f \
               and 'ttnn' not in f \
               and 'ttnnp' not in f \
               and '_f' not in f:
                file_path = os.path.join(root, f)
                fout_dir = os.path.join(root, f"_f{tnp_name}")
                os.makedirs(fout_dir)

                fout_path = os.path.join(fout_dir, f"{ntpath.basename(file_path)[:-4]}_f{tnp_name}.csv")
                fout = open(fout_path, 'w+')

                with open(file_path, 'r') as fin:
                    curr_line = 1
                    for line in fin:
                        if 'First_Sample' in line:
                            fout.write(line)
                        else:
                            line_arr = line.split(',')
                            first = line_arr[0].split(swp)[0]
                            second = line_arr[1].split(swp)[0]
                            if (first in tumor_arr or first in normal_arr):
                                if (second in tumor_arr or second in normal_arr):
                                    line = ','.join(line_arr[0:-1])
                                    for t, n in zip(tumor_arr, normal_arr):
                                        if ((first == t and second == n) or (first == n and second == t)):
                                            line = f"{line},true\n"
                                            break
                                    if 'true' not in line:
                                        line = f"{line},false\n"
                                    print(line)
                                    fout.write(line)
                        curr_line += 1

                subprocess.call(['python3', os.path.join(script_dir, 'Analysis_Scripts', 'python', 'processing_graphing2.py'), \
                                 fout_path, tnp_key])

if __name__ == '__main__':
    main()
