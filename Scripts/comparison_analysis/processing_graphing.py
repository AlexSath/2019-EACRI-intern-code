#!/usr/bin/python

import os
import sys
import ntpath
import subprocess

script_path = os.path.abspath(__file__)
curr_dir = os.path.dirname(__file__)
input_csv = sys.argv[1]
tnp_csv_key = sys.argv[2]
out_dir = os.path.dirname(input_csv)
new_out_dir = os.path.join(out_dir, f"_f{ntpath.basename(tnp_csv_key[:-4])}")

def main():
    os.makedirs(new_out_dir)
    subprocess.call(['python3', os.path.join(curr_dir, 'python', 'add_pair_column.py'), input_csv, tnp_csv_key])
    subprocess.call(['Rscript', os.path.join(curr_dir, 'R', 'violinplot.R'), \
                     f"{os.path.dirname(input_csv)}/_f{ntpath.basename(tnp_csv_key[:-4])}/{ntpath.basename(input_csv)[:-4]}.pc.csv", new_out_dir])
    subprocess.call(['Rscript', os.path.join(curr_dir, 'R', 'comptype.violin.R'), \
                     f"{os.path.dirname(input_csv)}/_f{ntpath.basename(tnp_csv_key[:-4])}/{ntpath.basename(input_csv)[:-4]}.ttnnp.csv", new_out_dir])

if __name__ == '__main__':
    main()
