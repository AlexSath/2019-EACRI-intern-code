#!/usr/bin/python
import fingerprint_processing_tools as fpt
import itertools
import ntpath
import csv
import sys
import os

current_file_path = os.path.abspath(__file__)
current_dir_path = os.path.dirname(current_file_path)
dir_with_comparisons = sys.argv[1]
all_comparisons_path = os.path.join(dir_with_comparisons, "all_fingerprint_comparisons.csv")
all_comparisons_dir = os.path.dirname(all_comparisons_path)

def main():
    #if the all_fingerprint_comparisons file exists:
    for fname in os.listdir(current_dir_path):
        if '.csv' in fname:
            sys.exit(f"A .csv file was found in the directory")
    else:
        fout = open(all_comparisons_path, "w+")
        #search through Comparisons directory:
        for root, dir, file in os.walk(dir_with_comparisons):
            num_files = len(file)
            file_num = 1
            for file_name in file:
                file_comparison = file_name.replace(".tsv", "").split("~")
                current_file_path = os.path.join(root, file_name)

                if ".tsv" in file_name and (file_name != "all_fingerprint_comparisons.csv"):
                    print(f"Currently processing {file_name}... [{file_num}/{num_files}]")
                    fin = open(current_file_path, 'r')

                    if (file_num == 1):
                        current_line = 1
                        for line in fin:
                            if (current_line == 1):
                                fout.write("First_Sample,Second_Sample," + line.replace("\t", ",").replace("=","_"))
                                current_line += 1
                            else:
                                fout.write(f"{file_comparison[0]},{file_comparison[1]}," + line.replace("\t", ","))
                                current_line += 1
                    else:
                        current_line = 1
                        for line in fin:
                            if (current_line == 1):
                                current_line += 1
                            else:
                                fout.write(f"{file_comparison[0]},{file_comparison[1]}," + line.replace("\t", ","))
                                current_line += 1
                file_num += 1


if __name__ == "__main__":
    main()
