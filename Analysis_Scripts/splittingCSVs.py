#!/urs/bin/python3

import os
import sys
import ntpath

input_path = sys.argv[1]
input_file_name = ntpath.basename(input_path)
input_dir = os.path.dirname(input_path)

def main():
    output_paths = []
    inputobj = os.stat(input_path)
    file_size = inputobj.st_size
    print(f"{input_file_name} is {file_size} bytes...")

    num_splits = float(file_size / 10000000)
    remainder = get_remainder_float(num_splits)

    if (remainder != 0):
        if (remainder >= 0.5):
            num_splits = round(num_splits)
        if (remainder < 0.5):
            num_splits = round(num_splits)
            num_splits += 1
    print(f"Will split {input_file_name} into {num_splits} sub-files")

    output_array = []
    current_output = 1
    while (current_output <= num_splits):
        csv_name = f"output{current_output}.csv"
        output_array.append(os.path.join(input_dir, csv_name))
        current_output += 1

    num_lines = get_numlines_or_header(input_path, 0)
    file_header = get_numlines_or_header(input_path, 1)

    lines_per_file = float(num_lines / num_splits)
    remainder = get_remainder_float(lines_per_file)
    if (remainder != 0):
        if (remainder >= 0.5):
            lines_per_file = round(lines_per_file)
        if (remainder < 0.5):
            lines_per_file = round(lines_per_file)
            lines_per_file += 1


    min_line = 1
    max_line = lines_per_file
    for file in output_array:
        fout = open(file, 'a+')
        fout.write(file_header)

        fin = open(input_path, 'r')
        current_line = 0
        for line in fin:
            if (current_line >= min_line) and (current_line <= max_line):
                fout.write(line)
            current_line +=1

        min_line += lines_per_file
        max_line += lines_per_file


def get_remainder_float(float_num):
    remainder = str(float_num).split(".")[1]
    true_remainder = ""
    current_char = 0
    for char in remainder:
        if current_char == 0:
            if char == "0":
                true_remainder += "0.0"
        else:
            true_remainder += char
        current_char += 1
    remainder = 0.0
    remainder = float(true_remainder)
    return remainder


def get_numlines_or_header(fin, num):
    finput = open(fin, 'r')
    current_line = 0
    header = ""
    for line in finput:
        if current_line == 0:
            header = line
        current_line += 1

    if num == 0:
        return current_line
    elif num == 1:
        return header

if __name__ == "__main__":
    main()
