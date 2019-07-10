#!/usr/bin/python3

import sys
import itertools

file_to_edit = sys.argv[1]
key_file = sys.argv[2]

def main():
    first_sample_location = 0
    second_sample_location = 0
    first_sample_array = []
    second_sample_array = []
    fin = open(key_file, 'r')

    current_line = 1
    for line in fin:
        line = line.split(",")
        if (current_line == 1):
            if "First_Sample" not in line and "Second_Sample" not in line:
                print("Was not able to find columns 'First_Sample' and/or 'Second_Sample' in the first line of the file")
            else:
                current_column = 0
                for column in line:
                    if (column == "First_Sample"):
                        first_sample_location = current_column
                    elif (column == "Second_Sample"):
                        second_sample_location = current_column
                    current_column += 1
        else:
            first_sample_array.append(line[first_sample_location])
            second_sample_array.append(line[second_sample_location])
        current_line += 1
    fin.close()

    fin = open(file_to_edit, 'r')
    current_line = 1
    input_array = []
    for line in fin:
        pair_found = False
        line = line.replace("\n", "")
        line_array = line.split(",")
        if (current_line == 1):
            line += ",is_pair\n"
        else:
            for f, s in zip(first_sample_array, second_sample_array):
                if ((line_array[0] == f) and (line_array[1] == s)) or \
                   ((line_array[1] == f) and (line_array[0] == s)):
                    pair_found = True
            if pair_found:
                line += ",true\n"
            elif not pair_found:
                line += ",false\n"
        input_array.append(line)
        current_line += 1
    fin.close()

    fout = open(file_to_edit, 'w')
    for line in input_array:
        fout.write(line)
    fout.close()


if __name__ == "__main__":
    main()
