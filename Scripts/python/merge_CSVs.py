import os
import itertools
import csv
import sys
import fingerprint_processing_tools as fpt

current_file_path = os.path.abspath(__file__)
current_dir_path = os.path.dirname(current_file_path)
dir_with_comparisons = sys.argv[1]
all_comparisons_path = os.path.join(dir_with_comparisons, "all_fingerprint_comparisons.tsv")
all_comparisons_dir = os.path.dirname(all_comparisons_path)
create_all_comparisons = False

def main():
    #if the all_fingerprint_comparisons file exists:
    if os.path.isfile(all_comparisons_path):
        print("all_fingerprint_comparisons.tsv exists; pass...")
        create_all_comparisons = False
    else:
        #If the file doesn't exist, create it:
        all_comparisons_file = open(all_comparisons_path, "w+")
        all_comparisons_file.close()
        create_all_comparisons = True

    #search through Comparisons directory:
    for root, dir, file in os.walk(dir_with_comparisons):
        for file_name in file:

            #path of the current file in the loop:
            current_file_path = os.path.join(root, file_name)
            if ".tsv" in file_name and (file_name != "all_fingerprint_comparisons.tsv"):
                if (create_all_comparisons == True):
                    name_added_tsv_path = add_file_name_as_row(current_file_path, file_name)
                    create_all_comparisons = False
                else:
                    row_removed_tsv_path = remove_top_row(current_file_path, file_name)
                    name_added_tsv_path = add_file_name_as_row(row_removed_tsv_path, file_name)
                append_to_tsv_file(name_added_tsv_path, all_comparisons_path)

    create_all_comparisons_csv()
    return 0


'''A function that adds the file name of the comparison .tsv to the left of the
file as a column. This is necessary in order to identify the data once it is
added to the all_comparisons .tsv.'''
def add_file_name_as_row(path_to_tsv, tsv_file_name):
    name_added_dir_path = fpt.create_dir_if_absent(os.path.join(dir_with_comparisons, "tsvs_with_name_rows"))
    name_added_file_path = os.path.join(name_added_dir_path, tsv_file_name[:-4] + ".rowadd")

    fout = open(name_added_file_path, "w+")

    with open(path_to_tsv) as open_tsv_file:
        for line in open_tsv_file:
            line = tsv_file_name[:-4] + "\t" + line
            fout.write(line)

    fout.close()
    open_tsv_file.close()
    return name_added_file_path


'''A function that removes the top row of the given .tsv file so it can be appened
easily to the bottom of the all_comparisons .tsv'''
def remove_top_row(tsv_file_path, tsv_file_name):
    row_removed_dir_path = fpt.create_dir_if_absent(os.path.join(dir_with_comparisons, "row_removed_tsvs"))
    row_removed_file_path = os.path.join(row_removed_dir_path, tsv_file_name[:-4] + ".rmvd")

    fout = open(row_removed_file_path, "w+")

    current_readline = 1
    with open(tsv_file_path) as open_tsv_file:
        for line in open_tsv_file:
            if (current_readline == 1):
                current_readline += 1
            else:
                current_readline += 1
                fout.write(line)

    fout.close()
    open_tsv_file.close()
    return row_removed_file_path

'''A function that takes a smaller .tsv file and appends it to the end of the
larger file_output'''
def append_to_tsv_file(tsv_input_path, tsv_ouput_path):
    fout = open(tsv_ouput_path, "a")

    with open(tsv_input_path) as append_file:
        for line in append_file:
            line_array = line.split("\t")
            #print("Writing '" + line_array[0] + "' to all_fingerprint_comparisons.tsv")
            fout.write(line)

    append_file.close()
    fout.close()

'''A function that finds the all_comparisons.tsv and converts it into all_comparisons.csv.
In the process, it spits the first column identifying the data by comparison into two
columns, one for each sample being compared.'''
def create_all_comparisons_csv():
    fin = open(all_comparisons_path)
    fout = open(os.path.join(all_comparisons_dir, "all_comparisons.csv"), "w+")
    first_line = True

    for line in fin:
        line.replace("\n", "")
        line_array = line.split("\t")
        line_length = len(line_array)
        num_in_line = 0

        if "1.000" in line:
            continue
        else:
            while (num_in_line != line_length):

                if (num_in_line == 0): #if this is the first term in the array, we know it is the sample names
                    if (first_line == True): #if this is the first line, then we want to write headers to fout.
                        fout.write("First Sample,Second Sample,")
                        first_line = False
                    else: #if this is not the first line, split the sample names into separate columns
                        split_names = line_array[num_in_line].split("~")
                        fout.write(split_names[0] + "," + split_names[1] + ",")

                elif (num_in_line != 1 and num_in_line != 2):
                    if "\n" in line_array[num_in_line]: #if there is a newline character, do not add a comma
                        fout.write(line_array[num_in_line])
                    else:
                        fout.write(line_array[num_in_line] + ",")
                num_in_line += 1 #increment the array tracker by 1

            num_in_line = 0 #reset the array tracker to 0

    print("Created and populated all_comparisons.csv")


if __name__ == "__main__":
    main()
