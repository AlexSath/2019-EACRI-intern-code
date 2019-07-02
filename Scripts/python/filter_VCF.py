import os
import gzip
import shutil
import itertools
import linecache
import sys
import fingerprint_processing_tools as fpt
import create_fingerprints_from_VCFs as cffV

current_file_path = os.path.abspath(__file__)
dir_to_search = sys.argv[1]
dir_to_copy_to = sys.argv[2]
min_allele_depth = int(sys.argv[3])
min_allele_freq = float(sys.argv[4])
min_time_genes_read = int(sys.argv[5])

class Header: #The header class contains all of the information necessary to process the VCF
    file_path = ""
    line_number = 0
    columns = []
    format_position = 0
    algorithm_types = []


def main():
    vcf_file_paths = cffV.retrieve_vcf_files(dir_to_search, 1)
    vcf_file_names = cffV.retrieve_vcf_files(dir_to_search, 0)

    for path, name in zip(vcf_file_paths, vcf_file_names):
        file_path_array = path.split("/")
        copy_dir_array = dir_to_copy_to.split("/")
        make_path_arr = []
        current_dir_level = 0
        while (current_dir_level <= len(copy_dir_array) - 2):
            make_path_arr.append(copy_dir_array[current_dir_level])
            current_dir_level += 1
        while (current_dir_level <= len(file_path_array) - 2):
            make_path_arr.append(file_path_array[current_dir_level])
            current_dir_level += 1

        make_dir_path = "/" + os.path.join(*make_path_arr)
        os.makedirs(make_dir_path, exist_ok=True)

        gunzip_path = os.path.join(make_dir_path, name + ".vcf")

        print(f"Unzipping {path} to {gunzip_path}.vcf")
        with gzip.open(path, "rb") as fin:
            with open(gunzip_path, "wb+") as fout:
                shutil.copyfileobj(fin, fout)

        fin.close()
        fout.close()

        vcf_header = get_header_info(gunzip_path) #Create a header using the desired .vcf file
        file_to_gzip = parse_vcf(vcf_header)

        print(f"Zipping {name}.vcf to {name}.vcf.gz")
        with open(gunzip_path, "rb") as fin:
            with gzip.open(file_to_gzip + ".gz", "wb") as fout:
                shutil.copyfileobj(fin, fout)

        gzip_target_array = file_to_gzip.split("/")
        print(f"Deleting {gzip_target_array[len(gzip_target_array) - 1]}...")
        os.remove(file_to_gzip)
        os.remove(gunzip_path)
        fin.close()
        fout.close()


def parse_vcf(Header):
    fin = open(Header.file_path, "r") #input is the desired .vcf file
    fout_path = f"{Header.file_path[:-4]}_filtered_{str(min_allele_depth)}_{str(min_allele_freq)}_{str(min_time_genes_read)}.vcf"
    fout = open(fout_path, "w+") #output is path/<name>_parsed.vcf; open if necessary
    passed_header = False #to know if the script is still in the header portion of the VCF file

    for line in fin:
        allele_depth_ok = False
        allele_freq_ok = False
        num_calls_ok = False
        line_allele_depth = 0
        line_allele_freq = 0.0
        time_gene_read = 0

        line_array = line.replace("\n", "").split("\t") #parse the line by tabs into an array; remove the newline character

        if (line_array == Header.columns): #if the line array is the same as the header's columns, then this line is the header
            passed_header = True #the current line is the header, so passed_header can be set to True
            fout.write(line)
        else: #if the current line is not the column header...
            if (passed_header == False):
                fout.write(line) #If the header hasn't been passed, do write to the file

            else:
                '''If the header has been passed, then the formats can be processed
                 - The line_ref_array variable is the reference for the data produced by the algorithms
                 - It is in the same tab number as the "FORMAT" in the column_array
                 - The format must be split into references by colons, according to .vcf standards'''
                line_ref_array = line_array[Header.format_position].split(":")
                line_depth_freq = determine_data_format(line_ref_array)
                location_depth_freq = [[], []]

                for r in line_depth_freq[0]: #find locations of allele depth data
                    location_index = str_num_in_array(line_ref_array, r)
                    location_depth_freq[0].append(location_index)

                for r in line_depth_freq[1]: #find locations of allele frequency data
                    location_index = str_num_in_array(line_ref_array, r)
                    location_depth_freq[1].append(location_index)

                line_alg_data = get_line_alg_data(line_array, len(line_array), Header.format_position)
                allele_depth = []
                allele_freq = []

                times_gene_absent = 0
                for alg_data in line_alg_data:
                    alg_data = alg_data.split(":")

                    if "./." in alg_data:
                        times_gene_absent += 1

                    now_depth = True
                    for string, index in zip(line_depth_freq, location_depth_freq):
                        while (now_depth == True): #When you are dealing with allele depth:
                            for s, i in zip(string, index):
                                if (s == "DP") and (alg_data[i] != "."):
                                    allele_depth.append(int(alg_data[i]))
                            now_depth = False #Tell the script you are no long in allele depth
                        if "FREQ" in string:
                            for s, i in zip(string, index):
                                if (s == "FREQ") and (alg_data[i] != "."):
                                    freq_string = alg_data[i]
                                    allele_freq.append(float(freq_string.replace("%", "")))
                        if "AD" in string and "VD" in string:
                            for s, i in zip(string, index):
                                try:
                                    if (s == "RD") and (alg_data[i] != "."):
                                        regular_reads = alg_data[i]
                                    elif (s == "AD") and (alg_data[i] != "."):
                                        alternate_reads = alg_data[i]
                                    allele_freq.append(float( (int(alternate_reads) / (int(regular_reads)+int(alternate_reads)) )*100 ))
                                except:
                                    continue
                        if "AD" in string:
                            for s, i in zip(string, index):
                                try:
                                    if (s == 'AD') and (alg_data[i] != "."):
                                        reads = alg_data[i].split(",")
                                        regular_reads = reads[0]
                                        alternate_reads = reads[1]
                                        allele_freq.append(float( (int(alternate_reads) / (int(regular_reads)+int(alternate_reads)) )*100 ))
                                except:
                                    continue
                        if "VD" in string:
                            for s, i in zip(string, index):
                                try:
                                    if (s == "VD") and (alg_data[i] != "."):
                                        allele_freq.append(float( (int(alg_data[i]) / (int(alg_data[i])+allele_depth) )*100 ))
                                except:
                                    continue
                        else:
                            continue
                    #End <string, index> for loop
                #End <alg_data> for loop
                time_gene_read = len(alg_data) - times_gene_absent
                if allele_freq is not None and allele_depth is not None:
                    line_allele_freq = average_array_value(allele_freq)
                    line_allele_depth = average_array_value(allele_depth)
            #End <passed_header> if statement
        #End <Header.columns> if statement

        if (time_gene_read >= 2) and (line_allele_freq >= min_allele_freq) and (line_allele_depth >= min_allele_depth):
            fout.write(line)
    #End <line> for loop

    fin.close()
    fout.close()
    return fout_path


def average_array_value(array):
    average = 0.0

    running_total = 0.0
    for float in array:
        running_total += float

    if (len(array) != 0):
        average = running_total/len(array)

    return average


def determine_data_format(ref_array):
    allele_depth_type = []
    allele_freq_type = []

    for ref in ref_array:
        if (len(ref) == 2):
            if (ref == "DP"):
                allele_depth_type.append(ref)
            if (ref == "AD"):
                allele_freq_type.append(ref)
            if (ref == "RD"):
                allele_freq_type.append(ref)
            if (ref == "VD"):
                allele_freq_type.append(ref)
        if (len(ref) == 4):
            if (ref == "FREQ"):
                allele_freq_type.append(ref)

    depth_freq_type = [allele_depth_type, allele_freq_type]
    return depth_freq_type


def str_num_in_array(array, str_to_find):
    str_location = 0

    current_location = 0
    for str in array:
        if (str == str_to_find):
            str_location = current_location
        current_location += 1

    return str_location


def get_line_alg_data(line_array, array_length, column_num):
    alg_data = []

    current_column = column_num + 1
    while (current_column != array_length):
        alg_data.append(line_array[current_column])
        current_column += 1

    return alg_data


def get_header_info(file_input):
    def get_header_line_number(file_input):
        line_number = 0
        fin = open(file_input, "r")
        current_line_number = 1
        for line in fin:
            if "#CHROM" in line:
                line_number = current_line_number
                break
            current_line_number += 1
        fin.close()
        return line_number

    def get_column_array(file_input, line_number):
        column_array = []
        header_line = linecache.getline(file_input, line_number)
        column_array = header_line.replace("\n", "").split("\t")
        linecache.clearcache()
        return column_array

    def get_format_position(file_input, column_array):
        format_position = 0
        curr_column = 0
        for c in column_array:
            if (c == "FORMAT"):
                format_position = curr_column
            curr_column += 1
        return format_position

    def get_header_algorithms(column_array, format_position):
        algorithm_array = []
        num_columns = len(column_array)
        current_algorithm = format_position + 1
        while (current_algorithm != num_columns):
            try:
                algorithm_array.append(column_array[current_algorithm])
            except Exception as e:
                print("Could not append to algorithm_array because {}".format(e))
                break
            current_algorithm += 1
        return algorithm_array

    new_header = Header()
    new_header.file_path = file_input
    new_header.line_number = get_header_line_number(new_header.file_path)
    new_header.columns = get_column_array(new_header.file_path, new_header.line_number)
    new_header.format_position = get_format_position(new_header.file_path, new_header.columns)
    new_header.algorithm_types = get_header_algorithms(new_header.columns, new_header.format_position)

    return new_header

if __name__ == "__main__":
    main()
