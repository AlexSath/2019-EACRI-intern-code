#!/usr/bin/python
#Created by Alex R. Sathler for the EACRI
#Perl Script for computing genome fingerprints found at /gglusman/genome-fingerprints/ on GitHub
import os
import gzip
import shutil
import itertools
import linecache
import sys
import fingerprint_processing_tools as fpt
import create_fingerprints_from_VCFs as cffV

current_file_path = os.path.abspath(__file__)
dir_to_search = sys.argv[1] #The folder where the raw .vcf.gz files are found. Nested folders O.K.
dir_to_copy_to = sys.argv[2]  #The folder in which the filtered .vcf.gz files will output.
# - - - - - - - - - - - - - - #File/Folder hierarchies will be copied from the first folder
min_allele_depth = int(sys.argv[3]) #Minimum number of reads for a line to be copied to the new file
min_allele_freq = float(sys.argv[4]) #Minimum percentage of variant reads for a line to be copied to the new file
try:
    min_time_genes_read = int(sys.argv[5]) #Minimum number of reads at the locus for a line to be copied to the new file
except:
    min_time_genes_read = 0
try:
    max_population_frequency = float(sys.argv[6])
except:
    max_population_frequency = 1.0

class Header: #The header class contains all of the information necessary to process the VCF
    file_path = ""
    line_number = 0 #The file line number at which the header is located
    columns = []
    format_position = 0 #The position of the 'FORMAT' column
    algorithm_types = [] #An array of the different algorithms that were used to read the genome.


def main():
    print('Retrieving VCF file paths...') #Retrieve all .vcf.gz paths from the desired folder
    vcf_file_paths = cffV.retrieve_vcf_files(dir_to_search, 0)
    print('Retrieving VCF file names...') #Retrieve the names of the .vcf.gz files
    vcf_file_names = cffV.retrieve_vcf_files(dir_to_search, 2)

    #Main loop of this file - iterates through filepath/filename pairs and
    #finds the desired output location for the new .vcf file. It then creates
    #a header for the input file, filters the input, outputs to the desired
    #directory with the same folder hierarchy, and gzips the output file.
    #NOTE: THIS ASSUMES THE INPUT DIRECTORY AND OUTPUT DIRECTORY ARE AT THE SAME LEVEL IN THE HIERARCHY
    current_VCF = 1
    for path, name in zip(vcf_file_paths, vcf_file_names):
        file_path_array = path.split("/")
        copy_dir_array = dir_to_copy_to.split("/")
        make_path_arr = [] #Declaration of array with the new path
        current_dir_level = 0
        #Use the path of the ouput directory as a parent for the other dirs
        while (current_dir_level <= len(copy_dir_array) - 1):
            make_path_arr.append(copy_dir_array[current_dir_level])
            current_dir_level += 1
        #Once you have reached the output directory level of the hierarchy,
        #append the dir structure of the input
        while (current_dir_level <= len(file_path_array) - 2):
            make_path_arr.append(file_path_array[current_dir_level])
            current_dir_level += 1

        #Creating the full dir path of the file output:
        make_dir_path = "/" + os.path.join(*make_path_arr)
        os.makedirs(make_dir_path, exist_ok=True)

        #Defining the .vcf output path of the gzipped raw vcf:
        gunzip_output_path = os.path.join(make_dir_path, name + ".vcf")

        #Opening the raw .vcf.gz into the desired output path
        with gzip.open(path, "rb") as fin:
            with open(gunzip_output_path, "wb+") as fout:
                shutil.copyfileobj(fin, fout)
        fin.close()
        fout.close()

        #Console info, declaring the Header for the raw .vcf to be filtered,
        #filtering the raw .vcf and returning the path of the filter process output
        print(f"Filtering {name}.vcf... [{current_VCF}/{len(vcf_file_names)}]")
        vcf_header = get_header_info(gunzip_output_path) #Create a header using the desired .vcf file
        file_to_gzip = parse_vcf(vcf_header)

        #Gzipping the new filtered .vcf file
        with open(file_to_gzip, "rb") as fin:
            with gzip.open(file_to_gzip + ".gz", "wb") as fout:
                shutil.copyfileobj(fin, fout)
        fin.close()
        fout.close()

        #Removing unnecessary files and iterating to the next version of the loop
        os.remove(gunzip_output_path)
        os.remove(file_to_gzip)
        current_VCF += 1


'''A complex .vcf filter function that reads every line (which contains one read
at a specified locus) and parses the line into the information needed to decide
if the line should be outputted into the new .vcf file. The allele depth, allele
frequency, and number of allele read cutoffs are arguments during the script call.

Header: A header-class object who's path is used to determine the output .vcf
file path
'''
def parse_vcf(Header):
    fin = open(Header.file_path, "r") #input is the desired .vcf file
    fout_path = f"{Header.file_path[:-4]}_filtered_{str(min_allele_depth)}" + \
                f"_{str(min_allele_freq)}" + \
                f"_{str(min_time_genes_read)}" + \
                f"_{str(max_population_frequency)}.vcf"
    fout = open(fout_path, "w+") #output is path/<name>_parsed.vcf; open if necessary
    passed_header = False #to know if the script is still in the header portion of the VCF file

    for line in fin:
        #Declaring necessary variables:
        allele_depth_ok = False
        allele_freq_ok = False
        num_calls_ok = False
        pop_freq_ok = False
        line_allele_depth = 0
        line_allele_freq = 0.0
        time_gene_read = 0
        line_alg_data = []

        #Split tab-separated line into an array and remove newline characters:
        line_array = line.replace("\n", "").split("\t")

        #if the line array is the same as the header's columns, then this line is the header:
        if (line_array == Header.columns):
            passed_header = True #the current line is the header, so passed_header can be set to True
            fout.write(line)
        else: #if the current line is not the column header...
            if (passed_header == False):
                #If the header hasn't been passed, do write to the file:
                fout.write(line)

            else:
                # If the header has been passed, then the formats can be processed
                # - The line_ref_array variable is the reference for the data produced by the algorithms
                # - It is in the same tab number as the "FORMAT" in the file's column_array
                # - The format must be split into references by colons, according to .vcf standards
                line_ref_array = line_array[Header.format_position].split(":")
                line_depth_freq = determine_data_format(line_ref_array)
                location_depth_freq = [[], []]

                #find locations of allele depth data:
                for r in line_depth_freq[0]:
                    location_index = str_num_in_array(line_ref_array, r)
                    location_depth_freq[0].append(location_index)

                #find locations of allele frequency data:
                for r in line_depth_freq[1]:
                    location_index = str_num_in_array(line_ref_array, r)
                    location_depth_freq[1].append(location_index)

                #Getting read algorithm data and declaring depth and freq arrays
                line_alg_data = get_line_alg_data(line_array, len(line_array), Header.format_position)
                line_alg_data = line_alg_data[1:4]
                allele_depth = []
                allele_freq = []


                #Loop that iterates through each piece of algorithm data and parses
                #it to find what each algorithm says about the line read in terms of
                #depth and freqeuency
                times_gene_absent = 0
                for alg_data in line_alg_data:
                    alg_data = alg_data.split(":")

                    #If "./." is in the data, then log that the read is absent
                    if "./." in alg_data:
                        times_gene_absent += 1

                    now_depth = True
                    #Iterate through string/string location pairs within their
                    #respective arrays and search for key strings. When a key
                    #string is found, use its location index to extract its
                    #corresponding data
                    for string, index in zip(line_depth_freq, location_depth_freq):
                        #When you are dealing with allele depth:
                        while (now_depth == True):
                            for s, i in zip(string, index):
                                if (s == "DP") and (alg_data[i] != "."):
                                    allele_depth.append(int(alg_data[i]))
                            now_depth = False #Tell the script you are no longer in allele depth
                        #When moving on to allele frequency from allele depth:
                        if "FREQ" in string:
                            for s, i in zip(string, index):
                                if (s == "FREQ") and (alg_data[i] != "."):
                                    freq_string = alg_data[i]
                                    #FREQ uses % format - get rid of the special character:
                                    allele_freq.append(float(freq_string.replace("%", "")))
                        #When both 'AD' and 'VD' are in the string, one is for regular reads,
                        #the other for variant reads
                        if "AD" in string and "VD" in string:
                            for s, i in zip(string, index):
                                try:
                                    #Assign the necessary variables:
                                    if (s == "RD") and (alg_data[i] != "."):
                                        regular_reads = alg_data[i]
                                    elif (s == "AD") and (alg_data[i] != "."):
                                        alternate_reads = alg_data[i]
                                    #Calculate the allele freq (=alt/(norm+alt)):
                                    allele_freq.append(float( (int(alternate_reads) / (int(regular_reads)+int(alternate_reads)) )*100 ))
                                except:
                                    continue
                        #If only 'AD' is in the string, normal and variant reads are separated by a comma:
                        if "AD" in string:
                            for s, i in zip(string, index):
                                try:
                                    if (s == 'AD') and (alg_data[i] != "."):
                                        #Split by the comma and assign normal and variant data to variables
                                        reads = alg_data[i].split(",")
                                        regular_reads = reads[0]
                                        alternate_reads = reads[1]
                                        #Calculate the allele frequency as in the above case:
                                        allele_freq.append(float( (int(alternate_reads) / (int(regular_reads)+int(alternate_reads)) )*100 ))
                                except:
                                    continue
                        #If only 'VD' is found, then it is expressing the number of variant reads
                        #and total allele depth will need to be leveraged for allele frequency calculation
                        if "VD" in string:
                            for s, i in zip(string, index):
                                try:
                                    if (s == "VD") and (alg_data[i] != "."):
                                        #Allele freq = alt/(depth+alt)
                                        allele_freq.append(float( (int(alg_data[i]) / (int(alg_data[i])+allele_depth) )*100 ))
                                except:
                                    continue
                        else:
                            continue
                    #End <string, index> for loop
                #End <alg_data> for loop

                #Calculate the number of reads at the locus:
                time_gene_read = len(alg_data) - times_gene_absent
                #Average the deth and frequency values stored in their respective arrays:
                if allele_freq is not None and allele_depth is not None:
                    line_allele_freq = average_array_value(allele_freq)
                    line_allele_depth = average_array_value(allele_depth)

                #The following assumes the INFO column comes before the FORMAT column
                exac_all_ok = False
                thousand_g_ok = False

                line_info_array = line_array[Header.format_position - 1].split(";")
                for info in line_info_array:
                    info = info.lower().strip()
                    if "exac_all" in info:
                        info_array = info.split("=")
                        try:
                            if (float(info_array[1]) <= max_population_frequency):
                                exac_all_ok = True
                        except:
                            if (info_array[1] == "."):
                                exac_all_ok = True
                    elif '1000g2014oct_all' in info:
                        info_array = info.split("=")
                        try:
                            if (float(info_array[1]) <= max_population_frequency):
                                thousand_g_ok = True
                        except:
                            if (info_array[1] == "."):
                                thousand_g_ok = True

                if (thousand_g_ok == True) and (exac_all_ok == True):
                    pop_freq_ok = True

            #End <passed_header> if statement
        #End <Header.columns> if statement

        #If the line meets the necessary conditions, write it to the output file:
        if (time_gene_read >= min_time_genes_read) and \
           (line_allele_freq >= min_allele_freq) and \
           (line_allele_depth >= min_allele_depth) and \
           (pop_freq_ok == True):
            fout.write(line)

    #End <line> for loop

    fin.close()
    fout.close()
    return fout_path

'''Function that takes the integers or floats in an array and reterns their average'''
def average_array_value(array):
    average = 0.0

    running_total = 0.0
    for float in array:
        running_total += float

    if (len(array) != 0):
        average = running_total/len(array)

    return average

'''Function that sorts through a data format array and looks for key two letter
codes. These codes are stored in a two-dimentional array. First dimention for
letters to do with read depth, second dimention to do with allele frequency'''
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

'''A function that takes an array and find the number location of the string
in the array. Will return 0 if the string is not found.'''
def str_num_in_array(array, str_to_find):
    str_location = 0

    current_location = 0
    for str in array:
        if (str == str_to_find):
            str_location = current_location
        current_location += 1

    return str_location

'''A function that uses the column array and the format column number to return
all of the algorithm data after the format column'''
def get_line_alg_data(line_array, array_length, column_num):
    alg_data = []

    current_column = column_num + 1
    while (current_column != array_length):
        alg_data.append(line_array[current_column])
        current_column += 1

    return alg_data

'''This function initializes a Header object by filling its variables with data
from a file provided with the file_input variable'''
def get_header_info(file_input):
    #Finds the line number at which the header is found
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

    #Finds header line and separates the individual columns by tabs and stores them in an array
    def get_column_array(file_input, line_number):
        column_array = []
        header_line = linecache.getline(file_input, line_number)
        column_array = header_line.replace("\n", "").split("\t")
        linecache.clearcache()
        return column_array

    #Finds the position of the 'Format' column (what number is it in the header?)
    def get_format_position(file_input, column_array):
        format_position = 0
        curr_column = 0
        for c in column_array:
            if (c == "FORMAT"):
                format_position = curr_column
            curr_column += 1
        return format_position

    #Finds the names of the .vcf read algorithms, which are found after the 'Format column'
    #Stores this information in an array
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

    new_header = Header() #Initialize the new header object
    #Populating the new Header object with variables:
    new_header.file_path = file_input
    new_header.line_number = get_header_line_number(new_header.file_path)
    new_header.columns = get_column_array(new_header.file_path, new_header.line_number)
    new_header.format_position = get_format_position(new_header.file_path, new_header.columns)
    new_header.algorithm_types = get_header_algorithms(new_header.columns, new_header.format_position)

    #Function returns the new Header object, fully populated
    return new_header

if __name__ == "__main__":
    main()
