import os
import subprocess
import itertools
import fingerprint_processing_tools as fpt

current_file_path = os.path.abspath(__file__)
current_dir_path = os.path.dirname(current_file_path)
twice_parent_folder_path = os.path.join(current_dir_path, "..", "..")

def main():
    fingerprint_location = fpt.ask_for_path('where the fingerprints are found')

    fingerprint_paths = [] #array containing the path of all of the fingerprints that will be compared
    fingerprint_paths = find_fingerprint_paths(fingerprint_location, 1)
    fingerprint_names = [] #array where the names of the fingerprints that will be compared
    fingerprint_names = find_fingerprint_paths(fingerprint_location, 0)

    key_location = input('Please input the path where the .csv key is found:')
    compare_fingerprints_from_key(fingerprint_paths, fingerprint_names, key_location)
    #compare_fingerprints(fingerprint_paths, fingerprint_names)
    mk_fngp_name_file(fingerprint_location, fingerprint_names)


'''This function finds the paths of the desired fingerprints and the names of
the fingerprint, which are the names of the outn.gz files without the extension'''
def find_fingerprint_paths(fingerprint_dir_path, is_path):
    file_paths = []
    # r=root, d=direcotries, f=outFiles
    for r, d, f in os.walk(fingerprint_dir_path):
        if (is_path == 1):
            for file in f:
                if '.outn.gz' in file:
                    file_paths.append(os.path.join(r, file))
        else:
            for file in f:
                if '.outn.gz' in file:
                    file = file[:-8]
                    file_paths.append(file)
    return file_paths


def compare_fingerprints(fingerprint_path_array, fingerprint_name_array):
    comparisons_folder_path = fpt.create_dir_if_absent(os.path.join(twice_parent_folder_path, "Key_Comparisons"))
    script_path = fpt.find_file_path(current_dir_path, "compareDMFs.pl", 2)
    misc_tsvs_path = fpt.create_dir_if_absent(os.path.join(comparisons_folder_path, "misc_tsvs"))

    '''this loop goes through both name and path arrays and uses compareDMFs.pl
    to compare each fingerprint with each other fingerprint if it doesn't exist'''
    #these nested loops use Itertools to associate each path with its associated fingerprint name
    for (p, n) in zip(fingerprint_path_array, fingerprint_name_array):
        for (p2, n2) in zip(fingerprint_path_array, fingerprint_name_array):

            print('Comparing {} with {}...'.format(n, n2)) #so you know what is going on...

            output_file_path = os.path.join(misc_tsvs_path, n + "~" + n2 + ".tsv")
            reverse_file_path = os.path.join(misc_tsvs_path, n2 + "~" + n + ".tsv")
            #this if statement checks if the file already exists
            if os.path.isfile(output_file_path):
                print('File already exists, pass...')
            elif os.path.isfile(reverse_file_path):
                print('Comparison already exists, pass...')
            else: #if the comparison doesn't exist, create the csv and fill it with comparison data.
                file_output = open(output_file_path, "w+")
                subprocess.call(["perl", script_path, p, p2], stdout=file_output)


'''Creates comparisons based on a .csv file tumor_normal key.'''
def compare_fingerprints_from_key(fingerprint_path_array, fingerprint_name_array, key_path):
    comparisons_folder_path = fpt.create_dir_if_absent(os.path.join(twice_parent_folder_path, "Key_Comparisons"))
    script_path = fpt.find_file_path(current_dir_path, "compareDMFs.pl", 2)
    misc_tsvs_path = fpt.create_dir_if_absent(os.path.join(comparisons_folder_path, "misc_tsvs"))

    fin = open(key_path)
    current_line = 2
    for line in fin:
        line = line.replace("\n", "")
        name_array = line.split(',')

        first_comp_path = fpt.find_file_path(current_dir_path, name_array[0] + ".outn.gz", 2)
        scnd_comp_path = fpt.find_file_path(current_dir_path, name_array[1] + ".outn.gz", 2)

        try:
            first_path_exists = os.path.isfile(first_comp_path)
            second_path_exists = os.path.isfile(scnd_comp_path)
            if (first_path_exists and specification):
                print('Comparing {} with {}...'.format(name_array[0], name_array[1])) #so you know what is going on...

                output_file_path = os.path.join(misc_tsvs_path, n + "~" + n2 + ".tsv")
                reverse_file_path = os.path.join(misc_tsvs_path, n2 + "~" + n + ".tsv")
                #this if statement checks if the file already exists
                if os.path.isfile(output_file_path):
                    print('File already exists, pass...')
                elif os.path.isfile(reverse_file_path):
                    print('Comparison already exists, pass...')
                else: #if the comparison doesn't exist, create the csv and fill it with comparison data.
                    #file_output = open(output_file_path, "w+")
                    #subprocess.call(["perl", script_path, name_array[1], name_array[2]], stdout=file_output)
                    continue
        except:
            print("Could not compare {} with {}".format(first_path_exists, second_path_exists))
        current_line += 1


def mk_fngp_name_file(output_path, name_array):
    fout = open(os.path.join(output_path, 'fingerprint_names.csv'), "w+")

    array_length = len(name_array)
    current_name = 0
    fout.write('sample_id,sample_number\n')
    while (array_length != current_name):
        cn_string = str(current_name + 1)
        fout.write(name_array[current_name] + ',' + cn_string + '\n')
        current_name += 1

    fout.close()


if __name__ == "__main__":
    main()
