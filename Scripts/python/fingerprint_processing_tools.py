import os
import itertools
import csv

current_file_path = os.path.abspath(__file__)
current_dir_path = os.path.dirname(current_file_path)

def create_dir_if_absent(dir_path):
    if (os.path.isdir(dir_path)):
        pass
    else:
        os.makedirs(dir_path)
        print("Making directory {}".format(dir_path))
    return dir_path


'''This function will find a file with a desired fill name if it is in any
directory that is x levels above the provided file path.

current_path: the path from which all others will be file_found'
name_of_desired_file: the name of the file you want to find - includes extension
levels_of_search: how many levels of directories you want to search above the current_path'''
def find_file_path(current_path, name_of_desired_file, levels_of_search):
    number_of_jumpbacks = 0 #how many times the script has added ..
    string_to_join = "" #the string that will be joined to current_path
    file_found = False #will turn to true if the desired file is found

    while (number_of_jumpbacks != levels_of_search):
        if (number_of_jumpbacks == 0): #if there have been 0 jumbacks, only add ..
            string_to_join = string_to_join + ".."
            number_of_jumpbacks += 1
        else:
            string_to_join = string_to_join + "/.." #if there has been 1 or more jumpbacks, add /..
            number_of_jumpbacks += 1

    directory_to_search = os.path.join(current_path, string_to_join)

    for r, d, f in os.walk(directory_to_search):
        for file_name in f:
            if (file_name == name_of_desired_file):
                return os.path.join(r, file_name) #join the file name with its root
                file_found = True

    if (file_found == False): #if the file is not found, notify the user
        print('Desired file {} not found in dirs {} levels above current'.format(name_of_desired_file, levels_of_search))


'''
This is a function that will determine if the string that it is given is a valid
path to a directory in the computer that it is running within. If the path points
to a directory, then the function will return 1. Otherwise, it will notify the
user and return 0.

path [string]: For the function to properly check if the path is valid, this
variable must be in single quotes.
'''
def check_if_directory_exists(path):
    if os.path.isdir(os.path.abspath(path)):
        return 1
    else:
        print("Directory either does not exist or is not readable.")
        return 0


'''
This function repeatedly asks the user for a string until the string is proven
to be a valid directory using check_if_directory_exists()

specification [string]: This parameter is added to the end of the phrase 'please
input the path,' and serves to specify which path should be inputed by the user.
'''
def ask_for_path(specification):
    directory_path = input('Please input the path {}: '.format(specification))
    path_exists = check_if_directory_exists(directory_path)
    while (path_exists != 1):
        directory_path = input('Please try again: ')
        path_exists = check_if_directory_exists(directory_path)
    print(directory_path)
    return directory_path
