'''rebase_terms.py
Takes a csv with paired person name(s) and test names and standardizes the
names of the people by eliminating possessives using hyphens to separate
multiple names.

input file : 323 chemistry eponyms - split - edits - utf8.csv
output file : terms_re-base.csv

'''

__title__ = 'rebase_terms.py'
__version__ = '1.0.1'
__author__ = 'Toby C. Cornish, MD, PhD'
__license__ = 'GPLv3'
__copyright__ = 'Copyright 2020'

import os
import sys
import csv
import re

input_file_path = os.path.abspath(r'gi_eponyms_split.csv')
output_file_path = os.path.abspath(r'terms_re-base.csv') 

def main():
    '''Main function'''
    # if the output files exists, delete them
    delete_if_exists(output_file_path)

    # read the csv of terms
    terms = read_file(input_file_path)

    terms_out = []

    print('='*40)
    for term in terms:
        name = term['Name(s)']
        test = term['Term']
        print(f"{name}, {test}")
        term_out = alter_name(term)
        print(f"  {name} -> {term_out['Name(s)']} {test}")
        terms_out.append(term_out)
        print('-'*40)

    # remove exact duplicates of Name(s) and Test
    terms_out = remove_duplicates(terms_out)
    
    for term in terms_out:
        print(term)

    write_file(output_file_path,terms_out)

def remove_duplicates(l):
    '''Takes a list of dicts and returns a list of unique elements'''
    return [i for n, i in enumerate(l) if i not in l[n + 1:]] 

def alter_name(term):
    '''Standardize the names in the terms'''
    name = term['Name(s)']

    # remove possesives
    altered_name = name.replace("'s","")
    altered_name = altered_name.replace("s'","s")

    # convert "and" to hyphen
    altered_name = altered_name.replace(" and ","-")

    # convert ", " to hyphen
    altered_name = altered_name.replace(", ","-")


    term['Name(s)'] = altered_name

    return term

def write_file(file_path, terms):
    '''Write the terms dict to the output file'''
    fieldnames = ['Name(s)','Test']

    with open(file_path,'w', newline='') as f:
        writer = csv.DictWriter(f,fieldnames=fieldnames,dialect='excel')
        writer.writeheader()
        for term in terms:
            writer.writerow(term)

def read_file(file_path):
    '''Read the terms from csv into a dict'''
    terms = []
    with open(file_path,'r') as f:
        terms = list(csv.DictReader(f))
    return terms

def delete_if_exists(file_path):
    '''Delete a file if it exists'''
    if os.path.exists(file_path):
        os.remove(file_path)

if __name__ == '__main__':
    main()