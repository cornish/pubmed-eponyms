'''remove_pmid_dupes.py
Remove duplicates within a base term

input file : pmid_results.csv
output files : pmid_results - dupes removed.csv

'''

__title__ = 'remove_pmid_dupes.py'
__version__ = '1.0.0'
__author__ = 'Toby C. Cornish, MD, PhD'
__license__ = 'GPLv3'
__copyright__ = 'Copyright 2020'


import os
import sys
import csv
from collections import defaultdict

input_file_path = os.path.abspath(r'pmid_results.csv')
output_file_path = os.path.abspath(r'pmid_results - dupes removed.csv') 

def main():
    # if the output file exists, delete them
    delete_if_exists(output_file_path)

    results = read_file(input_file_path)

    write_pmid_results(results)

def read_file(file_path):
    results = defaultdict(dict)
    with open(file_path,'r') as f:
        for row in csv.DictReader(f):
            # if the PMID hasn't already been added, add it
            if not row['PMID'] in results[row['base_term']]:
                results[row['base_term']][row['PMID']] = row
            else: 
                print(f"Duplicate PMID {row['PMID']} suppressed for {row['base_term']}")
    return results

def write_pmid_results(results):
    # PMIDS in results should now be unique per base_term; so we will just strip a couple
    # layers off the dict and use the original rows we read in

    fieldnames = ['base_term','term','PMID','Year','DP','TA','JT']

    with open(output_file_path,'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames, dialect='excel')
        writer.writeheader()

        for pmids in results.values():
            for pmid in pmids.values():
                writer.writerow(pmid)

def delete_if_exists(file_path):
    if os.path.exists(file_path):
        os.remove(file_path)

if __name__ == '__main__':
    main()