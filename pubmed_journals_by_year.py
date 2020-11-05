'''pubmed_journals_by_year.py
This script will use the input file to create a list of journal abbreviations
and a range of years from start_year to end_year. It will then search pubmed
for each of those journals over the range of years to determine how many
papers each journal published each year in the range. 

Please note that this script is known to be unstable for unknown reasons, but
it was designed to be restarted. It will automatically skip any journals that
are already present in the output file. Apologies for that; it is really not 
clear why the script occassionally locks up completely. This may be a windows
issue.

input_file: 'pmid_results.csv'
output_file: 'journal_counts.csv' 

'''

__title__ = 'pubmed_journals_by_year.py'
__version__ = '1.0.0'
__author__ = 'Toby C. Cornish, MD, PhD'
__license__ = 'GPLv3'
__copyright__ = 'Copyright 2020'

import os
import sys
import csv
import time
from collections import defaultdict
from configparser import ConfigParser

from Bio import Entrez
from Bio import Medline

input_file_path = os.path.abspath(r'pmid_results.csv')
journal_output_file_path = os.path.abspath(r'journal_counts.csv') 

def main():
    '''Main function'''
    # read our config file
    config = read_config_file('config.ini')
    if not config['use_api_key']:
        config['api_key'] = None
    print(config)
    
    # Read the journals to search from the input file
    start_year, end_year, journals = read_file(input_file_path)

    # Read the journals we have finished with from the output file
    finished = read_progress_from_file(journal_output_file_path)

    # remove the finished files from the journals list 
    journals = [j for j in journals if j not in finished]

    print('='*40)
    for journal in journals:
        print(f'{journal}')
        result = []
        chunk_size = 20
        # 99% of the journals will be sparsely populated over the range, so
        # we'll check a range of years first, then get specific years if needed
        print('='*40)
        for batch in chunks(range(start_year,end_year+1),chunk_size):
            years = list(batch)
            start = min(years)
            end = max(years)
            print('-'*40)
            search_string = f'''"{journal}" [Journal] AND "{start}"[PPDAT]:"{end}"[PPDAT]'''
            print(f'   {start} - {end}: {search_string}')           
            count = search(search_string,config['email'],config['api_key'])
            print(f'    Found {count} publications.')
            if count > 0:
                for year in years:
                    search_string = f'''"{journal}" [Journal] AND "{year}"[PPDAT]'''
                    print(f'  {year} : {search_string}')
                    count = search(search_string,config['email'],config['api_key'])
                    print(f'    Found {count} publications.')
                    result.append(count)
            else:
                result.extend([0]*len(years))
        write_journal_result(journal,result,start_year,end_year)

def search(query, email, api_key=None):
    '''Use Entrez.esearch to identify all publications in a journal, return count'''

    Entrez.email = email
    if api_key:
        Entrez.api_key = api_key

    # see https://www.ncbi.nlm.nih.gov/books/NBK25499/#chapter4.ESearch for information on parameters
    print("  Retrieving PMIDs ...")

    while(True):
        try:
            retstart = 0
            handle = Entrez.esearch(db='pubmed', 
                                    sort='relevance', 
                                    retmax='100000',
                                    retstart=retstart,
                                    retmode='xml',
                                    term=query
            )
            results = Entrez.read(handle)
            break
        except KeyboardInterrupt:
            sys.exit()
        except:
            # this will probably work itself out
            print('  Error occurred, retrying...')
            time.sleep(1)

    count = int(results['Count'])

    print("    Done.")
    return count

def chunks(l, n): 
    '''Generator for list l as chunks of size n'''
    for i in range(0, len(l), n):  
        yield l[i:i + n] 

def read_file(file_path):
    '''Read the input file and return a list of unique journals, sorted''' 
    journals = []
    start_year = 3000
    end_year = 0
    with open(file_path,'r') as f:
        reader = csv.DictReader(f)
        for row in reader:
            if row['TA'] != '':
                journals.append(row['TA'])
            if is_int(row['Year']):
                start_year = min(start_year,int(row['Year']))
                end_year = max(end_year,int(row['Year']))

    # make the list unique and sort alphabetically (case insentive)
    journals = sorted(list(set(journals)), key=str.casefold)
    return start_year, end_year, journals

def is_int(val):
    '''Return True if the string can be convereted to an int'''
    try:
        num = int(val)
    except ValueError:
        return False
    return True

def read_progress_from_file(file_path):
    '''Read the output file and return a list of journals that are done'''
    if os.path.isfile(file_path):
        with open(file_path,'r') as f:
            reader = csv.DictReader(f)
            return [row['journal'] for row in reader] 
    else: 
        return []

def write_journal_result(journal,result,start_year,end_year):
    '''Append our results to the output file'''
    file_exists = os.path.isfile(journal_output_file_path) # if it doesn't exist, we are creating it and need to write a header

    fieldnames = ['journal',]
    fieldnames.extend(list(range(start_year,end_year+1)))

    with open(journal_output_file_path,'a', newline='') as f: # mode is 'append'
        csv_writer = csv.writer(f, dialect='excel')

        if not file_exists:
            csv_writer.writerow(fieldnames)

        data = [journal,]
        data.extend(result)
        csv_writer.writerow(data)

def delete_if_exists(file_path):
    '''Delete the file if it exists'''
    if os.path.exists(file_path):
        os.remove(file_path)

def read_config_file(config_file_path):
    '''Read the INI style configuration file into a dict'''
    config = {}
    parser = ConfigParser()
    parser.read(config_file_path)

    try:
        # read Entrez section
        config['email'] = parser.get('Entrez', 'email')
        config['api_key'] = parser.get('Entrez', 'api_key')
        config['use_api_key'] = parser.getboolean('Entrez', 'use_api_key')
    except Exception as e:
        print("An error occurred reading the config file.")
        print("  error: %s" % e)  
        sys.exit()

    return config

if __name__ == '__main__':
    main()