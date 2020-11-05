'''pubmed_search_to_csv.py
Searchs pubmed for the EXACT terms in the input file in the title & abstract
Returns all the PMIDs identified categorized by base term and permuted term
Does not remove duplicates within a base term

input file : terms_permuted.csv
output files : term_results.csv, pmid_results.csv
'''

__title__ = 'pubmed_search_to_csv.py'
__version__ = '1.0.0'
__author__ = 'Toby C. Cornish, MD, PhD'
__license__ = 'GPLv3'
__copyright__ = 'Copyright 2020'

import os
import sys
import csv
import time
from configparser import ConfigParser
from collections import defaultdict

from Bio import Entrez
from Bio import Medline

input_file_path = os.path.abspath(r'terms_permuted.csv')
term_output_file_path = os.path.abspath(r'term_results.csv') 
pmid_output_file_path = os.path.abspath(r'pmid_results.csv')

# this version has been updated to accomodate the terms_permuted update

def main():
    '''Main function'''
    # read our config file
    config = read_config_file('config.ini')
    if not config['use_api_key']:
        config['api_key'] = None
    print(config)

    # if the output files exists, delete them
    delete_if_exists(term_output_file_path)
    delete_if_exists(pmid_output_file_path)

    search_terms = read_file(input_file_path)

    print(search_terms)

    print('='*40)
    for term in search_terms:
        search_string = term['term'] #use the permuted term
        base_term = term['base_term']
        print(f'{base_term} : {search_string}')
        results = search(term,'title/abstract',config['email'],config['api_key'])
        print(results)
        id_list = results['pmids']
        print(f'  Found {len(id_list)} PMIDs')
        papers = list(fetch_details(id_list,config['email'],config['api_key'])) #convert generator to list
        for i, paper in enumerate(papers):
            PMID = paper.get('PMID') # get returns None if the key doesnt exist
            TA = paper.get('TA')
            JT = paper.get('JT')
            print(f"      {i}: {PMID}, {TA}, {JT}")
        write_pmid_results(term, papers)
        write_term_results(results)
        print('='*40)
        time.sleep(1) # required if too many nulls returns in a row

def search(term, field, email, api_key=None):
    '''Use Entrez.esearch to identify PMIDs that match our phrase exactly'''

    Entrez.email = email
    if api_key:
        Entrez.api_key = api_key
    base_term = term['base_term']
    query = '"' + term['term'] + '"' # add quotes to phrase
    
    # see https://www.ncbi.nlm.nih.gov/books/NBK25499/#chapter4.ESearch for information on parameters
    # of note: To retrieve more than 100,000 UIDs, submit multiple esearch requests while incrementing the value of retstart
    print("  Retrieving PMIDs ...")
    print(f"    Search term : {query}")
    query_result = {
        'base_term' : base_term,
        'term' : query,
        'pmids' : [],
        'query_translation' : None,
        'quoted_phrase_found' : False,
        'count' : 0,
        'field' : field
    }
    pmids = []
    retstart = 0
    while True:    
        handle = Entrez.esearch(db='pubmed', 
                                sort='relevance', 
                                retmax='100000',
                                retstart=retstart,
                                retmode='xml',
                                term=query,
                                field=field)
        results = Entrez.read(handle)

        # populate part of the query_result
        query_result['query_translation'] = results['QueryTranslation']
        query_result['count'] = int(results['Count'])

        # if the quoted phrase isn't found we get whacko results
        if quoted_phrase_found(results):
            pmids.extend(results['IdList'])
            count = int(results['Count'])
            query_result['quoted_phrase_found'] = True
        else:
            query_result['quoted_phrase_found'] = False
            print("      The quoted phrase wasn't found (0 results).")
            break

        print(f'      Retrieved {len(pmids)} of {count} PMIDs.')
        # if we haven't retrieved all the PMIDs, increment the retstart and loop
        if count > len(pmids):
            retstart = len(pmids)
        else:
            break

    print("      Done.")
    query_result['pmids'] = pmids
    return query_result

def quoted_phrase_found(results):
    '''If the quoted phrase was found return True'''
    if 'WarningList' in results:
        if 'QuotedPhraseNotFound' in results['WarningList']:
            return False
        else:
            return True
    else:
        return True

def fetch_details(id_list, email, api_key=None):
    '''Retrieve the PMID details for the list of PMIDs'''
    num_ids = len(id_list) # a single call only returns up to 10,000
    Entrez.email = email
    if api_key:
        Entrez.api_key = api_key

    # see https://www.ncbi.nlm.nih.gov/books/NBK25499/#chapter4.EFetch for information on parameters
    # see https://www.nlm.nih.gov/bsd/mms/medlineelements.html for the medline elements
    # efetch can only return a max of 10000 results per request, so we'll work in chunks
    results = []
    start = 0
    end = 0
    print(f"  Retrieving details for {num_ids} PMIDs ...")
    # we'll use a the chunks generator to get batches of 10,000 PMIDs
    for batch in chunks(id_list,10000):
        start = end
        end = start + len(batch) 
        print(f"    Retrieving details for PMIDs {start} to {end} ...")
        pmids = ','.join(batch)
        handle = Entrez.efetch( db="pubmed", 
                                id=pmids, 
                                rettype="medline",
                                retmode="text")
        results.extend(list(Medline.parse(handle)))
    print("    Done.")
    return results

def chunks(l, n): 
    '''Generator for list l as chunks of size n'''
    for i in range(0, len(l), n):  
        yield l[i:i + n] 

def read_file(file_path):
    '''Read the CSV of terms and return a dict'''
    terms = []
    with open(file_path,'r') as f:
        terms = list(csv.DictReader(f))
    return terms

def write_term_results(results):
    '''Write CSV with one term per row'''
    file_exists = os.path.isfile(term_output_file_path) # if it doesn't exist, we are creating it and need to write a header
    fieldnames = ['base_term','term','field','count','quoted_count','quoted_phrase_found','query_translation']

    with open(term_output_file_path,'a', newline='') as f: # mode is 'append'
        csv_writer = csv.writer(f, dialect='excel')

        if results['quoted_phrase_found']:
            quoted_count = results['count']
        else:
            quoted_count = 0

        if not file_exists:
            csv_writer.writerow(fieldnames)

        csv_writer.writerow([   results['base_term'],
                                results['term'], 
                                results['field'],
                                results['count'],
                                quoted_count,
                                results['quoted_phrase_found'],
                                results['query_translation']
        ])

def write_pmid_results(term, papers):
    '''Write CSV with one PMID per row'''
    file_exists = os.path.isfile(pmid_output_file_path) # if it doesn't exist, we are creating it and need to write a header
    fieldnames = ['base_term','term','PMID','Year','DP','TA','JT']

    base_term = term['base_term']
    search_string = '"' + term['term'] + '"' # add quotes to phrase

    with open(pmid_output_file_path,'a', newline='') as f: # mode is 'append'
        csv_writer = csv.writer(f, dialect='excel')

        if not file_exists:
            csv_writer.writerow(fieldnames)

        for paper in papers:
            PMID = paper.get('PMID') # get returns None if the key doesn't exist
            TA = paper.get('TA')
            JT = paper.get('JT')
            DP = paper.get('DP')
            if DP:
                Year = DP[:4]
            else:
                Year = ''
            csv_writer.writerow([base_term, search_string, PMID, Year, DP, TA, JT])

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