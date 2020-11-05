'''permute_terms.py
Takes a csv with paired standardized (re-based) person name(s) and test names 
and creates a series of variations (permutations) that attempt to include all 
common equivalents for how eponyms might be written in the literature.

input file : terms_re-base.csv
output file : terms_permuted.csv

'''

__title__ = 'permute_terms.py'
__version__ = '1.0.0'
__author__ = 'Toby C. Cornish, MD, PhD'
__license__ = 'GPLv3'
__copyright__ = 'Copyright 2020'

import os
import sys
import csv
import re
import itertools

input_file_path = os.path.abspath(r'terms_re-base.csv')
output_file_path = os.path.abspath(r'terms_permuted.csv') 

def main():
    '''Main function'''
    # if the output files exists, delete them
    delete_if_exists(output_file_path)

    terms = read_file(input_file_path)

    terms_out = []

    print('='*40)
    for term in terms:
        name = term['Name(s)']
        test = term['Test']
        print(f"{name}, {test}")
        terms_out.extend(permute(term,terms_out))
        print('-'*40)

    for term in terms_out:
        print(term)
        pass
    write_file(output_file_path,terms_out)

def permute(term,terms_out):
    '''Create extensive variations on our base terms'''
    test = term['Test']
    name = term['Name(s)']

    terms_out = []

    # split name into list
    if name != '':
        name_list = split_names(term['Name(s)'])
    else:
        name_list = []
    print(name_list)

    if len(name_list) == 0:
        # no permutation
        terms_out.append({ 'base_term' : f'{test}', 'term' : f'{test}' })

    if len(name_list) == 1:
        # no permutation
        terms_out.append({ 'base_term' : f'{name} {test}', 'term' : f'{name} {test}' })

        # possessives
        if name[-1] == 's':
            terms_out.append({ 'base_term' : f'{name} {test}', 'term' : f"{name}' {test}" })

        terms_out.append({ 'base_term' : f'{name} {test}', 'term' : f"{name}'s {test}" })
        
        # inversion
        terms_out.append({ 'base_term' : f'{name} {test}', 'term' : f"{test} of {name}" })  

    if len(name_list) == 2:
        base_term = f"{'-'.join(name_list)} {test}"

        # no permutation
        terms_out.append({ 'base_term' : base_term, 
                           'term' : f"{'-'.join(name_list)} {test}" })

        # join variants        
        terms_out.append({ 'base_term' : base_term, 
                           'term' : f"{' '.join(name_list)} {test}" })

        terms_out.append({ 'base_term' : base_term, 
                           'term' : f"{' and '.join(name_list)} {test}" })

        # all double possessives
        possessive_names_list = list(map(get_possessives, name_list))
        for pn_pairs in itertools.product(possessive_names_list[0],
                                          possessive_names_list[1]):
            terms_out.append({ 'base_term' : base_term, 
                               'term' : f"{' and '.join(pn_pairs)} {test}" })

        # all final possessives
        possessive_names_list = get_possessives(name_list[-1])
        for pn in possessive_names_list:
            terms_out.append({ 'base_term' : base_term, 
                               'term' : f"{name_list[0]} and {pn} {test}" })
            terms_out.append({ 'base_term' : base_term, 
                               'term' : f"{name_list[0]} {pn} {test}" })
            terms_out.append({ 'base_term' : base_term, 
                               'term' : f"{name_list[0]}-{pn} {test}" })
                                                              
        # inversions
        terms_out.append({ 'base_term' : f'{name} {test}', 
                           'term' : f"{test} of {' and '.join(name_list)}" })
        terms_out.append({ 'base_term' : base_term, 
                           'term' : f"{test} of {' '.join(name_list)}" })
        terms_out.append({ 'base_term' : base_term, 
                           'term' : f"{test} of {'-'.join(name_list)}" })                           

    if len(name_list) == 3:
        base_term = f"{'-'.join(name_list)} {test}"

        # no permutation
        terms_out.append({ 'base_term' : base_term, 
                           'term' : f"{'-'.join(name_list)} {test}" })

        # join variants        
        terms_out.append({ 'base_term' : base_term, 
                           'term' : f"{' '.join(name_list)} {test}" })

        terms_out.append({ 'base_term' : base_term, 
                            'term' : f"{name_list[0]}, {name_list[1]} and {name_list[2]} {test}" })         
        terms_out.append({ 'base_term' : base_term, 
                            'term' : f"{name_list[0]}, {name_list[1]}, and {name_list[2]} {test}" })

        # all triple possessives
        possessive_names_list = list(map(get_possessives, name_list))
        for pn_triples in itertools.product(possessive_names_list[0],
                                            possessive_names_list[1], 
                                            possessive_names_list[2]):
            terms_out.append({ 'base_term' : base_term, 
                               'term' : f"{pn_triples[0]}, {pn_triples[1]} and {pn_triples[2]} {test}" })
            terms_out.append({ 'base_term' : base_term, 
                               'term' : f"{pn_triples[0]}, {pn_triples[1]}, and {pn_triples[2]} {test}" })

        # all final possessives
        possessive_names_list = get_possessives(name_list[-1])
        for pn in possessive_names_list:
            terms_out.append({ 'base_term' : base_term, 
                               'term' : f"{name_list[0]}, {name_list[1]} and {pn} {test}" })         
            terms_out.append({ 'base_term' : base_term, 
                               'term' : f"{name_list[0]}, {name_list[1]}, and {pn} {test}" }) 
            terms_out.append({ 'base_term' : base_term, 
                               'term' : f"{name_list[0]} {name_list[1]} {pn} {test}" })                    
            terms_out.append({ 'base_term' : base_term, 
                               'term' : f"{name_list[0]}-{name_list[1]}-{pn} {test}" })                    

        # inversions
        terms_out.append({ 'base_term' : f'{name} {test}', 
                           'term' : f"{test} of {' and '.join(name_list)}" }) 
        terms_out.append({ 'base_term' : f'{name} {test}', 
                           'term' : f"{test} of {' '.join(name_list)}" })
        terms_out.append({ 'base_term' : f'{name} {test}', 
                           'term' : f"{test} of {'-'.join(name_list)}" })
        terms_out.append({ 'base_term' : f'{name} {test}', 
                           'term' : f"{test} of {name_list[0]}, {name_list[1]} and {name_list[2]}" }) 
        terms_out.append({ 'base_term' : f'{name} {test}', 
                           'term' : f"{test} of {name_list[0]}, {name_list[1]}, and {name_list[2]}" }) 

    return terms_out

def get_possessives(name): 
    '''Return a list of possible possessives for a name'''
    possessives = []
    if name[-1] == 's':
        possessives.append(f"{name}'")
    possessives.append(f"{name}'s")
    return possessives

def split_names(name):
    '''Split name string on - and , return list of names'''
    name = name.replace(', ',',')
    return re.split(r'[-,]',name)

def write_file(file_path, terms_out):
    '''Write the file of terms with permuted names'''
    fieldnames = ['base_term','term']

    with open(file_path,'w', newline='') as f:
        writer = csv.DictWriter(f,fieldnames=fieldnames,dialect='excel')
        writer.writeheader()
        for term in terms_out:
            writer.writerow(term)

def read_file(file_path):
    '''Write file of re-based (uniform) terms'''
    terms = []
    with open(file_path,'r') as f:
        terms = list(csv.DictReader(f))
    return terms

def delete_if_exists(file_path):
    '''Delete file if it exists'''
    if os.path.exists(file_path):
        os.remove(file_path)

if __name__ == '__main__':
    main()