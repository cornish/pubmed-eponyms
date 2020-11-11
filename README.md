# pubmed-eponyms

This is a collection of Python scripts for searching pubmed using [BioPython](https://pypi.python.org/pypi/biopython) and working with [eponymous](https://en.wikipedia.org/wiki/Eponym) terms.

This project serves to deposit the code used in the manuscript indicated below. In addition to being of interest to those studying medical [eponyms](https://en.wikipedia.org/wiki/Eponym), it should also be of general use for anyone looking to develop software for automatically searching Pubmed / Medline. The script `pubmed_search_to_csv.py` provides a good example of how to use [BioPython's](https://pypi.python.org/pypi/biopython) [Entrez.esearch](https://biopython.org/DIST/docs/api/Bio.Entrez-module.html) and [Entrez.efetch](https://biopython.org/DIST/docs/api/Bio.Entrez-module.html) to search Pubmed and return search results even when they exceed the [NCBI's Entrez eutils](https://www.ncbi.nlm.nih.gov/books/NBK25500/) built-in limits. While the [eutils API](https://www.ncbi.nlm.nih.gov/books/NBK25500/) can be used directly, [BioPython](https://pypi.python.org/pypi/biopython) greatly simplifies the tedious aspects of making http requests (such as throttling, re-attempting and error handling) and is highly recommended for this task.

## Citing
In addition to citing this GitHub repository (https://github.com/cornish/pubmed-eponyms), please cite the following paper:

- A Biopython-based method for comprehensively searching for eponyms in Pubmed. Toby C. Cornish, Larry J. Kricka, and Jason Y. Park, In preparation


## License
Gnu Public License v3, see text of the full license in project.


## Dependencies:
- [Python 3.6](https://www.python.org/downloads/) and up
- [BioPython](https://pypi.python.org/pypi/biopython)


## Script files
1. `rebase_terms.py`
2. `permute_terms.py`
3. `pubmed_search_to_csv.py`
4. `remove_pmid_dupes.py`
5. `pubmed_journals_by_year.py`

### A diagrammatic representation of data flow indicating the scripts used in the process. Please see individual scripts for usage and details.

![](https://github.com/cornish/pubmed-eponyms/blob/main/data%20flow%20diagram/data_flow.png)

### An example of the term permutations created by permute_terms.py for terms with zero, one, two and three separate names.

<img src="https://github.com/cornish/pubmed-eponyms/blob/main/data%20flow%20diagram/permutations.png" height="100" />

![](https://github.com/cornish/pubmed-eponyms/blob/main/data%20flow%20diagram/permutations.png)

## Configuration file:
1. `config.ini`
   - This is an INI-style configuration file where the scripts will look for Entrez-related credentials including your email and API key
   - An API key for the e-utilities is not required at the time this was written, but may be in the future; currently it permits more requests per second to Entrez
   - See [here](https://ncbiinsights.ncbi.nlm.nih.gov/2017/11/02/new-api-keys-for-the-e-utilities/) for more information about API keys for NCBI's E-utilities

## Data files:
1. `gastrointestinal eponyms.txt`
   - This is the original list of terms collected from review articles:
     + [Kanne JP, Rohrmann CA, Lichtenstein JE. Eponyms in radiology of the digestive tract: historical perspectives and imaging appearances. Part I. Pharynx, esophagus, stomach, and intestine. Radiographics. 26(1) (2006) 129-42.](https://doi.org/10.1148/rg.261055084) 
     + [Kanne JP, Rohrmann CA, Lichtenstein JE. Eponyms in radiology of the digestive tract: historical perspectives and imaging appearances. Part 2. Liver, biliary system, pancreas, peritoneum, and systemic disease. Radiographics. 26(2) (2006) 465-80.](https://doi.org/10.1148/rg.262055130)

2. `gi_eponyms_split.csv`
   - This is the original list with terms split into Name(s) and Term fields; multiple name eponyms should be separated by by hyphens to distinguish them from last names with internal spaces (i.e. "Van Slyke")
   - Input to `rebase_terms.py`
3. `data/terms_re-base.csv`
   - Output of the `rebase_terms.py` script
   - Input to `permute_terms.py`
   - Standardized version of base names including removal of possessives and use of hyphens for multiple names
   - Version of the data from the paper
4. `data/terms_permuted.csv`
   - Output of the `permute_terms.py` script
   - Input to `pubmed_search_to_csv.py`
   - Permutations of terms to include possesives, various forms of joining multiple names, and inversions
   - See examples above
   - Version of the data from the paper
5. `data/term_results.csv`
   - Output of the `pubmed_search_to_csv.py` script
   - Summarizes pubmed search results for all terms (including terms with no results)
   - One row per term permutation
   - Version of the data from the paper
5. `data/pmid_results.csv`
   - Output of the `pubmed_search_to_csv.py` script
   - Input to `remove_pmid_dupes.py` 
   - Input to `pubmed_journals_by_year.py`
   - Pubmed search results for all terms with hits
   - One row per PMID 
   - Version of the data from the paper
6. `data/pmid_results - dupes removed.csv`
   - Output of the `remove_pmid_dupes.py` script
   - Duplicate PMIDs removed within base terms
   - Version of the data from the paper
7. `data/journal_counts.csv`
   - Output of the `pubmed_journals_by_year.py` script
   - A matrix of total publications per year for all journals for which we have hits 
   - Ranges from the earliest year with hits to the latest year with hits
   - Version of the data from the paper
