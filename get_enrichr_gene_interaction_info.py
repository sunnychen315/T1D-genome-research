"""
INPUT:      A enrichr_gene_names.csv  file with a list of genes of interest in first column,
            the tissue of interest in the second column, and the gene set library to use in the
            third column.
RETURNS:    A list of interacting genes for each gene and specified tissue.

PARAMETERS:

"""

from math import nan
import pandas as pd
import json
import requests
import re
import csv 

default_gene_set_libraries = [
    'GO_Biological_Process_2015',
    "ChEA_2015",
    "KEGG_2016",
    "ESCAPE",
    "Epigenomics_Roadmap_HM_ChIP-seq",
    "ENCODE_TF_ChIP-seq_2015",
    "ENCODE_Histone_Modifications_2015",
    "OMIM_Expanded",
    "TF-LOF_Expression_from_GEO",
    "Single_Gene_Perturbations_from_GEO_down",
    "Single_Gene_Perturbations_from_GEO_up",
    "Disease_Perturbations_from_GEO_down",
    "Disease_Perturbations_from_GEO_up",
    "Drug_Perturbations_from_GEO_down",
    "Drug_Perturbations_from_GEO_up",
    "WikiPathways_2016",
    "Reactome_2016",
    "BioCarta_2016",
    "NCI-Nature_2016"
]

#####
# HELPER FUNS
#####
def get_gene_ids(genes):
    # ENRICHR URL
    ENRICHR_URL = 'https://maayanlab.cloud/Enrichr/addList'

    # INIT GENE IDE LIST
    gene_ids = []

    # NOW GET ID FOR EACH GENE
    for gene in genes:
        genes_str = '\n'.join(list(gene))
        description = 'Query gene list'
        payload = {
            'list': (None, genes_str),
            'description': (None, description)
        }

        # GET ID
        response = requests.post(ENRICHR_URL, files=payload)
        if not response.ok:
            raise Exception('Error analyzing gene list for gene: ', gene)

        data = json.loads(response.text)
        gene_ids.append(data['userListId'])

    return gene_ids

###
# FIND A TISSUE TYPE (I.E., SUBSTRING) WITHIN THE
# LARGER DESCRIPTION OF A RETURNED RESULT
def is_substring(sub_string, larger_string, chk_lower_upper_cap = True):
    # If chk_lower_upper_cap is False, then only the provided sub_string, without modification
    # is checked for being part of the larger_string
    # First make all sub_string characters lower case, uppercase, and capitalized
    # then check for being paer of the larger string.
    lc_substring = sub_string.lower()
    uc_substring = sub_string.upper()
    cap_substring = lc_substring.capitalize()

    # Now check if any of the lower, upper, or capitalized sub strings
    # are within the larger string.
    lc_res = re.search("("+lc_substring+")", larger_string)
    uc_res = re.search("("+uc_substring+")", larger_string)
    cap_res = re.search("("+cap_substring+")", larger_string)
    unchanged_res = re.search("("+sub_string+")", larger_string)

    if chk_lower_upper_cap:
        if unchanged_res:
            return unchanged_res[1]
        elif lc_res:
            return lc_res[1]
        elif uc_res:
            return uc_res[1]
        elif cap_res:
            return cap_res[1]
        else:
            return False
    else: # We checked only the original sub string
        if unchanged_res:
            return unchanged_res[1]
        else:
            return False

###
def chk_and_remove_nans(list_of_names):
   new_name_list = []
   for name in list_of_names:
       if name != nan:
           new_name_list.append(name)

   return new_name_list

#####
# START OF MAIN CODE
#####
# READ THE  enrichr_gene_tissue_lib_names.csv FILE TO GET THE GENES TO PROCESS.
# READ GENE NAME, TISSUE, AND LIBRARY INFO
gene_tissue_library_names = pd.read_csv("enrichr_gene_names.csv")

# EXTRACT GENE NAMES INTO A LIST
gene_names = gene_tissue_library_names.loc[:,'GENE']
if len(gene_names) > 1:
    gene_names = gene_names.squeeze().tolist()
    gene_names = chk_and_remove_nans(gene_names)
elif len(gene_names) == 1:
    gene_names = [gene_names.squeeze()]
else:
    import sys
    print("gen_enrich_gene_interaction_info: Problem with reading gene names.")
    sys.exit(0)

# EXTRACT TISSUE TYPES INTO A LIST
tissue_names = gene_tissue_library_names.loc[:,'TISSUE']
if len(tissue_names) > 1:
    tissue_names = tissue_names.squeeze().tolist()
    tissue_names = chk_and_remove_nans(tissue_names)
elif len(tissue_names) == 1:
    tissue_names = [tissue_names.squeeze()]
else:
    import sys
    print("gen_enrich_gene_interaction_info: Problem with reading tissue names.")
    sys.exit(0)

# EXTRACT LIBRARY TYPES INTO A LIST
library_names = gene_tissue_library_names.loc[:,'LIBRARY']
if len(library_names) > 1:
    library_names = library_names.squeeze().tolist()
    library_names = chk_and_remove_nans(library_names)
elif len(library_names) == 1:
    library_names = [library_names.squeeze()]
else:
    import sys
    print("gen_enrich_gene_interaction_info: Problem with reading library names.")
    sys.exit(0)

# GET GENE IDs
gene_ids = get_gene_ids(gene_names)

print("gene_names: ", gene_names)
print("gene_ids: ", gene_ids)
print("tissue_names: ", tissue_names)
print("library names: ", library_names)
print()

# SPECIFY URL OF ENRICHR DB TO USE AND CREATE QUERY STRING
ENRICHR_URL = 'https://maayanlab.cloud/Enrichr/enrich'
query_string = '?userListId=%s&backgroundType=%s'

# CREATE A RESULT PANDAS OBJECT
results_col_names= ["Gene Target", "Interacting Gene", "Tissue", "Enrichr Library"]
results = pd.DataFrame(columns=results_col_names)

# Loop through each gene of interest to find interacting genes
for gene_name, gene_id, tissue_name, lib_name in zip(gene_names, gene_ids, tissue_names, library_names):
    print("Retrieving interacting genes for "+gene_name+" with id: "+str(gene_id)+" for tissue: "+tissue_name+" and enrichment library: "+lib_name)
    response = requests.get(ENRICHR_URL + query_string % (gene_id, lib_name))
    # MAKE SURE WE HAVE VALID RESPONSE
    if not response.ok:
        raise Exception('Error fetching enrichment results for gene: '+gene_name)

    # EXTRACT RESULTS
    data = json.loads(response.text)

    # RETRIEVE INFO FOR THE DESIRED TISSUE AND FROM THE DESIRED LIBRARY
    if lib_name in data.keys():
        tissue_info = data[lib_name]    # Extract the entry for the specified library
        entry_found = False
        for entry in tissue_info:
            tissue_match_res = is_substring(tissue_name, entry[1])
            if tissue_match_res:
                interacting_gene = entry[1].split()[0]
                print("      The Gene Name that interacts with "+gene_name+" is "+interacting_gene+" in "+tissue_match_res+" tissue.")
                ["Gene Target", "Interacting Gene", "Tissue", "Enrichr Library"]
                new_result_row = {"Gene Target" : gene_name, "Interacting Gene" : interacting_gene,
                                  "Tissue" : tissue_match_res, "Enrichr Library" : lib_name}
                results = results.append(new_result_row, ignore_index = True)
                entry_found = True

        if not entry_found:
            print("      Interacting Gene(s) for "+gene_name+" in tissue "+tissue_name+" in library "+lib_name+" not found")
            entry_found = False

        print("==========\n")

# WRITE RESULTS TO .csv  FILE
results.to_csv("enrichr_interacting_gene_results.csv", index=False, encoding='utf-8')

print("#####################################")
print("Result Data Frame")
print(results)
print("#####################################")
print("Completed writing enrichr results to the file named 'enrichr_interacting_gene_results.csv'")
print()
print(">>>>>>>>>> DONE <<<<<<<<<<")

# MAKE SURE THIS MODULE IS EXECUTED AND NOT IMPORTED
if(__name__ != '__main__'):
    print("Module get_enrichr_gene_interaction_info.py  is intended to be executed and not imported.")