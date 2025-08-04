###############
### Imports ###

# General imports
import os
import argparse
import copy
import multiprocessing as mp
from pathlib import Path
from typing import List, Set
from collections import Counter
import pandas as pd
import numpy as np
from typing import List, Set, Optional
from pydantic import BaseModel

# Pheval phenopacket utils, and semsimian algorithm
from pheval.utils.phenopacket_utils import phenopacket_reader
from pheval.utils.phenopacket_utils import PhenopacketUtil
from semsimian import Semsimian




# For divying up data into batches for parallel processing
def divide_workload(data_list, num_proc: int=1) -> list:
    """
    Meant to divide up the elements in data_list into num_proc equal portions
    by iteratively adding each element to a basket never repeating the same basket until all baskets have an equal amount
    If num_proc == 1 then the original input list will be returned nested in a top layer list i.e. [data_list]
    """

    # Deal with our edge case at the very begginning which then is used as input into the second potential edge case
    ndata_elements = len(data_list)
    if ndata_elements < num_proc:
        num_proc = ndata_elements

    # Edge case
    if num_proc <= 1:
        return [data_list]
    else:
        baskets = [[] for i in range(0, num_proc)]
        index_count = 0
        for d in data_list:
            baskets[index_count].append(d)
            if index_count == (num_proc-1):
                index_count = 0
            else:
                index_count += 1

        #print("- Workload divided into {} portions with each portion recieving {} elements respectively...".format(num_proc, [format(len(b), ',') for b in baskets]))
        return baskets


####################################################
### Preprocessing and processing data structures ###
class Patient(BaseModel):
    sample_name: str
    phenopacket_path: str

    # Pulled from phenopacket
    phenotype_ids: list
    phenotype_count: int

    disease_name: str
    disease_id: str

    gene_symbol: str
    gene_id: str

    # For sssom mappings
    disease_mapped: Optional[str] = None
    gene_mapped: Optional[str] = None

        
def read_sssom_to_lookup(fpath, exact_match=True):
    """
    Assumes monarch initiative header style of # and then first line is header
    """

    sssom_map = {}
    with open(fpath, 'r') as infile:
        header = 0
        cols = {}
        for line in infile:


            line = line.strip('\r').strip('\n')

            # Weird thing happening with the monarch sssom files (new line character most likely as the final line)
            if len(line) == 0:
                break

            if line[0] == "#":
                continue

            header += 1
            cols = line.split('\t')
            if header == 1:
                col_inds = {v:i for i,v in enumerate(cols)}
                continue

            # Our actual data here
            map_key = cols[col_inds["object_id"]]
            map_val = cols[col_inds["subject_id"]]
            map_type = cols[col_inds["predicate_id"]]

            # Only deal with exact matches
            if map_type != "skos:exactMatch":
                continue

            sssom_map.update({map_key:map_val})

    val_count = len(set(list(sssom_map.values())))
    print("- {} unique terms mapping to {} unique terms...".format(format(len(sssom_map), ','), format(val_count, ',')))

    # Now add self referencing keys
    curr_keys = list(sssom_map.keys())
    for k in curr_keys:
        v = sssom_map[k]
        sssom_map.update({v:v})

    return sssom_map


def gather_input_output_info(input_path, output_path, results_suffix="_results.tsv"):
    """
    Input is allowed to be a directory containing .json phenopacket files
    Or input is allowed to be a filepath specifying a .json phenopacket
    """

    if os.path.isdir(input_path):
        process_data = [os.path.join(input_path, fname) for fname in os.listdir(input_path) if fname.endswith(".json")]

    elif os.path.isfile(input_path) and input_path.endswith(".json"):
        process_data = [input_path]

    else:
        return None, None

    # Pre-format output file names
    out_data = [os.path.join(output_path, pname.split('/')[-1].replace(".json", results_suffix)) for pname in process_data]

    # Convert to dictionaries
    process_data = {f.split("/")[-1].replace(".json", ""):f for f in process_data}
    out_data = {f.split("/")[-1].replace(results_suffix, ""):f for f in out_data}

    return process_data, out_data


def phenopacket_paths_to_data(sample_path_dict):

    # Return data structure and function variables
    pobjs = {}
    multi_gene, multi_dis = {}, {}
    processed, total_samples = 0, len(sample_path_dict)

    # Loop through all phenopackets, extract / map data and return Patient objects
    for sname, spath in sample_path_dict.items():
        
        phenopacket_util = PhenopacketUtil(phenopacket_reader(Path(spath))) # Requires Path object
        #phenopacket_util = PhenopacketUtil(phenopacket_reader(spath))
        observed_phenotypes = phenopacket_util.observed_phenotypic_features()
        phenotype_ids = [observed_phenotype.type.id for observed_phenotype in observed_phenotypes]
        phenotype_count = len(phenotype_ids)

        # Default is to take first term 
        # Will display how many multi terms there are at end... should be few if not none)
        dis_obj = phenopacket_util.diseases()[0]
        dis_name, dis_id = dis_obj.term.label, dis_obj.term.id

        gene_obj = phenopacket_util.diagnosed_genes()[0]
        gene_symbol, gene_id = gene_obj.gene_symbol, gene_obj.gene_identifier

        if len(phenopacket_util.diagnosed_genes()) > 1:
            multi_gene.update({sname:''})

        if len(phenopacket_util.diseases()) > 1:
            multi_dis.update({sname:''})

        #print(sname)
        #print(dis_name, dis_id)
        #print(gene_symbol, gene_id)

        pobjs.update({sname:Patient.model_validate({"sample_name":sname,
                                                    "phenopacket_path":spath,
                                                    "phenotype_ids":phenotype_ids,
                                                    "phenotype_count":phenotype_count,
                                                    "disease_name":dis_name,
                                                    "disease_id":dis_id,
                                                    "gene_symbol":gene_symbol,
                                                    "gene_id":gene_id})})

        processed += 1
        if processed % 1_000 == 0:
            print("- {}/{} phenopackets read into memory...".format(format(processed, ','),
                                                                    format(total_samples, ',')))


    print("- Multi gene diagnosis phenopackets found {}...".format(len(multi_gene)))
    print("- Multi disease diagnosis phenopackets found {}...".format(len(multi_dis)))
    print("- {} Phenopackets information read into memory...".format(format(len(pobjs), ',')))
    return pobjs


def filter_zero_data(phen_data, sub_sample=False):
    
    # Remove zero phenotype count samples
    removed = 0
    for k in list(phen_data.keys()):
        if phen_data[k].phenotype_count == 0:
            del phen_data[k]
            removed += 1       
    print("- {} samples removed with zero phenotypes...".format(removed))

    # Default is no subsampleing
    if sub_sample != False:

        # Subsample (to make full pipeline connection easier)
        sub_samp = 0
        removed = 0
        for k in list(phen_data.keys()):
            if sub_samp >= sub_sample:
                del phen_data[k]
                removed += 1

            sub_samp += 1
        print("- {} samples removed via subsampling...".format(removed))
        print("- {} samples remaining...".format(len(phen_data)))
    
    # Means we are left with zero samples to process so we should error out
    if len(phen_data) == 0:
        raise ValueError("- Error, No remaining samples after filtering for zero phenotype counts...")
    
    return phen_data


def parse_input_arguments_to_process_data(in_arg, out_arg, hp_map, mode="disease"):
    
    # Hard coded (for now, but easy enough to add)
    # Note, only used for bulk processing... otherwise user specified out_arg is used as filepath/name
    results_suffix="_{}_results.tsv".format(mode)

    #####################################
    ### Multiple phenopacket scenario ###
    if os.path.isdir(in_arg):

        # Make sure our output directory exists
        out_dir = out_arg
        if not os.path.isdir(out_dir):
            print("- Creating output directory {} as it was not found to exist yet...".format(out_dir))
            os.makedirs(out_dir, exist_ok=True)
        
        # Gather filepaths, and produce pre-formated output filenames
        process_data = [os.path.join(in_arg, fname) for fname in os.listdir(in_arg) if fname.endswith(".json")]
        out_data = [os.path.join(out_dir, pname.split('/')[-1].replace(".json", results_suffix)) for pname in process_data]
        
        # Make sure we actually have sampels to process
        if len(process_data) == 0:
            raise FileNotFoundError("- Error, No phenopacket (.json) files found in {}... Exiting".format(in_arg))
        
        # Convert to dictionaries
        in_paths = {f.split("/")[-1].replace(".json", ""):f for f in process_data}
        out_paths = {f.split("/")[-1].replace(results_suffix, ""):f for f in out_data}
    
        # Now read all phenopacket(s) into memory and remove phenopackets with zero phenotypes
        phen_data = phenopacket_paths_to_data(in_paths)
        phen_data = filter_zero_data(phen_data, sub_sample=False)
        

    ######################################################
    ### Single file input scenario (.json phenopacket) ###
    elif in_arg.endswith(".json"):
        
        # Make sure we have a valid input file, and that an output file is specified (cannot be a directory for single sample)
        if not os.path.isfile(args.in_arg):
            raise FileNotFoundError("- Input file {} does not exist... Exiting".format(args.in_arg))

        out_dir = os.path.dirname(out_arg)
        if not os.path.isdir(out_dir):
            print("- Creating output directory {} as it was not found to exist yet...".format(out_dir))
            os.makedirs(out_dir, exist_ok=True)
        
        if out_dir == args.out_arg.rstrip('/'):
            raise ValueError("- Error, output file name must be specified in out_arg {}".format(args.out_arg))
        else:
            outfile = out_arg

        sname = args.in_arg.split("/")[-1].replace(".json", "")
        in_paths = {sname:args.in_arg}
        out_paths = {sname:outfile}
        
        phen_data = phenopacket_paths_to_data(in_paths)
        phen_data = filter_zero_data(phen_data, sub_sample=False)
    

    #####################################################
    ### Single sample with HP:XYZ,HP:ZYX,... as input ###
    elif args.in_arg.startswith("HP:"):

        # Split out phenotype terms
        phenotype_terms = args.in_arg.split(",")
        
        # Perfrom some initial sanity checks / warnings
        if len(phenotype_terms) != len(set(phenotype_terms)):
            print("- Warning, repeat input phenotype terms found...")
        
        if len(phenotype_terms) != len([v for v in phenotype_terms if v.startswith("HP:")]):
            raise ValueError("- Error, all input phenotype terms must begin with HP:   Only hpo encoded phenotypes are considered")
        
        # Make sure all phenotype terms are in the hp_map
        invalid_terms = [v for v in phenotype_terms if v not in hp_map]
        if len(invalid_terms) > 0:
            raise ValueError("- Error, input phenotype terms {} not found in HP mappings... Exiting".format(invalid_terms))

        # Make sure an output file is specified (cannot be a directory for single sample)
        out_dir = os.path.dirname(out_arg)
        if not os.path.isdir(out_dir):
            print("- Creating output directory {} as it was not found to exist yet...".format(out_dir))
            os.makedirs(out_dir, exist_ok=True)
        
        if out_dir == args.out_arg.rstrip('/'):
            raise ValueError("- Error, output file name must be specified in out_arg {}".format(args.out_arg))
        else:
            outfile = out_arg
        
        # Now mimic the output file structure
        out_paths = {"single_sample":outfile}
        phen_data = {"single_sample":Patient.model_validate({"sample_name":"single_sample",
                                                             "phenopacket_path":'',
                                                             "phenotype_ids":phenotype_terms,
                                                             "phenotype_count":len(phenotype_terms),
                                                             "disease_name":'',
                                                             "disease_id":'',
                                                             "gene_symbol":'',
                                                             "gene_id":''})}
    

    ######################
    ### Last condition ###
    else:
        raise ValueError("- Error, Invalid input arguments... Note, HP terms must start with HP: if using terms as input instead of phenopacket(s). Exiting")
    
    # Note, handling of no phenotype containing samples found is handled in the filter_zero_data function
    return phen_data, out_paths


##########################
### Semsimain wrappers ###
def get_phenotype_associations(semsim, phenotype_ids, outfile, symbol_map, name_map, mode="disease"):
    """
    This algorithm leverages Semsimian + Monarchs phenio ontology to find the disease(s) 
    that are most associated with a patients phenotypes. A single .json phenopacket file can be passed in or
    a directory containing multiple .json phenopacket files. The patients observed phenotype terms are pulled
    from the data and are used as input to semsimian. Disease information is returned with the top associated
    ids appearing first in the list. 

    Subject prefixes within the phenio db begginning with "MONDO:" are compared
    """

    # Define our db prefix term
    if mode == "disease":
        subject_prefix = "MONDO:"
    elif mode == "gene":
        subject_prefix = "HGNC:"

    # Ensure phenotype_ids are of set type
    if type(phenotype_ids) != type(set()):
        phenotype_ids = set(phenotype_ids)

    # Perform search (results are sorted in order of best ranking to worst ranking)
    results =  semsim.associations_search(object_closure_predicate_terms={"biolink:has_phenotype"},
                                            object_terms=phenotype_ids, # Must be set
                                            include_similarity_object=False,
                                            subject_terms=None,
                                            subject_prefixes=[subject_prefix],

                                            search_type="full",
                                            #score_metric="ancestor_information_content",
                                            score_metric="phenodigm",
                                            limit=10000,
                                            direction="object_to_subject")

    ###results = [[0,0,"A"], [1,1,"B"], [2,2,"C"]] Testing purposes

    # Results are originally in form of [[score, details, mondo_id], ...]
    results = np.asarray(results).T

    # Convert our non-human readable names to human readable names
    symbols, names = [], []
    for v in results[2]:
        symbols.append(symbol_map[v])
        names.append(name_map[v])

    # Convert to df
    results_df = pd.DataFrame({"{}_id".format(mode):results[2],
                               "{}_symbol".format(mode):symbols,
                               "{}_name".format(mode):names,
                               "score":np.round(results[0].astype(float), decimals=4)})

    # Must set back to float otherwise results are not sorted properly
    results_df['score'] = results_df['score'].astype(float)

    # Now "rank" our results based on the score (pre sorted by best first)
    results_df['rank'] = results_df['score'].rank(method='dense', ascending=False)

    # Write out our data
    results_df.to_csv(outfile, sep='\t', header=True, index=False)

    return results_df


# Allows us to bulk process multiple samples without having to instantiate a new semsimian object for each sample
# This is a sub function of the main function, designed to be called in parallel (or single core)	
def proccess_samples(phenio_path, mode, mode_symbols, mode_names, patients_to_process):

    # Load necessary data into memory for semsimian processing
    semsim = Semsimian(predicates=["rdfs:subClassOf"], 
                       spo=None, 
                       resource_path=phenio_path)
    print("- Semsimian object loaded...")

    tt = 0
    for p in patients_to_process:
        outpath = p[1]
        phenotype_ids = p[0].phenotype_ids

        # Perform search (results are sorted in order of best ranking to worst ranking)
        results = get_phenotype_associations(semsim, 
                                             phenotype_ids, 
                                             outpath,
                                             mode_symbols,
                                             mode_names, 
                                             mode=mode)

        tt += 1
        # Progress statement
        if tt % 100 == 0:
            print("- Processed {}/{}".format(format(tt, ','), format(len(patients_to_process), ',')))

    return None



if __name__ == '__main__':
    ################
	## ARG PARSE ###
    def parse_input_command():
        parser = argparse.ArgumentParser(description='Runs semphen algorithm on phenopacket(s) or HP term(s) as input')
        parser.add_argument("-i", "--in_arg", help="Path to directory containing phenopackets, or path to single phenopacket JSON file, or comma separated list of HP terms", required=True, type=str, default=None)
        parser.add_argument("-o", "--out_arg", help="Path to output directory in case of multiple samples. Or path to output filename for single sample. Output directory will be created if doesn't exist already.", required=True, type=str, default=None)
        parser.add_argument("-d", "--data_dir", help="Directory containing necessary data files for semphen to run", required=True, type=str)
        parser.add_argument("-m", "--mode", help="Prioritization mode... disease or gene are allowed", required=True, choices=["disease", "gene"], type=str)
        parser.add_argument("-c", "--num_proc", help="Number of cores to use for parallel processing", required=False, type=int, default=1)
        return parser.parse_args()
    
    args = parse_input_command()
    ############################
    
    ###########################
    ### Initial file checks ###

    # Make sure all necessary download data is present before doing anything
    req_files = ["phenio.db", "monarch-kg_nodes.tsv", "gene_mappings.sssom.tsv"]
    if not os.path.isdir(args.data_dir):
        raise FileNotFoundError("- Error, Input directory to necessary data not found {}...".format(args.data_dir))

    missing_files = [f for f in req_files if not os.path.isfile(os.path.join(args.data_dir, f))]
    if len(missing_files) > 0:
        raise FileNotFoundError("- Error, Necessary files {} not found in {}...".format(missing_files, args.data_dir))
    
    ################################################
    ### Monarch kg hgnc, and mondo term mappings ###

    # Read kg nodes into memory, and pull out hgnc, and mondo relevant data, and HP term mappings
    kg_nodes = pd.read_csv(os.path.join(args.data_dir, "monarch-kg_nodes.tsv"), sep="\t", header=0, low_memory=False)
    hgnc_symbols = {k:v for k,v in zip(kg_nodes["id"], kg_nodes["name"]) if k.startswith("HGNC:")}
    hgnc_names = {k:v for k,v in zip(kg_nodes["id"], kg_nodes["full_name"]) if k.startswith("HGNC:")}
    print("- Gene mappings loaded | {} HGNC ids found...".format(format(len(hgnc_symbols), ',')))

    mondo_names = {k:v for k,v in zip(kg_nodes["id"], kg_nodes["name"]) if k.startswith("MONDO:")}
    mondo_descriptions = {k:v for k,v in zip(kg_nodes["id"], kg_nodes["description"]) if k.startswith("MONDO:")}
    print("- Disease mappings loaded | {} MONDO ids found...".format(format(len(mondo_names), ',')))

    hp_map = {k:v for k,v in zip(kg_nodes["id"], kg_nodes["name"]) if k.startswith("HP:")}
    print("- HP term mappings loaded | {} HP ids found...".format(format(len(hp_map), ',')))


    ########################################
    ### Input and output data formatting ###

    # Set our input namespace mappings for semsimian function / results
    # Note, gene and disease are only valid options imposed by higher level argparse choices option
    if args.mode == "gene":
        entity_symbols = hgnc_symbols
        entity_names = hgnc_names

    elif args.mode == "disease":
        entity_symbols = mondo_names
        entity_names = mondo_descriptions

    # Deal with input / output arguments for sample data and pull phenotypes from phenopackets (if necessary)
    phen_data, outpaths = parse_input_arguments_to_process_data(args.in_arg,
                                                                args.out_arg,
                                                                hp_map,
                                                                args.mode)

    # Copy our patient information / output paths for parallel processing (or single core processing)
    phen_base_data = [[copy.copy(v),copy.copy(outpaths[k])] for k,v in phen_data.items()] ###[0:10] # For testing
    print("- Input output formatting complete...")

    ##################################
    ### Parallel processing option ###
    if (args.num_proc > 1) and (len(phen_base_data) > 1):

        # Divy up necessary input data for semsimian into parallel processing chunks
        div_entity_symbols = [copy.copy(entity_symbols) for i in range(0, args.num_proc)]
        div_entity_names = [copy.copy(entity_names) for i in range(0, args.num_proc)]
        div_phenio_paths = [os.path.join(args.data_dir, "phenio.db") for i in range(0, args.num_proc)]
        div_modes = [args.mode for i in range(0, args.num_proc)]
        div_phen_data = divide_workload(phen_base_data, num_proc=args.num_proc)
        print("- Parallel processing with {} cores...".format(args.num_proc))

        # Setup parallel processing overhead, kick off jobs via asynchronous processing, and retrieve results
        output = mp.Queue()
        pool = mp.Pool(processes=args.num_proc)
        results = [pool.apply_async(proccess_samples, args=(ph, m, msy, mn, pdata)) for ph, m, msy, mn, pdata in zip(div_phenio_paths, 
                                                                                                                     div_modes, 
                                                                                                                     div_entity_symbols, 
                                                                                                                     div_entity_names, 
                                                                                                                     div_phen_data)]
        output = [p.get() for p in results]
        pool.close()
        pool.join()

    ##############################
    ### Single core processing ###
    else:
        print("- Single core processing...")
        proccess_samples(os.path.join(args.data_dir, "phenio.db"), 
                         args.mode, 
                         entity_symbols, 
                         entity_names, 
                         phen_base_data)

    print("- Done processing {} samples".format(format(len(phen_base_data), ',')))