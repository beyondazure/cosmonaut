# Written by beyondazure (Daniel Kabanovsky) at Boston College

import argparse
import bisect
import csv
import download
import math
import plotly.express as px
import re
import sys
import os
from tqdm import tqdm

# the information for each mutation is stored in an instance of the Mutation class
class Mutation:

    def __init__(self, name, mut_type, mutation, histology, identifier):
        self.name = name
        self.mut_type = mut_type
        self.mutation = mutation
        self.identifier = identifier
        self.sample_count = 0
        self.total_count = 0
        self.histology = histology
        self.gsf = 0
        self.rsf = 0
        self.tgf = 0
        self.tsf = 0
        self.escore_gsf = 0
        self.escore_rsf = 0
        self.escore_tgf = 0
        self.escore_tsf = 0
        self.rank = 0

# for efficiency, TSV files are treated as plain text and regex is used to extract
# useful data from each line
    
def re_capture(line, amino_acid, mutation_type):
    
    # defining capture patterns
    name_pattern =  re.compile(r'^(.*?)\t')
    mut_pattern = re.compile(r'p\.(.*?)\t')
    phen_pattern = re.compile(r'COSO(.*?)\t')
    id_pattern = re.compile(r'COSG(.*?)\t')
    
    # matching patterns in text
    name_match = name_pattern.search(line)
    mut_match = mut_pattern.search(line)
    phen_match = phen_pattern.search(line)
    id_match = id_pattern.search(line)
    
    # extracting data into variables
    name = name_match.group(1) if name_match != None else "?"
    mut = mut_match.group(1) if mut_match != None else "?"
    phen = phen_match.group(1) if phen_match != None else "?"
    identifier = id_match.group(1) if id_match != None else "?"

    # searching for the gene name of the current mutation in gene_mutation_names
    target_index = bisect.bisect_left(gene_mutation_names, name)
    # gene name is encountered for the first time
    if target_index >= len(gene_mutation_names) or gene_mutation_names[target_index] != name:
        bisect.insort(gene_mutation_names, name) # add name into correct location to maintain sorting 
        gene_mutation_list[name] = []
        gene_total_dict[name] = 1
    # gene name has been encountered already
    else:
        gene_total_dict[name] += 1
    # creating a histology search string and resolving it to tissue type
    search_string = "COSO" + phen
    hist = hist_dict[search_string] if search_string in hist_dict else "unknown"
    # saving the number of mutations associated with a given gene AND histology
    if hist in hist_total_dict:
        if name in hist_total_dict[hist][1]:
            hist_total_dict[hist][1][name] += 1
            hist_total_dict[hist][0] += 1 
            # first variable stores total number of mutations for given tissue type
        else:
            hist_total_dict[hist][1][name] = 1
            hist_total_dict[hist][0] += 1
    else:
        print("Unknown histology type. ")
        sys.exit()

    # using MOAR regex to determine if a mutation is a loss/gain of target amino acid
    if mut[0] == amino_acid or mut[-1] == amino_acid:
        
        global all_mutation_types
        all_mutation_types += 1 # mutation matches amino acid but not necessarily mutation type
        # establishing mutation type using regex 
        if re.search("\w*del\w*", mut) != None or re.search("\w*ins\w*", mut) != None:
            mut_type = "indel"
        elif re.search("\w*fs\w*", mut) != None:
            mut_type = "frameshift"
        elif re.search("\w*dup\w*", mut) != None:
            mut_type = "duplication"
        elif re.search("\w*[\?\=\*]\w*", mut) != None:
            mut_type = "non-specified"
        else:
            mut_type = "point"
            if mut[-1] == amino_acid and choice != "gain":
                return None # we're looking for a loss of target amino acid but found a gain
        if mut_type != mutation_type:
            return None # not the target mutation type
        # instantiating a new Mutation object
        new_mutation = Mutation(name = name, mut_type = mut_type, mutation = mut, histology = hist, identifier = identifier)

        if not args.disable_collapse:
        # collapsing two amino acids into one if they have the same tissue type and residue
            for other_mutation in gene_mutation_list[name]:
                if other_mutation.mutation[:-1] == mut[:-1] and other_mutation.histology == hist:
                    other_mutation.sample_count += 1
                    global amino_acid_collapse
                    amino_acid_collapse += 1
                    return None # "merged" new mutation into the old one, nothing more to do
        # if the mutation is different from all previous ones, save it into gene_mutation_list
        # since it's encountered for the first time, its sample count is 1
        gene_mutation_list[name].append(new_mutation)
        new_mutation.sample_count = 1
        
def sample_parser(line, dict1, dict2):
    
    # defining capture patterns
    id_pattern =  re.compile(r'COSO(.*?)\t')
    hist_pattern = re.compile(r'^[^\t]+\t([^\t]+)\t.*') 
    # NOTE: this searches histology based on position within a line; 
    # could break if COSMIC changes its file structure

    # matching patterns in text
    id_match = id_pattern.search(line)
    hist_match = hist_pattern.search(line)
    
    # extracting data into variables and populating dicts
    ids = "COSO" + id_match.group(1) if id_match != None else "?"
    hist = hist_match.group(1) if id_match != None else "?"
    dict1[ids] = hist
    dict2[hist] = [0, {}]

# pretty banner
print("""
_________                                                   __   
\_   ___ \  ____  ______ _____   ____   ____ _____   __ ___/  |_ 
/    \  \/ /  _ \/  ___//     \ /  _ \ /    \\\\__  \ |  |  \   __\\
\     \___(  <_> )___ \|  Y Y  (  <_> )   |  \/ __ \|  |  /|  |  
 \______  /\____/____  >__|_|  /\____/|___|  (____  /____/ |__|  
        \/           \/      \/            \/     \/             
      """)
print("* * * Command-line parser for bulk mutation data from COSMIC * * *\n")

# parsing command-line arguments
my_args = argparse.ArgumentParser(description="Command-line parser for bulk mutation data from COSMIC")
my_args.add_argument('input_files', nargs='*', type=str, default=None, help="Input files")
my_args.add_argument('-d', '--disable_collapse', action='store_true', help="Disable mutation collapse feature")
my_args.add_argument('-s', '--script', action='store_true', help="Download files from COSMIC automatically")
my_args.add_argument('-f', '--classify', type=str, default=None, help="COSMIC histology classification file (REQUIRED)")
my_args.add_argument('-t', '--tissue', type=str, default=None, help="Search by tissue type")
my_args.add_argument('-o', '--output', type=str, default=None, help="Save output to file")
my_args.add_argument('-c', '--cutoff', type=int, default=None, help="Filter out mutations with sample count below (strictly less than) the cutoff")
args = my_args.parse_args()

# handling invalid input
if len(args.input_files) != 2 and not args.script:
    print("Usage: python cosmonaut.py [/path/to/genomic_screen_data.tsv] [/path/to/targeted_screen.tsv] -f [/path/to/classification_file.tsv] -o [output_file] -t [tissue type] -c [cutoff value] [-s -d]")
    sys.exit()

if args.classify == None and not args.script:
    print("Error: COSMIC histology classification file is required.")
    sys.exit()

if args.input_files == None and not args.script:
    print("Error: either provide COSMIC input files or use scripted download (-s flag).")
    sys.exit()

# saving command-line input
if args.script:
    file_obj = download.scripted_download() # launching scripted download
    filepath_genomic = file_obj[0]
    filepath_targeted = file_obj[1]
    classification_file = file_obj[2]
else:
    filepath_genomic = args.input_files[0]
    filepath_targeted = args.input_files[1]
    classification_file = args.classify

# defining possible amino acid/mutation choices
amino_acids = ["alanine", "arginine", "asparagine", "aspartate", "cysteine", "glycine", "glutamine", "glutamate", "histidine", "isoleucine", "leucine", "lysine", "methionine", "phenylalanine", "proline", "serine", "threonine", "tryptophan", "tyrosine", "valine", "non-specific"]
three_letter_codes = ["ala", "arg", "asn", "asp", "cys", "gly", "gln", "glu", "his", "ile", "leu", "lys", "met", "phe", "pro", "ser", "thr", "trp", "tyr", "val", "ns"]
one_letter_codes = ["A", "R", "N", "D", "C", "G", "Q", "E", "H", "I", "L", "K", "M", "F", "P", "S", "T", "W", "Y", "V", "NS"]
ns_mutations = ["indel", "frameshift", "duplication", "non-specified", "point"]

# acquiring/validating user input
amino_acid = input("Enter the amino acid to search for (letter/code/name): ")
if amino_acid.lower() not in amino_acids and amino_acid.lower() not in three_letter_codes and amino_acid.upper() not in one_letter_codes:
    while amino_acid.lower() not in amino_acids and amino_acid.lower() not in three_letter_codes and amino_acid.upper() not in one_letter_codes:
        amino_acid = input("Invalid amino acid designation. Please enter a valid amino acid (letter/code/name): ")

mutation_type = input("Enter the mutation type to search for (indel/frameshift/duplication/point/non-specified): ")
if mutation_type.lower() not in ns_mutations:
    while mutation_type.lower() not in ns_mutations:
        mutation_type = input("Invalid mutation type specified. Please provide a valid mutation type (indel/frameshift/duplication/point/non-specified): ")
mutation_type = mutation_type.lower()

if mutation_type == "point":
    choice = input("Would you like to search for gain or loss mutations? (gain/loss): ")
    if choice.lower() not in ["gain", "loss"]:
        while choice.lower() not in ["gain", "loss"]:
            choice = input("Invalid option: please enter \"gain\" or \"loss\": ")
    choice = choice.lower()
else:
    choice = None

if amino_acid.lower() in amino_acids:
    amino_acid = one_letter_codes[amino_acids.index(amino_acid.lower())]
elif amino_acid.lower() in three_letter_codes:
    amino_acid = one_letter_codes[three_letter_codes.index(amino_acid.lower())]
else:
    amino_acid = amino_acid.upper()

# defining data structures
gene_mutation_names = [] # stores dictionaries 
hist_dict = {} # stores a mapping of COSMIC histology identifiers to tissue type strings
gene_mutation_list = {} # stores a list of mutations (Mutation objects)(value) occurring in each gene (key)
gene_total_dict = {} # stores the total number of samples for each gene
hist_total_dict = {} # stores the total number of samples for each tissue type

print("Starting...")

# populating hist_dict and hist_total_dict
with open(classification_file, 'r') as class_file:
    header = class_file.readline()
    for line in class_file:
        sample_parser(line, hist_dict, hist_total_dict)

all_mutation_types = 0 # incremented in re_capture() every time a mutation of correct amino acid (but not necessarily mutation type) is encountered
amino_acid_collapse = 0 # incremented in re_capture() every time two mutations are collapsed into one based on similarity

# finding the number of mutations in file 1
with open(filepath_genomic, 'rb') as mut_file1:
    FILE_LENGTH1 = -1 # subtracting the header
    for line in mut_file1:
        FILE_LENGTH1 += 1

# parsing every line in file 1 using re_capture()
with open(filepath_genomic, 'r') as mut_file1:
    progress_bar = tqdm(total=FILE_LENGTH1, desc="  Parsing genomic mutation screen... ", unit=" mutations")
    header = mut_file1.readline()
    for line in mut_file1:
        re_capture(line, amino_acid, mutation_type)
        progress_bar.update(1)

# closing progress bar and flush to stdout
progress_bar.close()
sys.stdout.flush()
print("Scanning " + filepath_genomic + " done.")

# finding the number of mutations in file 2
with open(filepath_targeted, 'r') as mut_file2:
    FILE_LENGTH2 = -1 # subtracting the header
    for line in mut_file2:
        FILE_LENGTH2 += 1

# parsing every line in file 2 using re_capture()
with open(filepath_targeted, 'r') as mut_file2:
    progress_bar = tqdm(total=FILE_LENGTH2, desc="  Parsing targeted mutation screen... ", unit=" mutations")
    header = mut_file2.readline()
    for line in mut_file2:
        re_capture(line, amino_acid, mutation_type)
        progress_bar.update(1)

# closing progress bar and flushing to stdout
progress_bar.close()
sys.stdout.flush()
print("Scanning " + filepath_targeted + " done.")

# creating another progress bar
progress_bar = tqdm(total=len(gene_mutation_names), desc="  Calculating mutation frequencies... ", unit=" genes")

# obtaining mutation frequencies
for gene in gene_mutation_names:
    gene_mutation_count = 0 # total number of samples in the current gene
    sorting_dict = {} # used for ranking mutations based on the number of samples

    for mutation in gene_mutation_list[gene]:
        # find total number of mutations associated with the current gene
        mutation.total_count = len(gene_mutation_list[gene])
        gene_mutation_count += mutation.sample_count # add sample count of current mutation

        # inserting mutation into the corresponding list in sorting_dict with number of samples as the key
        if mutation.sample_count not in sorting_dict:
            sorting_dict[mutation.sample_count] = []
            sorting_dict[mutation.sample_count].append(mutation)
        else:
            sorting_dict[mutation.sample_count].append(mutation)
    
    # sorting the keys in sorting dict based on descending order and turning it into a list
    sorting_list = sorted(sorting_dict.items(), reverse=True)
    # traversing the resulting list and initializing the rank, frequency, and E-score values
    for i in range(len(sorting_list)):
        for mutation in sorting_list[i][1]:
            mutation.rank = i + 1 # i.e., mutation at index 0 would be assigned rank 1
            # calculating frequencies and E-scores
            mutation.gsf = (mutation.sample_count / gene_total_dict[mutation.name]) * 100
            mutation.rsf = (mutation.sample_count / gene_mutation_count) * 100
            mutation.tgf = (mutation.sample_count / hist_total_dict[mutation.histology][1][mutation.name]) * 100
            mutation.tsf = (mutation.sample_count / hist_total_dict[mutation.histology][0]) * 100
            mutation.escore_gsf = (mutation.sample_count * (mutation.gsf / 100))
            mutation.escore_rsf = (mutation.sample_count * (mutation.rsf / 100))
            mutation.escore_tgf = (mutation.sample_count * (mutation.tgf / 100))
            mutation.escore_tsf = (mutation.sample_count * (mutation.tsf / 100))
    progress_bar.update(1)

# closing progress bar and flushing to stdout
progress_bar.close()
sys.stdout.flush()
print("Calculating mutation frequencies done.")

# transferring mutations from gene_mutation_list (list of lists) 
# to a linear list for easier traversal/analysis
mutations_list = []
for gene in gene_mutation_list:
    for mutation in gene_mutation_list[gene]:
        mutations_list.append(mutation)

# sorting mutations based on E-score (TGF)
mutations_list_sorted = sorted(mutations_list, key = lambda x : x.escore_tgf)
after_collapsing = len(mutations_list_sorted) # final number of filtered mutations
total_amino_acid = after_collapsing + amino_acid_collapse # number of mutations before collapsing

if after_collapsing == 0:
    print("No results :(") # T.T
    sys.exit()

# formatting all the data for printing and writing to output file line by line
if args.output:
    output_file = args.output
    out_file = open(args.output, "w")
else:
    # trying to create an output file name that won't overwrite a file already on the system
    output_file = "out"
    if os.path.exists(output_file + ".tsv"):
            postfix = 1
            while os.path.exists(output_file + "_" + str(postfix) + ".tsv"):
                postfix += 1
    output_file = output_file + "_" + str(postfix) + ".tsv"
    out_file = open(output_file, "w") 

# creating a data structure for output data
data = []
data.append(['GENE', 'COSG IDENTIFIER', 'TYPE', 'MUTATION', 'HISTOLOGY', 'SAMPLE COUNT', 'GENE-SPECIFIC FREQUENCY (%)', 'RESIDUE-SPECIFIC FREQUENCY (%)', 'TISSUE-GENE FREQUENCY (%)', 'TISSUE-SPECIFIC FREQUENCY', 'E-SCORE (GENE)', 'E-SCORE (RESIDUE)', 'E-SCORE (TISSUE-GENE)', 'E-SCORE (TISSUE)', 'RANK'])
for mutation in mutations_list_sorted:
    # if tissue type and cutoff  conditions are passed, then we write the mutation to output
    if (args.tissue == None or mutation.histology == args.tissue) and (args.cutoff == None or mutation.sample_count >= args.cutoff):
        # formatting a list of attributes of an individual mutation
        data_line = [mutation.name, "COSG" + mutation.identifier, mutation.mut_type, mutation.mutation[:-1], mutation.histology, mutation.sample_count, round(mutation.gsf, 3), round(mutation.rsf, 3), round(mutation.tgf, 3), round(mutation.tsf, 5), round(mutation.escore_gsf, 3), round(mutation.escore_rsf, 3), round(mutation.escore_tgf, 3), round(mutation.escore_tsf, 5), str(mutation.rank) + " / " + str(mutation.total_count)]
        # adding formatted list to data structure
        data.append(data_line)

# setting up csv.writer() to write a TSV file
writer = csv.writer(out_file, delimiter='\t')
for line in data:
    writer.writerow(line) # writing each line to output file
print("Output written to " + output_file + ".")

out_file.close()

total_mutations = FILE_LENGTH1 + FILE_LENGTH2
target_tissue = len(data) - 1 # subtracting the header
# specifying tissue type for graphics
if args.tissue == None:
    tissue = "all"
else:
    tissue = args.tissue

# setting up data structure for plotly
funnel_data = dict(
    # numerical values
    number = [round(total_mutations, 3), round(all_mutation_types, 3), round(total_amino_acid, 3), round(after_collapsing, 3), round(target_tissue, 3)],
    # side labels
    stage = ["Initial mutations", "Target amino acids (" + amino_acid + ")", "Target mutation type (" + mutation_type + ")", "After collapsing ", "Target tissue type (" + tissue + ") + over cutoff"]
)

if args.disable_collapse:
    funnel_data["number"].pop(3)
    funnel_data["stage"].pop(3)

if tissue == "all":
    funnel_data["number"].pop(-1)
    funnel_data["stage"].pop(-1)

# plot figure
fig = px.funnel(funnel_data, x='number', y='stage')
# convert to a PNG
image_bytes = fig.to_image(format="png")
# write to an output file
with open("funnel_plot_" + tissue + ".png", "wb") as file:
    file.write(image_bytes)

# repeating the same but with log values for easier visualization
funnel_data = dict(
    number = [round(math.log(total_mutations, 10), 3), round(math.log(all_mutation_types, 10), 3), round(math.log(total_amino_acid, 10), 3), round(math.log(after_collapsing, 10), 3), round(math.log(target_tissue, 10), 3)],
    stage = ["Initial mutations", "Target amino acids (" + amino_acid + ")", "Target mutation type (" + mutation_type + ")", "After collapsing ", "Target tissue type (" + tissue + ") + over cutoff"]
)

if args.disable_collapse:
    funnel_data["number"].pop(3)
    funnel_data["stage"].pop(3)

if tissue == "all":
    funnel_data["number"].pop(-1)
    funnel_data["stage"].pop(-1)

fig = px.funnel(funnel_data, x='number', y='stage')
image_bytes = fig.to_image(format="png")
with open("funnel_plot_" + tissue + "_log" + ".png", "wb") as file:
    file.write(image_bytes)

print("Saved funnel plot output.")
print(" *** REPORT ***") # printing report
print("Total mutations: " + str(total_mutations) + "\nAfter filtering non-target amino acids: " + str(all_mutation_types) + "\nAfter filtering non-target mutation types: " + str(total_amino_acid))
if not args.disable_collapse:
    print("After collapsing similar mutations: " + str(after_collapsing))
if tissue != "all":
    print("After filtering non-target tissue types: " + str(target_tissue))
# done!
