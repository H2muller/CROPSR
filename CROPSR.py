#!/usr/bin/env python3

### Importing required libraries
import argparse
import cropsr_functions
import re
import sys
from multiprocessing import cpu_count, Pool
import pandas as pd
import numpy as np
from numpy import vectorize
from time import gmtime, strftime
import random
import string
import csv
import array
import time

### CROPSR Version
__version__ = '1.11b'


### Defining the arguments
parser = argparse.ArgumentParser()

parser.add_argument('-f', '--fasta', metavar='', required=True, dest='f', 
                    help='[required] path to input file in FASTA format'
                    )
parser.add_argument('-g', '--gff', metavar='', dest='g', 
                    help='path to input file in GFF format'
                    )
parser.add_argument('-p', '--phytozome', metavar='', dest='p', default=None,
                    help='path to input annotation info file in TXT format, default = None'
                    )
parser.add_argument('-o', '--output', metavar='', dest='o', default='data.csv',
                    help='path to output file, default = data.csv'
                    )
parser.add_argument('-l', '--length', metavar='', dest='l', type=int, default=20,
                    help='length of the gRNA se3quence, default = 20'
                    )
parser.add_argument('-L', '--flanking', metavar='', dest='L', type=int, default=200,
                    help='length of flanking region for verification, default = 200'
                    )
parser.add_argument('--cas9', action='store_true',
                    help='specifies that design will be made for the Cas9 CRISPR system'
                    )
parser.add_argument('-v', '--verbose', action='store_true',
                    help='prints visual indicators for each iteration'
                    )

args = parser.parse_args()


def import_fasta_file(fasta):
    '''
    imports and formats a genome file in the fasta format for use
    '''
    with open(fasta, 'r') as f:
        myFASTA = f.read()
        if args.verbose:
            print(f'Genome file {fasta} successfully imported')
        linecount = myFASTA.count('\n')
        if 2 * myFASTA.count('>') != linecount + 1:
            if args.verbose:
                print('formatting genome')
            from cropsr_functions import formatted
            myFASTA = formatted(myFASTA)
            if args.verbose:
                print(f'Genome file {fasta} successfully formatted')
    from cropsr_functions import generate_dictionary as gendict
    genome_dictionary = gendict(myFASTA)
    if args.verbose:
        print(f'The genome was successfully converted to a dictionary')
    return genome_dictionary


def import_gff_file(gff):
    import pandas as pd
    '''
    imports and formats a genome annotation file in the GFF format for use
    '''
    start_index = 0
    with open(gff,"r") as raw_gff:
        if args.verbose:
            print(f'Annotation file {gff} successfully imported')
        gff_lines = raw_gff.readlines()
        for index in range(len(gff_lines)):
            if ("##" not in gff_lines[index]):
                start_index = index
                break
    col_names = ["chromosome", "source", "feature", "start", "end", "score", "strand", "phase", "attributes"]
    gff_df = pd.read_csv(gff, sep='\t', skiprows = start_index, header = None, names = col_names)
    if args.verbose:
        print(f'Annotation database successfully generated')
    return gff_df


def find_PAM_site(target,input_sequence):
    import re
    '''
    locates the target PAM motif in input sequence
    '''
    PAM_site = [match.span() for match in re.finditer(target,input_sequence)]
    return PAM_site


def convert_seq_to_int(input_sequence):
    '''
    converts a DNA sequence to the coded integers
    '''
    new_sequence = str(input_sequence).replace('A',1).replace('T',2).replace('C',3).replace('G',4)
    int_seq = int(new_sequence)
    return int_seq


def get_reverse_complement(input_sequence):
    '''
    converts a DNA sequence to its reverse complement
    '''
    bases = input_sequence.replace('A','U').replace('C','Z').replace('G','C').replace('Z','G').replace('T','A').replace('U','T')[::-1]
    return bases


def get_gRNA_sequence(input_sequence):
    '''
    converts a DNA sequence to its complimentary RNA sequence
    '''
    RNA = input_sequence.replace('A','U').replace('C','Z').replace('G','C').replace('Z','G').replace('T','A')[::-1]
    return RNA


def create_dataframe():
    import pandas as pd
    '''
    creates a dataframe to store sgRNA information
    '''
    df_cols = [
                'crispr_id',        # STR
                'crispr_sys',       # CAT
                'sequence',         # STR
                'long_sequence',    # STR
                'chromosome',       # CAT
                'start_pos',        # INT
                'end_pos',          # INT
                'cutsite',          # INT
                'strand',           # CAT
                'on_site_score',    # FLOAT
                'features',         # LIST
                'status'               # STR
                ]
    df = pd.DataFrame(columns=df_cols)
    return df


def apply_cutsite(start_pos, end_pos, crispr_sys):
    if crispr_sys == 'cas9':
        cutsite = end_pos-3
    return cutsite


intersect = 0.59763615
low_gc = -0.2026259
high_gc = -0.1665878

first_matrix = np.array([ 0.        ,  0.        ,  0.        ,  0.        ,  0.        ,
        0.        ,  0.        , -0.2753771 , -0.3238875 ,  0.        ,
        0.17212887,  0.        ,  0.        ,  0.        , -0.1006662 ,
        0.        ,  0.        ,  0.        , -0.2018029 ,  0.24595663,
        0.03644004,  0.        ,  0.09837684,  0.        ,  0.        ,
        0.        , -0.7411813 , -0.3932644 ,  0.        ,  0.        ,
        0.        ,  0.        ,  0.        ,  0.        ,  0.        ,
        0.        ,  0.        ,  0.        ,  0.        ,  0.        ,
        0.        ,  0.        ,  0.        ,  0.        , -0.466099  ,
        0.        ,  0.        ,  0.        ,  0.        ,  0.        ,
        0.        ,  0.        ,  0.        ,  0.        ,  0.        ,
        0.        ,  0.08537695,  0.        , -0.013814  ,  0.        ,
        0.27262051, -0.2859442 ,  0.1190226 ,  0.        ,  0.09745459,
        0.        ,  0.        , -0.1755462 ,  0.        ,  0.        ,
       -0.3457955 , -0.6780964 ,  0.22508903,  0.        , -0.5077941 ,
        0.        ,  0.        , -0.054307  ,  0.        , -0.4173736 ,
        0.        , -0.0907126 ,  0.        ,  0.37989937,  0.        ,
       -0.5305673 ,  0.05782332,  0.        ,  0.        , -0.8770074 ,
        0.        ,  0.        ,  0.        , -0.4031022 , -0.8762358 ,
        0.27891626, -0.0773007 , -0.2216372 ,  0.28793562,  0.        ,
        0.        ,  0.        ,  0.        ,  0.        ,  0.        ,
        0.        ,  0.        ,  0.        ,  0.        ,  0.11787758,
        0.        , -0.6890167 ,  0.        ,  0.        , -0.1604453,
        0.        ,  0.        ,  0.        ,  0.        ,  0.38634258 ])

second_matrix = np.array([  
        0.        ,  0.        ,  0.        ,  0.        ,  0.        ,
        0.        ,  0.        ,  0.        ,  0.        ,  0.        ,
        0.        ,  0.        ,  0.        ,  0.        ,  0.        ,
        0.        ,  0.        ,  0.        ,  0.        ,  0.        ,
        0.        ,  0.        ,  0.        ,  0.        ,  0.        ,
        0.        ,  0.        ,  0.        ,  0.        , -0.6257787 ,
        0.        ,  0.        ,  0.        ,  0.        ,  0.        ,
        0.        ,  0.        ,  0.        ,  0.        ,  0.        ,
        0.        ,  0.        ,  0.        ,  0.        ,  0.        ,
        0.        ,  0.        ,  0.        ,  0.        ,  0.        ,
        0.        ,  0.        ,  0.        ,  0.        ,  0.        ,
        0.        ,  0.        ,  0.        ,  0.        ,  0.        ,
        0.        ,  0.        ,  0.        ,  0.        ,  0.        ,
        0.        ,  0.        ,  0.        ,  0.        ,  0.        ,
        0.        ,  0.        ,  0.        ,  0.        ,  0.        ,
        0.        ,  0.        ,  0.        ,  0.30004332,  0.        ,
       -0.8348362 ,  0.        ,  0.        ,  0.        ,  0.76062777,
        0.        ,  0.        ,  0.        ,  0.        ,  0.        ,
        0.        ,  0.        ,  0.        ,  0.        ,  0.        ,
        0.        ,  0.        ,  0.        ,  0.        ,  0.        ,
        0.        ,  0.        ,  0.        ,  0.        ,  0.        ,
        0.        ,  0.        ,  0.        ,  0.        ,  0.        ,
        0.        , -0.4908167 ,  0.        ,  0.        ,  0.        ,
        0.        ,  0.        ,  0.        ,  0.        ,  0.        ,
        0.        ,  0.        ,  0.        ,  0.        ,  0.        ,
        0.        ,  0.        ,  0.        ,  0.        ,  0.        ,
        0.        ,  0.        ,  0.        ,  0.        ,  0.        ,
        0.        ,  0.        ,  0.        ,  0.        ,  0.        ,
        0.        ,  0.        ,  0.        ,  0.        ,  0.        ,
        0.        ,  0.        ,  0.        ,  0.        ,  0.        ,
        0.        ,  0.        ,  0.        ,  0.        ,  0.        ,
        0.        ,  0.        ,  0.        ,  0.        ,  0.        ,
        0.        ,  0.        ,  0.        ,  0.        ,  0.        ,
        0.        ,  0.        ,  0.        ,  0.        ,  0.        ,
        0.        ,  0.        ,  0.        ,  0.        ,  0.        ,
        0.        ,  0.        ,  0.        ,  0.        ,  0.        ,
        0.7092612 , -0.5868739 ,  0.49629861,  0.        ,  0.        ,
        0.        ,  0.        ,  0.        ,  0.        ,  0.        ,
        0.        , -1.5169074 ,  0.        ,  0.        ,  0.        ,
        0.        ,  0.        ,  0.        ,  0.        ,  0.        ,
        0.        ,  0.        ,  0.        ,  0.        ,  0.        ,
        0.        ,  0.        , -0.3345637 ,  0.        ,  0.        ,
        0.        ,  0.        ,  0.        ,  0.        ,  0.        ,
        0.        ,  0.        ,  0.        ,  0.        ,  0.        ,
        0.76384993,  0.        , -0.5370252 ,  0.        ,  0.        ,
        0.        ,  0.        ,  0.        ,  0.        ,  0.        ,
        0.        ,  0.        ,  0.        ,  0.        ,  0.        ,
        0.        ,  0.        ,  0.        ,  0.        ,  0.        ,
        0.        ,  0.        ,  0.        ,  0.        ,  0.        ,
        0.        ,  0.        ,  0.        ,  0.        ,  0.        ,
        0.        ,  0.        ,  0.        ,  0.        ,  0.        ,
        0.        ,  0.        ,  0.        ,  0.        ,  0.        ,
        0.        ,  0.        ,  0.        , -0.7981461 ,  0.        ,
        0.        ,  0.        ,  0.        ,  0.        ,  0.        ,
        0.        ,  0.        ,  0.        ,  0.        ,  0.        ,
        0.        ,  0.        ,  0.        ,  0.        ,  0.        ,
        0.        ,  0.        ,  0.        ,  0.        ,  0.        ,
        0.        ,  0.        ,  0.        ,  0.        ,  0.        ,
        0.        ,  0.        ,  0.        ,  0.        ,  0.35318325,
        0.        ,  0.        ,  0.        ,  0.        ,  0.        ,
        0.        ,  0.        ,  0.        , -0.6668087 ,  0.        ,
        0.        ,  0.        ,  0.        ,  0.        ,  0.        ,
        0.        , -0.3672668 ,  0.        ,  0.        ,  0.74807209,
        0.        ,  0.        ,  0.        ,  0.        ,  0.        ,
        0.        ,  0.        ,  0.56820913,  0.        ,  0.        ,
        0.        ,  0.        ,  0.        ,  0.        ,  0.        ,
        0.        ,  0.32907207, -0.8364568 ,  0.        ,  0.        ,
       -0.7822076 ,  0.        ,  0.        ,  0.        ,  0.        ,
        0.        ,  0.        , -1.029693  ,  0.        ,  0.        ,
        0.        ,  0.        ,  0.        ,  0.        ,  0.        ,
        0.        ,  0.        ,  0.        ,  0.        ,  0.        ,
        0.        ,  0.        ,  0.        ,  0.        ,  0.        ,
        0.        , -0.4632077 ,  0.        ,  0.85619782,  0.        ,
        0.        ,  0.        ,  0.        , -0.5794924 ,  0.        ,
        0.        ,  0.64907554,  0.        ,  0.        ,  0.        ,
        0.        ,  0.        ,  0.        ,  0.        ,  0.        ,
        0.        ,  0.        ,  0.        ,  0.        ,  0.        ,
        0.        ,  0.        , -0.0773007 ,  0.        ,  0.        ,
        0.        , -0.2216372 ,  0.        ,  0.        ,  0.        ,
        0.28793562,  0.        ,  0.        ,  0.        ,  0.        ,
        0.        ,  0.        ,  0.        ,  0.        ,  0.        ,
        0.        ,  0.        ,  0.        ,  0.        ,  0.        ,
        0.        ,  0.        ,  0.        ,  0.        ,  0.        ,
        0.        ,  0.        ,  0.        ,  0.        ,  0.        ,
        0.        ,  0.        ,  0.        ,  0.        ,  0.        ,
        0.        ,  0.        ,  0.        ,  0.        ,  0.11787758,
        0.        ,  0.        ,  0.        ,  0.        ,  0.        ,
        0.        ,  0.        ,  0.        ,  0.        ,  0.        ,
        0.        ,  0.        ,  0.        ,  0.        ,  0.        ,
        0.        ,  0.        ,  0.        ,  0.        ,  0.        ,
        0.        ,  0.        ,  0.        ,  0.        ,  0.        ,
        0.        ,  0.        ,  0.        ,  0.        ,  0.        ,
        0.        ,  0.        ,  0.        , -0.69774   ])

def rs1_score(sequences):

    size = len(sequences)
    seq1 = sequences[:,0:29]
    seq2 = sequences[:,1:30]

    sequences = np.repeat(sequences,4,axis=1)
    seq1 = np.repeat(seq1,16,axis=1)
    seq2 = np.repeat(seq2,16,axis=1)

    matrix1 = np.empty((size,120))
    matrix2 = np.empty((size,464))
    temp1 = np.empty((size,464))
    temp2 = np.empty((size,464))

    matrix_1_cmp = np.array([[65, 84, 67, 71]*30])
    matrix_2_cmp1 = np.array([[65,65,65,65,84,84,84,84,67,67,67,67,71,71,71,71]*29])
    matrix_2_cmp2 = np.array([[65,84,67,71,65,84,67,71,65,84,67,71,65,84,67,71]*29])

    np.equal(sequences,matrix_1_cmp,out=matrix1)
    score_first = np.matmul(matrix1,first_matrix)

    np.equal(seq1,matrix_2_cmp1,out=temp1)
    np.equal(seq2,matrix_2_cmp2,out=temp2)
    np.logical_and(temp1,temp2,out=matrix2)

    score_second = np.matmul(matrix2,second_matrix)
    score = (score_first + score_second + intersect + low_gc) * -1
    return 1/(1 + np.exp(score))


alphanum = np.array(list('ABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789'), dtype="|U1")
def get_id(num_to_gen):
    return np.random.choice(alphanum, [num_to_gen, 7])


def preprocess_PAM_sites(DF):
    DF['crispr_sys'] = DF['raw'][5]
    DF['sequence'] = DF['raw'][3]
    DF['long_sequence'] = DF['raw'][4]
    DF['chromosome'] = DF['raw'][2]
    DF['start_pos'] = DF['raw'][0]
    DF['end_pos'] = DF['raw'][1]
    DF['strand'] = DF['raw'][6]
    DF['status'] = 'completed'
    return DF


def main():
    begin = time.time()
    if not args.cas9:
        sys.exit('Please select at least one CRISPR system: Cas9') 
    
    if args.verbose:
        print(f"""
################################################################################
##                                                                            ##
##                                                                            ##
##          .o88b.   d8888b.    .d88b.    d8888b.   .d8888.   d8888b.         ##
##         d8P  Y8   88  `8D   .8P  Y8.   88  `8D   88'  YP   88  `8D         ##
##         8P        88oobY'   88    88   88oodD'   `8bo.     88oobY'         ##
##         8b        88`8b     88    88   88ººº       `Y8b.   88`8b           ##
##         Y8b  d8   88 `88.   `8b  d8'   88        db   8D   88 `88.         ##
##          `Y88P'   88   YD    `Y88P'    88        `8888Y'   88   YD         ##
##                                                                            ##
##                                                                            ##
################################################################################
U.S. Dept. of Energy's Center for Advanced Bioenergy and Bioproducts Innovation
University of Illinois at Urbana-Champaign

        You are currently utilizing the following settings:

        CROPSR version:                                 {__version__}
        Path to genome file in FASTA format:            {args.f}
        Path to output file:                            {args.o}
        Length of the gRNA sequence:                    {args.l}
        Length of flanking region for verification:     {args.L}
        Number of available CPUs:                       {cpu_count()}
        Path to annotation file in GFF format:          {args.g}
        Path to annotation_info file in TXT format:     {args.p}
        Designing for CRISPR system:
            Streptococcus pyogenes Cas9                 {args.cas9}
        """)


    # Set up timing
    file1 = open("time.txt","w")

    ### Import genome files
    fasta_file = import_fasta_file(args.f)
    gff_df = import_gff_file(args.g)

    ### Locate PAMs by nuclease type
    if args.verbose:
        print(f'''
            Initiating PAM site detection.
            
            Please wait, this may take a while...
            ''')

    ### Create Dataframe containing all PAM site information
    data = [
                'crispr_id',        # STR
                'crispr_sys',       # CAT
                'sequence',         # STR
                'long_sequence',    # STR
                'chromosome',       # CAT
                'start_pos',        # INT
                'end_pos',          # INT
                'cutsite',          # INT
                'strand',           # CAT
                'on_site_score',    # FLOAT
                'features',          # LIST
                'status', # STR
                ]

    # Set up output CSV file
    with open(args.o, 'w') as file:
        writer = csv.writer(file)
        writer.writerow(data)
    file.close()
    
    Complete_dataset = []

    for chromosome,sequence in fasta_file.items():
        print("Searching on Chromosome: ", chromosome[:25])
        print("With start of sequence: ", sequence[:25])

        if args.cas9:
            # + strand
            motif = re.compile(r'(?=.GG)')
            cas9_target_list = find_PAM_site(motif,sequence)
            for target in cas9_target_list:
                pam_location = (target[0]-args.l, target[0]) # Updated Regex for proper PAM Site alignment
                if pam_location[0] >= 5 and pam_location[0]+5 <= len(sequence)+10 and pam_location[1] >= 5 and pam_location[1] <= len(sequence)+10:
                    shortseq = get_gRNA_sequence(sequence[pam_location[0]:pam_location[1]])
                    longseq = get_gRNA_sequence(sequence[pam_location[0]-5:pam_location[1]+5])
                    crispr_guide = [pam_location[0],pam_location[1],chromosome[1::],shortseq,longseq,'cas9','+']
                    Complete_dataset.append(crispr_guide)

            # - strand
            motif = re.compile(r'(?=CC.)')
            cas9_target_list2 = find_PAM_site(motif,sequence)
            for target in cas9_target_list2:
                pam_location = (target[0]+3, target[0]+3+args.l) # Updated Regex for proper PAM Site alignment
                if pam_location[0] >= 5 and pam_location[0]+5 <= len(sequence)+10 and pam_location[1] >= 5 and pam_location[1] <= len(sequence)+10:
                    shortseq = get_gRNA_sequence(get_reverse_complement(sequence[pam_location[0]:pam_location[1]]))
                    longseq = get_gRNA_sequence(get_reverse_complement(sequence[pam_location[0]-5:pam_location[1]+5]))
                    crispr_guide = [pam_location[1],pam_location[0],chromosome[1::],shortseq,longseq,'cas9','-']
                    Complete_dataset.append(crispr_guide)

            if args.verbose:
                print (f'''
                {len(cas9_target_list + cas9_target_list2):n} Cas9 PAM sites were found on {chromosome[1::]}
                ''')


        size = len(Complete_dataset)
        count = 0

        # Score sequences, fill rows, and manually write to CSV
        with open(args.o, 'a') as file:
            writer = csv.writer(file)
            ids = get_id(size)
            ids= [array.array('B', map(ord,z)).tobytes().decode("utf-8") for z in ids.tolist()]
            counter = 0
            for i in range(size): 
                count += 1
                if ((count == 1000000 and i < size-1) or (count < 1000000 and i == size-1)):
                    index_range = count*counter 

                    lesser_list = Complete_dataset[index_range:index_range+count]

                    sequences = [np.frombuffer(bytes(str(item[4].replace('U','T')).upper(),"ascii"), 'uint8') if len(item[4]) == 30
                    else np.empty(30,) for item in lesser_list ]

                    score = rs1_score(np.array(sequences))

                    write_csv = [ (ids[index_range-index-1], lesser_list[index][5], lesser_list[index][3], lesser_list[index][4], 
                    lesser_list[index][2], lesser_list[index][0], lesser_list[index][1], 
                    apply_cutsite(lesser_list[index][0],lesser_list[index][1],lesser_list[index][5]), lesser_list[index][6], 
                    score[index],'','completed') if len(lesser_list[index][4]) == 30 else 
                    (ids[index_range-index-1], lesser_list[index][5], lesser_list[index][3], lesser_list[index][4], 
                    lesser_list[index][2], lesser_list[index][0], lesser_list[index][1], lesser_list[index][6], -1,'','completed') 
                    for index in range(len(lesser_list)) ]

                    count = 0
                    counter += 1

                    writer.writerows(write_csv)

        end = time.time()
        file1.write("Total runtime of the program is " + str(end-begin))
        time.sleep(5)
        #file1.close()

        file.close()

                
    ### CONFIRMATION MESSAGE
    if args.verbose:
        print(f'The output file has been generated at {args.o}')


if __name__ == '__main__':
    main()
