#!/usr/bin/env python3

### Importing required libraries
import argparse
import cropsr_functions
import re
import sys
from multiprocessing import cpu_count, Pool
import pandas as pd
from numpy import vectorize
from time import gmtime, strftime

### Defining the arguments
parser = argparse.ArgumentParser()

parser.add_argument('-f', '--fasta', required=True, dest='f', 
                    help='[required] path to input file in FASTA format'
                    )
parser.add_argument('-g', '--gff', metavar='', dest='g', 
                    help='path to input file in GFF format'
                    )
parser.add_argument('-p', '--phytozome', metavar='', dest='p', default=None,
                    help='path to input annotation info file in TXT format, default = None'
                    )
parser.add_argument('-o', '--output', dest='o', default='data.csv',
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
# parser.add_argument('--cpf1', action='store_true',
#                     help='specifies that design will be made for the Cpf1 CRISPR system'
#                     )
# parser.add_argument('--CUDA', action='store_true',
#                     help='runs processing steps utilizing the GPU instead of CPU where possible'
#                     )
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
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
    bases = list(input_sequence)
    bases = reversed([complement.get(base, base) for base in bases])
    bases = ''.join(bases)
    return bases


def get_gRNA_sequence(input_sequence):
    '''
    converts a DNA sequence to its complimentary RNA sequence
    '''
    complement = {'A': 'U', 'C': 'G', 'G': 'C', 'T': 'A'}
    RNA = list(input_sequence)
    RNA = reversed([complement.get(base, base) for base in RNA])
    RNA = ''.join(RNA)
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
                'features'          # LIST
                ]
    df = pd.DataFrame(columns=df_cols)
    return df


def apply_cutsite(start_pos, end_pos, crispr_sys):
    if crispr_sys == 'cas9':
        cutsite = end_pos-3
    elif crispr_sys == 'cpf1':
        cutsite = end_pos-5
    return cutsite


def rs1_score(sequence):
    """
    Generates a binary matrix for DNA/RNA sequence, where each column is a possible base
    and each row is a position along the sequence. Matrix column order is A, T/U, C, G
    """
    import math
    import numpy as np
    seq = str(sequence).upper()
    seq = list(seq)
    matrix1  = np.zeros([len(sequence),4], dtype=int)
    for i,item in enumerate(sequence):
        if item == 'A':
            matrix1[i,0] = 1
        if item == 'T':
            matrix1[i,1] = 1
        if item == 'U':
            matrix1[i,1] = 1
        if item == 'C':
            matrix1[i,2] = 1
        if item == 'G':
            matrix1[i,3] = 1


    """
    Generates a binary matrix for DNA/RNA sequence, where each column is a possible
    pair of adjacent bases, and each row is a position along the sequence.
    Matrix column order is AA, AT, AC, AG, TA, TT, TC, TG, CA, CT, CC, CG, GA, GT, GC, GG
    """
    sequence = sequence.replace('U','T')
    pairwise_sequence = []
    for i in range(len(sequence)):
        if i < len(sequence)-1:
            basepair = sequence[i]+sequence[i+1]
            pairwise_sequence.append(basepair)
    matrix2 = np.zeros([len(pairwise_sequence),16], dtype=int)
    for i,item in enumerate(pairwise_sequence):
        if item == 'AA':
            matrix2[i,0] = 1
        if item == 'AT':
            matrix2[i,1] = 1
        if item == 'AC':
            matrix2[i,2] = 1
        if item == 'AG':
            matrix2[i,3] = 1
        if item == 'TA':
            matrix2[i,4] = 1
        if item == 'TT':
            matrix2[i,5] = 1
        if item == 'TC':
            matrix2[i,6] = 1
        if item == 'TG':
            matrix2[i,7] = 1
        if item == 'CA':
            matrix2[i,8] = 1
        if item == 'CT':
            matrix2[i,9] = 1
        if item == 'CC':
            matrix2[i,10] = 1
        if item == 'CG':
            matrix2[i,11] = 1
        if item == 'GA':
            matrix2[i,12] = 1
        if item == 'GT':
            matrix2[i,13] = 1
        if item == 'GC':
            matrix2[i,14] = 1
        if item == 'GG':
            matrix2[i,15] = 1


    """
    Scoring matrix
    """
    intersect = 0.59763615
    low_gc = -0.2026259
    high_gc = -0.1665878

    first_order = ['G02','A03','C03','C04','C05',
                    'G05','A06','C06','C07','G07',
                    'A12','A15','C15','A16','C16',
                    'T16','A17','G17','C18','G18',
                    'A19','C19','G20','T20','G21',
                    'T21','C22','T22','T23','C24',
                    'G24','T24','A25','C25','T25',
                    'G28','T28','C29','G30']
    first_scores = [-0.2753771,-0.3238875,0.17212887,-0.1006662,-0.2018029,
                    0.24595663,0.03644004,0.09837684,-0.7411813,-0.3932644,
                    -0.466099,0.08537695,-0.013814,0.27262051,0.1190226,
                    -0.2859442,0.09745459,-0.1755462,-0.3457955,-0.6780964,
                    0.22508903,-0.5077941,-0.4173736,-0.054307,0.37989937,
                    -0.0907126,0.05782332,-0.5305673,-0.8770074,-0.8762358,
                    0.27891626,-0.4031022,-0.0773007,0.28793562,-0.2216372,
                    -0.6890167,0.11787758,-0.1604453,0.38634258]
    first_order_scores = dict(zip(first_order,first_scores))

    second_order = ['GT02','GC05','AA06','TA06','GG07',
                    'GG12','TA12','TC12','TT12','GG13',
                    'GA14','GC14','TG17','GG19','TC19',
                    'CC20','TG20','AC21','CG21','GA21',
                    'GG21','TC22','CG23','CT23','AA24',
                    'AG24','AG25','CG25','TG25','GT27',
                    'GG29']
    second_scores = [-0.6257787,0.30004332,-0.8348362,0.76062777,-0.4908167,
                    -1.5169074,0.7092612,0.49629861,-0.5868739,-0.3345637,
                    0.76384993,-0.5370252,-0.7981461,-0.6668087,0.35318325,
                    0.74807209,-0.3672668,0.56820913,0.32907207,-0.8364568,
                    -0.7822076,-1.029693,0.85619782,-0.4632077,-0.5794924,
                    0.64907554,-0.0773007,0.28793562,-0.2216372,0.11787758,
                    -0.69774]
    second_order_scores = dict(zip(second_order,second_scores))


    # order 1 score matrix
    """ row order == A T/U C G """
    first_matrix = np.zeros([4,30], dtype=float)
    def posit(key):
        return int(key[1:])-1
    for k,v in first_order_scores.items():
        if k[0] == 'A':
            first_matrix[0,posit(k)] = v
        elif k[0] == 'T':
            first_matrix[1,posit(k)] = v
        elif k[0] == 'C':
            first_matrix[2,posit(k)] = v
        elif k[0] == 'G':
            first_matrix[3,posit(k)] = v


    # order 2 score matrix
    """ row order == AA AT AC AG TA TT TC TG CA CT CC CG GA GT GC GG """
    second_matrix = np.zeros([16,29], dtype=float)
    for k,v in second_order_scores.items():
        if k[0:2] == 'AA':
            second_matrix[0,int(k[2:])-1] = v
        if k[0:2] == 'AT':
            second_matrix[1,int(k[2:])-1] = v
        if k[0:2] == 'AC':
            second_matrix[2,int(k[2:])-1] = v
        if k[0:2] == 'AG':
            second_matrix[3,int(k[2:])-1] = v
        if k[0:2] == 'TA':
            second_matrix[4,int(k[2:])-1] = v
        if k[0:2] == 'TT':
            second_matrix[5,int(k[2:])-1] = v
        if k[0:2] == 'TC':
            second_matrix[6,int(k[2:])-1] = v
        if k[0:2] == 'TG':
            second_matrix[7,int(k[2:])-1] = v
        if k[0:2] == 'CA':
            second_matrix[8,int(k[2:])-1] = v
        if k[0:2] == 'CT':
            second_matrix[9,int(k[2:])-1] = v
        if k[0:2] == 'CC':
            second_matrix[10,int(k[2:])-1] = v
        if k[0:2] == 'CG':
            second_matrix[11,int(k[2:])-1] = v
        if k[0:2] == 'GA':
            second_matrix[12,int(k[2:])-1] = v
        if k[0:2] == 'GT':
            second_matrix[13,int(k[2:])-1] = v
        if k[0:2] == 'GC':
            second_matrix[14,int(k[2:])-1] = v
        if k[0:2] == 'GG':
            second_matrix[15,int(k[2:])-1] = v

    item_gc = sequence[0][5:-5]
    gc_count = item_gc.count('G') + item_gc.count('C')
    if gc_count < 10:
        gc_score = low_gc
    else:
        gc_score = high_gc
    first_first = np.matmul(first_matrix,matrix1)
    score_first = np.trace(first_first)
    score_second = np.trace(np.matmul(second_matrix,matrix2))
    score = (1/(1 + math.exp(-(intersect + gc_score + score_first + score_second))))
    return score


def get_id(sys_type, chr):
    'generates a unique ID for each target site'
    import random
    import string
    #  defining first digit
    if sys_type == 'cas9':
        first = 'A'
    elif sys_type == 'cpf1':
        first = 'B'
    else:
        first = '0'
    # defining two following digits
    second = str(chr[-2::])
    #defining remaining sequence
    remaining = ''.join(random.choices(string.ascii_letters + string.digits,k=7)).upper()
    return ''.join([first, second, remaining])


def preprocess_PAM_sites(DF):
    DF['crispr_sys'] = DF['raw'][5]
    DF['sequence'] = DF['raw'][3]
    DF['long_sequence'] = DF['raw'][4]
    DF['chromosome'] = DF['raw'][2]
    DF['start_pos'] = DF['raw'][0]
    DF['end_pos'] = DF['raw'][1]
    DF['strand'] = DF['raw'][6]
    DF['raw'] = 'completed'
    return DF


def fill_row(DF):
    from numpy import vectorize
    DF['crispr_id'] = vectorize(get_id)(DF['crispr_sys'],DF['chromosome'])
    if len(DF['long_sequence'])==30:
        DF['cutsite'] = vectorize(apply_cutsite)(DF['start_pos'],DF['end_pos'],DF['crispr_sys'])
        DF['on_site_score'] = vectorize(rs1_score)(DF['long_sequence'])
    else:
        DF['on_site_score'] = -1
        start_pos = DF['start_pos']
        end_pos = DF['end_pos']
        chrom = DF['chromosome']
        problem_seq = DF['long_sequence']
        print(f'length error occurred at guide in position {start_pos} - {end_pos} of {chrom}, sequence: {problem_seq}')
    return DF


def main():
    if not args.cas9 and not args.cpf1:
        sys.exit('Please select at least one CRISPR system: Cas9 or Cpf1')
    
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
##                                               Version: 1.01.1-beta         ##
##                                                                            ##
################################################################################
U.S. Dept. of Energy's Center for Advanced Bioenergy and Bioproducts Innovation
University of Illinois at Urbana-Champaign

        You are currently utilizing the following settings:

        Path to genome file in FASTA format:            {args.f}
        Path to output file:                            {args.o}
        Length of the gRNA sequence:                    {args.l}
        Length of flanking region for verification:     {args.L}
        Number of available CPUs:                       {cpu_count()}
        Path to annotation file in GFF format:          {args.g}
        Designing for CRISPR system:
            Streptococcus pyogenes Cas9                 {args.cas9}
            Prevotella/Francisella Cpf1:                {args.cpf1}
        """)


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
    CRIPSR_dataframe = create_dataframe()
    CRIPSR_dataframe.to_csv(args.o, header=True, index=False, mode='w')

    for chromosome,sequence in fasta_file.items():
        CRIPSR_dataframe = create_dataframe()
        Complete_dataset = []

        if args.cas9:
            # + strand
            motif = re.compile(r'(?=.GG)')
            cas9_target_list = find_PAM_site(motif,sequence)
            for target in cas9_target_list:
                pam_location = (target[0]-(args.l+1),target[0]-1)
                if pam_location[0] >= 5 and pam_location[0]+5 <= len(sequence)+10 and pam_location[1] >= 5 and pam_location[1] <= len(sequence)+10:
                    shortseq = get_gRNA_sequence(sequence[pam_location[0]:pam_location[1]])
                    longseq = get_gRNA_sequence(sequence[pam_location[0]-5:pam_location[1]+5])
                    crispr_guide = [pam_location[0],pam_location[1],chromosome[1::],shortseq,longseq,'cas9','+']
                    Complete_dataset.append(crispr_guide)

            # - strand
            motif = re.compile(r'(?=CC.)')
            cas9_target_list2 = find_PAM_site(motif,sequence)
            for target in cas9_target_list2:
                pam_location = (target[0]+1,target[0]+(args.l+1))
                if pam_location[0] >= 5 and pam_location[0]+5 <= len(sequence)+10 and pam_location[1] >= 5 and pam_location[1] <= len(sequence)+10:
                    shortseq = get_gRNA_sequence(get_reverse_complement(sequence[pam_location[0]:pam_location[1]]))
                    longseq = get_gRNA_sequence(get_reverse_complement(sequence[pam_location[0]-5:pam_location[1]+5]))
                    crispr_guide = [pam_location[1],pam_location[0],chromosome[1::],shortseq,longseq,'cas9','-']
                    Complete_dataset.append(crispr_guide)

            if args.verbose:
                print (f'''
                {len(cas9_target_list + cas9_target_list2):n} Cas9 PAM sites were found on {chromosome[1::]}
                ''')

        if args.cpf1:
        # + strand
            motif = re.compile(r'(?=TTT.)')
            cpf1_target_list = find_PAM_site(motif,sequence)
            for target in cpf1_target_list:
                pam_location = (target[0]+4,target[0]+(args.l+4))
                if pam_location[0] >= 0 and pam_location[0] <= len(sequence) and pam_location[1] >= 0 and pam_location[1] <= len(sequence):
                    shortseq = get_gRNA_sequence(sequence[pam_location[0]:pam_location[1]])
                    longseq = get_gRNA_sequence(sequence[pam_location[0]-5:pam_location[1]+5])
                    crispr_guide = pd.Series((pam_location[0],pam_location[1],chromosome[1::],shortseq,longseq,'cpf1','+'))
                    Complete_dataset.append(crispr_guide)

        # - strand
            motif = re.compile(r'(?=.AAA)')
            cpf1_target_list2 = find_PAM_site(motif,sequence)
            for target in cpf1_target_list2:
                pam_location = (target[0]-(args.l+1),target[0]-1)
                if pam_location[0] >= 0 and pam_location[0] <= len(sequence) and pam_location[1] >= 0 and pam_location[1] <= len(sequence):
                    shortseq = get_gRNA_sequence(sequence[pam_location[0]:pam_location[1]])
                    longseq = get_gRNA_sequence(sequence[pam_location[0]-5:pam_location[1]+5])
                    crispr_guide = pd.Series((pam_location[1],pam_location[0],chromosome[1::],shortseq,longseq,'cpf1','-'))
                    Complete_dataset.append(crispr_guide)

            if args.verbose:
                print (f'''
                {len(cpf1_target_list + cpf1_target_list2):n} Cpf1 PAM sites were found on {chromosome[1::]}
                ''')

        lesser_list = []
        # print(Complete_dataset)
        for i,item in enumerate(Complete_dataset):
            lesser_list.append((item))
            if len(lesser_list) == 5000 or i == len(Complete_dataset)-1:
                CRIPSR_dataframe['raw'] = lesser_list
                CRIPSR_dataframe.apply(preprocess_PAM_sites,axis=1)
                CRIPSR_dataframe.apply(fill_row,axis=1)
                CRIPSR_dataframe.to_csv(args.o, header=False, index=False, mode='a+')
                CRIPSR_dataframe = create_dataframe()
                lesser_list = []

                
    ### CONFIRMATION MESSAGE
    if args.verbose:
        print(f'The output file has been generated at {args.o}')


if __name__ == '__main__':
    main()
