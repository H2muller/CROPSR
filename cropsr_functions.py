#!/usr/bin/env python3


def generate_dictionary(input):
    '''
    Written by: Hans Müller Paul
    '''
    import itertools
    dictionary = input.split()
    dictionary = dict(itertools.zip_longest(*[iter(dictionary)] * 2, fillvalue=""))
    return dictionary


def formatted(input_genome):
    '''
    Written by: Hans Müller Paul and Joao Paulo Gomes Viana
    '''
    import re
    formatted = re.sub('\n','',input_genome)
    formatted = re.sub('>', '\n>',formatted)
    formatted = formatted[1:]
    formatted = re.sub('([0-9]+)','\\1 \n',formatted)
    return formatted


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
            myFASTA = formatted(myFASTA)
            if args.verbose:
                print(f'Genome file {fasta} successfully formatted')
    genome_dictionary = generate_dictionary(myFASTA)
    if args.verbose:
        print(f'The genome was successfully converted to a dictionary')
    return genome_dictionary


def import_gff_file(gff):
    from pandas import read_csv
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
    gff_df = read_csv(gff, sep='\t', skiprows = start_index, header = None, names = col_names)
    if args.verbose:
        print(f'Annotation database successfully generated')
    return gff_df


def find_PAM_site(target,input_sequence):
    '''
    locates the target PAM motif in input sequence
    '''
    from re import finditer
    PAM_site = [match.span() for match in finditer(target,input_sequence)]
    return PAM_site


def create_dataframe():
    '''
    creates a dataframe to store sgRNA information
    '''
    from pandas import DataFrame
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
    df = DataFrame(columns=df_cols)
    return df


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


def apply_cutsite(start_pos, end_pos, crispr_sys):
    if crispr_sys == 'cas9':
        cutsite = end_pos-3
    elif crispr_sys == 'cpf1':
        cutsite = end_pos-5 # this number is currently a placeholder
    return cutsite


def rs1_score(sequence):
    """
    Generates a binary matrix for DNA/RNA sequence, where each column is a possible base
    and each row is a position along the sequence. Matrix column order is A, T/U, C, G
    """
    import math
    from numpy import zeros, matmul, trace
    seq = str(sequence).upper()
    seq = list(seq)
    matrix1  = zeros([len(sequence),4], dtype=int)
    for i, item in enumerate(sequence):
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
    for i in range(len(sequence)-1):
        basepair = sequence[i]+sequence[i+1]
        pairwise_sequence.append(basepair)
    matrix2 = zeros([len(pairwise_sequence),16], dtype=int)
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
    Scoring algorithm
    """
    intersect = 0.59763615
    low_gc = -0.2026259
    high_gc = -0.1665878

    first_order = ['G02','A03','C03','C04','C05','G05','A06','C06','C07','G07','A12','A15','C15','A16',
                    'C16','T16','A17','G17','C18','G18','A19','C19','G20','T20','G21','T21','C22','T22','T23',
                    'C24','G24','T24','A25','C25','T25','G28','T28','C29','G30']
    first_scores = [-0.2753771,-0.3238875,0.17212887,-0.1006662,-0.2018029,
                    0.24595663,0.03644004,0.09837684,-0.7411813,-0.3932644,
                    -0.466099,0.08537695,-0.013814,0.27262051,0.1190226,
                    -0.2859442,0.09745459,-0.1755462,-0.3457955,-0.6780964,
                    0.22508903,-0.5077941,-0.4173736,-0.054307,0.37989937,
                    -0.0907126,0.05782332,-0.5305673,-0.8770074,-0.8762358,
                    0.27891626,-0.4031022,-0.0773007,0.28793562,-0.2216372,
                    -0.6890167,0.11787758,-0.1604453,0.38634258]
    first_order_scores = dict(zip(first_order,first_scores))

    second_order = ['GT02','GC05','AA06','TA06','GG07','GG12','TA12','TC12','TT12','GG13','GA14','GC14',
                    'TG17','GG19','TC19','CC20','TG20','AC21','CG21','GA21','GG21','TC22','CG23','CT23',
                    'AA24','AG24','AG25','CG25','TG25','GT27','GG29']
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
    first_matrix = zeros([4,30], dtype=int)
    for k,v in first_order_scores.items():
        if k[0] == 'A':
            first_matrix[0,int(k[1:])-1] = v
        if k[0] == 'T':
            first_matrix[1,int(k[1:])-1] = v
        if k[0] == 'C':
            first_matrix[2,int(k[1:])-1] = v
        if k[0] == 'G':
            first_matrix[3,int(k[1:])-1] = v

    # order 2 score matrix
    """ row order == AA AT AC AG TA TT TC TG CA CT CC CG GA GT GC GG """
    second_matrix = zeros([16,29], dtype=int)
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
    score_first = trace(matmul(first_matrix,matrix1))
    score_second = trace(matmul(second_matrix,matrix2))
    score = (1/(1 + math.exp(-(intersect + gc_score + score_first + score_second))))
    return round(score,5)


# Helper function for retrieve features by position, will be used for every line.
def retrieve_gff_helper(chrom, cutsite, gff_dataframe_csv):
    import pandas as pd
    gff_dataframe = pd.read_csv(gff_dataframe_csv)
    # Filter gff by input_chrom
    filtered_gff_df_1 = gff_dataframe[gff_dataframe['chromosome'] == chrom]

    # Check for each row of gff_dataframe, if cutsite is between start and end, list features
    # If source of GFF is phytozome, look for additional file

    # Use pandas filtering
    start_end_filter = (filtered_gff_df_1['start'] <= cutsite) & (filtered_gff_df_1['end'] <= cutsite)
    filtered_gff_df_2 = filtered_gff_df_1[start_end_filter]

    feature_list = filtered_gff_df_2['attributes'].tolist()

    feature_string = ''.join(feature_list)
    return feature_string


def retrieve_features_by_position(gff_dataframe_csv, DF, vectorized_helper_function):
    '''
    searches a gff dataframe for all the genomic features at a given chromosome
    and position
    '''

    # For each cutsite, if it's between gff_dataframe['start'] and gff_dataframe['end'], add DF['feature']
    feature_array = vectorized_helper_function(DF['chromosome'], DF['cutsite'], gff_dataframe_csv)

    print(feature_array)
    return


def get_id(sys_type, chr):
    'generates a unique ID for each target site'
    import random
    import string
    #  defining first digit
    if sys_type == 'cas9':
        first = 'A'
    elif sys_type == 'cpf1':
        first = 'B'
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
    # pd.to_numeric(DF, errors='coerce')
    return DF


def fill_row(DF):
    from numpy import vectorize
    DF['crispr_id'] = vectorize(get_id)(DF['crispr_sys'],DF['chromosome'])
    DF['cutsite'] = vectorize(apply_cutsite)(DF['start_pos'],DF['end_pos'],DF['crispr_sys'])
    DF['on_site_score'] = vectorize(rs1_score)(DF['sequence'])
    # DF['features'] = vectorize(retrieve_features_by_position)(gff_dataframe,DF['chromosome'],DF['cutsite'])
    # pd.to_numeric(DF, errors='coerce')
    return DF


def optimize_dataframe(DF):
    '''
    Currently not in use
    '''
    import pandas as pd
    # int
    DF['start_pos'] = pd.to_numeric(DF['start_pos'],downcast='unsigned')
    DF['end_pos'] = pd.to_numeric(DF['end_pos'],downcast='unsigned')
    # DF['cutsite'] = pd.to_numeric(DF['cutsite'])
    #float
    DF['on_site_score'] = pd.to_numeric(DF['on_site_score'],downcast='float')
    # category
    DF['crispr_sys'] = DF['crispr_sys'].astype('category')
    DF['chromosome'] = DF['chromosome'].astype('category')
    DF['strand'] = DF['strand'].astype('category')
    return DF


def parallelize(data, func):
    '''
    Currently not in use
    '''
    from multiprocessing import cpu_count, Pool
    if cpu_count() > 2:
        cores = cpu_count()-1    # Runs in all cores except for one
    else:
        cores = 1                   # Runs in a single core
    pool = mp.Pool(cores)
    pool.imap(func,data,chunksize = 500)
    pool.close()
    pool.join()
    return 