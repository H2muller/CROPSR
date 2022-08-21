#!/usr/bin/env python3


def one_base_matrix(sequence):
    """
    Generates a binary matrix for DNA/RNA sequence, where each column is a possible base
    and each row is a position along the sequence. Matrix column order is A, T/U, C, G
    """
    # Import libraries
    from numpy import zeros
    # Function
    seq = str(sequence).upper()
    seq = list(seq)
    matrix = zeros([len(sequence),4], dtype=int)
    baseIndex = {
        'A': 0,
        'T': 1,
        'U': 1,
        'C': 2,
        'G': 3
    }
    for i,item in enumerate(sequence):
        index = baseIndex[item.upper()]
        matrix[i,index] = 1
    return matrix


def pairwise_matrix(sequence):
    """
    Generates a binary matrix for DNA/RNA sequence, where each column is a possible
    pair of adjacent bases, and each row is a position along the sequence.
    Matrix column order is AA, AT, AC, AG, TA, TT, TC, TG, CA, CT, CC, CG, GA, GT, GC, GG
    """
    # Import libraries
    from numpy import zeros
    # Function
    sequence = sequence.replace('U','T')
    pairwise_sequence = []
    for i in range(len(sequence)):
        if i < len(sequence)-1:
            basepair = sequence[i]+sequence[i+1]
            pairwise_sequence.append(basepair)
    matrix = zeros([len(pairwise_sequence),16], dtype=int)
    pairIndex = {
        'AA': 0,
        'AT': 1,
        'AC': 2,
        'AG': 3,
        'TA': 4,
        'TT': 5,
        'TC': 6,
        'TG': 7,
        'CA': 8,
        'CT': 9,
        'CC': 10,
        'CG': 11,
        'GA': 12,
        'GT': 13,
        'GC': 14,
        'GG': 15,
    }
    for i,item in enumerate(pairwise_sequence):
        index = pairIndex[item.upper()]
        matrix[i,index] = 1
    return matrix



def rs1_score(sequence):
    """
    Generates a binary matrix for DNA/RNA sequence, where each column is a possible base
    and each row is a position along the sequence. Matrix column order is A, T/U, C, G
    """
    # Import Libraries
    import math
    from numpy import zeros, sum
    # Function
    """
    Scoring algorithm
    """
    intersect = 0.59763615
    low_gc = -0.2026259
    high_gc = -0.1665878

    """
    Weight matrixes derived from Doench 2014
    """
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
    baseIndex = {
        'A': 0,
        'T': 1,
        'U': 1,
        'C': 2,
        'G': 3
    }
    for k,v in first_order_scores.items():
        index = baseIndex[k[0].upper()]
        first_matrix[index, int(k[1:])-1] = v

    # order 2 score matrix
    """ row order == AA AT AC AG TA TT TC TG CA CT CC CG GA GT GC GG """
    second_matrix = zeros([16,29], dtype=int)
    pairIndex = {
        'AA': 0,
        'AT': 1,
        'AC': 2,
        'AG': 3,
        'TA': 4,
        'TT': 5,
        'TC': 6,
        'TG': 7,
        'CA': 8,
        'CT': 9,
        'CC': 10,
        'CG': 11,
        'GA': 12,
        'GT': 13,
        'GC': 14,
        'GG': 15,
    }
    for k,v in second_order_scores.items():
        index = pairIndex[k[0:2].upper()]
        second_matrix[index, int(k[2:])-1] = v

    item_gc = sequence[0][5:-5]
    gc_count = item_gc.count('G') + item_gc.count('C')
    if gc_count < 10:
        gc_score = low_gc
    else:
        gc_score = high_gc
    matrix1 = one_base_matrix(sequence)
    score_first = sum(matrix1.dot(first_matrix))
    matrix2 = pairwise_matrix(sequence)
    score_second = sum(matrix2.dot(second_matrix))
    score_sum = score_first + score_second + gc_score

    score = math.pow((1 - math.exp(-(intersect + score_sum))),-1)
    return score


def generate_dictionary(input):
    """

    """
    # import libraries
    import itertools
    # function
    dictionary = input.split()
    dictionary = dict(itertools.zip_longest(*[iter(dictionary)] * 2, fillvalue=""))
    return dictionary


def location(primer, genome):
    """
    Written by: Hans Müller Paul and Zhiwen Jiang
    """
    a = True
    list_of_beginning = []
    list_of_end = []
    start = 0
    primer_location = []
    while a:
        beginning = genome.find(primer, start)+1
        if beginning + len(primer)-1 >= len(genome) or beginning == 0:
            a = False
        else:
            end = beginning + len(primer)-1
            list_of_beginning.append(beginning)
            list_of_end.append(end)
            start = beginning
        primer_location = list(zip(list_of_beginning, list_of_end))
    return primer_location


def formatted(input_genome):
    """
    Written by: Hans Müller Paul and Joao Paulo Gomes Viana
    """
    # import libraries
    import re
    # function
    formatted = re.sub('\n','',input_genome)
    formatted = re.sub('>', '\n>',formatted)
    formatted = formatted[1:]
    formatted = re.sub('([0-9]+)','\\1 \n',formatted)
    return formatted


def printProgressBar (iteration, total, prefix = '', suffix = '', decimals = 1, length = 100, fill = '█'):
    """
    Call in a loop to create terminal progress bar
    Args:
        iteration   - Required  : current iteration (Int)
        total       - Required  : total iterations (Int)
        prefix      - Optional  : prefix string (Str)
        suffix      - Optional  : suffix string (Str)
        decimals    - Optional  : positive number of decimals in percent complete (Int)
        length      - Optional  : character length of bar (Int)
        fill        - Optional  : bar fill character (Str)
    """
    percent = ("{0:." + str(decimals) + "f}").format(100 * (iteration / float(total)))
    filledLength = int(length * iteration // total)
    bar = fill * filledLength + '-' * (length - filledLength)
    print('\r%s |%s| %s%% %s' % (prefix, bar, percent, suffix), end = '\r')
    # Print New Line on Complete
    if iteration == total: 
        print()
    """ Paste the following section within loop """
    # """ Update Progress Bar """
    # printProgressBar(i + 1, l, prefix = 'Analyzing guide sequences:', suffix = 'Complete', length = 50)


def parallelize(data, func):
    """
    """
    # Import Libraries
    import multiprocessing as mp
    import numpy as np
    # Function
    if mp.cpu_count() > 2:
        cores = mp.cpu_count()-1    # Runs in all cores except for one
    else:
        cores = 1                   # Runs in a single core
    
    data_split = np.array_split(data, cores)
    pool = mp.Pool(cores)
    pool.map(func,data_split)
    pool.close()
    pool.join()
    return pool


def create_dataframe():
    """
    creates a dataframe to store information
    """
    # Import Libraries
    import pandas as pd
    # Function
    df_cols = [
                'sequence',         # STR
                'on_site_score'    # FLOAT
                ]
    df = pd.DataFrame(columns=df_cols)
    """
    implement memory optimization by assigning appropriate dtype
    """
    return df


def save_dataframe_to_tmp(data):
    """ save one h_matrix and one permutation in temorary files with sequence_number appended names.
    Args:
        data: list that will have a function applied to it
        unique_id: temporary file name suffix. 
    """
    # Import Libraries
    import os
    import pandas as pd
    import numpy as np
    # Function
    tmp_dir = os.getcwd().join("/tmp_directory")
    os.makedirs(tmp_dir, mode=0o755, exist_ok=True)
    tmp_file_name = os.path.join(tmp_dir, f'tmp_{os.getpid()}')
    dataframe = create_dataframe()
    for i,item in enumerate(data):
        score = rs1_score(item)
        dataline = pd.Series([item,score],index=dataframe.columns)
        # print(dataline)
        dataframe = dataframe.append(dataline,ignore_index=True)
    with open(tmp_file_name, 'wb') as temp_path:
        return dataframe.to_csv(temp_path)
