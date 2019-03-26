#!/usr/bin/env python3
# Written by: Hans Müller Paul
#                           NOTES:

import random
import sys
import numpy as np
import pandas as pd
import math


def sequence_generator(number,length=20,seqtype='DNA'):
    list_of_sequences = []
    if number > 0:
        for i in range(number):
            if seqtype is not 'DNA':
                seqtype = 'RNA'
                list_of_sequences.append(''.join(random.choices(['A','U','C','G'], k=length)))
            elif seqtype is 'DNA':
                list_of_sequences.append(''.join(random.choices(['A','T','C','G'], k=length)))
    else:
        sys.exit('please provide a non-zero number of sequences to be generated')
    return list_of_sequences


def one_base_matrix(sequence):
    '''
    Generates a binary matrix for DNA/RNA sequence, where each column is a possible base
    and each row is a position along the sequence. Matrix column order is A, T/U, C, G
    '''
    seq = str(sequence).upper()
    seq = list(seq)
    matrix = np.zeros([len(sequence),4], dtype=int)
    for i,item in enumerate(sequence):
        if item == 'A':
            matrix[i,0] = 1
        if item == 'T':
            matrix[i,1] = 1
        if item == 'U':
            matrix[i,1] = 1
        if item == 'C':
            matrix[i,2] = 1
        if item == 'G':
            matrix[i,3] = 1
    return matrix


def pairwise_matrix(sequence):
    '''
    Generates a binary matrix for DNA/RNA sequence, where each column is a possible
    pair of adjacent bases, and each row is a position along the sequence.
    Matrix column order is AA, AT, AC, AG, TA, TT, TC, TG, CA, CT, CC, CG, GA, GT, GC, GG
    '''
    sequence = sequence.replace('U','T')
    pairwise_sequence = []
    for i in range(len(sequence)):
        if i < len(sequence)-1:
            basepair = sequence[i]+sequence[i+1]
            pairwise_sequence.append(basepair)
    matrix = np.zeros([len(pairwise_sequence),16], dtype=int)
    for i,item in enumerate(pairwise_sequence):
        if item == 'AA':
            matrix[i,0] = 1
        if item == 'AT':
            matrix[i,1] = 1
        if item == 'AC':
            matrix[i,2] = 1
        if item == 'AG':
            matrix[i,3] = 1
        if item == 'TA':
            matrix[i,4] = 1
        if item == 'TT':
            matrix[i,5] = 1
        if item == 'TC':
            matrix[i,6] = 1
        if item == 'TG':
            matrix[i,7] = 1
        if item == 'CA':
            matrix[i,8] = 1
        if item == 'CT':
            matrix[i,9] = 1
        if item == 'CC':
            matrix[i,10] = 1
        if item == 'CG':
            matrix[i,11] = 1
        if item == 'GA':
            matrix[i,12] = 1
        if item == 'GT':
            matrix[i,13] = 1
        if item == 'GC':
            matrix[i,14] = 1
        if item == 'GG':
            matrix[i,15] = 1
    return matrix


def printProgressBar (iteration, total, prefix = '', suffix = '', decimals = 1, length = 100, fill = '█'):
    '''
    Call in a loop to create terminal progress bar
    @params:
        iteration   - Required  : current iteration (Int)
        total       - Required  : total iterations (Int)
        prefix      - Optional  : prefix string (Str)
        suffix      - Optional  : suffix string (Str)
        decimals    - Optional  : positive number of decimals in percent complete (Int)
        length      - Optional  : character length of bar (Int)
        fill        - Optional  : bar fill character (Str)
    '''
    percent = ("{0:." + str(decimals) + "f}").format(100 * (iteration / float(total)))
    filledLength = int(length * iteration // total)
    bar = fill * filledLength + '-' * (length - filledLength)
    print('\r%s |%s| %s%% %s' % (prefix, bar, percent, suffix), end = '\r')
    # Print New Line on Complete
    if iteration == total: 
        print()


def rs1_score(sequence):

    '''
    Generates a binary matrix for DNA/RNA sequence, where each column is a possible base
    and each row is a position along the sequence. Matrix column order is A, T/U, C, G
    '''
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

    '''
    Generates a binary matrix for DNA/RNA sequence, where each column is a possible
    pair of adjacent bases, and each row is a position along the sequence.
    Matrix column order is AA, AT, AC, AG, TA, TT, TC, TG, CA, CT, CC, CG, GA, GT, GC, GG
    '''
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

    '''
    Scoring algorithm
    '''
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
    ''' row order == A T/U C G '''
    first_matrix = np.zeros([4,30], dtype=int)
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
    ''' row order == AA AT AC AG TA TT TC TG CA CT CC CG GA GT GC GG '''
    second_matrix = np.zeros([16,29], dtype=int)
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
    score_first = np.sum(matrix1.dot(first_matrix))
    score_second = np.sum(matrix2.dot(second_matrix))
    score_sum = score_first + score_second + gc_score

    score = math.pow((1 - math.exp(-(intersect + score_sum))),-1)
    return score


def main():
    print('Generating random sequences, please wait')
    lista = sequence_generator(10000,30)
    l = len(lista)
    nova_lista = []
    ''' Initial call to print 0% progress '''
    printProgressBar(0, l, prefix = 'Analyzing guide sequences:', suffix = 'Complete', length = 50)
    for i,item in enumerate(lista):
        nova_lista.append([item,rs1_score(item)])
        # nova_lista.append([item, one_base_matrix(item), pairwise_matrix(item)])
        ''' Update Progress Bar '''
        printProgressBar(i + 1, l, prefix = 'Analyzing guide sequences:', suffix = 'Complete', length = 50)
    df = pd.DataFrame(nova_lista, columns=['sequence', 'score'])
    # df = pd.DataFrame(nova_lista, columns=['sequence','base matrix','pairwise matrix', 'score'])
    print(df)


if __name__ == '__main__':
    main()