#!/usr/bin/env python3
# Written by: Hans Müller Paul
#                           NOTES:


### Importing required libraries
import argparse
import datetime
import gffutils
import json
import math
import numpy as np
import pandas as pd
import re
import subprocess
import sys
from os import path,remove

### Defining the arguments
parser = argparse.ArgumentParser()

parser.add_argument('-f', '--fasta', required=True, dest='f', 
                    help='[required] path to input file in FASTA format'
                    )
# parser.add_argument('-g', '--gff', required=True, dest='g', 
#                     help='[required] path to input file in GFF format'
#                     )
parser.add_argument('-o', '--output', dest='o', default='output.txt',
                    help='path to output file'
                    )
parser.add_argument('-l', '--length', metavar='', dest='l', type=int, default=23,
                    help='length of the gRNA sequence, default = 23'
                    )
parser.add_argument('-L', '--flanking', metavar='', dest='L', type=int, default=500,
                    help='length of flanking region for verification, default = 500'
                    )
parser.add_argument('--cas9', action='store_true',
                    help='specifies that design will be made for the Cas9 CRISPR system'
                    )
parser.add_argument('--cpf1', action='store_true',
                    help='specifies that design will be made for the Cpf1 CRISPR system'
                    )
# parser.add_argument('--full', action='store_true',
#                     help='generates an output file containing all sgRNAs. WARNING: OUTPUT FILE MAY BE LARGE'
#                     )
parser.add_argument('-v', '--verbose', action='store_true',
                    help='prints each step of each iteration (for debugging)'
                    )

args = parser.parse_args()


### Defining class functions
class PAM:
    def __init__(self, location, chromosome, guide_type, mainstrand=True):
        if mainstrand:
            strand = '+'
        else:
            strand = '-'
        self.strand = strand
        self.location = location
        self.chr = chromosome[1::]
        self.type = guide_type
        self.id = self.__get_id()
    
    def __get_id(self):
        return f'[{self.type}]{self.chr}.({self.location[0]}:{self.location[1]}).strand({self.strand})'


class CRISPR:
    def __init__(self, sequence, full_sequence, location, cutsite, chromosome, guide_type, strand, off_site_score):
        self.strand = strand
        self.sequence = sequence
        self.long = full_sequence
        self.location = location
        self.cutsite = cutsite
        self.chr = chromosome
        self.type = guide_type
        self.id = self.__get_id()
        self.ds = self.__doench_score(sequence)
        # self.flank = self.__flanking_region_for_sequencing()
        self.off = off_site_score

    # def __flanking_region_for_sequencing(self):
    #     offset = args.L
    #     start = self-offset
    #     end = self+len(self.sgrna)+offset
    #     if start <=0:
    #         start == 0
    #     if end >= len(self.sequence):
    #         end == len(self.sequence)
    #     genomic_region = self.sequence[start:end]
    #     return genomic_region

    def __get_id(self):
        return f'[{self.type}]{self.chr}.({self.location[0]}:{self.location[1]})'

    def __doench_score(self,full_sequence):
        '''
        Generates a binary matrix for DNA/RNA sequence, where each column is a possible base
        and each row is a position along the sequence. Matrix column order is A, T/U, C, G
        '''
        seq = str(full_sequence).upper()
        seq = list(seq)
        matrix1  = np.zeros([len(full_sequence),4], dtype=int)
        for i,item in enumerate(full_sequence):
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
        full_sequence = full_sequence.replace('U','T')
        pairwise_sequence = []
        for i in range(len(full_sequence)):
            if i < len(full_sequence)-1:
                basepair = full_sequence[i]+full_sequence[i+1]
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

        item_gc = full_sequence[0][5:-5]
        gc_count = item_gc.count('G') + item_gc.count('C')
        if gc_count < 10:
            gc_score = low_gc
        else:
            gc_score = high_gc
        score_first = np.sum(matrix1.dot(first_matrix))
        score_second = np.sum(matrix2.dot(second_matrix))
        score_sum = score_first + score_second + gc_score
        '''
        Logistic regression module to calculate the score for a particular guide sequence
        '''
        doench_score = math.pow(1 - math.exp(-(intersect + score_sum)),-1)
        return doench_score


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
            from functions import formatted
            myFASTA = formatted(myFASTA)
            if args.verbose:
                print(f'Genome file {fasta} successfully formatted')
    from functions import generate_dictionary as gendict
    genome_dictionary = gendict(myFASTA)
    if args.verbose:
        print(f'The genome was successfully converted to a dictionary')
    return genome_dictionary


def import_gff_file(gff):
    '''
    imports and formats a genome annotation file in gff3 format for feature identification
    '''
    if path.exists(str(f'{gff[::-5]}.db')):
        featuredb = gffutils.FeatureDB(str(f'{gff[::-5]}.db'))
        if args.verbose:
            print(f'Feature database successfully imported for {gff}')
        return featuredb
    else:
        featuredb = gffutils.create_db(gff,
                                        str(f'{gff[::-5]}.db'),
                                        merge_strategy="create_unique",
                                        keep_order=True)
        featuredb = gffutils.helpers.sanitize_gff_db(featuredb)
        gffutils.FeatureDB(featuredb)
        if args.verbose:
            print(f'Feature database successfully generated from {gff}')
        return featuredb


def find_PAM_site(target,input_sequence):
    '''
    locates the target PAM motif in input sequence
    '''
    PAM_site = [match.span() for match in re.finditer(target,input_sequence)]
    return PAM_site


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


def off_site_score(sequence, id, genome):
    off_site = 'b'
    return off_site


def create_dataframe():
    '''
    creates a dataframe to store sgRNA information
    '''
    df_cols = [
                'crispr_id',        # STR
                'crispr_sys',       # CAT
                'sequence',         # STR
                'chromosome',       # STR (CAT?)
                'start_pos',        # INT
                'end_pos',          # INT
                'cutsite',          # INT
                'strand',           # CAT
                'on_site_score',    # FLOAT
                'off_site_score'    # FLOAT
                ]
    df = pd.DataFrame(columns=df_cols)
    df.set_index('crispr_id', inplace=True)
    '''
    implement memory optimization by assigning appropriate dtype
    '''
    return df


def add_to_dataframe(dataframe, guide_id, system, sequence, chromosome, start_location,
                    end_location, cutsite, strand, on_site_score, off_site_score):
    '''
    adds formatted data to fill the rows of the dataframe containing
    sgRNA information for each potential target site
    '''
    data_line = pd.Series([guide_id,
                        system,
                        sequence,
                        chromosome,
                        start_location,
                        end_location,
                        cutsite,
                        strand,
                        on_site_score,
                        off_site_score]) #,
                        # dtype=['object',
                        # 'category',
                        # 'object',
                        # 'category',
                        # 'int64',
                        # 'int64',
                        # 'int64',
                        # 'category',
                        # 'float64',
                        # 'float64'])
    return dataframe.append(data_line, verify_integrity=True, ignore_index=True)


# def get_feature_by_location(sgrna,feature_db):
#     feat = feature_db.region(f'{sgrna.chr}:{sgrna.location[0]}-{sgrna.location[1]}')
#     for item in feat:
#         feat_list=(item.feature.seqid,
#         item.feature.featuretype,
#         item.feature.frame)
#         return feat_list


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


def main():
    if not args.cas9 and not args.cpf1:
        sys.exit('Please select at least one CRISPR system: Cas9 or Cpf1')
    if args.verbose:
        print(args)


    ### Import genome files
    fasta_file = import_fasta_file(args.f)
    # gff_file = import_gff_file(args.g)
    # gff_db = gffutils.FeatureDB(gff_file)


    ### Locate PAMs by nuclease type
    if args.verbose:
        print(f'''
            Initiating PAM site detection.
            
            Please wait, this may take a while...
            ''')


    cas9_count = 0
    cpf1_count = 0


    CRISPR_dataframe = create_dataframe()

    if args.cas9:
        invalid_cas9_targets = []
        for k,v in fasta_file.items():
            # + strand
            motif = re.compile(r'(?=.GG)')
            cas9_target_list = find_PAM_site(motif,v)
            for target in cas9_target_list:
                ''' designing guide from PAM site '''
                pam_location = (target[0]+1,target[0]+3)
                pam_site = PAM(pam_location,k,'cas9')
                sgrna_position = (pam_site.location[0]-(args.l+1),pam_site.location[0]-1)
                if sgrna_position[0] >= 0 and sgrna_position[0] <= len(v) and sgrna_position[1] >= 0 and sgrna_position[1] <= len(v):
                    sgrna_sequence = v[sgrna_position[0]:sgrna_position[1]]
                    sgrna_sequence = get_reverse_complement(sgrna_sequence)
                    sgrna_full_sequence = get_gRNA_sequence(v[sgrna_position[0]-5:sgrna_position[1]+5])
                    sgrna_sequence = get_gRNA_sequence(sgrna_sequence)
                    sgrna_oss = off_site_score(sgrna_sequence,pam_site.id,args.f)
                    sgrna = CRISPR(sgrna_sequence,sgrna_full_sequence,sgrna_position,sgrna_position[1]-3,pam_site.chr,pam_site.type,pam_site.strand, sgrna_oss)
                    cas9_count = cas9_count + 1
                    # feat = gff_file.region(seqid=sgrna.chr, start=sgrna.location[0], end=sgrna.location[1])
                    # feat = get_feature_by_location(sgrna,gff_db)
                    # print(f'''{sgrna.id}:{sgrna.sequence}
                    
                    # {feat}
                    
                    # ''')
                    ''' inclusion of guides in dataframe '''
                    add_to_dataframe(CRISPR_dataframe,
                                    sgrna.id,
                                    sgrna.type,
                                    sgrna.sequence,
                                    sgrna.chr,
                                    sgrna.location[0],
                                    sgrna.location[1],
                                    sgrna.cutsite,
                                    sgrna.strand,
                                    sgrna.ds,
                                    sgrna.off)
                else:
                    invalid_cas9_targets.append(pam_site)
            # - strand
            motif = re.compile(r'(?=CC.)')
            cas9_target_list2 = find_PAM_site(motif,v)
            for target in cas9_target_list2:
                pam_location = (target[0]+1,target[0]+3)
                pam_site = PAM(pam_location[::-1],k,'cas9',False)
                sgrna_position = (pam_site.location[0]+(args.l+1),pam_site.location[0]+1)
                if sgrna_position[0] >= 0 and sgrna_position[0] <= len(v) and sgrna_position[1] >= 0 and sgrna_position[1] <= len(v):
                    sgrna_sequence = v[sgrna_position[1]:sgrna_position[0]]
                    sgrna_sequence = get_gRNA_sequence(sgrna_sequence)
                    sgrna_full_sequence = get_gRNA_sequence(v[sgrna_position[1]-5:sgrna_position[0]+5])
                    sgrna_oss = off_site_score(sgrna_sequence,pam_site.id,args.f)
                    sgrna = CRISPR(sgrna_sequence,sgrna_full_sequence,sgrna_position,sgrna_position[1]-3,pam_site.chr,pam_site.type,pam_site.strand, sgrna_oss)
                    cas9_count = cas9_count + 1
                    print(sgrna.id,':',sgrna.sequence)
                    add_to_dataframe(CRISPR_dataframe,
                                    sgrna.id,
                                    sgrna.type,
                                    sgrna.sequence,
                                    sgrna.chr,
                                    sgrna.location[0],
                                    sgrna.location[1],
                                    sgrna.cutsite,
                                    sgrna.strand,
                                    sgrna.ds,
                                    sgrna.off)
                else:
                    invalid_cas9_targets.append(pam_site)
            if args.verbose:
                print (f'''
                {len(cas9_target_list + cas9_target_list2)} Cas9 PAM sites were found on {k[1::]}
                ''')

    if args.cpf1:
        invalid_cpf1_targets = []
        for k,v in fasta_file.items():
        # + strand
            motif = re.compile(r'(?=TTT.)')
            cpf1_target_list = find_PAM_site(motif,v)
            for target in cpf1_target_list:
                pam_location = (target[0]+1,target[0]+4)
                pam_site = PAM(pam_location,k,'cpf1')
                del pam_location
                sgrna_position = (pam_site.location[1]+1,pam_site.location[1]+(args.l+1))
                if sgrna_position[0] >= 0 and sgrna_position[0] <= len(v) and sgrna_position[1] >= 0 and sgrna_position[1] <= len(v):
                    sgrna_sequence = v[sgrna_position[0]:sgrna_position[1]]
                    sgrna_sequence = get_gRNA_sequence(sgrna_sequence)
                    sgrna_full_sequence = get_gRNA_sequence(v[sgrna_position[0]-5:sgrna_position[1]+5])
                    sgrna_oss = off_site_score(sgrna_sequence,pam_site.id,args.f)
                    sgrna = CRISPR(sgrna_sequence,sgrna_full_sequence,sgrna_position,sgrna_position[1]-3,pam_site.chr,pam_site.type,pam_site.strand,sgrna_oss)
                    cpf1_count = cpf1_count + 1
                    print(sgrna.id,':',sgrna.sequence)
                    add_to_dataframe(CRISPR_dataframe,
                                    sgrna.id,
                                    sgrna.type,
                                    sgrna.sequence,
                                    sgrna.chr,
                                    sgrna.location[0],
                                    sgrna.location[1],
                                    sgrna.cutsite,
                                    sgrna.strand,
                                    sgrna.ds,
                                    sgrna.off)
                else:
                    invalid_cpf1_targets.append(pam_site)

        # - strand
            motif = re.compile(r'(?=.AAA)')
            cpf1_target_list2 = find_PAM_site(motif,v)
            for target in cpf1_target_list2:
                pam_location = (target[0]+1,target[0]+4)
                pam_site = PAM(pam_location[::-1],k,'cpf1',False)
                del pam_location
                sgrna_position = (pam_site.location[1]-1,pam_site.location[1]-(args.l+1))
                if sgrna_position[0] >= 0 and sgrna_position[0] <= len(v) and sgrna_position[1] >= 0 and sgrna_position[1] <= len(v):
                    sgrna_sequence = v[sgrna_position[1]:sgrna_position[0]]
                    sgrna_sequence = get_gRNA_sequence(sgrna_sequence)
                    sgrna_full_sequence = get_gRNA_sequence(v[sgrna_position[1]-5:sgrna_position[0]+5])
                    sgrna_oss = off_site_score(sgrna_sequence,pam_site.id,args.f)
                    sgrna = CRISPR(sgrna_sequence,sgrna_full_sequence,sgrna_position,sgrna_position[1]-3,pam_site.chr,pam_site.type,pam_site.strand,sgrna_oss)
                    cpf1_count = cpf1_count + 1
                    print(sgrna.id,':',sgrna.sequence)
                    add_to_dataframe(CRISPR_dataframe,
                                    sgrna.id,
                                    sgrna.type,
                                    sgrna.sequence,
                                    sgrna.chr,
                                    sgrna.location[0],
                                    sgrna.location[1],
                                    sgrna.cutsite,
                                    sgrna.strand,
                                    sgrna.ds,
                                    sgrna.off)
                else:
                    invalid_cpf1_targets.append(pam_site)
            if args.verbose:
                print (f'''
                {len(cpf1_target_list + cpf1_target_list2)} Cpf1 PAM sites were found on {k[1::]}
                ''')

    ### Gather info on each target motif

    if args.cas9:
        print(f'a total of {cas9_count} gRNAs was designed for Cas9')
    if args.cpf1:
        print(f'a total of {cpf1_count} gRNAs was designed for Cpf1')


    """
    ### WRITE TO OUTPUT FILE (add to head)

    output_dict = {
    "Target genome": 'Sorghum Bicolor',
    "Total Cas9 gRNAs": len(cas9_guides),
    # "Total Cpf1 gRNAs": len(cpf1_guides),
    "Generated with": "$software name",
    "Date": str(datetime.datetime.now()),
    "Cas9 sgRNA Info": Cas9_guide_info,
    # "Cpf1 sgRNA Info": Cpf1_guide_info
    }
    with open(args.o, 'w') as f:
    json.dump(output_dict, f, indent=2)
    """

### CONFIRMATION MESSAGE
    print(f'The output file has been generated at {args.o}')

if __name__ == '__main__':
    main()