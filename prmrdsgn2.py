#!/usr/bin/env python3

#  Written by: Hans Müller Paul and Jacob Heldenbrand

# IMPORT LIBRARIES
import argparse
import datetime
import itertools
import json
import math
import multiprocessing as mp
import pandas as pd
import re
import subprocess
import sys
from os import path,remove


# DEFINING VARIABLES
parser = argparse.ArgumentParser()

parser.add_argument('-i', '--input', dest='i', required=True,
                    help='path to input file in FASTA format'
                    )
parser.add_argument('-g', '--genome', dest='g', required=True, 
                    help='path to input genome file in FASTA format'
                    )
# parser.add_argument('-n', '--number', metavar='', dest='n', type=int, default=5,
#                     help='number of candidate primer pairs to pick, default = 5'
#                     )
parser.add_argument('-e', '--extension', metavar='', dest='e', type=int, default=100,
                    help='number of bases from the start and end of the sequence to look for primers in, default = 100'
                    )
parser.add_argument('-a', '--amplicon', metavar='', dest='a', type=int, default=1000,
                    help='expected length of the amplicon (for alignment), default = 1000'
                    )
parser.add_argument('-s', '--short', metavar='', dest='s', type=int, default=20,
                    help='shortest acceptable primer, default = 20'
                    )
parser.add_argument('-l', '--long', metavar='', dest='l', type=int, default=30,
                    help='longest acceptable primer, default = 30'
                    )
parser.add_argument('-m', '--mintemp', metavar='', dest='m', type=float, default=50,
                    help='min Tm in celsius, default = 50'
                    )
parser.add_argument('-x', '--maxtemp', metavar='', dest='x', type=float, default=65,
                    help='max Tm in celsius, default = 65'
                    )
parser.add_argument('-M', '--mingc', metavar='', dest='M', type=float, default=35, 
                    help='min GC percentage, default = 35'
                    )
parser.add_argument('-X', '--maxgc', metavar='', dest='X', type=float, default=65, 
                    help='max GC percentage, default = 65'
                    )
parser.add_argument('-D', '--tmdiff', metavar='', dest='D', type=float, default=0.5,
                    help='accepted TM difference to form primer pair, default = 0.5'
                    )
parser.add_argument('-o', '--output', metavar='', dest='o', default='output.txt',
                    help='path to output file in FASTA format'
                    )
parser.add_argument('-v', '--verbose', action='store_true',
                    help='prints each step of each iteration (for debugging)'
                    )
parser.add_argument('-F', '--full', action='store_true',
                    help='keeps all temporary and log files'
                    )

args = parser.parse_args()

class Primer:
    """
    Creates an object of the Primer class, containing the sequence and
    relevant metadata, such as GC content and TM
    """
    def __init__(self, sequence):
        self.sequence = sequence
        self.GC_percentage = self.__calculate_GC()
        self.Tm = self.__calculate_Tm()
        self.Anneal = self.__calculate_Annealing()

    def __calculate_GC(self):
        upper_seq = self.sequence.upper()
        gc_count = upper_seq.count('G') + upper_seq.count('C')
        gc_fraction = float(gc_count) / len(self.sequence)
        return 100 * gc_fraction

    def __calculate_Annealing(self):
        N = len(self.sequence)
        anneal_temp = 81.5 + 0.41 * self.GC_percentage - (675 / N)
        return anneal_temp

    def __calculate_Tm(self):
        N = len(self.sequence)
        upper_seq = self.sequence.upper()
        if N < 13:
            melt_temp = ((upper_seq.count('A') + upper_seq.count('T')) * 2 +
            (upper_seq.count('C') + upper_seq.count('G')) * 4)
        else:
            melt_temp = 64.9 + 41 * (upper_seq.count('G') + upper_seq.count('C') - 16.4) / N
        return melt_temp

    def get_primer_elements(self):
        return {"sequence": self.sequence,
                "GC": round(self.GC_percentage, 2),
                "Tm": round(self.Tm, 2)  # added round to display only 2 decimal points
                }


def create_reverse_complement(input_sequence):
    """
    Given an input sequence, returns its reverse complement.
    """
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
    bases = list(input_sequence)
    bases = reversed([complement.get(base, base) for base in bases])
    bases = ''.join(bases)
    return bases


def get_primers(input_sequence):
    """
    Given an input sequence, return a list of all possible Primer objects.
    """
    length = args.e
    primer_sequence_list = [input_sequence[i:j + 1] for i in range(length) for j in range(i + args.s, i + args.l)]
    # Turn primer sequences into Primer objects
    primer_list = [Primer(primer) for primer in primer_sequence_list]
    return primer_list


def filter_primers(list_of_primers):
    """
    For each Primer object, verify that its parameters (TM and GC content)
    are within specified ranges. Objects out of defined ranges are removed
    from the list.
    """
    filtered_list = list_of_primers.copy()
    for primer in list_of_primers:
        if primer.GC_percentage < args.M or primer.GC_percentage > args.X or primer.Tm < args.m or primer.Tm > args.x:
            filtered_list.remove(primer)
    return filtered_list


def bowtie2_alignment(forward_primer_list,reverse_primer_list,genome):
    """
    For a given list of Primer objects and a given genome, generates temporary files
    containing the primer sequences in FASTA format, and submits these files for alignment
    with the genome. The resulting SAM file is returned as a variable.
    """
    if not path.exists(f'{genome}_index.1.bt2'):
        build_index = f"bowtie2-build -f {genome} {genome}_index"
        subprocess.run(build_index,stdout=subprocess.PIPE, shell=True)
    quiet = '--quiet'
    if args.verbose:
        quiet = ''
    make_sam = ''
    if args.full:
        make_sam = '-S alignment_output.sam'
    align = f'bowtie2 -p {mp.cpu_count()} -x {genome}_index -X {args.a + 2*(args.l)} -f -1 {forward_primer_list} -2 {reverse_primer_list} {make_sam} -D 50 -N 1 --local --no-overlap --no-contain {quiet} 2> markers.txt'
    if args.verbose:
        print(f'Aligning with Bowtie2 – {mp.cpu_count()} threads will be utilized')
    alignment = subprocess.run(align,stdout=subprocess.PIPE, shell=True)
    if not args.full:
        remove({make_sam[2:]})
    return alignment

def create_sam_dataframe():
    '''
    creates a dataframe to store alignment results
    '''
    df_cols = [
                "Primer ID",        # STR
                "Primer Sequence",  # STR
                "5' Position",      # INT
                "Amplicon Size",    # INT
                "Alignment Score"   # STR
                ]
    dataframe = pd.DataFrame(columns=df_cols)
    dataframe.set_index('Primer ID', inplace=True)
    return dataframe

def parse_sam_file(input_line,dataframe):
    '''
    takes a line from a sam file generated by Bowtie2
    and stores alignment results in a given dataframe
    '''
    if not input_line.startswith('@'):
        line = input_line.split('\t')
        if 'fwd' in line[0]:
            five_prime = line[2]+'.'+line[3]
        elif 'rvs' in line[0]:
            fiveprimepos = str(int(line[3])+len(line[9]))
            five_prime = line[2]+'.'+fiveprimepos
        else: five_prime = None
        if five_prime is not None:
            data = {
                    "Primer ID":line[0],
                    "Primer Sequence":line[9],
                    "5' Position":five_prime,
                    "Amplicon Size":abs(int(line[8])),
                    "Alignment Score":line[5]
                    }
        else:
            sys.exit(f'Error while parsing Primer {line[0]}')
        df = pd.DataFrame(data=data,index=[0])
        df.set_index('Primer ID', inplace=True)
        dataframe.append(df,ignore_index=False,verify_integrity=False,sort=False)
    return


def main():
    if args.verbose:
        print(f'''
        Running script with the following parameters:
        {args}
        ''')

    # Importing target fragment fasta file
    if path.exists(args.i):
        with open(args.i, 'r') as f:
            contents = f.read()
            seqlist = contents.split()
        name_tag = seqlist[0]
        seq = ''.join(seqlist[1:])
        if args.verbose:
            print(f'file name is {name_tag} and sequence is {seq}')
    else:
        sys.exit('Input file does not exist at the specified path')

    # Importing genome file
    if path.exists(args.g):
        with open(args.g, 'r') as g:
            genome = g.read()
            if args.verbose:
                print(f'Genome file {args.g} successfully imported')
            linecount = genome.count('\n')
            if 2 * genome.count('>') != linecount + 1:
                if args.verbose:
                    print('formatting genome')
                from functions import formatted
                genome = formatted(genome)
                if args.verbose:
                    print(f'Genome file {args.g} successfully formatted')
        # from functions import generate_dictionary as gendict
        # genome_dictionary = gendict(genome)
        # if args.verbose:
        #     print(f'The genome was successfully converted to a dictionary')
    else:
        sys.exit('Genome file does not exist at specified path')

    ### FORWARD PRIMERS ###
    # Listing the Primer objects
    forward_primers = get_primers(seq)
    if args.verbose:
        print(f'{len(forward_primers)} forward primers successfully generated')
    # Filtering list of Priemr objects for TM / GC content
    filtered_forward_primers = filter_primers(forward_primers)
    if args.verbose:
        print(f'{len(forward_primers) - len(filtered_forward_primers)} forward primers discarded based on GC / TM')

    ### REVERSE PRIMERS ###
    # Listing the Primer objects
    reverse_primers = get_primers(create_reverse_complement(seq))
    if args.verbose:
        print(f'{len(reverse_primers)} reverse primers successfully generated')
    # Filtering list of Priemr objects for TM / GC content
    filtered_reverse_primers = filter_primers(reverse_primers)
    if args.verbose:
        print(f'{len(reverse_primers) - len(filtered_reverse_primers)} reverse primers discarded based on GC / TM')

    ### GENERATE PRIMER PAIRS ###
    # Generating all possible primer combinations
    primer_pairs = list(itertools.product(filtered_forward_primers, filtered_reverse_primers))
    # Filtering out all pairs that have Tm's that are too far apart
    filtered_primer_pairs = [pair for pair in primer_pairs if math.isclose(pair[0].Tm, pair[1].Tm, abs_tol=args.D)]
    if args.verbose:
        print(f'{len(filtered_primer_pairs)} primer pairs successfully created')

    ### ALIGNMENTS ###
    n1 = 1
    n2 = 1
    for primer_pair in filtered_primer_pairs:
        with open ('forward_list.fa','a') as for_list:    
            for_list.write(f'>{name_tag}fwd_primer_{n1}\n{primer_pair[0].sequence}\n\n')
        n1 = n1+1
        with open ('reverse_list.fa','a') as rev_list:    
            rev_list.write(f'>{name_tag}rvs_primer_{n2}\n{primer_pair[1].sequence}\n\n')
        n2 = n2+1
    alignment_result = bowtie2_alignment('forward_list.fa','reverse_list.fa',args.g)
    if not args.full:
        remove(f'forward_list.fa')
        remove(f'reverse_list.fa')
    if args.verbose:
        print(alignment_result)
    
    # with open(args.o, 'w') as outfile:
    #     outfile.write(alignment_result)


    # Creating elements directly stored in output dictionary
    # flat_primer_pair_info = []
    # for i in filtered_primer_pairs:
    #     flat_primer_pair_info.append(
    #         {
    #             "Forward primer": i[0].get_primer_elements(),
    #             "Reverse primer": i[1].get_primer_elements()
    #         }
    #     )

    # output_dict = {
    #     "Target gene": name_tag,
    #     "Genetic sequence": seq,
    #     "Total primer pairs": len(flat_primer_pair_info),
    #     "Generated with": "$software name",
    #     "Date": str(datetime.datetime.now()),
    #     "Primer Pair Info": flat_primer_pair_info
    # }

    # with open(args.o, 'w') as outfile:
    #     json.dump(output_dict, outfile, indent=2)
    # print(f'The output file has been generated at {args.o}')

if __name__ == '__main__':
    main()
