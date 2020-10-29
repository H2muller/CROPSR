#!/usr/bin/env python3
# Written by: Hans Müller Paul and Dave Istanto
#                           NOTES: adding switch-case statements replacing if clauses in functions

### Importing required libraries
import argparse
import cropsr_functions
import re
import sys
from multiprocessing import cpu_count
from pandas import Series

### Defining the arguments
parser = argparse.ArgumentParser()

parser.add_argument('-f', '--fasta', required=True, dest='f', 
                    help='[required] path to input file in FASTA format'
                    )
parser.add_argument('-g', '--gff', required=True, dest='g', 
                    help='[required] path to input file in GFF format'
                    )
parser.add_argument('-o', '--output', dest='o', default='data.csv',
                    help='path to output file, default = data.csv'
                    )
parser.add_argument('-m','--model', metavar ='', dest='m', type=str, default='Doench',
                    help='selects model that will be run to evaluate gRNAs for on-site activity, default = Doench'
                    )
parser.add_argument('-l', '--length', metavar='', dest='l', type=int, default=20,
                    help='length of the gRNA sequence, default = 20'
                    )
parser.add_argument('-L', '--flanking', metavar='', dest='L', type=int, default=200,
                    help='length of flanking region for verification, default = 200'
                    )
parser.add_argument('--cas9', action='store_true',
                    help='specifies that design will be made for the Cas9 CRISPR system'
                    )
parser.add_argument('--cpf1', action='store_true',
                    help='specifies that design will be made for the Cpf1 CRISPR system'
                    )
parser.add_argument('--CUDA', action='store_true',
                    help='runs processing steps utilizing the GPU instead of CPU where possible'
                    )
parser.add_argument('-v', '--verbose', action='store_true',
                    help='prints visual indicators for each iteration'
                    )

args = parser.parse_args()


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
##                                               Version: 1.01.1              ##
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
        Enable CUDA (requires nvidia GPU):              {args.CUDA}
        Designing for CRISPR system:
            Streptococcus pyogenes Cas9                 {args.cas9}
            Prevotella spp. Cpf1:                       {args.cpf1}
        """)

        # Path to annotation file in GFF format:          {args.g}
    ### Import genome files
    fasta_file = cropsr_functions.import_fasta_file(args.f)
    gff_df = cropsr_functions.import_gff_file(args.g)

    ### Locate PAMs by nuclease type
    if args.verbose:
        print(f'''
            Initiating PAM site detection.
            
            Please wait, this may take a while...
            ''')

    ### Create Dataframe containing all PAM site information

    CRIPSR_dataframe = cropsr_functions.create_dataframe()
    CRIPSR_dataframe.to_csv(args.o, header=True, index=False, mode='w')

    for chromosome,sequence in fasta_file.items():
        CRIPSR_dataframe = cropsr_functions.create_dataframe()
        Complete_dataset = []

        if args.cas9:
            # + strand
            motif = re.compile(r'(?=.GG)')
            cas9_target_list = cropsr_functions.find_PAM_site(motif,sequence)
            for target in cas9_target_list:
                pam_location = (target[0]-(args.l+1),target[0]-1)
                if pam_location[0] >= 0 and pam_location[0] <= len(sequence) and pam_location[1] >= 0 and pam_location[1] <= len(sequence):
                    shortseq = cropsr_functions.get_gRNA_sequence(sequence[pam_location[0]:pam_location[1]])
                    longseq = cropsr_functions.get_gRNA_sequence(sequence[pam_location[0]-5:pam_location[1]+5])
                    crispr_guide = [pam_location[0],pam_location[1],chromosome[1::],shortseq,longseq,'cas9','+']
                    Complete_dataset.append(crispr_guide)

            # - strand
            motif = re.compile(r'(?=CC.)')
            cas9_target_list2 = cropsr_functions.find_PAM_site(motif,sequence)
            for target in cas9_target_list2:
                pam_location = (target[0]+1,target[0]+(args.l+1))
                if pam_location[0] >= 0 and pam_location[0] <= len(sequence) and pam_location[1] >= 0 and pam_location[1] <= len(sequence):
                    shortseq = cropsr_functions.get_gRNA_sequence(cropsr_functions.get_reverse_complement(sequence[pam_location[0]:pam_location[1]]))
                    longseq = cropsr_functions.get_gRNA_sequence(cropsr_functions.get_reverse_complement(sequence[pam_location[0]-5:pam_location[1]+5]))
                    crispr_guide = [pam_location[1],pam_location[0],chromosome[1::],shortseq,longseq,'cas9','-']
                    Complete_dataset.append(crispr_guide)

            if args.verbose:
                print (f'''
                {len(cas9_target_list + cas9_target_list2)} Cas9 PAM sites were found on {chromosome[1::]}
                ''')

        if args.cpf1:
        # + strand
            motif = re.compile(r'(?=TTT.)')
            cpf1_target_list = cropsr_functions.find_PAM_site(motif,sequence)
            for target in cpf1_target_list:
                pam_location = (target[0]+4,target[0]+(args.l+4))
                if pam_location[0] >= 0 and pam_location[0] <= len(sequence) and pam_location[1] >= 0 and pam_location[1] <= len(sequence):
                    shortseq = cropsr_functions.get_gRNA_sequence(sequence[pam_location[0]:pam_location[1]])
                    longseq = cropsr_functions.get_gRNA_sequence(sequence[pam_location[0]-5:pam_location[1]+5])
                    crispr_guide = Series((pam_location[0],pam_location[1],chromosome[1::],shortseq,longseq,'cpf1','+'))
                    Complete_dataset.append(crispr_guide)

        # - strand
            motif = re.compile(r'(?=.AAA)')
            cpf1_target_list2 = cropsr_functions.find_PAM_site(motif,sequence)
            for target in cpf1_target_list2:
                pam_location = (target[0]-(args.l+1),target[0]-1)
                if pam_location[0] >= 0 and pam_location[0] <= len(sequence) and pam_location[1] >= 0 and pam_location[1] <= len(sequence):
                    shortseq = cropsr_functions.get_gRNA_sequence(cropsr_functions.get_reverse_complement(sequence[pam_location[0]:pam_location[1]]))
                    longseq = cropsr_functions.get_gRNA_sequence(cropsr_functions.get_reverse_complement(sequence[pam_location[0]-5:pam_location[1]+5]))
                    crispr_guide = Series((pam_location[1],pam_location[0],chromosome[1::],shortseq,longseq,'cpf1','-'))
                    Complete_dataset.append(crispr_guide)

            if args.verbose:
                print (f'''
                {len(cpf1_target_list + cpf1_target_list2)} Cpf1 PAM sites were found on {chromosome[1::]}
                ''')

        lesser_list = []
        for i,item in enumerate(Complete_dataset):
            lesser_list.append((item))
            if len(lesser_list) == 5000 or i == len(Complete_dataset)-1:
                CRIPSR_dataframe['raw'] = lesser_list
                CRIPSR_dataframe.apply(cropsr_functions.preprocess_PAM_sites,axis=1)
                CRIPSR_dataframe.apply(cropsr_functions.fill_row,axis=1)
                CRIPSR_dataframe.to_csv(args.o, header=False, index=False, mode='a+')
                CRIPSR_dataframe = cropsr_functions.create_dataframe()
                lesser_list = []

    ### CONFIRMATION MESSAGE
    if args.verbose:
        print(f'The output file has been generated at {args.o}')


if __name__ == '__main__':
    main()