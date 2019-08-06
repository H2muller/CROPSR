#!/usr/bin/env python3

import itertools
import re
import sys

def generate_dictionary(input):
    dictionary = input.split()
    dictionary = dict(itertools.zip_longest(*[iter(dictionary)] * 2, fillvalue=""))
    return dictionary


def location(primer, genome):
    '''
    Written by: Hans Müller Paul and Zhiwen Jiang
    '''
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
    '''
    Written by: Hans Müller Paul and Joao Paulo Gomes Viana
    '''
    formatted = re.sub('\n','',input_genome)
    formatted = re.sub('>', '\n>',formatted)
    formatted = formatted[1:]
    formatted = re.sub('([0-9]+)','\\1 \n',formatted)
    return formatted

