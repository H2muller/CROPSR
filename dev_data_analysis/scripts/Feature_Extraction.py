#!/usr/bin/env python
# coding: utf-8

# In[2]:


# Module for sgRNA / site feature extraction
import numpy as np


# In[3]:


# 1st Order, one hot encoding: 30 sites * 4 states (bases) = 120 length
# Hierarchy -> Site index (1 before 2)
#           -> Base: A,T,G,C

def get_first_order_vec(sequence):
    first_order_vec = np.array([])
    for base in sequence:
        one_base_vec = np.zeros((4))
        if base == "A":
            one_base_vec[0] += 1
        elif base == "T":
            one_base_vec[1] += 1
        elif base == "G":
            one_base_vec[2] += 1
        elif base == "C":
            one_base_vec[3] += 1
        first_order_vec = np.concatenate([first_order_vec, one_base_vec])
    return first_order_vec


# In[4]:


# 2nd Order, one hot encoding: 29 sites * 16 states (base pairs possible) = 464 length
# Hierarchy -> Site index (1 before 2)
#           -> first base: A,T,G,C
#           -> second base: A,T,G,C
# Example: (AA, AT, AG, AC, TA, TT, ...)
def get_second_order_vec(sequence):
    second_order_vec = np.array([])
    for i in range(0,len(sequence)-1):
        two_base_vec = np.zeros((16))
        two_base = sequence[i:i+2]
        if two_base == "AA":
            two_base_vec[0] += 1
        elif two_base == "AT":
            two_base_vec[1] += 1
        elif two_base == "AG":
            two_base_vec[2] += 1
        elif two_base == "AC":
            two_base_vec[3] += 1    
        elif two_base == "TA":
            two_base_vec[4] += 1
        elif two_base == "TT":
            two_base_vec[5] += 1
        elif two_base == "TG":
            two_base_vec[6] += 1
        elif two_base == "TC":
            two_base_vec[7] += 1  
        elif two_base == "GA":
            two_base_vec[8] += 1
        elif two_base == "GT":
            two_base_vec[9] += 1
        elif two_base == "GG":
            two_base_vec[10] += 1
        elif two_base == "GC":
            two_base_vec[11] += 1 
        elif two_base == "CA":
            two_base_vec[12] += 1
        elif two_base == "CT":
            two_base_vec[13] += 1
        elif two_base == "CG":
            two_base_vec[14] += 1
        elif two_base == "CC":
            two_base_vec[15] += 1 
        second_order_vec = np.concatenate([second_order_vec, two_base_vec])
    return second_order_vec


# In[5]:


# Calculate melting temperature given sequence
def calculate_Tm(sequence):
    N = len(sequence)
    if N < 13:
        melt_temp = ((sequence.count('A') + sequence.count('T')) * 2 +
        (sequence.count('C') + sequence.count('G')) * 4)
    else:
        melt_temp = 64.9 + 41 * (sequence.count('G') + sequence.count('C') - 16.4) / N
    return melt_temp


# In[6]:


# Gives the feature vector of 30 bp site, based on:
# 1. First order position (120 length)
# 2. Second order position (464 length)
# 3. G/C content (2 vectors, 1 length each)
# 4. Melting temperature (4 length)
# Total: 590 length (features)

def ext_sgRNA_feat(site_30_bp):
    site_30_upper = site_30_bp.upper()
    
    # Get first order vector
    first_order_vec = get_first_order_vec(site_30_bp)
    
    # Get second order vector
    second_order_vec = get_second_order_vec(site_30_bp)
    
    # Getting G/C related features, from 20 bp guide, specific to 30 bp inputs
    twenty_bp = site_30_bp[5:25]
    
    # 1. G+C Count (length 1, 20 possible states)
    GC_count_vec = np.zeros((1))
    GC_count = twenty_bp.count("G") + twenty_bp.count("C")
    GC_count_vec[0] += GC_count
    
    # 2. G+C > or < 10 (length 1, 2 possible states)
    GC_fifty_vec = np.zeros((1))
    if GC_count >= len(twenty_bp)/2:
        GC_fifty_vec[0] += 1
    
    # Vectors melting temperature related features:
    # Melting Temperature Vector (probably in celcius)
    # Tm0: Melting Temp of all 30 bp site
    # Tm1: Melting Temp of first 5 bp of 20 bp site
    # Tm2: Melting Temp of the next 8 bp of 20 bp site
    # Tm3: Melting Temp of the next 5 bp (end of 20 bp sgRNA and NGG PAM site)
    
    tm_vec = np.zeros((4))
    
    tm_vec[0] = calculate_Tm(site_30_bp)
    tm_vec[1] = calculate_Tm(site_30_bp[5:10])
    tm_vec[2] = calculate_Tm(site_30_bp[10:18])
    tm_vec[3] = calculate_Tm(site_30_bp[18:23])
    
    #Extra info about biophysics
    
    full_feature_vector = np.concatenate([first_order_vec, second_order_vec, GC_count_vec, GC_fifty_vec, tm_vec])
    
    return full_feature_vector
    

