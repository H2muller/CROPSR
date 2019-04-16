#!/usr/bin/env python3 

import pandas as pd 
import numpy as np
import math 

#Read in data 
supp_3_csv = "/Users/daveistanto/Dropbox/UIUCGraduateSchool/Researches/CROPSR/Papers:Supplementary_Materials/Supp_3.csv"
supp_3_df = pd.read_csv(supp_3_csv)

#Sort experimental condition by name
feat_df = supp_3_df.iloc[:,3:]
feat_df = feat_df.reindex(sorted(feat_df.columns), axis=1)

supp_3_df = supp_3_df.drop(columns=(supp_3_df.columns.values[3:]))
supp_3_df = pd.concat([supp_3_df, feat_df], axis = 1)

#Filter row by given gene
target_gene = "ENST00000300060" # Can be changed later, as an input

def filter_by_gene(input_df, input_gene):
    filter_bool = input_df["Target"] == input_gene
    input_df = input_df[filter_bool]
    return input_df


#Get RPM and average RPM of given experimental condition
#Returns dataframe with transformed values
exp_cond = "MOLM13 CD15"

def get_processed_df(input_df, input_cond):
    #Filter by experimental condition
    drop_column_list = []
    for col_name in input_df.columns:
        if input_cond not in col_name and col_name != "Sequence" and col_name != "Target" and col_name != "Position":
            drop_column_list.append(col_name)
    input_df = input_df.drop(columns = drop_column_list)

    #Get RPM and transform it
    for col_index in range(3,len(input_df.columns)):
        col_sum = input_df[input_df.columns[col_index]].sum()
        input_df.iloc[:,col_index] /= col_sum
    
    #input_df.iloc[:,3:] = np.log2((input_df.iloc[:,3:] * 1000000) + 1)
    input_df.iloc[:,3:] = input_df.iloc[:,3:] * 1000000
    #Get average transformed RPM per row
    input_df["Mean Transformed RPM " + input_cond] = input_df.iloc[:,3:].mean(axis = 1)
    return input_df

def rank_sgRNA(input_processed_df):
    #Sort by last row (Mean RPM)
    input_processed_df = input_processed_df.sort_values(by = input_processed_df.columns[len(input_processed_df.columns) - 1], ascending = False)

    #Get percent rank and top 20% label
    rank_array = np.zeros(len(input_processed_df))
    quant_label = [0]*len(input_processed_df)
    for row_index in range(len(input_processed_df)):
        rank = row_index
        rank_array[row_index] += 1 - (rank/(len(input_processed_df) - 1))
        if rank_array[row_index] >= 0.8:
            quant_label[row_index] = True
        else:
            quant_label[row_index] = False
        
    input_processed_df["Percent Rank"] = rank_array
    input_processed_df["Top 20%"] = quant_label




    return input_processed_df.iloc[:,3:]



sample_filtered = filter_by_gene(supp_3_df, target_gene)
proc_df = get_processed_df(sample_filtered, exp_cond)
# print(rank_sgRNA(proc_df))
print(proc_df)





