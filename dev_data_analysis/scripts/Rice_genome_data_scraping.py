#!/usr/bin/env python
# coding: utf-8

# In[1]:


# Module for getting 30 bp location near cut site, then adding it to the dataset as extra column
import pandas as pd


# In[2]:


# Import original dataset
rice_sgRNA_data_path = "/Users/daveistanto/Dropbox/UIUCGraduateSchool/Researches/CROPSR/data_files/rice_data.csv"
rice_df = pd.read_csv(rice_sgRNA_data_path, skiprows=1)


# In[3]:


# Get rice genome fasta file, separated per locus, and make dictionary out of it

rice_loc_seq = "/Users/daveistanto/Dropbox/UIUCGraduateSchool/Researches/CROPSR/data_files/Rice_genome/all.seq"
with open(rice_loc_seq) as loc_file:
    loc_data = loc_file.readlines()
    loc_seq_string = ""
    loc_seq_dict = dict()
    for line_index in range(len(loc_data)):
        if loc_data[line_index][0] == ">":
            curr_locus = loc_data[line_index][1:15]
        else:
            loc_seq_string += loc_data[line_index][:-1]
            if line_index == len(loc_data) - 1 or loc_data[line_index + 1][0] == ">":
                loc_seq_dict[curr_locus] = loc_seq_string
                loc_seq_string = ""


# In[4]:


# Get complement strand

def get_complement(input_strand):
    rev_strand = input_strand[::-1]
    com_strand = []
    for base in rev_strand:
        if base == "A":
            com_strand.append("T")
        elif base == "T":
            com_strand.append("A")
        elif base == "G":
            com_strand.append("C")
        elif base == "C":
            com_strand.append("G")
    com_strand = ''.join(com_strand)
    return com_strand
    


# In[5]:


# Make complement strand dictionary
comp_dict = dict()

for k, v in loc_seq_dict.items():
    comp_dict[k] = get_complement(v)


# In[6]:


# Iterate over rice_df rows, and add actual relative position and 30 bp site
rice_df["Actual_relative_position_in_respective_strand"] = "NA"
rice_df["Site_30_bp"] = "NA"
non_found_index = []

for i, row in rice_df.iterrows():
    print(i/len(rice_df))
    if (row["Strand"] == "+"):
        try:
            # Get index of sequence in the locus, 5' - 3' strand
            act_loc_pos = loc_seq_dict[row["Locus"]].index(row["SG20_PAM"])
            site_30_bp = loc_seq_dict[row["Locus"]][act_loc_pos - 5: act_loc_pos + 25]
            rice_df["Actual_relative_position_in_respective_strand"][i] = act_loc_pos
            rice_df["Site_30_bp"][i] = site_30_bp
        except:
            non_found_index.append(i)
    elif (row["Strand"] == "-"):
        try:
            # Get index of sequence in the locus, 3' - 5' strand
            act_loc_pos = comp_dict[row["Locus"]].index(row["SG20_PAM"])
            site_30_bp = comp_dict[row["Locus"]][act_loc_pos - 5: act_loc_pos + 25]
            rice_df["Actual_relative_position_in_respective_strand"][i] = act_loc_pos
            rice_df["Site_30_bp"][i] = site_30_bp
        except:
            non_found_index.append(i)


# In[10]:


# Write new csv to secondary memory (Storage) and print non_found_index to stdout
print(non_found_index)
rice_df.to_csv("/Users/daveistanto/Dropbox/UIUCGraduateSchool/Researches/CROPSR/data_files/complete_rice_data.csv", index=False)


# In[ ]:





# In[ ]:




