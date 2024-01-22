#!/usr/bin/env python
# coding: utf-8

# In[1]:


import os
def get_last_folder_name(folder_path):

    #This function takes a folder path as input and returns the name of the last folder in the path.

    return os.path.basename(folder_path)

def get_subfolders(folder_path):

    #This function takes a folder path as input and returns a list of all subfolder paths.
    subfolders = [os.path.join(folder_path, f) for f in os.listdir(folder_path) if os.path.isdir(os.path.join(folder_path, f))]
    folder_dict = {}
    for subfolders_path in subfolders:
        names = get_last_folder_name(subfolders_path)
        folder_dict[names] = subfolders_path
    return folder_dict

#Test the function with a sample path (won't execute here as it's a code cell)



# In[2]:


def read_fasta(file_path):
    sequences = {}
    with open(file_path, 'r') as fasta_file:
        sequence_id = None
        sequence = ''
        for line in fasta_file:
            line = line.strip()
            if line.startswith('>'):
                if sequence_id:
                    sequences[sequence_id] = sequence
                sequence_id = line[1:]
                sequence = ''
            else:
                sequence += line
        if sequence_id:
            sequences[sequence_id] = sequence
    return sequences


# In[3]:


def read_csv_check_columns(csv_path, column_names):
    """
    This function reads a CSV file into a pandas DataFrame and checks if specific columns exist.
    """
    # Read the CSV file
    df = pd.read_csv(csv_path)
    
    # Check if the specified columns are in the DataFrame
    columns_exist = all(column in df.columns for column in column_names)
    
    return df, columns_exist


# In[4]:


def check_and_create_folder(folder_path):
    """
    This function checks if a folder exists at the given path. 
    If it doesn't exist, the function creates it.
    """
    if not os.path.exists(folder_path):
        os.makedirs(folder_path)
        print(f"Folder created at: {folder_path}")
    else:
        print(f"Folder exists at: {folder_path}")


# In[5]:


def get_all_file_paths(directory):
    file_paths = []  # List to store file paths

    # Walking top-down from the root
    for root, directories, files in os.walk(directory):
        for filename in files:
            # Join the two strings to form the full file path.
            filepath = os.path.join(root, filename)
            file_paths.append(filepath)  # Add it to the list.

    return file_paths


# In[6]:


import pandas as pd
def check_key_in_dict(dictionary, key_to_check):
    """
    This function checks if a given key exists in a dictionary.

    """
    return key_to_check in dictionary
def process_dataframe_rows(df,path_dict,fasta_sequences):
    """
    This function iterates through a DataFrame row by row, checks each row using a provided function,
    and separates the rows into two DataFrames based on the check result.

    """

    NGS_id,Sample_id,Locus,NGS_id_fail,Sample_id_fail,Locus_fail,reason = [],[],[],[],[],[],[]
  
    for _, row in df.iterrows():
        if check_key_in_dict(path_dict,row['NGS_id']) == False or check_key_in_dict(fasta_sequences,row['Locus']) == False:
            if check_key_in_dict(path_dict,row['NGS_id']) == False:
                reason.append('no sequence file found')
                NGS_id_fail.append(row['NGS_id'])
                Sample_id_fail.append(row['Sample_id'])
                Locus_fail.append(row['Locus'])
            else:
                reason.append('no reference seuqnce found')

                NGS_id_fail.append(row['NGS_id'])
                Sample_id_fail.append(row['Sample_id'])
                Locus_fail.append(row['Locus'])
        else:
                NGS_id.append(row['NGS_id'])
                Sample_id.append(row['Sample_id'])
                Locus.append(row['Locus'])
    dict_pass = {'NGS_id':NGS_id,'Sample_id':Sample_id,'Locus':Locus}
    dict_fail = {'NGS_id':NGS_id_fail,'Sample_id':Sample_id_fail,'Locus':Locus_fail,'reason':reason}
    df_pass = pd.DataFrame.from_dict(dict_pass)
    df_fail = pd.DataFrame.from_dict(dict_fail)
    return df_pass , df_fail

import zipfile
import os
def unzip(zip_file_path):
    output_folder = os.path.dirname(zip_file_path)

    # Extract the contents of the zip file
    with zipfile.ZipFile(zip_file_path, 'r') as zip_ref:
        zip_ref.extractall(output_folder)


# In[9]:


file_path = 'retron_reference.fasta'
fasta_sequences = read_fasta(file_path)


# In[10]:


sample_path = './data'
path_dict = get_subfolders(sample_path)


# In[11]:


csv_path="retron_batch9_test.csv"
working_path = './work'
out_path = './out'
column_names = ['NGS_id','Sample_id','Locus']


# In[85]:


import pandas as pd
import subprocess
import os
import shutil

dataframe, columns_exist = read_csv_check_columns(csv_path, column_names)
if columns_exist == False:
    print("Columns not right")
else:
    df_pass , df_fail = process_dataframe_rows(dataframe,path_dict,fasta_sequences)
    df_fail.to_csv('fail_list.csv',sep = ',')
check_and_create_folder(working_path)
check_and_create_folder(out_path)
for num in range(len(df_pass)):
    site_ = df_pass['Locus'][num]
    reference_seq = fasta_sequences[site_.replace(" ", "")]
    NGS_id_ = df_pass['NGS_id'][num]
    sub_seq_folder = path_dict[NGS_id_]
    seq_paths = get_all_file_paths(sub_seq_folder)
    r1 , r2 = 'none','none'
    for path in seq_paths:
        if '_1.fq' in path:
            r1 = path
            #print(path)
        elif '_2.fq' in path:
            r2 = path
    command = 'CRISPResso --fastq_r1 ' + r1 + ' --fastq_r2 '+r2 + ' --amplicon_seq '+ reference_seq + ' -o '+ working_path
    print("processing: "+df_pass['NGS_id'][num])
    process = subprocess.run(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    print(process)
out_paths = get_all_file_paths(working_path)
for num in range(len(df_pass)):
    _NGS_id_ = '_'+df_pass['NGS_id'][num]+'_'
    Sample_id_ = df_pass['Sample_id'][num]
    for path_ in out_paths:
        if _NGS_id_ in path_ and 'Alleles_frequency_table.zip' in path_:
            unzip(path_)
            unzip_path = path_.replace(".zip", ".txt")
            new_file_path = os.path.join(out_path,Sample_id_+'.txt' ) 
            shutil.copy2(unzip_path, new_file_path)
        else:
            pass


# In[ ]:




