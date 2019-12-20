#!/usr/bin/env python

# TO BE EXECUTED BEFORE THIS SCRIPT: module load anaconda/2019.10

# usage: postpro_qc.py <inputfile_directory> <outputfile_directory>
# EXAMPLE: python postpro_qc.py /nobackup/ummz/analysis_nov19/results/1_quality_control/ /nobackup/ummz/analysis_nov19/results/1_quality_control/postprocessed

import sys
import os
import zipfile
import pandas as pd
import numpy as np

# running from terminal
dir_in = sys.argv[1]
dir_out = sys.argv[2]

# running in consol 
#dir_in = '/nobackup/ummz/analysis_nov19/results/1_quality_control'
#dir_out = '/nobackup/ummz/analysis_nov19/results/1_quality_control/postprocessed' 

print('Input directory: ', dir_in)
print('Output directory: ', dir_out)

# check if output directory exists, otherwise create it
if not os.path.exists(dir_out):
    print('Output directory does not exist. Creating: ', dir_out)
    os.makedirs(dir_out)

# create a directory for unzipped files, then unzip all files
print(dir_out + '/unzipped')
if not os.path.exists(dir_out + '/unzipped'):
    #print('Creating new folder for unzipped files: ', dir_out + '/unzipped')
    os.makedirs(dir_out + '/unzipped')

# unzip all files and place them in new directory
#print('Unzipping files ...')
for files in os.listdir(dir_in):
    if files.endswith(".zip"):
        name = os.path.join(dir_in, files)
        zipfile.ZipFile(name).extractall(dir_out + '/unzipped')

dir_unzip = dir_out + '/unzipped'

with open(dir_out + '/all_merged.txt', 'wb') as outfile:
    for dirs in sorted(os.listdir(dir_out + '/unzipped')):
        with open(dir_out + '/unzipped/' + dirs + '/summary.txt', 'rb') as infile: 
            outfile.write(infile.read())

# read all_merged.txt file as dataframe and then filter and count how many rows meet given conditions
df = pd.read_table(dir_out + '/all_merged.txt', header=None, names=('res', 'metric', 'file'))


# subset df to create one dataframe par metric 
df_1 = df.loc[df['metric'] == 'Basic Statistics']
df_2 = df.loc[df['metric'] == 'Per base sequence quality']
df_3 = df.loc[df['metric'] == 'Per tile sequence quality']
df_4 = df.loc[df['metric'] == 'Per sequence quality scores']
df_5 = df.loc[df['metric'] == 'Per base sequence content']
df_6 = df.loc[df['metric'] == 'Per sequence GC content']
df_7 = df.loc[df['metric'] == 'Per base N content']
df_8 = df.loc[df['metric'] == 'Sequence Length Distribution']
df_9 = df.loc[df['metric'] == 'Sequence Duplication Levels']
df_10 = df.loc[df['metric'] == 'Overrepresented sequences']
df_11 = df.loc[df['metric'] == 'Adapter Content']

# create a list with all metrics together (separately for R1 and R2)
lista_r1 = [df_1.iloc[::2, :], df_2.iloc[::2, :], df_3.iloc[::2, :], df_4.iloc[::2, :], df_5.iloc[::2, :], 
            df_6.iloc[::2, :], df_7.iloc[::2, :], df_8.iloc[::2, :], df_9.iloc[::2, :], df_10.iloc[::2, :], 
            df_11.iloc[::2, :]]

lista_r2 = [df_1.iloc[1::2, :], df_2.iloc[1::2, :], df_3.iloc[1::2, :], df_4.iloc[1::2, :], df_5.iloc[1::2, :], 
            df_6.iloc[1::2, :], df_7.iloc[1::2, :], df_8.iloc[1::2, :], df_9.iloc[1::2, :], df_10.iloc[1::2, :], 
            df_11.iloc[1::2, :]]

print()
warned_r1 = []
idx_w_r1 = []
#list_out_warned_r1 = list()

print('Metrics flagged with WARN in Read 1:')
for i in range(len(lista_r1)):
	warned_r1 = lista_r1[i]['res'] == 'WARN'
	if any(lista_r1[i]['res'] == 'WARN') == True:
		print(str(np.sum(warned_r1)) + ' samples were flagged with WARN => ' + np.unique(lista_r1[i]['metric']))
		idx_w_r1 = np.where(warned_r1)[0]
		#list_out_warned_r1.append(lista_r1[i].iloc[idx_w_r1,2])	
		print(lista_r1[i].iloc[idx_w_r1,2].to_string(index=False))

print()
failed_r1 = []
idx_f_r1 = []
#list_out_failed_r1 = list()

print('Metrics flagged with FAIL in Read 1:')
for i in range(len(lista_r1)):
	failed_r1 = lista_r1[i]['res'] == 'FAIL'
	if any(lista_r1[i]['res'] == 'FAIL') == True:
		print(str(np.sum(failed_r1)) + ' samples were flagged with FAIL => ' + np.unique(lista_r1[i]['metric']))
		idx_f_r1 = np.where(failed_r1)[0]
		#list_out_failed_r1.append(lista_r1[i].iloc[idx_f_r1,2])
		print(lista_r1[i].iloc[idx_f_r1,2].to_string(index=False))

print()
warned_r2 = []
idx_w_r2 = []
#list_out_warned_r2 = list()

print('Metrics flagged with WARN in Read 2:')
for i in range(len(lista_r2)):
	warned_r2 = lista_r2[i]['res'] == 'WARN'
	if any(lista_r2[i]['res'] == 'WARN') == True:
		print(str(np.sum(warned_r2)) + ' samples were flagged with WARN => ' + np.unique(lista_r2[i]['metric']))
		idx_w_r2 = np.where(warned_r2)[0]
		#list_out_warned_r2.append(lista_r2[i].iloc[idx_w_r2,2])	
		print(lista_r2[i].iloc[idx_w_r2,2].to_string(index=False))

print()
failed_r2 = []
idx_f_r2 = []
#list_out_failed_r2 = list()

print('Metrics flagged with FAIL in Read 2:')
for i in range(len(lista_r2)):
	failed_r2 = lista_r2[i]['res'] == 'FAIL'
	if any(lista_r2[i]['res'] == 'FAIL') == True:
		print(str(np.sum(failed_r2)) + ' samples were flagged with FAIL => ' + np.unique(lista_r2[i]['metric']))
		idx_f_r2 = np.where(failed_r2)[0]
		#list_out_failed_r2.append(lista_r2[i].iloc[idx_f_r2,2])
		print(lista_r2[i].iloc[idx_f_r2,2].to_string(index=False))








