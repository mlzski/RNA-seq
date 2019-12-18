#!/usr/bin/env python

# usage: postpro_qc.py <inputfile_directory> <outputfile_directory>
# EXAMPLE: python postpro_qc.py /nobackup/ummz/analysis_nov19/results/1_quality_control/ /nobackup/ummz/analysis_nov19/results/1_quality_control/postprocessed

import sys
import os
import zipfile
import pandas

#import numpy as np

dir_in = sys.argv[1]
dir_out = sys.argv[2]

print('Input directory: ', dir_in)
print('Output directory: ', dir_out)

# check if output directory exists, otherwise create it
if not os.path.exists(dir_out):
	print('Output directory does not exist. Creating: ', dir_out)
	os.makedirs(dir_out)

# create a directory for unzipped files, then unzip all files
print(dir_out + '/unzipped')
if not os.path.exists(dir_out + '/unzipped'):
	print('Creating new folder for unzipped files: ', dir_out + '/unzipped')
	os.makedirs(dir_out + '/unzipped')

# unzip all files and place them in new directory
print('Unzipping files ...')
for files in os.listdir(dir_in):
    if files.endswith(".zip"): 
	name = os.path.join(dir_in, files)
	zipfile.ZipFile(name).extractall(dir_out + '/unzipped')

dir_unzip = dir_out + '/unzipped'

# merge all summary.txt files
#for dirs in os.listdir(dir_out + '/unzipped'): 
#	f = open(dir_out + '/unzipped/' + dirs + '/summary.txt', 'r')
#	file_contents = f.read()
#	print(file_contents)
#	f.close

with open(dir_out + '/all_merged.txt', 'wb') as outfile:
	for dirs in sorted(os.listdir(dir_out + '/unzipped')):
		with open(dir_out + '/unzipped/' + dirs + '/summary.txt', 'r') as infile: 
       			outfile.write(infile.read())

# read all_merged.txt file as dataframe and then filter and count how many rows meet given conditions

	
# extract these lines that either WARN or FAIL
#for line in open(dir_out + '/all_merged.txt'):
#	if line.startswith('WARN') or line.startswith('FAIL'):
#		print(line)




