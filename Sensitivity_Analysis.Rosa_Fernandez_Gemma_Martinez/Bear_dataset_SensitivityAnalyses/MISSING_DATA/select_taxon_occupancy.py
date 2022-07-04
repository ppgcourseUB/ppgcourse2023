#!/usr/bin/python

from Bio import SeqIO
import glob
import os

j = int(input('Choose a minimum taxon occupancy: '))

all_files = glob.glob('*.fa*')

dirname = 'orthologs_min_' + str(j) + '_taxa'

os.mkdir(dirname)
	
for i in all_files:
	parsed = SeqIO.to_dict(SeqIO.parse(i, 'fasta'))
	if len(parsed) >= j:
		select = True
	elif len(parsed) > 4:
		if 'Gnat' in parsed:
			select = True
		if 'Lepa' in parsed:
			select = True
		else: 
			select = False
	else:
		select = False
	if select:
		os.chdir(dirname)
		for seq in parsed:
			out = open(i, 'a')
			out.write('>' + seq + '\n' + str(parsed[seq].seq) + '\n')
			out.close()
		os.chdir('..')
