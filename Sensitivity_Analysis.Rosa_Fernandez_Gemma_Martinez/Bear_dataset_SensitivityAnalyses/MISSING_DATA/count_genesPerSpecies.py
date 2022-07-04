#!/bin/python

import glob
from cogent3 import load_aligned_seqs, PROTEIN

orthogroups = glob.glob('*.fa')

global_counts = {}
count = 0

for i in orthogroups:
	names = []
	aln = load_aligned_seqs(i, moltype=PROTEIN)
	names = aln.names
	for j in names:
		if j not in global_counts.keys():
			global_counts[j] = 1
		elif j in global_counts.keys():
			global_counts[j] += 1
	count += 1

for sp in global_counts:
	print(sp, '\t', global_counts[sp], '\t', count, '\t', round(global_counts[sp]/float(count), 2))
