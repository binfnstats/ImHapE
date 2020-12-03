#!/usr/bin/env python
''' Description: Takes a directory of EMBOSS stretcher aligned sequences and generates a pseudo-VCF file and an identifier file'''

import pandas as pd
import numpy as np

def processStretched(alignment_dir, reference_path):
	import pandas as pd
	import os
	reference = ''.join([seq.rstrip() for seq in open(reference_path,'r')])
	reference = list(reference)
	positions = [i for i in range(1,len(reference)+1)]
	vcf = pd.DataFrame({'positions': positions, 'reference': reference})
	check, count = [], 1
	for file in os.listdir(alignment_dir):
		print(count)
		count += 1
		try:
		    alignment = [seqs.strip() for seqs in open(alignment_dir + file, 'r')]
		    end = alignment.index('>EMBOSS_002')
		    alignment = ''.join(alignment[end+1:])
		    alignment = list(alignment)
		    if len(alignment) == len(reference):
		        vcf[file] = alignment
		    else:
		        check.append(file)
		except:
			pass
	names = vcf.columns[:1]
	country = [i.split('+')[0] for i in names]
	identifier = [i.split('+')[1] if len(i.split('+')) > 1 else None for i in names]
	date = [i.split('+')[2].split('.')[0] if len(i.split('+')) > 2 else None for i in names]
	ids = pd.DataFrame({'filename':names, 'country':country, 'id':identifier, 'date':date})

	return (check, vcf, ids)


check, vcf, ids = processStretched(alignment_dir = '/.mounts/labs/awadallalab/private/touellette/covid/alignments/', reference_path = '/.mounts/labs/awadallalab/private/touellette/covid/data/wuhan_1_reference.txt')

# Save files
output = '/.mounts/labs/awadallalab/private/touellette/covid/processed/'
vcf.to_csv(output + 'covid_vcf.csv', index=False)

