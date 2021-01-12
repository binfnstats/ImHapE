#!/usr/bin/env python3
''' Description: This script trains the deep learning (CNN) model on haplotype data from simulations '''

WUHAN_REFERENCE = '/WUHAN_PATH/data/wuhan_1_reference.txt'
ENGLAND_REFERENCE = '/COG_REFERENCE_PATH/cog_sequences/England-LOND-D51C5-2020.txt'
PROBLEMATIC_SITES = '/PROBLEMATIC_SITES_PATH/problematic_sites_sarsCov2.v2.caution.vcf'

GISAID_DIRECTORY = '/GISAID_PATH/gisaid_sequences/'
COG_DIRECTORY = '/COG_PATH/cog_sequences/'

import numpy as np 
import os
import pandas as pd

wuhan_ref = ''.join([i.rstrip() for i in open(WUHAN_REFERENCE, 'r')])
england_ref = [i.rstrip() for i in open(ENGLAND_REFERENCE, 'r')][0]
problematic = pd.read_csv(PROBLEMATIC_SITES,sep='\t', comment = '#', header=0)
problematic = problematic.reset_index().T.reset_index().T
problematic.iloc[:,2] = [int(i) for i in problematic.iloc[:,2].tolist()]
problem_sites = [i for i in problematic.iloc[:,2].tolist()]

def convertToBinary(reference, sequence):
  converted = [int(sequence[i] != reference[i] and sequence[i] != 'N' and sequence[i] != '-' and i not in problem_sites) for i in range(len(reference))]
  converted[21764] = 1 if sequence[21764] == '-' else 0  
  converted[21765] = 1 if sequence[21765] == '-' else 0
  converted[21766] = 1 if sequence[21766] == '-' else 0   # Check for deletion between 21765-21770; not confident in other gaps
  converted[21767] = 1 if sequence[21767] == '-' else 0
  converted[21768] = 1 if sequence[21768] == '-' else 0
  converted[21769] = 1 if sequence[21769] == '-' else 0   
  tostring = ''.join([str(i) for i in converted])
  return (tostring)

# Process GISAID sequences with Wuhan reference
# --------------------------------------------------
GISAID_WUHAN_OUTPUT = '/.mounts/labs/awadallalab/private/touellette/covid/haplotypes_gisaid_wuhan/'
os.chdir(GISAID_DIRECTORY)
for file in os.listdir():
  seq = [i.rstrip() for i in open(file, 'r')][0]
  out = convertToBinary(reference = wuhan_ref, sequence = seq)
  with open(GISAID_WUHAN_OUTPUT + file, "w") as text:
    text.write(out)

# Process COG England sequences with Wuhan reference
# --------------------------------------------------
COG_WUHAN_OUTPUT = '/.mounts/labs/awadallalab/private/touellette/covid/haplotypes_cog_wuhan/'
os.chdir(COG_DIRECTORY)
for file in os.listdir():
  seq = [i.rstrip() for i in open(file, 'r')][0]
  out = convertToBinary(reference = wuhan_ref, sequence = seq)
  with open(COG_WUHAN_OUTPUT + file, "w") as text:
    text.write(out)

# Process COG England sequences with England reference
# ----------------------------------------------------
COG_ENGLAND_OUTPUT = '/.mounts/labs/awadallalab/private/touellette/covid/haplotypes_cog_england/'
os.chdir(COG_DIRECTORY)
for file in os.listdir():
  seq = [i.rstrip() for i in open(file, 'r')][0]
  out = convertToBinary(reference = england_ref, sequence = seq)
  with open(COG_ENGLAND_OUTPUT + file, "w") as text:
    text.write(out)