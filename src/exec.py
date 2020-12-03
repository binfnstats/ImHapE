#!/usr/bin/env python3
''' Description: This script performs simulations of viral evolution using specified parameters such as replication rate, death rate, fitness, etc. '''

from simulation import simulateViralEvolution
import model as mod
import numpy as np
import argparse
import random
import time
from virus import Virus, Population
from copy import deepcopy,copy

# Arguments
parser = argparse.ArgumentParser(description = 'Simulation parameters')
parser.add_argument('-r', '--rep_rate', default=2.02, type = float, help='Replication rate (r)')
parser.add_argument('-w', '--fitness', default=1, type = float, help='Fitness (1 + s)')
parser.add_argument('-x', '--death_rate', default=1, type = float, help='Extinction/death rate (x)')
parser.add_argument('-p', '--probBen', default=0.01, type = float, help='Probability of positively selected mutation')
parser.add_argument('-u', '--mutRate', default=1e-2, type = float, help='Mutation rate')
parser.add_argument('-i', '--initSize', default=110, type = int, help='Initial population size')
parser.add_argument('-gs', '--genomeSize', default=1000, type = int, help='Genome size (bp)')
parser.add_argument('-g', '--gens', default=250, type = int, help='Generations')
parser.add_argument('-ms', '--maxPopSize', default=1e5, type = float, help='Generations')
parser.add_argument('-n', '--numberSims', default=500, type = int, help='Number of simulations to perform.')
parser.add_argument('-ims', '--imageSize', default=200, type = int, help='Number of aligned haplotypes per image.')
parser.add_argument('-out', '--output', help='Output directory')
parser.add_argument('-max', '--max_mutations', type = int, default = 1000, help='Maximum number of mutations per aligned haplotype block.')
parser.add_argument('-min', '--min_mutations', type = int, default = 10, help='Minimum number of mutations per aligned haplotype block.')
parser.add_argument('-sort', '--sort', default = 'none', help='Sort rows or columns')

# Example argument
# - python3 exec.py -r 2.02 -w 1.1 -x 1 -p 0.01 -u 1e-2 -i 110 -gs 1000 -g 250 -ms 1e5 -n 1 -out test

args = parser.parse_args()
r = args.rep_rate
w = args.fitness
x = args.death_rate
p = args.probBen
u = args.mutRate
gs = args.genomeSize
g = args.gens
i = args.initSize
ims = args.imageSize
ms = int(args.maxPopSize)
n = args.numberSims
out = str(args.output)
max_muts = args.max_mutations
min_muts = args.min_mutations
sort = str(args.sort)

if sort == 'row':
	sortrow, sortcol = (True, False)
elif sort == 'col':
	sortrow, sortcol = (False, True)
elif sort == 'row_col':
	sortrow, sortcol = (True, True)
else:
	sortrow, sortcol = (False, False)

randid = str(random.randint(1, 1e7)) # Assign random id to for saving simulations

for val in range(n):

	# Run simulations
	print('[1] Simulation started \n')
	sim = simulateViralEvolution(r = r, w = w, x = x, probBen = p, mutRate = u, genomeSize = gs, initSize = i, gens = g, maxPopSize = ms, print_per_gen = 10)
	
	# Only accept simulations with sufficient mutations across each aligned haplotypes
	print('[2] Extracting haplotype sample \n')
	nmuts, count = 0, 0
	while nmuts < min_muts or nmuts > max_muts:
		haplotypes, modes = mod.sampleData(sims = sim, size = ims, reps = 1, sort_row = sortrow, sort_col = sortcol)
		nmuts = np.sum(haplotypes != 0)
		count += 1
		if count > 20:
			break

	# Simulation names
	if w > 1:
		name = w
	else:
		name = 'neutral'
	
	np.save(out + str(name) + '_' + randid + '_' + str(val), haplotypes)

	sim = None
	haplotypes = None
	modes = None
