#!/usr/bin/env python3
''' Description: This script contains the viral evolution simulation script and is the basis for evolutionary inferences. It works as a random birth-death process which is parametrized by birth/death rates, fitness, mutation rates, and probability of beneficial mutations '''

def simulateViralEvolution(genomeSize, initSize, gens, mutRate, probBen, r, w, x, maxPopSize = 1e7, maxBenSize = 1e7, print_per_gen = 1):

    import random
    import numpy as np
    import time
    from virus import Virus, Population
    from copy import deepcopy,copy

    start_time = time.time()

    # Initialize 
    reference = [random.choice([1, 2, 3, 4]) for i in range(0, genomeSize)]
    neutral = [Virus(genomeSize) for i in range(0, initSize)]
    beneficial = []
    track_ben = []
    track_neut = []
    if print_per_gen != -1:
        # Print parameters
        print("¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬")
        print("Viral simulation")
        print("¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬")
        print('')
        print('1/3 Simulation parameters')
        print('-------------------------')
        print(f'Genome size (bp): {genomeSize}')
        print(f'Generations (divisions): {gens}')
        print(f'Mutation rate (µ): {mutRate}')
        print(f'Probability of a beneficial mutation: {probBen}')
        print(f'Viral replication rate (r): {r}')
        print(f'Viral death rate (x): {x}')
        print(f'Fitness (1 + s): {w}')
        print('')
        print('Replication coefficient neutral viruses : %s' %(r - x))
        print('Replication coefficient beneficial viruses: %s' %(r - x/w))
        print('')

    for j in range(0, gens+1):
        v1_prev,v2_prev = len(neutral),len(beneficial)

        if v1_prev > 0:
            # Update neutral population size
            neutral_rep = np.int(v1_prev * (r - x))
            neutral = [deepcopy(neutral[i]) for i in random.choices(range(len(neutral)), k = neutral_rep)]
            v1_new = len(neutral)

            # Add neutral mutations based on mutRate
            neutral_mut_ind = random.choices(range(v1_new), k = np.random.poisson(mutRate * v1_new))
            for i in neutral_mut_ind:
                neutral[i].neutralMutation(track_neut = track_neut)
            
            # Conversion to beneficial class at rate of mutRate * probBen
            new_ben_ind = random.choices(range(v1_new), k = np.random.poisson(mutRate * probBen * v1_new))
            new_ben = [(neutral[i]) for i in new_ben_ind]
            for new in new_ben:
                new.beneficialMutation(track_ben = track_ben)
            beneficial = beneficial + new_ben
            v2_new = len(beneficial)
            neutral = [neutral[i] for i in range(v1_new) if i not in new_ben_ind]

        if v2_new > 0:
            # Update beneficial population size
            beneficial_rep = np.int(len(beneficial) * (r - ( (1/w)*x ))) # Beneficial eliminated slower by a factor of 1/w
            beneficial = [deepcopy(beneficial[i]) for i in random.choices(range(len(beneficial)), k = beneficial_rep)]
            v2_new = len(beneficial)

            # Add neutral mutations to beneficial viruses
            neut_ben_mut_ind = random.choices(range(v2_new), k = np.random.poisson(mutRate * v2_new))
            for i in neut_ben_mut_ind:
                beneficial[i].neutralMutation(track_neut = track_neut)

            # Add beneficial mutations to beneficial viruses
            ben_mut_ind = random.choices(range(v2_new), k = np.random.poisson(mutRate * probBen * v2_new))
            for i in ben_mut_ind:
                beneficial[i].beneficialMutation(track_ben = track_ben)
        
        if print_per_gen != -1:
            if j == 0:
                print('2/3 Growth and rates per generation')
                print('-----------------------------------')
                print('{:6s} {:6s} {:6s} {:6s} {:6s}'.format('Gen', 'NP', 'BP', 'rN', 'rB'))
            else:
                if j % print_per_gen == 0:
                    print('{:2d} {:6d} {:5d} {:7.2f} {:6.2f}'.format(j, int(len(neutral)), int(len(beneficial)), v1_new/(v1_prev+1), v2_new/(v2_prev+1)))

        if (len(beneficial) + len(neutral)) > maxPopSize:
            break 
        if len(beneficial) > maxBenSize:
            break
    
    sample_matrix = []
    for i in range(len(neutral)):
        row = deepcopy(reference)
        for (pos,mut) in neutral[i].neutral_d.items():
            row[pos] = mut
        for (pos,mut) in neutral[i].beneficial_d.items(): # Can remove this; neutral don't have beneficial mutations
            row[pos] = mut
        sample_matrix.append(row)

    for i in range(len(beneficial)):
        row = deepcopy(reference)
        for (pos,mut) in beneficial[i].neutral_d.items():
            row[pos] = mut
        for (pos,mut) in beneficial[i].beneficial_d.items():
            row[pos] = mut
        sample_matrix.append(row)
  
    all_muts = np.array(sample_matrix)
                
    #final_beneficial = [[beneficial[j].neutral[i] if beneficial[j].beneficial[i] == 0 else beneficial[j].beneficial[i] for i in range(0, len(beneficial[j].beneficial))] for j in range(0, len(beneficial))]
    #final_neutral = [neutral[i].neutral for i in range(0, len(neutral))]
    #all_muts = np.array(final_beneficial + final_neutral) # (Rows = Samples, Columns = Genome Position)
    
    output = Population(r = r, w = w, x = x, probBen = probBen, mutRate = mutRate, genomeSize = genomeSize, initSize = initSize, gens = gens, numberBeneficial = len(beneficial), numberNeutral = len(neutral), mutations = all_muts, reference = reference, positiveLoci = track_ben, neutralLoci = [i for i in track_neut if i not in track_ben])
    
    if print_per_gen != -1:
        print('')
        print('3/3 Population description')
        print('--------------------------')
        print(f'Number of beneficial viruses: {output.numberBeneficial}')
        print(f'Number of neutral viruses: {output.numberNeutral}')
        elapsed_time = time.time() - start_time
        print('Time elapsed is ',elapsed_time, 'seconds')

    return(output)