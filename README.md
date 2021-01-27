# Image-based haplotype-guided evolutionary inference <img align="right" width="150" height="150" src="https://github.com/tomouellette/ImHapE/blob/master/icon.svg">

This repository contains code for *image-based haplotype-guided evolutionary inferences* for populations with phased genomes (**ImHapE**). 

### Overview
---
<img align ="left" src="https://github.com/tomouellette/ImHapE/blob/master/infograph.svg" width = "200" height = "400"> ImHapE treats evolutionary inference as an image recognition problem, harnessing the power of convolutional neural networks and utilizes an LSTM to account for potential linkage or positional information across full viral genomes. Images are formed by aligning haplotypes, or entire genomes together, to make evolutionary inferences. ImHapE is written in python and uses numpy to generate haplotype arrays. Any simulation software can be used to generate haplotypes for analysis, but we provide a simple evolutionary simulation framework that models an exponentially growing population under a a random birth-death process. It was developed for modeling SARS-CoV-2 evolution (allows for back mutations under a neutral model) but can be adapted to any non-recombining genome.

### Getting started
---

Running a simulation interactively is straight forward:

```python
from simulation import simulateViralEvolution
import model as mod

neutral = simulateViralEvolution(r = 2.02, x = 1, w = 1, probBen = 0, mutRate = 1e-4, initSize = 250, genomeSize = 5000)
haplotypes, modes = mod.sampleData(sims = neutral, size = 200, reps = 1, sort_row = False, sort_col = True)
```

If you would like to perform many simulations with `simulation.py` and automatically save the output in numpy format, you can run `exec.py` from the command line.

```bash
python3 exec.py -r 2.02 -w 1 -x 1 -p 0 -u 1e-2 -i 110 -gs 1000 -g 250 -ms 1e5 -n 250 -out ./output_folder/
```

### Making inferences
---

If you would like to build your own CNN/RNN models, you must (i) simulate thousands of genomic windows of length N (e.g. 2500 bp) under a set of realistic parameters (simulation.py), (ii) train CNNs on the small genomic windows (cnn.py), (iii) run sliding window analysis on full 29903 bp genomes (genome.py), (iii) train RNN on sliding window estimates across simulated genomes. This pipeline is ideally run on a compute cluster.

Most functions can be found within model.py. As stated above, any simulation software can be used as long as you convert your simulated haplotypes/alignments to numpy format. 

Lots of potential avenues to improve inferences, feel free to use the code however you see fit. And always open to comments!

### Contact
---

If you have any questions, please email t.ouellette@mail.utoronto.ca.

---

Preprint is available at https://www.biorxiv.org/content/10.1101/2021.01.13.426571v1.

<b>A few references to other papers using CNN/RNN approaches in population genetics</b>

<sub>i. ImaGene: a convolutional neural network to quantify natural selection from genomic data. *BMC Bioinformatics*. Torada et al. (2019)</sub>

<sub>ii. The Unreasonable Effectiveness of Convolutional Neural Networks in Population Genetic Inference. *Molecular Biology and Evolution*. Flagel et al. (2018)</sub>

<sub>iii. Inferring the landscape of recombination using recurrent neural networks. *bioRxiv*. Adrion et al. (2019)</sub>
