# Image-based haplotype-guided evolutionary inference <img align="right" width="150" height="150" src="https://github.com/tomouellette/ImHapE/blob/master/icon.svg">

This repository contains code for *image-based haplotype-guided evolutionary inferences* for populations with phased genomes (**ImHapE**). 

### Overview
---

ImHapE treats evolutionary inference as an image recognition problem, harnessing the power of convolutional neural networks and utilizes LSTM to account for potential linkage or positional information. Images are formed by aligning haplotypes, or entire genomes together, to make evolutionary inferences. ImHapE is written in python and uses numpy to generate haplotype arrays. Any simulation software can be used to generate haplotypes for analysis, but we provide a simple evolutionary simulation framework that models an exponentially growing population under a a random birth-death process. It was developed for modeling SARS-CoV-2 evolution (allows for back mutations under a neutral model) but can be adapted to any non-recombining genome.

<img align="right" width="429" height="281" src="https://github.com/tomouellette/ImHapE/blob/master/infograph.pdf">

### Getting started
---

We provided a basic overview of some scripts in the code overview documentation.


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


### Contact
---

If you have any questions, please email t.ouellette@mail.utoronto.ca.
