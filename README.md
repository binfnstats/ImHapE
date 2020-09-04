# Image-based haplotype-guided evolutionary inference <img align="right" width="150" height="150" src="https://github.com/tomouellette/IHapE/blob/master/icon.svg">

This repository contains code for *image-based haplotype-guided evolutionary inferences* for populations with phased genomes. **IHapE** treats evolutionary inference as an image recognition problem, harnessing the power of convolutional neural networks. Images are formed by aligning haplotypes, or entire genomes together, to make evolutionary inferences.

### Overview
---

IHapE is written in python and uses numpy to generate haplotype arrays. Any simulation software can be used to generate haplotypes for analysis, but we provide a simple evolutionary simulation framework that models an exponentially growing population under a a random birth-death process. It was developed for modeling SARS-CoV-2 evolution (allows for back mutations under a neutral model), but can be adapted to any non-recombining genome.

### Getting started
---

Inferences require four basic steps: (1) simulations; (2) training the CNN; and (3) converting empirical haplotypes into numpy arrays; (4) predicting evolutionary model. 

#### (1) Simulations

Simulations require one script, simulation.py (Note: virus.py must be accessible to simulation.py as it defines classes used within the function)

```python
from virus import Virus, Population
from simulation import simulateViralEvolution

neutral = simulateViralEvolution(r = 2.02, x = 1, w = 1, probBen = 0, mutRate = 1e-4, initSize = 250, genomeSize = 5000)
```
<sub>The parameters are replication rate (r), death rate (x), fitness (w), mutation rate (mutRate), probability of a beneficial mutation (probBen; beneficial mutation rate equals mutRate \* probBen), initial population size (initSize), haplotype or genome size (genomeSize). Check `simulation.py` to see other default parameters. </sub>

<sub>**Note**: If you would like to perform many simulations with `simulation.py` and automatically save the output in numpy format, you can run `exec.py` from the command line. </sub>

<pre><code>python3 exec.py -r 2.02 -w 1 -x 1 -p 0 -u 1e-2 -i 110 -gs 1000 -g 250 -ms 1e5 -n 250 -out ./output_folder/
</code></pre>

<sub> The example above generates 250 unique simulations and saves output in the output folder (Note: / is required at end of output folder). `exec.py` will return each simulation formatted as [fitness]\_[id]\_[time].npy e.g. 1.05_7092001_18.npy. If you plan on using further scripts please do not change the [fitness] position. </sub>

#### (2) Training a CNN

- We implemented a CNN in tensorflow to analyze aligned haplotype data. We additionaly provide an option to sort positions by mutation frequency using the `flip = True` option in the `model.py -> sampleData` function. Sorting positions has previously been shown to improve alignment based CNN inferences. 

- To train the CNN, the simulated haplotypes must be converted into training and test data. This involves two steps using functions from `model.py`

```python
import model as mod

haplotypes, modes = mod.mergeData(positive_dir = ./positive_haplotypes/, neutral_dir = ./neutral_haplotypes, n = 100)
train_dataset, test_dataset = trainTestData(haplotypes = haplotypes, modes = modes, p = 0.2)
```
- The parameters for `mergeData` are the directory containing your positive simulations (positive_dir), the directory containing your neutral simulations (neutral_dir), and the number of samples per mode (n). The parameters for `trainTestData` are the haplotypes and modes generated from mergeData and the proportion of the data you are setting aside for testing/validation.




