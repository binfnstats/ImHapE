# Image-based haplotype-guided evolutionary inference <img align="right" width="150" height="150" src="https://github.com/tomouellette/IHapE/blob/master/icon.svg">

This repository contains code for *image-based haplotype-guided evolutionary inferences* for populations with phased genomes. **IHapE** treats evolutionary inference as an image recognition problem, harnessing the power of convolutional neural networks. Images are formed by aligning haplotypes, or entire genomes together, to make evolutionary inferences.

### Overview
---

IHapE is written in python and uses numpy to generate haplotype arrays. Any simulation software can be used to generate haplotypes for analysis, but we provide a simple evolutionary simulation framework that models an exponentially growing population under a a random birth-death process. It was developed for modeling SARS-CoV-2 evolution (allows for back mutations under a neutral model), but can be adapted to any non-recombining genome.

### Getting started
---

Evolutionary inference requires four basic steps: (1) simulations; (2) training the CNN; (3) converting empirical haplotypes into numpy arrays; (4) predicting evolutionary models. 

#### (1) Simulations

Simulations requires two scripts: `simulation.py` and `model.py` which define the evolutionary simulation and the sampling scheme to generate simulated haplotypes.  

```python
from simulation import simulateViralEvolution
import model as mod

neutral = simulateViralEvolution(r = 2.02, x = 1, w = 1, probBen = 0, mutRate = 1e-4, initSize = 250, genomeSize = 5000)
haplotypes, modes = mod.sampleData(sims = neutral, size = 200, reps = 1, flip = True)
```
<sub>*\** Parameter description: `simulateViralEvolution`, replication rate (r), death rate (x), fitness (w), mutation rate (mutRate), probability of a beneficial mutation (probBen; beneficial mutation rate equals mutRate \* probBen), initial population size (initSize), haplotype or genome size (genomeSize). Check `simulation.py` to see other default parameters. `sampleData`, the simulation (sims), the the height of the image or number of aligned haplotypes (size), the number of random samples from the population (reps), sort positions by frequency of mutations (flip)</sub>

If you would like to perform many simulations with `simulation.py` and automatically save the output in numpy format, you can run `exec.py` from the command line.

<pre><code>python3 exec.py -r 2.02 -w 1 -x 1 -p 0 -u 1e-2 -i 110 -gs 1000 -g 250 -ms 1e5 -n 250 -out ./output_folder/
</code></pre>

<sub>*\** The example above generates 250 unique simulations and saves output in the output folder (Note: / is required at end of output folder). `exec.py` will return each simulation formatted as [fitness]\_[id]\_[time].npy e.g. 1.05_7092001_18.npy. If you plan on using further scripts please do not change the [fitness] position. </sub>

#### (2) Training a CNN

We implemented a CNN with L2 regularization in tensorflow to analyze aligned haplotype data. To train the CNN, the simulated haplotypes must be converted into training and test data. This involves two steps using functions from `model.py`

```python
import model as mod

haplotypes, modes = mod.mergeData(positive_dir = ./positive_haplotypes/, neutral_dir = ./neutral_haplotypes, n = 100)
train_dataset, test_dataset = mod.trainTestData(haplotypes = haplotypes, modes = modes, p = 0.2)
```
<sub>*\** Parameter description: `mergeData`, the directory containing your positive simulations (positive_dir), the directory containing your neutral simulations (neutral_dir), and the number of samples of each evolutionary mode (n); `trainTestData`, haplotypes and modes generated from mergeData and the proportion of the data you are setting aside for testing/validation (p)</sub>

Once training and test data is generated, training the model is a simple function call.
```python
import model as mod

model, history, test_loss, test_acc = mod.trainCNN(train = train_dataset, test = test_dataset)
```
<sub>*\** If you would like to try out a different network architecture, please modify `trainCNN` in `model.py`</sub>


#### (3) Converting empirical haplotypes into numpy arrays

This task is dependent on your data, but the numpy array should share the same dimensions as your simulated data e.g (number of haplotypes, genome or haplotype size). The genome or haplotype size is defined by the `simulateViralEvolution(genomeSize = M)` call and the number of haplotypes is determined by `sampleData(size = N)` call as that defines how many aligned haplotypes per image.

#### (4) Predicting evolutionary models

Coming soon.


