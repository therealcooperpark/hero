# HERO
**H**ighways **E**numerated by **R**ecombination **O**bservations

HERO is a pipeline designed to parse and visualize highways of genome-wide recombination using the output of the recombination detection tool fastGEAR. Currently, it works in 3 steps:

1) For each recombined DNA fragment, uses BLASTN to find the most likely donor from the fastGEAR donor cluster.

2) Calculates which pairs have statistically high rates of recombination within the measured population. A *Highway* is defined as a pair whose recombination events are at least 1 Stdev above the population's mean of recombination events per pair.

3) Generates a file formatted for the connections dataset in iToL to visualize recombination highways over a phylogenetic tree of the population.

## INSTALLATION
HERO has the following dependencies:
- NCBI-blast+

`git clone https://github.com/therealcooperpark/hero.git`

## 
