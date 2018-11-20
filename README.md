# HERO
**H**ighways **E**numerated by **R**ecombination **O**bservations

HERO is a pipeline designed to parse and visualize highways of genome-wide recombination using the output of the recombination detection tool fastGEAR. Currently, it works in 3 steps:

1) For each recombined DNA fragment, uses BLASTN to find the most likely donor from the fastGEAR donor cluster.

2) Calculates which pairs have statistically high rates of recombination within the measured population. A *Highway* is defined as a pair whose recombination events are at least 1 Stdev above the population's mean of recombination events per pair.

3) Generates a file formatted for the connections dataset in iToL to visualize recombination highways over a phylogenetic tree of the population.

## INSTALLATION
HERO has the following dependencies:
- NCBI-blast+

### Download command:
`git clone https://github.com/therealcooperpark/hero.git`

## USAGE
```
usage: hero.py [options] fastgear

positional arguments:
  fastgear      Directory containing all fastGEAR runs for the population

optional arguments:
  -h, --help    show this help message and exit
  --log         Minimum log(BF) length to accept recombination event. [0]
  --length      Minimum fragment length for recombined sequence. [2]
  --pairs       Ignore direction of event when calculating recombination
                highways.
  --cleanup     Keep intermediate files.
  --cpus        Number of threads to use. [1]
  --nohighway   Hexadecimal code for non-highway color in output file
                [#D3D3D3]
```

## OUTPUT
- bad_genes.txt

List of genes that failed HERO parsing.

- highways.txt

List of donor:recipient pairs and their respective recombination counts

- highways_hero_itol.txt

iToL formatted file of recombination highways to be dragged into population's phylogenetic tree on iToL

- recombination_hero_itol.txt

iToL formatted file of non-highway recombination events to be dragged into population's phylogenetic tree on iToL


## QUESTIONS
Feel free to email questions or concerns to therealcooperpark@gmail.com

## CITATION
If you use this program, please cite:
**Place future paper here**

