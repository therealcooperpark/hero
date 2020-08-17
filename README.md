# HERO
**H**ighways **E**numerated by **R**ecombination **O**bservations


HERO is a pipeline, written in Python, designed to parse and visualize highways of genome-wide recombination using the output of the recombination detection tool [fastGEAR](https://mostowylab.com/news/fastgear?rq=fastgear). Currently, it works in 3 steps:

1) For each recombination event, we use BLASTN to find the most likely donor across a pool of potential genomes (determined by the clustering alogrithms found in fastGEAR). The most likely donor genome in each event earns a "donor potential" of 1. In the event of a tie, we split the value to each strain equally (probability = 1/#-of-strains). Once "donor potential" from every event is determined, the sum of these constitutes the overall "donating potential" of a genome. We iterate through the remaining events with undetermined donors and award the event to the donor with the highest total "donating potential".

2) Calculates which recombining pairs (a donor and its recipient) have statistically high rates of recombination within the measured population. A standard *Highway* is defined as a pair whose number of recombination events are at least 1 standard deviation above the population's mean of recombination events per recombining pair.

3) Generates a number of files related to individual rates of recombination, as well as two files that can be directly drag-and-dropped onto a phylogenetic tree of the population to view the Highways and other recombining pairs (see example below).


We also include the supporting script **sidekick.py** which can be used to run fastGEAR on each individual gene-alignment following a successful run of Roary.

1) Using Roary's *clustered_proteins* file to find core genes in the pan-genome and extracting individual gene alignments from *pan_genome_sequences* directory (requires the use of -z command in Roary). 

2) Uses Roary's original input gff files to rename FASTA headers in each alignment file to the original gff name rather than the gene ID (important for proper visualization in HERO)

3) Runs fastGEAR on each renamed gene alignment (with optional command). Resulting directory is immediately ready for HERO


## INSTALLATION
HERO has the following dependencies:
- [NCBI-blast+](https://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastDocs&DOC_TYPE=Download)
- [fastGEAR](https://mostowylab.com/news/fastgear)
- Python 3.5
	- [Biopython](https://biopython.org/wiki/Download)
	- [Matplotlib](https://matplotlib.org/)
	- [Pandas](https://pandas.pydata.org/)
	- [Plotnine](https://plotnine.readthedocs.io/en/stable/)

Optional Dependencies:
- Roary

### Download command:
`git clone https://github.com/therealcooperpark/hero.git`

## USAGE

### hero.py

```
usage: hero.py [options] fastgear

HERO - Highways Elucidated by Recombination Observations

positional arguments:
  fastgear         Directory containing all fastGEAR runs for the population

optional arguments:
  -h, --help       show this help message and exit
  --log            Minimum log(BF) length to accept recombination event. [0]
  --length         Minimum fragment length for recombined sequence. [2]
  --pairs          Ignore direction of event when calculating recombination
                   highways.
  --filter FILTER  Comma delimited list of strains and filtering categories.
                   Will match events within a category, and between
  --inter          Do not add inter-group events to box-plots when using
                   --filter [False]
  --population     Manually set size of population for statistical tests in
                   donor/recipient calculations, defaults to number of
                   recombining strains
  --cleanup        Keep intermediate files.
  --cpus           Number of threads to use. [1]
  --format         Output bar graph file formats [png]
  --highway        Hexadecimal code for highway color in output file [#0000FF]
  --nohighway      Hexadecimal code for non-highway color in output file
                   [#D3D3D3]
  --output         Output directory for HERO files
  
```

### sidekick.py
```
usage: sidekick.py [options] gffs

Use before HERO to extract individual alignment files for each core gene in
pan-genome. Will also run fastGEAR individually on genes.

positional arguments:
  gffs                Filepath to GFF files

optional arguments:
  -h, --help          show this help message and exit
  --core              Compile all core genes found in --cluster_file [True]
  --accessory         Compile all accessory genes found in --cluster_file
                      present in at least 2 isolates, but fewer than
                      'core_definition' [False]
  --core_definition   Fraction of genomes that a gene must be found in to be
                      'core' [0.99]
  --alignments        Filepath to individual gene alignments,
                      [./pan_genome_sequences]
  --cluster_file      Filepath to the 'clustered_proteins' output of Roary
                      [./clustered_proteins]
  --fastGEAR          Run fastGEAR on every gene collected by sidekick [False]
  --output            Filepath to output directory for renamed genome protein
                      files. [./]
		      
```

## OUTPUT

- **Discriminant_Highways_hero_itol.txt / Discriminant_Recombination_hero_itol.txt ("Indiscriminant" if using --pairs)**

Visualization files for iToL.

- **donations.png / recipients.png**

Bar graphs showing donation/receipt counts for each genome.

- **donor_frequency.txt / recipient_frequency.txt**

Tab-delimited files showing raw data from .png files above

- **external_recomb_fragments.fa**

FASTA format of all recombination fragments that originated outside of the population.

Header format: gene_GenomeID:fragment-length

- **failed_genes.txt**

List of genes that had an unsuccessful fastGEAR run

- **filtered_recombination_events/**

A collection of revised output files for each gene showing the donor and recipient pairs and size of recombination fragment for each event

- **gene_frequency.txt**

Tab-delimited file showing number of recombination events in each gene

- **hero_statistics.txt**

Details about number of recombination pairs/events

- **recombination_pairs.txt**

Number of recombination events between each pair

- **recombination_sizes.txt**

List of all fragment sizes for every recombination event

- **recombined_portions.png**

Graph detailing number of recombined bases in each genome

- **recombined_base_pairs.txt**

Raw data for recombined_portions.png


## QUESTIONS
Feel free to email questions, concerns, suggestions and scathing reviews to cjp1043@wildcats.unh.edu

## CITATION
Park CJ, Andam CP. *Hero* **Github** https://github.com/therealcooperpark/hero

