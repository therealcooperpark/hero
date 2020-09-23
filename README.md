# HERO
**H**ighways **E**numerated by **R**ecombination **O**bservations

HERO is a pipeline, written in Python, designed to parse and visualize highways of genome-wide homologous recombination between user-defined metadata groups using the output of the recombination detection tool [fastGEAR](https://pubmed.ncbi.nlm.nih.gov/28199698/). Currently, it works in X stages:

1) For each recombination event, HERO comapres the sequence similarity between the recombined DNA sequence and a pool of potential donor genomes (determined by the clustering algorithm found in fastGEAR) to identify the most likely donor metadata-group (defined by the user).

2) Calculates which recombining metadata-group pairs (a donor and its recipient) have statistically high rates of recombination within the measured population. A standard *Highway* is defined using a strict Interquartile fence which equals *3\*IQR + Q3* where *IQR* is the interquartile range and *Q3* is the third quartile of the events per pair.

3) Generates a number of files detailing specific information on recombination events and metadata-group pairs, as well as several files related to visualizing the network of recombination using [Circos](http://circos.ca/)
  

Because fastGEAR is most effective at predicting recombination on a gene-by-gene basis (as opposed to a concatenated gene alignment), a standard use of fastGEAR involves generating a pan-genome of the population using a program such as [Roary](http://github.com/sanger-pathogens/Roary) and then running fastGEAR on each gene alignment. To accomodate this workflow, we also provide **sidekick.py** to prepare the data for HERO. Briefly, sidekick will:

1) Use the original GFF files provided to Roary to replace the headers in the Roary gene alignments with proper genome IDs to unify fastGEAR results across genes.

2) Iteratively run fastGEAR on each modified gene alignment with optional multithreading for speed.

3) Prepare the primary input file necessary for HERO to process the new fastGEAR data.


## INSTALLATION
HERO has the following dependencies:
- Python3.6
    - [BioPython](https://github.com/biopython/biopython)
    - [Pandas](https://pandas.pydata.org/docs/getting_started/install.html)
- [fastGEAR](https://pubmed.ncbi.nlm.nih.gov/28199698/)
- [Circos](http://circos.ca/)

### Download command:
`git clone https://github.com/therealcooperpark/hero.git`

## USAGE

### hero.py

```
usage: hero.py --hero_table [table] --groups [groups_file] [options]

HERO - Highways Elucidated by Recombination Observations

optional arguments:
  -h, --help                  show this help message and exit
  --hero_table HERO_TABLE     HERO input table
  --groups GROUPS             Tab-deliminated file with genomes in 1st column and groups in 2nd
  -o OUTDIR, --outdir OUTDIR  Output directory [hero_results]
  -c CPUS, --cpus CPUS        CPUs to use [1]
  -l LENGTH, --length LENGTH  Minimum length required to process recomb event [0]
  -b BAYES, --bayes BAYES     Minimum bayes factor required to process recomb event [10]
```

### sidekick.py
```
usage: sidekick.py [options] gff_table

Use before HERO to convert Roary FASTA alignment headers to genome names

positional arguments:
  gff_table             Tab-delimited file of GFF file location and associated
                        genome for renaming

optional arguments:
  -h, --help             show this help message and exit
  --alns ALNS            Filepath to Roary pan_genome_sequences directory
                         (requires -z argument) [./pan_genome_sequences]
  --fastgear             Run fastGEAR on each gene alignment [False]
  --output OUTPUT        Output directory name
  --cpus CPUS            Number of cpus to use [1]

fastGEAR:
  --fgout FGOUT          Output directory name for fastgear runs
  --iters ITERS          Number of iterations [15]
  --bounds BOUNDS        Upper bound for number of clusters [10]
  --partition PARTITION  File containing a partition for strains [NA]
  --fg_output FG_OUTPUT  1=reduced output, 0=complete output
```

## OUTPUT
- summary_stats.txt

Basic information about amount of recombination detected.

- recombination_events.txt

Tab-delimited table of each recombination event including donor group, recipient group, recombined fragment length (bp), gene name, and a list of recipient genomes with evidence for the event.

- recombination_pairs.txt

Tab-delmited table of all unique metadata group pairs and the number of recombination events between them.

- circos.png/svg

PNG and SVG formatted circos networks visualizing the network of recombination measured by HERO. See below for details on interpretting the figure.

- highway_circos.png/svg

PNG and SVG formatted circos networks highlighting highways of recombination detected by HERO.

- circos.conf/circos_karyotype.txt/circos_links.txt/highway_circos.conf

Configuration files created by HERO to create the circos plots.


## FUTURE IMPROVEMENTS
- Prepare recombination stats on individual recipients (ETA: ~1 week)
- Provide a repo dedicated to an example workflow involving HERO (ETA: ~1 week)
- Provide a conda installation option (ETA: ~2 weeks)


## QUESTIONS
Please submit suggestions and bug reports to the [Issue Tracker](https://github.com/therealcooperpark/hero/issues)


## CITATION
If you use this program, please cite:
Park C, Andam C. *HERO* **Github** [https://github.com/therealcooperpark/hero](https://github.com/therealcooperpark/hero)

