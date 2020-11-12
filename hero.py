#! /usr/bin/env python3
''' 
HERO - Highways Enumerated by Recombination Observations
Author - Cooper Park
'''


from argparse import ArgumentParser
from Bio.SeqIO import parse as BioParse
from itertools import product
import math
import multiprocessing
import os
import pandas as pd
from plotnine import *
from random import randint
import subprocess
import time

start_time = time.time()

def get_args():
    parser = ArgumentParser(description='HERO - Highways Elucidated by Recombination Observations',
                            usage='hero.py --hero_table [table] --groups [groups_file] [options]')
    parser.add_argument('--hero_table', required=True, help='HERO input table')
    parser.add_argument('--groups', required=True, help='Tab-deliminated file with genomes in 1st column and groups in 2nd')
    parser.add_argument('-o', '--outdir', default='hero_results', type=str, help='Output directory [hero_results]')
    parser.add_argument('-c', '--cpus', default=1, type=int, help='CPUs to use [1]')
    parser.add_argument('-l', '--length', default=0, type=int, help='Minimum length required to process recomb event [0]')
    parser.add_argument('-b', '--bayes', default=1, type=float, help='Minimum bayes factor required to process recomb event [1]')

    return parser.parse_args()


def parse_metadata(metadata_file):
    ''' Parse metadata into a dictionary '''

    groups = {} # Key = Genome [1st column], Value = Group [2nd column]
    try:
        with open(metadata_file, 'r') as metafile:
            for line in metafile:
                line = line.strip().split()
                groups[line[0]] = line[1].lower()
    except FileNotFoundError:
        print('Groups file {0} could not be opened. Ensure filepath is correct'.format(metadata_file))
        exit(1)
    return groups


def parse_table(hero_table):
    ''' 
    Parse HERO table into list of arguments 
    Sanity check all paths
    '''

    genes = []
    with open(hero_table, 'r') as infile:
        for line in infile:
            line = line.strip().split()
            if os.path.exists(line[1]) and os.path.exists(line[2]):
                genes.append(line)
            else:
                print('Gene {0} has a bad filepath. Skipping.'.format(line[0]), flush = True)
    return genes


def unpack_arguments(arg_list):
    ''' Unpack arguments and parse recombination '''

    return parse_fastgear(arg_list[0], arg_list[1], arg_list[2])


def parse_fastgear(gene_name, fasta_path, fastgear_path):
    t0 = time.time()

    ''' Parse recent recombination events from fastgear run '''

    # Find FASTA file to parse sequence info
    if not any(BioParse(fasta_path, 'fasta')):
        print('{0} fasta file is bad. Removing from analysis.'.format(gene_name), flush=True)
        return gene_name
   
    # Parse FASTA file into dict
    seqs_dict = {}
    for record in BioParse(fasta_path, 'fasta'):
        seqs_dict[record.id] = record.seq

    # Setup genome class
    class Genome:
        strain_to_genome = {} # Key: Strain name, Value: Genome class ID
        lineages = {} # Key: Lineage, Value: Strain name

        def __init__(self, sequence, lineage, name):
            self.name = name
            self.lineage = lineage
            self.sequence = sequence

            # Update class dictionaries
            Genome.lineages.setdefault(lineage, [])
            Genome.lineages[lineage].append(name)
            Genome.strain_to_genome[name] = self


    # Parse lineage file and update Genome Class Dicts
    try:
        with open('{0}/output/lineage_information.txt'.format(fastgear_path), 'r') as fg_file:
            next(fg_file)
            for line in fg_file:
                line = line.strip().split()
                try:
                    seq = seqs_dict[line[3]]
                except KeyError:
                    print('{0} could not match a sequence to ID {1}. Removing from analysis.'.format(fastgear_path, line[3]), flush=True)
                    return gene_name
                
                # Add genome to Genome class
                Genome(seq, line[1], line[3])

    except FileNotFoundError:
        return gene_name

    # Parse recombination
    return parse_recombination(fastgear_path, Genome, gene_name)


def parse_recombination(fastgear_run, Genome, gene_name):
    ''' Parse recent recombination and filter events '''

    def add_event(d_lineage, s_idx, e_idx, recipient):
        ''' Update pair with new event and condense overlapping events '''

        # Make sure pair exists
        pair = d_lineage
        donor_lineages.setdefault(pair, [])

        # Append new index to list of events
        donor_lineages[pair].append([s_idx, e_idx, [recipient]])

        # Then condense events by index pairs
        donor_lineages[pair].sort(key = lambda x: x[0])

        merged_pairs = [] # final array to hold merged intervals
        merged_recipients = []
        start = -1
        end = -1
        for idx in range(len(donor_lineages[pair])):
            cur_event = donor_lineages[pair][idx]
            if cur_event[0] > end:
                if idx != 0:
                    merged_pairs.append([start, end, merged_recipients])
                    merged_recipients = []
                end = cur_event[1]
                start = cur_event[0]
                merged_recipients.extend(cur_event[2])
            elif cur_event[1] >= end:
                end = cur_event[1]
                merged_recipients.extend(cur_event[2])
            
        if end != -1 and [start,end] not in merged_pairs:
           merged_pairs.append([start, end, merged_recipients])

        donor_lineages[pair] = merged_pairs 


    # Open recent recomb file
    try:
        recomb_file = open('{0}/output/recombinations_recent.txt'.format(fastgear_run), 'r')
        next(recomb_file)
        next(recomb_file)
    except FileNotFoundError:
        print('{0} has an incomplete fastgear run. Removing from analysis.'.format(fastgear_run), flush=True)
        return gene_name

    # Find external donor lineage num for gene for filtering
    external_donor = str(max([int(x) for x in Genome.lineages]) + 1)
    # Logs all lineage pairs and tracks unique events
    donor_lineages = {} # Key:   donor lineage
                        # Value: List of unique events

    # Get event info
    for line in recomb_file:
        line = line.strip().split()
        s_idx, e_idx = int(line[0])-1, int(line[1]) # fastGEAR includes s_idx in the sequence, so subtract one for indexing
        d_lineage, strain_name = line[2], line[5]
        logbf = float(line[4])

        # If minimum length or bayes not met, move on (length/bayes are global vars)
        # If donor lineage is external donor, move on
        fragment_len = e_idx - s_idx # fastGEAR includes the start position in its len
        if fragment_len < length or (math.e**logbf) < bayes or d_lineage == external_donor:
            continue

        # Add event to lineage pair in dict
        add_event(d_lineage, s_idx, e_idx, strain_name)
    recomb_file.close() # Close recomb file

    # For each unique event, find the most likely donor(s).
    # Then for each unique metadata group in recipients, log an event
    events = set() # All recombination events

    for d_lineage in donor_lineages:
        for event in donor_lineages[d_lineage]:
            start, end = int(event[0]), int(event[1])
            sample_recipient = event[2][0]
            # All genome are expected to be roughly equal. So take the first genome
            recip_seq = Genome.strain_to_genome[sample_recipient].sequence[start:end]
            donor_group = find_pair(start, end, d_lineage, recip_seq, Genome)

            # Fit donor group to all unique recip groups
            if donor_group:
                for recipient in event[2]:
                    recip_group = metadata[recipient]
                    recip_strains = [strain for strain in event[2] if metadata[strain] == recip_group]
                    #final_info = (donor_group, recip_group, end-start, gene_name, ','.join(recip_strains))
                    final_info = (donor_group, recip_group, start+1, end, gene_name, ','.join(recip_strains))
                    events.add(final_info)

    return list(events)


def find_pair(s_idx, e_idx, d_lineage, recip_seq, Genome):
    ''' Try to find a metadata pair that is linked by this recombination event '''

    # Step 1: See if all donors in d_lineage are from same metadata group
    #         NOTE:Lots of checking metadata groups here.
    #         I always default to other in case genome wasn't established in
    #         metadata parsing.
    
    # Test donors for total consistency
    donors = Genome.lineages[d_lineage]
    metadata.setdefault(donors[0], 'other')
    metagroup = metadata[donors[0]] # metadata dict is a global var

    for donor in donors[1:]:
        metadata.setdefault(donor, 'other')
        cur_group = metadata[donor]
        if cur_group != metagroup: # Not all donors are the same
            break
    else: # All donors from same group! We can move on.
        return metagroup 

    # Step 2: Not all donors fit
    #         Get distance of recip seq to donor recomb fragments
    
    # Get distance of recip seq to all donor seqs
    shortest = None
    viable_donors = []
    for donor in donors:
        donor_frag = str(Genome.strain_to_genome[donor].sequence[s_idx:e_idx])
        # Calculate distance between donor and recip fragment
        dist = 0
        for idx, nuc in enumerate(donor_frag):
            if recip_seq[idx] != nuc:
                dist += 1

        # Compare dist to current best dist
        if not shortest: # This is the first comparison
            shortest = dist
            viable_donors.append(donor)
            continue

        # All other tests
        if dist < shortest:
            shortest = dist
            viable_donors = [donor]
        elif dist == shortest:
            viable_donors.append(donor)


    # Step 3 (2b?): If all likely donors from same metagroup, we win.
    #               Otherwise, discard the event.

    metagroup = metadata[viable_donors[0]]
    if len(viable_donors) > 1: # If multiple donors, check for consistency
        for donor in viable_donors[1:]:
            if metadata[donor] != metagroup:
                return None # If two metagroups exist, kill the search

    # We found a good metagroup! Send the event back
    return metagroup


def parse_events(recombination_events):
    ''' Parse events from multithreaded event finding '''

    good_events = []
    bad_genes = []
    for gene in recombination_events:
        if isinstance(gene, list): # Good events are lists, bad genes are str
            for event in gene:
                good_events.append(event)
        else:
            bad_genes.append(gene)
    return good_events, bad_genes


def calculate_highway(events, unique_groups):
    '''
    Calculate the theshold for highways of recombination
    highway = 3*IQR + Q3
    IQR = Interquartile range
    Q3  = Third quartile of the data
    '''

    recomb_events = {x:0 for x in unique_groups}
    # Get all unique combinations of group pairs
    for event in events:
        pair = (event[0], event[1])
        recomb_events.setdefault(pair, 0)
        recomb_events[pair] += 1

    # Calculate IQRs
    recomb_counts = list(recomb_events.values())
    recomb_df = pd.DataFrame({'Events': recomb_counts})
    q3 = recomb_df.quantile(q=0.75)['Events']
    q1 = recomb_df.quantile(q=0.25)['Events']
    IQR = q3 - q1
    significance_limit = q3 + (3*IQR)
    return recomb_events, significance_limit


class Metagroup:
    '''
    Each metadata group will be given an instance.
    Tracks recombination stats for each group
    '''

    metagroup_dict = {} # Key:   metagroup string name
                        # Value: Metagroup object instance

    def __init__(self, name):
        # Recombination variables
        self.name        = name
        self.donations   = 0    # Number of donations
        self.receipts    = 0    # Number of receipts
        self.group_stats = {}   # Key:   Other metagroup object string name
                                # Value: [donations_to, receipts_from]

        # Plotting variables
        self.d_pos       = 0    # Number of donations already plotted
        self.r_pos       = 0    # Number of receipts already plotted

    def total_events(self):
        return self.donations + self.receipts

    def add_event(event):
        '''
        Parse event to add donor and recipt credit
        to each metagroup in event
        '''

        donor, recipient = event[0], event[1]

        # Make object instance for each group if not already exists
        Metagroup.metagroup_dict.setdefault(donor, Metagroup(donor))
        Metagroup.metagroup_dict.setdefault(recipient, Metagroup(recipient))

        # Add donor/recipient credit to each group
        d_group = Metagroup.metagroup_dict[donor]
        r_group = Metagroup.metagroup_dict[recipient]
        d_group.donations += 1
        r_group.receipts  += 1
        d_group.group_stats.setdefault(recipient, [0, 0])
        r_group.group_stats.setdefault(donor, [0, 0])
        d_group.group_stats[recipient][0] += 1 # Add donor credit
        r_group.group_stats[donor][1] += 1     # Add recip credit


def make_circos(events, outdir):
    ''' Write circos files given events and list of genomes w/ metadata '''

    # Log all events in Metagroup class
    for event in events:
        Metagroup.add_event(event)

    # Write karyotype file for circos
    with open('{0}/circos_karyotype.txt'.format(outdir), 'w') as k_file:
        # Get random color for each group chunk
        rand_colors = random_colors(len(Metagroup.metagroup_dict.keys()))
        # Write color and group to karyotype file
        for idx, group in enumerate(Metagroup.metagroup_dict.values()):
            color = rand_colors[idx]
            k_file.write('chr - {0} {0} 0 {1} {0}\n'.format(group.name.lower(), group.total_events()))

    # Write link file
    with open('{0}/circos_links.txt'.format(outdir), 'w') as l_file:
        # Create links by the donor
        for d_group in Metagroup.metagroup_dict.values():
            donor = d_group.name

            # Get recipient from group_stats variable
            # If donor is in the list of recipients,
              # Put it on the end so it looks cleaner
            recipients = list(d_group.group_stats.keys())
            recipients.sort(key=donor.__eq__)
            for recipient in d_group.group_stats:
                donations = d_group.group_stats[recipient][0]
                r_group = Metagroup.metagroup_dict[recipient]

                ## Write link to file
                # Get donor plot range and update donor positions
                d_start = d_group.d_pos
                d_end   = d_start + donations
                d_group.d_pos += donations

                # Get recipient range and update recipient positions
                # All receipts should be plotted away from donations
                r_start = r_group.donations + r_group.r_pos
                r_end   = r_start + donations
                r_group.r_pos += donations

                # Write to file
                link  = donor + ' ' + str(d_start) + ' ' + str(d_end) + ' '
                link += recipient + ' ' + str(r_start) + ' ' + str(r_end) + '\n'
                l_file.write(link)

    # Write config_file
    # Tutorial to understanding circos config file can be found at:
    # circos.ca/documentation/tutorials/quick_start/
    with open('{0}/circos.conf'.format(outdir), 'w') as c_file:

        file_contents  = 'karyotype = {0}/circos_karyotype.txt\n'.format(outdir)

        # Global color scheme
        file_contents += '# Global color scheme\n'
        file_contents += '<colors>\n'
        for idx, name in enumerate(Metagroup.metagroup_dict.keys()):
            file_contents += '{0}* = {1}\n'.format(name, rand_colors[idx])
        file_contents += '</colors>\n'

        # Basic required content (karyotype file location, ideogram creation)
        file_contents += '<ideogram>\n\n<spacing>\n'
        file_contents += 'default = 0.005r # Spacing between out ring chunks\n'
        file_contents += '</spacing>\n\n'

        # Ideogram layout details
        file_contents += '# Ideogram layout details\n'
        file_contents += 'radius            = 0.9r # Size of radius for outer ring\n'
        file_contents += 'thickness         = 80p # Thickness of outer ring\n'
        file_contents += 'fill              = yes # Fill chunks with color?\n'
        file_contents += 'stroke_color      = dgrey # Color of chunk outline\n'
        file_contents += 'stroke_thickness  = 2p # Thickness of outline\n\n'

        # Ideogram label details
        file_contents += '# Ideogram label details\n'
        file_contents += 'show_label     = yes      # Show chunk labels?\n'
        file_contents += 'label_font     = default  # Font of the labels\n'
        file_contents += 'label_radius   = 1r + 75p # Where to place labels\n'
        file_contents += 'label_size     = 50       # Size of the label\n'
        file_contents += 'label_parallel = yes      # Set label parallel to chunks\n'
        file_contents += '</ideogram>\n\n'

        # Tick details
        # << SKIPPED FOR NOW >>

        # Link details
        file_contents += '# Links... The actual connections\n'
        file_contents += '<links>\n<link>\n'
        file_contents += 'file          = {0}/circos_links.txt # The file with links to draw\n'.format(outdir)
        file_contents += 'ribbon        = yes              # Turn links into fancy ribbons\n'
        file_contents += 'flat          = yes              # Flatten ribbons\n'
        file_contents += 'z             = 1                # importance for ribbon plotting\n'
        file_contents += 'radius1       = 0.8r             # Push donor end of ribbon inward\n'
        file_contents += 'color         = eval(var(chr2))  # Default link color\n'
        file_contents += 'radius        = 0.98r            # Where links will stop at\n'
        file_contents += 'bezier_radius = 0.1r             # How far from center the curves are drawn\n'
        file_contents += 'thickness     = 5                # Default thickness\n'

        # Establish rule to color links by donor chunk
        file_contents += '\n<rules>\n'
        file_contents += '\nflow = continue\n\n'
        file_contents += '<rule>\n'
        file_contents += 'condition = 1\n'
        file_contents += 'color     = eval(var(chr1))\n'
        file_contents += '</rule>\n<rule>\n'
        file_contents += 'condition = var(interchr)\n'
        file_contents += 'z         = 2\n'
        file_contents += '</rule>\n'
        file_contents += '</rules>\n\n'
        file_contents += '</link>\n</links>\n\n'

        # Default circos distributions to include
        file_contents += '# Default circos distributions to include\n'
        file_contents += '<image>\n<<include etc/image.conf>>\n</image>\n'
        file_contents += '<<include etc/colors_fonts_patterns.conf>>\n'
        file_contents += '<<include etc/housekeeping.conf>>\n'

        c_file.write(file_contents)


def make_highway_circos(highway, outdir):
    '''
    Create 2nd circos.conf file which filters the color
    of ribbons below the highway_definition threshold
    '''

    try:
        with open('{0}/circos.conf'.format(outdir), 'r') as circos_file, open('{0}/highway_circos.conf'.format(outdir), 'w') as outfile:
            for line in circos_file:
                if line == '</rules>\n':
                    outfile.write('<rule>\n')
                    outfile.write('condition = (var(end1) - var(start1)) < {0}\n'.format(highway))
                    outfile.write('color     = grey\n')
                    outfile.write('z         = 1\n')
                    outfile.write('</rule>\n')
                outfile.write(line)
    except IOError:
        print('Could not make highway circos file. Check circos.conf', flush=True)


def random_colors(num_colors):
    ''' Generate num_colors random colors '''
    # Current optimum maximum number of groups: 51 (255//5)

    colors = {k:[] for k in 'rgb'} # Dict of all R/G/B values
    for color in range(num_colors): # Make each color
        temp = {k: randint(0,255) for k in 'rgb'} # Get random RBG values
        for k in temp:
            # For each value, make sure it is at least 25 points
            # different from all other values in same position
            while True:
                c = temp[k]
                t = set(j for j in range(c-5, c+5) if 0 <= j <= 255)
                if t.intersection(colors[k]):
                    temp[k] = randint(0,255)
                else:
                    break
            colors[k].append(temp[k])
    # Format final colors
    final_colors = []
    for i in range(num_colors):
        final_colors.append( '{0},{1},{2}'.format(colors['r'][i], colors['g'][i], colors['b'][i]))

    return final_colors


def write_individual_stats(outdir, events):
    '''
    Write useful text files and plots for individual genome recomb data
    1) Histogram of recombination fragment sizes
    2) Histogram of recombination per gene
    3) Histogram of recombination per recipient
    '''

    # Step 1: Write out fragment data and collect gene/recipient data
    fragments  = open('{0}/fragment_sizes.txt'.format(outdir), 'w')
    recipient_counts = {}
    gene_counts      = {}
    fragments.write('Size\n')
    for event in events:
        # Write out fragment now
        fragments.write(str(event[3] - event[2])+'\n')

        # Add 1 to the count for the gene
        gene_counts.setdefault(event[4], 0)
        gene_counts[event[4]] += 1

        # Each genome gets 1 to its recipient count
        for genome in event[5].split(','):
            recipient_counts.setdefault(genome, 0)
            recipient_counts[genome] += 1
    fragments.close()

    # Write out recipient/gene data
    genes      = open('{0}/gene_counts.txt'.format(outdir), 'w')
    genes.write('Gene\tEvents\n')
    for gene, count in gene_counts.items():
        genes.write('{0}\t{1}\n'.format(str(gene), str(count)))
    genes.close()

    recipients = open('{0}/recipient_counts.txt'.format(outdir), 'w')
    recipients.write('Recipient\tEvents\n')
    for r, count in recipient_counts.items():
        recipients.write('{0}\t{1}\n'.format(str(r), str(count)))
    recipients.close()    

    # Step 2: Make each histogram
    make_histogram('{0}/gene_counts.txt'.format(outdir), 'gene', '{0}/gene_counts'.format(outdir))
    make_histogram('{0}/recipient_counts.txt'.format(outdir), 'recipient', '{0}/recipient_counts'.format(outdir))
    make_histogram('{0}/fragment_sizes.txt'.format(outdir), 'fragment', '{0}/fragment_sizes'.format(outdir))


def make_histogram(file_loc, plot_type, filename):
    '''
    Make a histogram given a file location and plot type
    '''

    # Load in each filetype properly
    if plot_type == 'gene':
        datas     = pd.read_csv(file_loc, header=0, sep='\t')
        x_lab     = '# of events per gene'
        histogram = (ggplot(datas, aes(x='Events')) 
                    + geom_histogram() + xlab(x_lab))
    elif plot_type == 'recipient':
        datas     = pd.read_csv(file_loc, header=0, sep='\t')
        x_lab     = '# of events per recipient'
        histogram = (ggplot(datas, aes(x='Events'))
                    + geom_histogram() + xlab(x_lab))
    elif plot_type == 'fragment':
        datas     = pd.read_csv(file_loc, header=0)
        x_lab     = 'Recombination fragment size (bp)'
        histogram = (ggplot(datas, aes(x='Size'))
                    + geom_histogram() + xlab(x_lab))
    else:
        datas = pd.read.csv(file_loc, header=1, sep='\t')
        x_lab = plot_type

    # Make histogram
    histogram += ylab('Frequency') 
    histogram += theme(panel_background=element_blank(),
                       axis_line=element_line(),
                       axis_text=element_text(size=20),
                       axis_title=element_text(size=25))
    # Save histogram
    histogram.save(filename=filename+'.svg')


#################
## Main Script ##
#################

# Get args and set some as global
args = get_args()
global length
global bayes
length = args.length
bayes = args.bayes

# Make new directory
try:
    os.mkdir(args.outdir)
except FileExistsError:
    dir_ext = str(time.time()).split('.')[0]
    new_dir = args.outdir+'_'+dir_ext
    args.outdir = new_dir
    os.mkdir(args.outdir)
    print('Directory already exists. Output directory will be {0}'.format(args.outdir), flush=True)

# Parse metadata from groups file
metadata = parse_metadata(args.groups)
print('Made it past metadata parsing!', flush = True)

## DEBUG MODE FOR FASTGEAR PARSING
#recombination_events = []
#for gene in args.fastgears:
#    print('Started {0}'.format(gene), flush=True)
#    recombination_events.append(parse_fastgear(gene))
#    print('Finished {0}'.format(gene), flush=True)

# Parse fastGEAR data from each run
genes = parse_table(args.hero_table)
pool = multiprocessing.Pool(processes=args.cpus)
recombination_events = pool.map(unpack_arguments, genes) 
pool.close()
pool.join()
print('Parsed all the genes!', flush = True)

# Filter good events out of all the data
events, genes = parse_events(recombination_events)
print('Got all the events!', flush=True)

# Plot results as circos plot
make_circos(events, args.outdir)
subprocess.run('circos --conf {0}/circos.conf -outputdir {0}'.format(args.outdir), shell=True)

# Find highways
unique_groups = set(product(Metagroup.metagroup_dict.keys(),repeat=2))
recomb_by_pair, highway_value = calculate_highway(events, unique_groups)
make_highway_circos(highway_value, args.outdir)
subprocess.run('circos --conf {0}/highway_circos.conf -outputdir {0} -outputfile highway_circos'.format(args.outdir), shell=True)


## Write some stats files ##
# Write out raw recombination event data
with open('{0}/recombination_events.txt'.format(args.outdir), 'w') as outfile:
    outfile.write('Donor_Group\tRecip_Group\tstart\tend\tGene\tRecipient_strains\n')
    for event in events:
        event = [str(x) for x in event]
        outfile.write('\t'.join(event) + '\n')

# Write recombinations per pair
with open('{0}/recombination_pairs.txt'.format(args.outdir), 'w') as outfile:
    outfile.write('Donor_Group\tRecipient_Group\tCount\n')
    for pair in recomb_by_pair:
        outfile.write('{0}\t{1}\t{2}\n'.format(pair[0], pair[1], recomb_by_pair[pair]))

# Write out summary statistics about recombination
with open('{0}/summary_stats.txt'.format(args.outdir), 'w') as outfile:
    outfile.write('Minimum recombinations for a highway: {0}\n'.format(highway_value))
    highways = [x for x in recomb_by_pair.values() if x > highway_value]
    outfile.write('Number of predicted events: {0}\n'.format(len(events)))
    outfile.write('Number of events from highways: {0}\n'.format(sum(highways)))
    outfile.write('Number of recombining pairs: {0}\n'.format(len(recomb_by_pair)))
    outfile.write('Number of highway pairs: {0}\n'.format(len(highways)))
    unique_genes = set()
    for event in events:
        unique_genes.add(event[4])
    outfile.write('Number of genes with recombination: {0}\n'.format(len(unique_genes)))

write_individual_stats(args.outdir, events)
print('Run time: ', time.time() - start_time, 'seconds')
