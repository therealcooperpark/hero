#! /usr/bin/env python3

'''
HERO: Highway Enumeration from Recombination Observations
'''


import warnings
warnings.filterwarnings('ignore')
import argparse
import Bio.SeqIO
import itertools
import math
import multiprocessing
from operator import itemgetter
import os
import pandas
from plotnine import *
import random
import shutil
import statistics
import string
import subprocess
import time

t0 = time.time()

# Args
parser = argparse.ArgumentParser(description = "HERO - Highways Elucidated by Recombination Observations", usage = "hero.py [options] fastgear")
parser.add_argument("fastgear", help = "Directory containing all fastGEAR runs for the population")
parser.add_argument("--log", metavar = '', default = 0, type = int, help = "Minimum log(BF) length to accept recombination event. [0]")
parser.add_argument("--length", metavar = '', default = 2, type = int, help = "Minimum fragment length for recombined sequence. [2]")
parser.add_argument("--pairs", action = "store_true", help = "Ignore direction of event when calculating recombination highways.")
parser.add_argument("--filter", help = "Comma delimited list of strains and filtering categories. Will match events within a category, and between")
parser.add_argument("--inter", action = "store_true", help = "Do not add inter-group events to box-plots when using --filter [False]")
parser.add_argument("--population", metavar = '', default = 0, type = int, help = "Manually set size of population for statistical tests in donor/recipient calculations, defaults to number of recombining strains")
parser.add_argument("--cleanup", action = "store_false", help = "Keep intermediate files.")
parser.add_argument("--cpus", metavar = '', default = 1, type = int, help = "Number of threads to use. [1]")
parser.add_argument("--format", metavar = '', default = "png", help = "Output bar graph file formats [png]")
parser.add_argument("--highway", default = "#0000FF", help = "Hexadecimal code for highway color in output file")
parser.add_argument("--nohighway", metavar = '', default = "#D3D3D3", help = "Hexadecimal code for non-highway color in output file [#D3D3D3]")
parser.add_argument("--output", metavar = '', default = "HERO", help = "Output directory for HERO files")
args = parser.parse_args()


def make_directory(name):
    # Make output directory for core genome alignments, rename if already in use
    try:
        os.mkdir(name)
        out = name
    except FileExistsError:
        ext = ''.join(random.choice(string.ascii_uppercase + string.digits) for _ in range(6))
        out = name + '_' + ext
        os.mkdir(out)
        print("Output directory already found, writing to {0}".format(out), flush=True)

    return out

def lineage_fastas(fastgear):
    # Check absolute vs relative path for appropriate pathways to gene fastGEARs
    if not fastgear.is_dir:
        print("Non-fastGEAR directory found. Skipping", flush=True)
        return fastgear.name

    try:
        with open("{0}/output/recombinations_recent.txt".format(fastgear.path), "r") as file:
            if len(file.readlines()) < 3:
                return fastgear.name
    except IOError:
        return fastgear.name

    # Setup main collections
    class Strains:
        alignment_length = 0 # Total alignment length will be set with each new gene processed

        def __init__(self, lineage, name):
            self.name = name
            self.lineage = lineage
            self.fasta = ""
            self.trim = 0
            self.unaligned_fasta = ""


            lineages.setdefault(lineage, [])
            query_lineages.setdefault(lineage, [])
            lineages[lineage].append(name)



    lineages = {} # Sort IDs into lineage
    fg_strains = {} # Collection of all strains in class Strains.
    query_lineages = {} # Sequences for query of blast, sorted by lineage
    # Parse and sort strain alignments by lineage, build blast DB for each file

    try:
        with open("{0}/output/lineage_information.txt".format(fastgear.path), "r") as file:
            next(file)
            for line in file:
                line = line.split()
                fg_strains[line[3]] = Strains(line[1], line[3])
    except FileNotFoundError:
        print("Gene {0} has an unsuccessful fastGEAR run. Skipping for results.".format(fastgear.name), flush=True)
        return fastgear.name

    # Identify alignment file
    alignment = ""
    for file in os.scandir(fastgear.path):
        if file.is_file() and ".mat" not in file.name:
            try:
                if any(Bio.SeqIO.parse(file.path, "fasta")):
                    alignment = file.path
            except UnicodeDecodeError:
                pass
    if alignment == "":
        print("Gene {0} has no alignment file present. Skipping for results.".format(fastgear.name), flush=True)
        return fastgear.name

    # Generate lineage based fasta files
    for record in Bio.SeqIO.parse(alignment, "fasta"):
        fg_strains[record.id].fasta = record
        fg_strains[record.id].unaligned_fasta = record
        fg_strains[record.id].unaligned_fasta.seq = unalignment(fg_strains[record.id].fasta.seq)

    # Write out reference lineage files
    os.mkdir("{0}/blast_files/{1}".format(args.output, fastgear.name))
    for num in lineages:
        records = []
        for strain in lineages[num]:
            records.append(fg_strains[strain].unaligned_fasta)

        Bio.SeqIO.write(records, "{0}/blast_files/{1}/lineage_{2}.fa".format(args.output, fastgear.name, num), "fasta")
    
        blastdb = subprocess.run("makeblastdb -in {0}/blast_files/{1}/lineage_{2}.fa -dbtype nucl".format(args.output, fastgear.name, num), shell = True, stdout = subprocess.PIPE, stderr = subprocess.PIPE)
        if blastdb.stderr:
            print(blastdb.stderr.decode('ascii'), flush=True)
            print("Revisit the gene alignment of {0} to troubleshoot. Skipping for results.".format(fastgear.name), flush=True)
            return fastgear.name

    return recomb_parser(fastgear, Strains, fg_strains, query_lineages, max(lineages.keys()))


def recomb_parser(fastgear, Strains, fg_strains, query_lineages, max_lineage):
    # Parse recombination events and build query blast files
    external = []
    overlap_check = {}
    overlap_removes = {}
    with open("{0}/output/recombinations_recent.txt".format(fastgear.path), "r") as recomb_file:
        next(recomb_file)
        next(recomb_file)

        # Check for base pair overlap in recipients from the same donor lineage
        skip_lines = set()
        line_num = 0
        possible_lines = {}
        for line in recomb_file:
            line = line.split()
            possible_lines.setdefault(line[2], [])

            # Find overlap
            for pos, idx in enumerate(possible_lines[line[2]]):
                new_idx = seq_overlap(int(idx[0]), int(idx[1]), int(line[0]), int(line[1]))
                if new_idx:
                    possible_lines[line[2]][pos] = new_idx
                    break
            else:
                possible_lines[line[2]].append((line[0], line[1]))
            line_num += 1

        # Identify overlap within the same recipient
        recomb_file.seek(0,0)
        next(recomb_file)
        next(recomb_file)
        line_num = 0
        for line in recomb_file:
            line = line.split()
            
            if float(line[4]) < args.log or line_num in skip_lines:
                continue
            line_num += 1
            donor, recipient, start_idx, end_idx  = line[2], fg_strains[line[5]], int(line[0]) - 1, int(line[1]) - 1
            # Start dictionary for donor lineage and recipient
            overlap_check.setdefault(donor, {})
            overlap_check[donor].setdefault(recipient, [])
            # Check for overlap event in each index
            for pos,pair in enumerate(list(overlap_check[donor][recipient])):
                new_idx = seq_overlap(start_idx, end_idx, pair[0], pair[1])
                if new_idx != 0: # If overlap, mark index of replacement
                    overlap_removes[donor][recipient][pos] = new_idx
                    break
            else: # Replace old index if overlap, otherwise add index as source for next
                overlap_check[donor][recipient].append((start_idx, end_idx))

    # For each unique recombination fragment, pull sequence, verify length meets minimum and set Bio.SeqIO record for event. Add this to query_lineages for blasting
    for donor in overlap_check:
        for recipient in overlap_check[donor]:
            for fragment in overlap_check[donor][recipient]:
                sequence = unalignment(recipient.fasta.seq[fragment[0]:fragment[1]])
                if len(sequence) <= args.length:
                    continue
                identity = "{0}:{1}".format(recipient.name, len(sequence)) # Recipient name and length of event

                if int(donor) > int(max_lineage): # Add event to external recombination or within population recombination
                    rec = Bio.SeqRecord.SeqRecord(id = "{0}_{1}".format(fastgear.name, identity), seq = sequence)
                    external.append(rec)
                else:
                    rec = Bio.SeqRecord.SeqRecord(id = identity, seq = sequence)
                    query_lineages[donor].append(rec)

    rec_indexes = {}
    recombination_events = []
    for query in query_lineages:
        filename = "{0}/blast_files/{1}/lineage_{2}_query.fa".format(args.output, fastgear.name, query)
        Bio.SeqIO.write(query_lineages[query], filename, "fasta")
        blast_output = subprocess.run("blastn -query {0} -subject {1}/blast_files/{2}/lineage_{3}.fa -evalue 1e-5 -outfmt 7 -max_hsps 1".format(filename, args.output, fastgear.name, query),\
            shell = True, stdout = subprocess.PIPE, stderr = subprocess.PIPE)

        if not args.cleanup: # Keeps blast output in the event of keeping intermediate files
            blast_file = open("{0}/blast_files/{1}/query_{2}.blast_output".format(args.output, fastgear.name, query), "w")

        # Parse blast_output
        for entry in blast_output.stdout.decode('ascii').split("\n"):

            if not args.cleanup:
                blast_file.write("{0}\n".format(entry))

            if "#" in entry: # Skip headers, but recognize that a new event is being processed
                potential_donors = []
                get_event, e_value, percent_identity  = 1, None, None
                continue

            if entry == "" or get_event == 0: # Checks for empty string (due to parsing) and if match has been found already
                continue

            entry = entry.split()
            
            if not e_value: # If first donor in list, add to list and move on
                e_value = entry[10]
                percent_identity = entry[2]
                potential_donors.append(entry[1])
            elif e_value == entry[10] and percent_identity == entry[2]: # Otherwise, check quality of next donor until they do not match
                potential_donors.append(entry[1])
            else: # When they don't match, use that list of potential donors for the event
                r = entry[0].split(":")
                recombination_events.append((fastgear.name, potential_donors, r[0], r[1]))
                get_event = 0                

        if not args.cleanup:
            blast_file.close()

    return [recombination_events, external]

def unalignment(seq):
    # Remove gap characters from sequence alignment
    new_seq = ""
    for nuc in seq:
        if nuc != "-":
            new_seq += nuc
    return Bio.Seq.Seq(new_seq)


def seq_overlap(start, end, a, b):
    # Check for overlap
    if a <= start and b >= start:
        if b > end:
            return (a,b)
        else:
            return (a,end)
    elif a >= start and a <= end:
        if b > end:
            return (start,b)
        else:
            return (start,end)
    else:
        overlap = 0
    
    return overlap
    

def cleanup(fastgear):
    # Clean intermediate files
    for file in os.scandir(fastgear.path):
        if "lineage" in file.name:
            os.remove(file.path)


def write_out_events(unique_recombination):
    # Write out predicted events as new fastGEAR styled files

    # Sort events by gene
    gene_sorted_dict = {}
    for pair in unique_recombination:
        p2 = pair.split(":")
        for event in unique_recombination[pair]:
            gene_sorted_dict.setdefault(event[1], [])
            gene_sorted_dict[event[1]].append((p2[0], p2[1], event[0]))

    # Write out events per gene
    for gene in gene_sorted_dict:
        with open("{0}/filtered_recombination_events_{1}bp_{2}logBF/{3}_events.txt".format(args.output, args.length, args.log, gene), "w") as output:
            output.write("Donor\tRecipient\tFragment Length\n")
            for event in gene_sorted_dict[gene]:
                output.write("{0}\t{1}\t{2}\n".format(event[0], event[1], event[2]))


def find_pairs(good_genes):
    # Get donor probability and assign pairs of recombination

    ### Output dictionaries
    unique_recombination = {} # Donor_recipient pairs and their recombination length/gene origin
    paired_recombination = {} # Donor_recipient pairs without direction
    recombination_fragments = {}

    ### Function dictionaries
    donors = {} # Donor probabilities 
    non_unique_events = [] # All events that have multiple possible donors

    for event in good_genes:
        if len(event[1]) == 1: # If only one donor, set pair as donor_recipient
            d,r = event[1][0], event[2]
            d_r = "{0}:{1}".format(event[1][0], event[2])
            paired_d_r = find_paired_events([d,r], d_r, paired_recombination)

            donors.setdefault(d, 0) # If there is a single donor, its probability for the event is one.
            donors[d] += 1

            unique_recombination.setdefault(d_r, [])
            unique_recombination[d_r].append((event[3], event[0])) # Append tuple of event length and gene source

            paired_recombination.setdefault(paired_d_r, [])
            paired_recombination[paired_d_r].append((event[3], event[0])) # Do the same for paired keys            

            recombination_fragments.setdefault(r, 0) # Set default for recipient and add fragment length to total bp recombined
            recombination_fragments[r] += int(event[3])
            continue

        probability = 1 / len(event[1])
        for donor in event[1]: # Multiple donors have probabilities added to their total
            donors.setdefault(donor, 0)
            donors[donor] += probability
        non_unique_events.append(event)

    # Sort donors by probability and collect their recombination events high to low
    prob_sorted_donors = sorted(list(donors.items()), key=lambda x: x[1], reverse = True)
    ancestral_events = 0

    next_idx = 0
    for idx, donor in enumerate(prob_sorted_donors):
        if idx < next_idx: # This will skip iterations when equally_likely donors have been found
            continue

        equally_likely = set()
        cur = idx
        while prob_sorted_donors[cur][1] == donor[1]: # Get all donors with equal likelihood
            equally_likely.add(donor[0])
            if cur + 2 > len(prob_sorted_donors): # +1 to counter index offset, +1 to go over list length appropriatelybz
                break
            cur += 1
        next_idx = cur # Push next_idx up to skip any equally likely donors

        for idx, event in enumerate(list(non_unique_events)): # Iterate every dictionary, check if donor present, act accordingly.
            overlap = equally_likely & set(event[1])
            if len(overlap) > 1:
                ancestral_events += 1
                non_unique_events.remove(event)
            elif len(overlap) == 1:
                potential_donor = list(overlap)[0]
                key = "{0}:{1}".format(potential_donor, event[2])
                paired_key = find_paired_events(key.split(":"), key, paired_recombination)

                unique_recombination.setdefault(key, [])
                unique_recombination[key].append((event[3], event[0])) # Append tuple of event length and gene source

                paired_recombination.setdefault(paired_key, [])
                paired_recombination[paired_key].append((event[3], event[0])) # Do the same for paired keys

                recombination_fragments.setdefault(event[2], 0) # Set default for recipient and add fragment length to total bp recombined
                recombination_fragments[event[2]] += int(event[3])

                non_unique_events.remove(event)
    print("{0} Leftover events".format(len(non_unique_events)), flush = True)
    return (unique_recombination, paired_recombination, recombination_fragments, ancestral_events)


def find_paired_events(strains, original_key, paired_recombination):
    # Find if iteration exists in dict already
    valid = 1
    if args.pairs:
        for iteration in itertools.permutations(strains, 2):
            if ":".join(iteration) in paired_recombination:
                paired_key = ":".join(iteration)
                valid = 0
                break

    if valid == 1:
        return original_key
    else:
        return paired_key


def get_filter_groups():
    # Make filtering dictionary for strains if user provides it
    strain_filters = {} # Strain level filtering. Key = strain, value = filter group
    with open(args.filter, "r") as filter_file:
        for line in filter_file:
            line = line.strip().split(",")
            strain_filters[line[0]] = line[1]
    return strain_filters


def filter_recombination(pairs, strain_filters):
    # Sort recombination events into user-defined categories
    filtered_recombination = {"inter": {}}

    for pair in pairs:
        split = pair.split(":")
        # Get category, set to "inter" if not found
        group1, group2 = strain_filters.setdefault(split[0], "inter"), strain_filters.setdefault(split[1], "inter")
        output_group = group1 if group1 == group2 else "inter"

        # Add to filtered_recombination
        filtered_recombination.setdefault(output_group, {})
        filtered_recombination[output_group].setdefault(pair, [])
        filtered_recombination[output_group][pair] += pairs[pair]

    return filtered_recombination


def write_out(group, recombination_dictionary):
    # Calculate Highway definition
    if len(group) == 0:
        print("{0} category has no recombination events. No output files to write".format(label), flush=True)
        return None
    # Sort recombination data
    recomb_nums = {}
    recomb_sizes = []
    recomb_genes = {}
    for pair in recombination_dictionary:
        recomb_nums.setdefault(pair, 0)
        for data in recombination_dictionary[pair]:
            recomb_nums[pair] += 1

            recomb_sizes.append(data[0])

            recomb_genes.setdefault(data[1], 0)
            recomb_genes[data[1]] += 1

    one_stdev_higher = round((sum(recomb_nums.values()) / len(recomb_nums.values())) + statistics.pstdev(recomb_nums.values()), 2)

    # Write output files
    names = name_collection(group)
    ## ^pair_file, stat_file, d_file, r_file, gene_file, size_file

    stats_writer(recomb_nums, names[0], names[1], one_stdev_higher)
    itol_out(recomb_nums, group, one_stdev_higher)
    gene_out(names[4], recomb_genes)
    size_out(names[5], recomb_sizes)
    if not args.pairs:
        filtered_stats[group] = frequent_finder(recomb_nums, names[2], names[3])
 

def name_collection(label):
    # Assign file names for each output function call
    if label == "inter":
        pair_file = "recombination_pairs.txt"
        stat_file = "hero_statistics.txt"
        d_file = "donor_frequency.txt"
        r_file = "recipient_frequency.txt"
        gene_file = "gene_frequency.txt"
        size_file = "recombination_sizes.txt"
    else:
        pair_file = "{0}_recombination_pairs.txt".format(label)
        stat_file = "{0}_hero_statistics.txt".format(label)
        d_file = "{0}_donor_frequency.txt".format(label)
        r_file = "{0}_recipient_frequency.txt".format(label)
        gene_file = "{0}_gene_frequency.txt".format(label)
        size_file = "{0}_recombination_sizes.txt".format(label)
    return (pair_file, stat_file, d_file, r_file, gene_file, size_file)


def stats_writer(recomb_dict, pair_file, stat_file, one_stdev_higher):
    # Write out recombination pairs and quantities + stats
    recombination = 0
    pair_num = 0
    highway_recomb = 0
    highway_num = 0

    print("Writing out results to {0}...".format(pair_file), flush=True)
    with open("{0}/{1}".format(args.output, pair_file), "w") as hero:
        if args.pairs:
            hero.write("Pair\tEvents\n")
        else:
            hero.write("Donor:Recipient\tEvents\n")

        for event in recomb_dict:
            recomb_num = int(recomb_dict[event])
            hero.write("{0}\t{1}\n".format(event, recomb_num))

            recombination += recomb_num
            pair_num += 1
            if recomb_num >= one_stdev_higher:
                highway_recomb += recomb_num
                highway_num += 1

    with open("{0}/{1}".format(args.output, stat_file), "w") as output:
        output.write("Total Recombination Events: {0}\n".format(recombination))
        output.write("Number of recombination pairs: {0}\n".format(pair_num))
        output.write("Recombination Events in Highways: {0}\n".format(highway_recomb))
        output.write("Number of Highways: {0}\n".format(highway_num))


def itol_out(recomb_dict, label, one_stdev_higher):
    # Write out highways/recombination in iToL format
    print("Total number of recombination pairs in {0}: {1}".format(label, len(recomb_dict)), flush=True)
    print("Writing iToL file...", flush=True)
    highways = {}
    non_highways = {}
    # Sort out highways and prepare for iToL writeout
    for i in recomb_dict:
        j = i.split(":")
        if recomb_dict[i] >= one_stdev_higher:
            highways[i] = (j[0], j[1], recomb_dict[i], args.highway)
        else:
            non_highways[i] = (j[0], j[1], recomb_dict[i], args.nohighway)

    if len(highways) > 0:
        collection = [highways, non_highways]
    else:
        collection = [non_highways]
        print("No Highways in {0}".format(label), flush=True)
        # Calculate percentiles for highway coverage
        # Write out files
    count = 0
    for file in collection:
        if len(collection) > 1 and count == 0:
            count += 1
            if args.pairs:
                if args.filter:
                    name = "{0}_Indiscriminant_Highways".format(label)
                else:
                    name = "Indiscriminant_Highways"
            else:
                if args.filter:
                    name = "{0}_Discriminant_Highways".format(label)
                else:
                    name = "Discriminant_Highways"
            color = args.highway
            legend_shapes = "1"
            legend_labels = "Highway"
        else:
            if args.pairs:
                if args.filter:
                    name = "{0}_Indiscriminant_Recombination".format(label)
                else:
                    name = "Indsicriminant_Recombination"
            else:
                if args.filter:
                    name = "{0}_Discriminant_Recombination".format(label)
                else:
                    name = "Discriminant_Recombination"
            color = args.nohighway
            legend_shapes = "1"
            legend_labels = "Recombination"

        with open("{0}/{1}_hero_itol.txt".format(args.output, name), "w") as itol:
            # Set mandatory information
            itol.write("DATASET_CONNECTION\n\nSEPARATOR COMMA\n\nDATASET_LABEL,{0}\n\nCOLOR,#ff0ff0\n\n".format(name))

            # Set Legend information
            itol.write("LEGEND_TITLE,{0}\n".format(name))
            itol.write("LEGEND_SHAPES,{0}\n".format(legend_shapes))
            itol.write("LEGEND_COLORS,{0}\n".format(color))
            itol.write("LEGEND_LABELS,{0}\n\n".format(legend_labels))

            # Set arrow information
            if args.pairs:
                itol.write("DRAW_ARROWS,0\n")
            else:
                itol.write("DRAW_ARROWS,1\n")
            itol.write("ARROW_SIZE,60\nMAXIMUM_LINE_WIDTH,5\nCURVE_ANGLE,0\nCENTER_CURVES,1\nALIGN_TO_LABELS,0\n\n")

            # Set Data for branches
            itol.write("DATA\n")
            for pair in file:
                itol.write("{0},{1},15,{2}\n".format(file[pair][0], file[pair][1], file[pair][3]))


def gene_out(gene_file, recomb_genes):
    # Write out file showing recombination frequency of genes
    print("Writing out gene recombination file...", flush=True)
    genes = []
    for gene in recomb_genes:
        genes.append((gene, recomb_genes[gene]))

    sorted(genes, key=itemgetter(1)) # Sort on count of events

    with open("{0}/{1}".format(args.output, gene_file), "w") as output:
        output.write("gene\tfrequency\n")
        for gene in genes:
            output.write("{0}\t{1}\n".format(gene[0],gene[1]))


def size_out(size_file, recomb_sizes):
    # Write out recombination sizes
    print("Writing out recombination size file...", flush=True)
    with open("{0}/{1}".format(args.output, size_file), "w") as output:
        for fragment in recomb_sizes:
            output.write("{0}\n".format(fragment))


def frequent_finder(recomb_dict, d_file, r_file):
    # Sort Isolates by most frequent donor and recipient
    donors = {}
    recipients = {}

    for pair in recomb_dict:
        split_pair = pair.split(":")
        d,r = split_pair[0], split_pair[1]

        donors.setdefault(d, 0)
        donors[d] += recomb_dict[pair]

        recipients.setdefault(r, 0)
        recipients[r] += recomb_dict[pair]

    sorted_donors = sorted([(d,n) for d,n in donors.items()], key=itemgetter(1))
    sorted_recipients = sorted([(r,n) for r,n in recipients.items()], key=itemgetter(1))

    # Define frequent donors/recipients
    d_count = [n[1] for n in sorted_donors]
    r_count = [n[1] for n in sorted_recipients]

    for collection in [d_count, r_count]:
        if args.population > len(collection):
            for x in range(args.population - len(collection)):
                collection.append(0)

    d_one_stdev_higher = (sum(d_count) / len(d_count)) + statistics.pstdev(d_count)
    r_one_stdev_higher = (sum(r_count) / len(r_count)) + statistics.pstdev(r_count)

    # Write out donor/recipient frequencies
    with open("{0}/{1}".format(args.output, d_file), "w") as output:
        output.write("Donor\tCount\t*** means significant\n")
        for donor in reversed(sorted_donors):
            if donor[1] >= d_one_stdev_higher:
                output.write("{0}\t{1}***\n".format(donor[0], donor[1]))
            else:
                output.write("{0}\t{1}\n".format(donor[0], donor[1]))

    with open("{0}/{1}".format(args.output, r_file), "w") as output:
        output.write("Recipient\tCount\t*** means significant\n")
        for recipient in reversed(sorted_recipients):
            if recipient[1] >= r_one_stdev_higher:
                output.write("{0}\t{1}***\n".format(recipient[0], recipient[1]))
            else:
                output.write("{0}\t{1}\n".format(recipient[0], recipient[1]))


    # Load donors and recipients into dictionaries for bar graphs    
    df_donors = [n[0] for n in sorted_donors if n[1] != 0]
    df_recipients = [n[0] for n in sorted_recipients if n[1] != 0]

    d_data = {"Donors": df_donors, "Events": [n for n in d_count if n != 0]}
    data_frame = pandas.DataFrame(d_data)
    r_data = {"Recipients": df_recipients, "Events": [n for n in r_count if n != 0]}
    rata_frame = pandas.DataFrame(r_data)

    # Build Bar Graphs
    d_bar = (ggplot(data_frame, aes(x="Donors")) + 
    geom_bar(aes(weight = "Events")) + 
    labs(title = "Donation Events", x = "Donor", y = "# of Donations") + 
    theme(axis_text_x = element_text(angle = 45, hjust = 1, size = 3)))
    d_bar.save(filename="{0}/donations.{1}".format(args.output, args.format), format=args.format)

    r_bar = (ggplot(rata_frame, aes(x="Recipients")) + 
    geom_bar(aes(weight = "Events")) + 
    labs(title = "Recipient Events", x = "Recipient", y = "# of Receipts") + 
    theme(axis_text_x = element_text(angle = 45, hjust = 1, size = 3)))
    r_bar.save(filename="{0}/recipients.{1}".format(args.output, args.format), format=args.format)   

    return (sorted_donors, sorted_recipients)


def base_pair_plot(recomb_proportions):
    # Build bar graph of recombined base pairs per recipient
    data = {"x": list(recomb_proportions.keys()), "y": list(recomb_proportions.values())}
    data_frame = pandas.DataFrame(data)
    graph = ggplot(data_frame, aes(x="x"))
    graph += geom_bar(aes(weight = "y"))
    graph += labs(title = "Total Recombination", x = "Genome", y = "# of Recombined Base Pairs")
    graph += theme(axis_text_x = element_text(angle = 45, hjust = 1, size = 3))
    graph.save(filename="{0}/recombined_portions.{1}".format(args.output, args.format), format=args.format)


def box_and_whisker(filtered_stats):
    # Make box and whisker plots for donation and reciept separated by filter groups
    bw_data_donor = {"Donation": [], "Group": []}
    bw_data_recipient = {"Receipt": [], "Group": []}
    for group in filtered_stats:
        if args.inter and group == "inter/inter":
            continue

        for idx,collection in enumerate(filtered_stats[group]):
            target = bw_data_donor if idx == 0 else bw_data_recipient
            keyword = "Donation" if idx == 0 else "Receipt"
            for strain in collection:
                target[keyword].append(strain[1])
                target["Group"].append(group.split("/")[0])

    donor_data_frame = pandas.DataFrame(bw_data_donor)
    recipient_data_frame = pandas.DataFrame(bw_data_recipient)

    d_bw = (ggplot(donor_data_frame, aes(x="Group", y="Donation", color = "Group")) +
    geom_boxplot())
    d_bw.save(filename="{0}/donor_box_and_whisker.{1}".format(args.output, args.format), format=args.format)    

    r_bw = (ggplot(recipient_data_frame, aes(x="Group", y="Receipt", color = "Group")) +
    geom_boxplot())
    r_bw.save(filename="{0}/recipient_box_and_whisker.{1}".format(args.output, args.format), format=args.format)    
    

####### Main script #######
print("Finding Donors for all recombination events...", flush=True)

class PsuedoDirEntry:
    # Allows multithreading of os.scandir
    def __init__(self, name, path, is_dir):
        self.name = name
        self.path = os.path.abspath(path)
        self.is_dir = is_dir

# Make a directory for the output files
args.output = make_directory(args.output)
os.mkdir("{0}/filtered_recombination_events_{1}bp_{2}logBF".format(args.output, args.length, args.log))
os.mkdir("{0}/blast_files".format(args.output))

# Iterate over every gene and add to global recombination dictionary *Multi-threaded section*
gene_dicts = [] # Each item in list is dictionary. Key - Recombination Pair. Value - Tuple of data (# of recomb events, length of event, gene recombined in)
for i in os.scandir(args.fastgear):
    gene_dicts.append(PsuedoDirEntry(i.name, i.path, i.is_dir()))

# Multithreading run through each gene to find recombination pairs
pool = multiprocessing.Pool(processes=args.cpus)
gene_recombs = pool.map(lineage_fastas, gene_dicts)
pool.close()
pool.join()

filtered_stats = {} # Donation/receipt data for each filtering group
external = [] # External recombination master list
bad_genes = [] # All failed genes
good_genes = []  # All good genes

# Sort recombination events from failed genes.
for gene in gene_recombs:
    if isinstance(gene, list):
        good_genes += gene[0]
        external += gene[1]
    else:
        bad_genes.append(gene)

unique_recombination, paired_recombination, recombination_fragments, ancestral_events  = find_pairs(good_genes)

write_out_events(unique_recombination)

# Filter if chosen, use paired results if chosen, and write output
if args.filter:
    # Get strain groups
    strain_filters = get_filter_groups()

    # Paired keys will always be correct regardless of user choice based on find_pairs logic
    filtered_recombination = filter_recombination(paired_recombination, strain_filters)

    for group in filtered_recombination:
        os.mkdir("{0}/{1}".format(args.output, group))
        write_out("{0}/{0}".format(group), filtered_recombination[group])

elif args.pairs:
    write_out('inter', paired_recombination)

else:
    write_out('inter', unique_recombination)

# Write out External recombination fragments
Bio.SeqIO.write(external, "{0}/external_recomb_fragments.fa".format(args.output), "fasta")

# Recomb_proportions plot
base_pair_plot(recombination_fragments)

# Make box-whisker plots of donors and recipients
if args.filter and not args.pairs:
    box_and_whisker(filtered_stats)

# Write out final files
with open("{0}/failed_genes.txt".format(args.output), "w") as bad: # Write out any failed genes for further evaluation
    for gene in bad_genes:
        bad.write("{0}\n".format(gene))

with open("{0}/recombined_base_pairs.txt".format(args.output), "w") as output: # Write out base pair total for recombination in each genome
    output.write("Genome\tBase Pairs\n")
    for genome in recombination_fragments:
        output.write("{0}\t{1}\n".format(genome, recombination_fragments[genome]))

if args.cleanup:
    shutil.rmtree("{0}/blast_files".format(args.output))

# Write out runtime
print("{0} Ancestral Events filtered out".format(ancestral_events), flush = True)
print("Finished in {0}".format(time.time() - t0), flush=True)


