#! /usr/bin/env python3

#### HERO: Highway Enumeration from Recombination Observations

import argparse
import Bio.SeqIO
import itertools
import math
import multiprocessing
import os
import random
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
parser.add_argument("--cleanup", action = "store_false", help = "Keep intermediate files.")
parser.add_argument("--cpus", metavar = '', default = 1, type = int, help = "Number of threads to use. [1]")
parser.add_argument("--nohighway", metavar = '', default = "#D3D3D3", help = "Hexadecimal code for non-highway color in output file [#D3D3D3]")
parser.add_argument("--output", metavar = '', default = "HERO", help = "Output directory for HERO files")
#parser.add_argument("--highway", default = "#0000FF", help = "Hexadecimal code for highway color in output file")
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

    with open("{0}/output/recombinations_recent.txt".format(fastgear.path), "r") as file:
        if len(file.readlines()) < 3:
            print("No Recombination in gene: {0}".format(fastgear.name), flush = True)
            return None

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
    for num in lineages:
        records = []
        for strain in lineages[num]:
            records.append(fg_strains[strain].unaligned_fasta)

        Bio.SeqIO.write(records, "{0}/lineage_{1}.fa".format(fastgear.path, num), "fasta")

        blastdb = subprocess.run("makeblastdb -in {0}/lineage_{1}.fa -dbtype nucl".format(fastgear.path, num), shell = True, stdout = subprocess.PIPE, stderr = subprocess.PIPE)
        if blastdb.stderr:
            print(blastdb.stderr.decode('ascii'), flush=True)
            print("Revisit the gene alignment of {0} to troubleshoot. Skipping for results.".format(fastgear.name), flush=True)
            return fastgear.name

    return recomb_parser(fastgear, Strains, fg_strains, query_lineages, max(lineages.keys()))


def recomb_parser(fastgear, Strains, fg_strains, query_lineages, max_lineage):
    # Parse recombination events and build query blast files
    external = []
    recomb_frags = {} # Recombination Event sizes
    with open("{0}/output/recombinations_recent.txt".format(fastgear.path), "r") as recomb_file:
        next(recomb_file)
        next(recomb_file)
        for line in recomb_file:
            line = line.split()
            
            if float(line[4]) < args.log:
                continue


            # Pull sequence fragment and set lineage to blast against
            recipient, start_idx, end_idx  = fg_strains[line[5]], int(line[0]) - 1, int(line[1]) - 1
            sequence = unalignment(recipient.fasta.seq[start_idx:end_idx])
            if len(sequence) <= args.length:
                continue
            identity = "{0}:{1}:{2}".format(line[5],start_idx,end_idx)
            recomb_frags[identity] = str(len(sequence)) # Log recombination fragment size

            # Prep fasta format for sequence fragment, then blast against lineage database
            if int(line[2]) > int(max_lineage):
                rec = Bio.SeqRecord.SeqRecord(id = "{0}:{1}".format(fastgear.name, identity), seq = sequence)
                external.append(rec)
            else:
                rec = Bio.SeqRecord.SeqRecord(id = identity, seq = sequence)
                query_lineages[line[2]].append(rec)


    rec_indexes = {}
    recombination_events = {}
    for query in query_lineages:
        filename = "{0}/lineage_{1}_query.fa".format(fastgear.path, query)
        Bio.SeqIO.write(query_lineages[query], filename, "fasta")
        blast_output = subprocess.run("blastn -query {0} -subject {1}/lineage_{2}.fa -evalue 1e-5 -outfmt 7 -max_hsps 1".format(filename, fastgear.path, query),\
            shell = True, stdout = subprocess.PIPE, stderr = subprocess.PIPE)

        # Parse blast_output
        get_event = 0
        for entry in blast_output.stdout.decode('ascii').split("\n"): # Skip Header lines
            if "BLASTN" in entry:
                get_event = 1
 
            if entry == "" or get_event == 0: # Checks for empty string (due to parsing) and if match has been found already
                continue

            if entry[0] == "#": # Skips if header line
                continue

            entry = entry.split()
            if entry[0] == entry[1]: # Read can align to itself in database. Skip this too. ***This might not be true. Need to investigate.***
                continue
            
            # Parse BLAST Results
            key = "{1}:{0}".format(entry[0].split(":")[0], entry[1])
            a = int(entry[8])
            b = int(entry[9])

            # Determine if overlap exists between any recombination events
            overlap = False
            if rec_indexes.setdefault(key, [(a,b)]) != [(a,b)]:
                for idx,coords in enumerate(rec_indexes[key]):
                    overlap = seq_overlap(coords[0], coords[1], a, b)
                    # Add event or extend previous range to list for that pair and edit fragment length of event
                    if overlap:
                        rec_indexes[key][idx] = overlap
                        old_range = recombination_events[key][idx][1] # Previous indices for recombination fragment
                        current_range = tuple(int(x) for x in entry[0].split(":")[2:]) # New indicies to extend the event
                        new_range = seq_overlap(old_range[0], old_range[1], current_range[0], current_range[1]) # Use overlap function to build a consensus range between the old and current
                        recombination_events[key][idx] = (1, new_range, fastgear.name)
                        break
                
            # Note recombination event if no overlap exists
            if not overlap:
                rec_indexes[key].append((a,b))
                recombination_events.setdefault(key, [])
                recomb_range = tuple(int(x) for x in entry[0].split(":")[1:])
                recombination_events[key].append((1, recomb_range, fastgear.name))
            get_event = 0

    with open("{0}/filtered_recombination_events_{1}bp_{2}logBF/{3}.txt".format(args.output, args.length, args.log, fastgear.name), "w") as file:
        file.write("Donor\tRecipient\tFragment Length\n")
        for key in recombination_events:
            for event in recombination_events[key]:
                donor, recipient = key.split(":")[0], key.split(":")[1]
                length = event[1][1] - event[1][0]
                file.write("{0}\t{1}\t{2}\n".format(donor, recipient, length))

    if args.cleanup == True:
        cleanup(fastgear)

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
    if a < start and b > start:
        if b > end:
#            new_length = (start - a) + (end - start) + (b - end)
            return (a,b)
        else:
#            new_length = (start - a) + (end - start)
            return (a,end)
    elif a > start and a < end:
        if b > end:
#            new_length = (end - start) + (b - end)
            return (start,b)
        else:
#            new_length = (end - start)
            return (start,end)
    else:
        overlap = 0
    
    return overlap
    

def cleanup(fastgear):
    # Clean intermediate files
    for file in os.scandir(fastgear.path):
        if "lineage" in file.name:
            os.remove(file.path)


def filtering(final_rec_events):
    # Sort recombination events by user list
    filters = {} # Sort strains by user categories
    dict_dict = {"other":{}} # Keeps labels and dictionaries together
    if args.filter:
        try:
            with open(args.filter, "r") as filter_list:
                # Sort strains by category into 'filters'
                for line in filter_list:
                    line = line.strip().split(",")
                    strain,cat = line[0],line[1]
                
                    dict_dict[cat] = {}
                    filters[strain] = cat
        except IOError:
            print("Filtering file not found! Skipping filteration",flush=True)

    for event in final_rec_events:
        i = event.split(":")
        # Test presence of strain in filter list
        for strain in i:
            try:
                filters[strain]
            except KeyError:
                filters[strain] = "other"

        # If categories are the same, add to category dictionary, otherwise add to 'other'
        t1, t2 = filters[i[0]], filters[i[1]]
        if t1 == t2:
            if filters[i[0]] == "other":
                dict_dict['other'][event] = final_rec_events[event]
            else:
                dict_dict[t1][event] = final_rec_events[event]
        else:
            cat = "{0}2{1}".format(t1, t2)
            dict_dict.setdefault(cat, {})
            dict_dict[cat][event] = final_rec_events[event]

    return dict_dict


def write_out(label, cat):
    # Calculate Highway definition
    if len(cat) == 0:
        print("{0} category has no recombination events. No output files to write".format(label), flush=True)
        return None
    # Sort recombination data
    recomb_nums = {}
    recomb_sizes = []
    recomb_genes = {}
    for pair in cat:
        for data in cat[pair]:
            recomb_nums.setdefault(pair, 0)
            recomb_nums[pair] += data[0]

            recomb_sizes.append(data[1])

            recomb_genes.setdefault(data[2], 0)
            recomb_genes[data[2]] += 1


    one_stdev_higher = (sum(recomb_nums.values()) / len(recomb_nums.values())) + statistics.pstdev(recomb_nums.values())

    # Write output files
    names = name_collection(label)
    ## ^pair_file, stat_file, d_file, r_file, gene_file, size_file

    stats_writer(recomb_nums, names[0], names[1], one_stdev_higher)
    itol_out(recomb_nums, label, one_stdev_higher)
    gene_out(names[4], recomb_genes)
    size_out(names[5], recomb_sizes)
    if not args.pairs:
        frequent_finder(names[0], names[2], names[3])


def name_collection(label):
    # Assign file names for each output function call
    if label == "other":
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
    non_highway_color = args.nohighway
    highway_color = ["#FFFF00","#FFA500","#D17B0F","#DF0B0B"]
    print("Writing iToL file...", flush=True)
    highways = {}
    non_highways = {}
    # Calculate stats to sort recombination pairs
    max_recomb = max(recomb_dict.values())
    quartile = (max_recomb - one_stdev_higher)/4
    twenty_fifth = "{0} - {1}".format(one_stdev_higher, one_stdev_higher + quartile)
    fiftieth = ">{0} - {1}".format(one_stdev_higher + quartile, one_stdev_higher + quartile*2)
    seventy_fifth = ">{0} - {1}".format(one_stdev_higher + quartile*2, one_stdev_higher + quartile*3)    
    hundredth = ">{0}".format(one_stdev_higher + quartile*3)
    # Sort out highways and prepare for iToL writeout
    for i in recomb_dict:
        j = i.split(":")
        if recomb_dict[i] >= one_stdev_higher:
            if recomb_dict[i] <= (one_stdev_higher + quartile):
                cur_color = highway_color[0]
            elif recomb_dict[i] <= (one_stdev_higher + quartile*2):
                cur_color = highway_color[1]
            elif recomb_dict[i] <= (one_stdev_higher + quartile*3):
                cur_color = highway_color[2]
            else:
                cur_color = highway_color[3]
            highways[i] = (j[0], j[1], recomb_dict[i], cur_color)
        else:
            non_highways[i] = (j[0], j[1], recomb_dict[i], non_highway_color)

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
            color = ",".join(highway_color)
            legend_shapes = "1,1,1,1"
            legend_labels = "{0},{1},{2},{3}".format(twenty_fifth, fiftieth, seventy_fifth, hundredth)
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
            color = non_highway_color
            legend_shapes = "1"
            legend_labels = name

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
    sorted_genes = [("placeholder",0)]
    for gene in recomb_genes:
        found = 0
        for pos,gene2 in enumerate(sorted_genes):
            if recomb_genes[gene] < gene2[1]:
                sorted_genes.insert(pos, (gene, recomb_genes[gene]))
                found = 1
                break
        if found == 0:
            sorted_genes.append((gene, recomb_genes[gene]))
    sorted_genes.remove(("placeholder",0))


    with open("{0}/{1}".format(args.output, gene_file), "w") as output:
        output.write("gene\tfrequency\n")
        for gene in sorted_genes:
            output.write("{0}\t{1}\n".format(gene[0],gene[1]))


def size_out(size_file, recomb_sizes):
    # Write out recombination sizes
    print("Writing out recombination size file...", flush=True)
    with open("{0}/{1}".format(args.output, size_file), "w") as output:
        for fragment in recomb_sizes:
            output.write("{0}\n".format(fragment[1] - fragment[0]))


def frequent_finder(pair_file, d_file, r_file):
    # Sort Isolates by most frequent donor and recipient
    donors = [("placeholder", 0)]
    recipients = [("placeholder", 0)]
    recomb_parse_dict = {}

    with open("{0}/{1}".format(args.output, pair_file), "r") as rec_file:
        next(rec_file)
        # Sort recombination_pairs into donor/recipient frequencies        
        for line in rec_file:
            line = line.split()

            pair,count = line[0].split(":"),int(line[1])
            donor,recipient = pair[0], pair[1]

            recomb_parse_dict.setdefault(donor, [0,0])
            recomb_parse_dict.setdefault(recipient, [0,0])

            recomb_parse_dict[donor][0] += count
            recomb_parse_dict[recipient][1] += count


    # Sort values into donor/recipient lists based on frequency
    d_count = []
    r_count = []
    for isolate in recomb_parse_dict:
        donation = recomb_parse_dict[isolate][0]
        receipt = recomb_parse_dict[isolate][1]

        d_trig = 0
        for pos,donor in enumerate(donors):
            if donation < donor[1]:
                donors.insert(pos, (isolate, donation))
                d_trig = 1
                break
        if d_trig == 0:
            donors.append((isolate, donation))
        d_count.append(donation)


        r_trig = 0
        for pos,recipient in enumerate(recipients):
            if receipt < recipient[1]:
                recipients.insert(pos, (isolate, receipt))
                r_trig = 1
                break
        if r_trig == 0:
            recipients.insert(pos, (isolate, receipt))
        r_count.append(receipt)

    # Remove placeholders
    donors.remove(("placeholder", 0))
    recipients.remove(("placeholder", 0))

    # Define frequent donors/recipients
    d_one_stdev_higher = (sum(d_count) / len(d_count)) + statistics.pstdev(d_count)
    r_one_stdev_higher = (sum(r_count) / len(r_count)) + statistics.pstdev(r_count)

    # Write out donor/recipient frequencies
    with open("{0}/{1}".format(args.output, d_file), "w") as output:
        output.write("Donor\tCount\t*** means significant\n")
        for donor in reversed(donors):
            if donor[1] >= d_one_stdev_higher:
                output.write("{0}\t{1}***\n".format(donor[0], donor[1]))
            else:
                output.write("{0}\t{1}\n".format(donor[0], donor[1]))

    with open("{0}/{1}".format(args.output, r_file), "w") as output:
        output.write("Recipient\tCount\t*** means significant\n")
        for recipient in reversed(recipients):
            if recipient[1] >= r_one_stdev_higher:
                output.write("{0}\t{1}***\n".format(recipient[0], recipient[1]))
            else:
                output.write("{0}\t{1}\n".format(recipient[0], recipient[1]))



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

# Iterate over every gene and add to global recombination dictionary *Multi-threaded section*
gene_dicts = [] # Each item in list is dictionary. Key - Recombination Pair. Value - Tuple of data (# of recomb events, length of event, gene recombined in)
external = []
recomb_proportions = {}
final_rec_events = {}
paired_final_rec_events = {}
bad_genes = []

for i in os.scandir(args.fastgear):
    gene_dicts.append(PsuedoDirEntry(i.name, i.path, i.is_dir()))

pool = multiprocessing.Pool(processes=args.cpus)
gene_dicts = pool.map(lineage_fastas, gene_dicts)
pool.close()
pool.join()

# Sort recombination events from failed genes.
for gene_dict in gene_dicts:
    if isinstance(gene_dict, list):
        external += gene_dict[1]
        gene_dict = gene_dict[0]
        for key in gene_dict:
            valid = 1
            if args.pairs:
                for iteration in itertools.permutations(key.split(":"), len(key.split(":"))):
                    if ":".join(iteration) in paired_final_rec_events:
                        paired_key = ":".join(iteration)
                        valid = 0                            
                        break

            if valid == 1:
                paired_key = key

            final_rec_events.setdefault(key, [])
            paired_final_rec_events.setdefault(paired_key, [])
            recomb_proportions.setdefault(key.split(":")[1], 0) # Add recombination fragment size to recipient in dictionary

            for event in gene_dict[key]:
                final_rec_events[key].append(event)
                paired_final_rec_events[paired_key].append(event)
                recomb_proportions[key.split(":")[1]] += int(event[1][1] - event[1][0])
    else:
        bad_genes.append(gene_dict)


# Filter if chosen, use paired results if chosen, and write output
if args.filter:
    if args.pairs:
        filtered_dicts = filtering(paired_final_rec_events)
    else:
        filtered_dicts = filtering(final_rec_events)
    for cat in filtered_dicts:
        os.mkdir("{0}/{1}".format(args.output, cat))
        write_out("{0}/{0}".format(cat), filtered_dicts[cat])
elif args.pairs:
    write_out('other', paired_final_rec_events)
else:
    write_out('other', final_rec_events)        

#if args.pairs:
#    paired_final_rec_events = filtering(paired_final_rec_events)
#else:
#    final_rec_events = filtering(final_rec_events)

#for cat in final_rec_events:
#    if args.filter:
#        os.mkdir("{0}/{1}".format(args.output, cat))
#        write_out("{0}/{0}".format(cat), final_rec_events[cat])
#    else:
#        write_out(cat, final_rec_events[cat])

with open("{0}/failed_genes.txt".format(args.output), "w") as bad: # Write out any failed genes for further evaluation
    for gene in bad_genes:
        bad.write("{0}\n".format(gene))

with open("{0}/recombined_base_pairs.txt".format(args.output), "w") as output: # Write out base pair total for recombination in each genome
    output.write("Genome\tBase Pairs\n")
    for genome in recomb_proportions:
        output.write("{0}\t{1}\n".format(genome, recomb_proportions[genome]))

Bio.SeqIO.write(external, "{0}/external_recomb_fragments.fa".format(args.output), "fasta")

# Write out runtime
print("Finished in {0}".format(time.time() - t0), flush=True)


