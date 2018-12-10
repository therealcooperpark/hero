#! /usr/bin/env python3

#### HERO: Highway Enumeration from Recombination Observations

import argparse
import Bio.SeqIO
import itertools
import math
from multiprocessing import Pool
import os
import random
import statistics
import subprocess
import time

t0 = time.time()

# Args
parser = argparse.ArgumentParser(description = "HERO - Highways Elucidated by Recombination Observations", usage = "hero.py [options] fastgear")
parser.add_argument("fastgear", help = "Directory containing all fastGEAR runs for the population")
parser.add_argument("--log", metavar = '', default = 0, type = int, help = "Minimum log(BF) length to accept recombination event. [0]")
parser.add_argument("--length", metavar = '', default = 2, type = int, help = "Minimum fragment length for recombined sequence. [2]")
parser.add_argument("--pairs", action = "store_true", help = "Ignore direction of event when calculating recombination highways.")
parser.add_argument("--cleanup", action = "store_false", help = "Keep intermediate files.")
parser.add_argument("--cpus", metavar = '', default = 1, type = int, help = "Number of threads to use. [1]")
parser.add_argument("--nohighway", metavar = '', default = "#D3D3D3", help = "Hexadecimal code for non-highway color in output file [#D3D3D3]")
parser.add_argument("--output", metavar = '', default = "HERO", help = "Output directory for HERO files")
#parser.add_argument("--highway", default = "#0000FF", help = "Hexadecimal code for highway color in output file")
args = parser.parse_args()



def lineage_fastas(fastgear):
    # Check absolute vs relative path for appropriate pathways to gene fastGEARs
    if not fastgear.is_dir:
        print("Non-fastGEAR directory found. Skipping")
        return fastgear.name

    # Setup main collections
    class Strains:
        alignment_length = 0 # Total alignment length will be set with each new gene processed

        def __init__(self, strain_num, lineage, name):
            self.name = name
            self.strain_num = strain_num
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
        with open("{0}/output/lineage_information.txt".format(os.path.abspath(fastgear.path)), "r") as file:
            next(file)
            for line in file:
                line = line.split()
                fg_strains[line[3]] = Strains(line[0], line[1], line[3])
    except FileNotFoundError:
        print("Gene {0} has an unsuccessful fastGEAR run. Skipping for results.".format(fastgear.name))
        return fastgear.name

    # Identify alignment file
    alignment = ""
    for file in os.scandir(os.path.abspath(fastgear.path)):
        if file.is_file() and ".mat" not in file.name:
            try:
                if any(Bio.SeqIO.parse(file.path, "fasta")):
                    alignment = file.path
            except UnicodeDecodeError:
                pass
    if alignment == "":
        print("Gene {0} has no alignment file present. Skipping for results.".format(fastgear.name))
        return fastgear.name

    # Generate lineage based fasta files
    for record in Bio.SeqIO.parse(alignment, "fasta"):
        fg_strains[record.id].fasta = record
        fg_strains[record.id].unaligned_fasta = record

    # Unalign sequence for blast
    for strain in fg_strains.values():	
        sequence = unalignment(strain.unaligned_fasta.seq)
        strain.unaligned_fasta.seq = sequence

    # Write out reference lineage files
    for num in lineages:
        records = []
        for strain in lineages[num]:
            records.append(fg_strains[strain].unaligned_fasta)

        Bio.SeqIO.write(records, "{0}/lineage_{1}.fa".format(os.path.abspath(fastgear.path), num), "fasta")

        blastdb = subprocess.run("makeblastdb -in {0}/lineage_{1}.fa -dbtype nucl".format(os.path.abspath(fastgear.path), num), shell = True, stdout = subprocess.PIPE, stderr = subprocess.PIPE)
        if blastdb.stderr:
            print(blastdb.stderr.decode('ascii'))
            print("Revisit the gene alignment of {0} to troubleshoot. Skipping for results.".format(fastgear.name))
            return fastgear.name

    return recomb_parser(fastgear, Strains, fg_strains, query_lineages, max(lineages.keys()))


def recomb_parser(fastgear, Strains, fg_strains, query_lineages, max_lineage):
    # Parse recombination events and build query blast files
    recombination_events = {}
    rec_indexes = {}
    with open("{0}/output/recombinations_recent.txt".format(os.path.abspath(fastgear.path)), "r") as recomb_file:
        next(recomb_file)
        next(recomb_file)
        for line in recomb_file:
            line = line.split()
            
            if float(line[4]) < args.log:
                continue

            if int(line[2]) > int(max_lineage):
                query_lineages.setdefault(line[2], [])

            # Pull sequence fragment and set lineage to blast against
            recipient = fg_strains[line[5]]
            start_idx = int(line[0]) - 1
            if start_idx < 1: # Recombined seqs can still have hyphens at the start. Trim value will force start index negative
                start_idx = 1 # This is the base case to avoid that, while letting the trim properly manage the end index.
            end_idx = int(line[1]) - 1
            if end_idx - start_idx <= args.length:
                print("Recombination event is too small. Skipping.")
                continue
            sequence = recipient.fasta.seq[start_idx:end_idx]
            sequence = unalignment(sequence)

            # Prep fasta format for sequence fragment, then blast against lineage database
            rec = Bio.SeqRecord.SeqRecord(id = line[5], seq = sequence)
            query_lineages[line[2]].append(rec)


    for query in query_lineages:
        # Print out external recombination events
        if query > max_lineage:
#            Bio.SeqIO.write(query_lineages[query], "{0}/external_recombinations.fa".format(os.path.abspath(fastgear.path)), "fasta")
            continue

        filename = "{0}/lineage_{1}_query.fa".format(os.path.abspath(fastgear.path), query)
        Bio.SeqIO.write(query_lineages[query], filename, "fasta")
        blast_output = subprocess.run("blastn -query {0} -subject {1}/lineage_{2}.fa -evalue 1e-5 -outfmt 7 -max_hsps 1".format(filename, os.path.abspath(fastgear.path), query),\
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
            key = "{0}:{1}".format(entry[0], entry[1])
            a = entry[8]
            b = entry[9]

            # Determine if overlap exists between any recombination events
            overlap = False
            if key not in rec_indexes:
                rec_indexes.setdefault(key, [(a,b)])
            else:
                for coords in rec_indexes[key]:
                    overlap = seq_overlap(coords[0], coords[1], a, b)
                    if overlap:
                        break
                # Add event range to list for that pair
                rec_indexes[key].append((a,b))

            # Note recombination event if no overlap exists
            recombination_events.setdefault(key, 0)
            if not overlap:
                recombination_events[key] += 1
            get_event = 0

    if args.cleanup == True:
        cleanup(fastgear)

    return recombination_events


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
        overlap = 1
    elif a > start and a < end:
        overlap = 1
    else:
        overlap = 0
    
    return overlap
    

def cleanup(fastgear):
    # Clean intermediate files
    for file in os.scandir(os.path.abspath(fastgear.path)):
        if "lineage" in file.name:
            os.remove(file.path)

def stats_writer(recomb_dict, one_stdev_higher):
    # Write out recombination pairs and quantities + stats
    recombination = 0
    pair_num = 0
    highway_recomb = 0
    highway_num = 0

    filename = "recombination_pairs.txt"
    print("Writing out results to {0}...".format(filename))
    with open("{0}/{1}".format(args.output, filename), "w") as hero:
        hero.write("Donor:Recipient\tEvents\n")

        for event in recomb_dict:
            recomb_num = recomb_dict[event]
            hero.write("{0}\t{1}\n".format(event, recomb_num))

            recombination += recomb_num
            pair_num += 1
            if recomb_num >= one_stdev_higher:
                highway_recomb += recomb_num
                highway_num += 1

    with open("{0}/HERO_statistics.txt".format(args.output), "w") as output:
        output.write("Total Recombination Events: {0}\n".format(recombination))
        output.write("Number of recombination pairs: {0}\n".format(pair_num))
        output.write("Recombination Events in Highways: {0}\n".format(highway_recomb))
        output.write("Number of Highways: {0}\n".format(highway_num))

    # Load other output files
    frequent_finder()
            

def frequent_finder():
    # Sort Isolates by most frequent donor and recipient
    donors = [("placeholder", 0)]
    recipients = [("placeholder", 0)]
    recomb_parse_dict = {}
    with open("{0}/recombination_pairs.txt".format(args.output), "r") as rec_file:
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
    for isolate in recomb_parse_dict:
        d_trig = 0
        for pos,donor in enumerate(donors):
            if recomb_parse_dict[isolate][0] < donor[1]:
                donors.insert(pos, (isolate, recomb_parse_dict[isolate][0]))
                d_trig = 1
                break
        if d_trig == 0:
            donors.append((isolate, recomb_parse_dict[isolate][0]))

        r_trig = 0
        for pos,recipient in enumerate(recipients):
            if recomb_parse_dict[isolate][1] < recipient[1]:
                recipients.insert(pos, (isolate, recomb_parse_dict[isolate][1]))
                r_trig = 1
                break
        if r_trig == 0:
            recipients.insert(pos, (isolate, recomb_parse_dict[isolate][1]))

    # Remove placeholders
    donors.remove(("placeholder", 0))
    recipients.remove(("placeholder", 0))

    # Write out donor/recipient frequencies
    with open("{0}/donor_frequency.txt".format(args.output), "w") as output:
        output.write("Donor\tCount\n")
        for donor in reversed(donors):
            output.write("{0}\t{1}\n".format(donor[0], donor[1]))

    with open("{0}/recipient_frequency.txt".format(args.output), "w") as output:
        output.write("Recipient\tCount\n")
        for recipient in reversed(recipients):
            output.write("{0}\t{1}\n".format(recipient[0], recipient[1]))


def itol_out(recomb_dict, one_stdev_higher):
    print("Total number of recombination pairs: {0}".format(len(recomb_dict)))
    non_highway_color = args.nohighway
    highway_color = ["#FFFF00","#FFA500","#D17B0F","#DF0B0B"]
    max_recomb = max(recomb_dict.values())
    print("Writing iToL file...")
    highways = {}
    non_highways = {}
    highway_recomb_values = []
    # Sort out highways and prepare for iToL writeout
    for i in recomb_dict:
        j = i.split(":")
        if recomb_dict[i] >= one_stdev_higher:
            if recomb_dict[i] <= math.ceil(max_recomb*0.25):
                cur_color = highway_color[0]
            elif recomb_dict[i] <= math.ceil(max_recomb*0.50):
                cur_color = highway_color[1]
            elif recomb_dict[i] <= math.ceil(max_recomb*0.75):
                cur_color = highway_color[2]
            else:
                cur_color = highway_color[3]
            highways[i] = (j[0], j[1], recomb_dict[i], cur_color)
            highway_recomb_values.append(recomb_dict[i])
        else:
            non_highways[i] = (j[0], j[1], recomb_dict[i], non_highway_color)

    # Calculate percentiles for highway coverage
    min_recomb = min(highway_recomb_values)
    quartile = math.ceil((max_recomb - min_recomb) / 4)
    twenty_fifth = "{0} - {1}".format(min_recomb, min_recomb + quartile)
    fiftieth = "{0} - {1}".format(min_recomb + quartile + 1, min_recomb + quartile*2)
    seventy_fifth = "{0} - {1}".format(min_recomb + quartile*2 + 1, min_recomb + quartile*3)    
    hundredth = ">{0}".format(min_recomb + quartile*3)
    # Write out files
    for num,file in enumerate([highways, non_highways]):
        if num == 0:
            if args.pairs:
                name = "Directionless_Highways"
            else:
                name = "Direction_based_Highways"
            color = ",".join(highway_color)
            legend_shapes = "1,1,1,1"
            legend_labels = "{0},{1},{2},{3}".format(twenty_fifth, fiftieth, seventy_fifth, hundredth)
        else:
            if args.pairs:
                name = "Directionless_recombination"
            else:
                name = "Direction_based_Recombination"
            color = non_highway_color
            legend_shapes = "1"
            legend_labels = name

        with open("{0}/{1}_hero_itol.txt".format(args.output, name), "w") as itol:
            # Set mandatory information
            itol.write("DATASET_CONNECTION\n\n")
            itol.write("SEPARATOR COMMA\n\n")
            itol.write("DATASET_LABEL,{0}\n\n".format(name))
            itol.write("COLOR,#ff0ff0\n\n")

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
            itol.write("ARROW_SIZE,60\n")
            itol.write("MAXIMUM_LINE_WIDTH,15\n")
            itol.write("CURVE_ANGLE,0\n")
            itol.write("CENTER_CURVES,1\n")
            itol.write("ALIGN_TO_LABELS,0\n\n")

            # Set Data for branches
            itol.write("DATA\n")
            for pair in file:
                itol.write("{0},{1},15,{2}\n".format(file[pair][0], file[pair][1], file[pair][3]))



### Main script ###

class PsuedoDirEntry:
    # Allows multithreading of os.scandir
    def __init__(self, name, path, is_dir):
        self.name = name
        self.path = path
        self.is_dir = is_dir

# Global collections
bad_genes = []
cwd = os.getcwd()
final_rec_events = {}

try:
    os.mkdir(args.output)
except IOError:
    pass

# Iterate over every gene and add to global recombination dictionary *Multi-threaded section*
gene_dicts = []
for i in os.scandir(args.fastgear):
    gene_dicts.append(PsuedoDirEntry(i.name, i.path, i.is_dir()))

if __name__ == '__main__':
    pool = Pool(processes=args.cpus)

    gene_dicts = pool.map(lineage_fastas, gene_dicts)

    pool.close()
    pool.join()


# Sort recombination events from failed genes.
for gene_dict in gene_dicts:
    if isinstance(gene_dict, dict):
        for pair in gene_dict:
            valid = 1
            # Ignores direction of recombination events when True
            if args.pairs:
                for iteration in itertools.permutations(pair.split(":"), len(pair.split(":"))):
                    if ":".join(iteration) in final_rec_events:
                        valid = 0
                        final_rec_events[":".join(iteration)] += gene_dict[pair]
                        break

            if valid == 1:
                final_rec_events.setdefault(pair, 0)
                final_rec_events[pair] += gene_dict[pair]
    else:
        bad_genes.append(gene_dict)




# Calculate Highway definition
one_stdev_higher = (sum(final_rec_events.values()) / len(final_rec_events.values())) + statistics.pstdev(final_rec_events.values())

# Write output files
stats_writer(final_rec_events, one_stdev_higher)
itol_out(final_rec_events, one_stdev_higher)
with open("{0}/HERO_failed_genes.txt".format(args.output), "w") as bad: # Write out any failed genes for further evaluation
    for gene in bad_genes:
        bad.write("{0}\n".format(gene))

# Write out time table
print("Finished in {0}".format(time.time() - t0))
