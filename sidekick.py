#! /usr/bin/env python3

import argparse
import Bio.SeqIO
import math
import os
import random
import shutil
import string
import subprocess

# Arguments
parser = argparse.ArgumentParser(description = "Use before HERO to extract individual alignment files for each core gene in pan-genome. Will also run fastGEAR individually on genes.", usage = "sidekick.py [options] gffs")
parser.add_argument("gffs", help = "Filepath to GFF files")
parser.add_argument("--core", action = 'store_false', help = "Compile all core genes found in --cluster_file [True]")
parser.add_argument("--accessory", action = 'store_true', help = "Compile all accessory genes found in --cluster_file present in at least 2 isolates, but fewer than 'core_definition' [False]")
parser.add_argument("--core_definition", metavar = '', default = "0.99", type = float, help = "Fraction of genomes that a gene must be found in to be 'core' [0.99]")
parser.add_argument("--alignments", metavar = '', default = "./pan_genome_sequences", help = "Filepath to individual gene alignments, [./pan_genome_sequences]")
parser.add_argument("--cluster_file", metavar = '', default = "./clustered_proteins", help = "Filepath to the 'clustered_proteins' output of Roary [./clustered_proteins]")
parser.add_argument("--fastGEAR", action = 'store_true', help = "Run fastGEAR on every gene collected by sidekick [False]")
parser.add_argument("--output", metavar = '', default = ".", help = "Filepath to output directory for renamed genome protein files. [./]")
args = parser.parse_args()


def pan_genome_extraction(cluster_file, genome_num):
    #### Find Core genes among clustered_proteins genes

    core_cluster = []
    accessory_cluster = []

    with open(args.cluster_file, "r") as cluster:
        for line in cluster:
            line = line.split(":")
        
            # Name of protein cluster
            group = line[0]
            if "/" in group:
                group = "_".join(group.split("/"))

            # Count number of occurances
            prots = line[1].split("\t")
            if args.core and len(prots) > math.floor(genome_num * args.core_definition):
                core_cluster.append(group)
            elif args.accessory and len(prots) > 1 and len(prots) <= math.floor(genome_num * 0.99):
                accessory_cluster.append(group)

    if args.accessory:
        print("{0} Groups found in Accessory Genome".format(len(accessory_cluster)), flush = True)
    print("{0} Groups found in core genome".format(len(core_cluster)), flush = True)


    return (core_cluster, accessory_cluster)


def file_finder(core_cluster, alignments, output):
    # Copy alignment files from pan_genome_sequences folder to new folder
    for group in core_cluster:
        try:
            shutil.copy("{0}/{1}.fa.aln".format(alignments, group), output)
        except FileNotFoundError:

            # Split by the error causing character
            if "-" in group:
                group = group.split("-")
            elif "'" in group:
                group = group.split("'")
            else:
                group = group.split(" ")

            # Join by an underscore to find file
            group = ("_").join(group)
            shutil.copy("{0}/{1}.fa.aln".format(alignments, group), output)


def make_directory(name):
    # Make output directory for core genome alignments
    try:
        os.mkdir(name)
        out = name
    except FileExistsError:
        ext = ''.join(random.choice(string.ascii_uppercase + string.digits) for _ in range(6))
        out = name + '_' + ext
        os.mkdir(out)
        print("Output directory already found, writing to {0}".format(out))
    
    return out

def gff_2_GeneID(gff):
    # Read GFF file for Gene ID
    with open("{0}/{1}".format(args.gffs, gff.name), "r") as file:
        for line in file:
            if "##" in line:
                continue
            if line.split()[2] == "CDS":
                try:
                    ID = line.split("=")[1].split("_")[0]
                except IndexError:
                    print(gff)
                break
    return ID


def alignment_parser(group, directory, gff_dict, final_out):
    # Read Alignment file (fasta) and sort into dictionary
    fasta_dict = {}

    try:
        file = open("{0}/{1}.fa.aln".format(directory, group))
    except FileNotFoundError:
            # Split by the error causing character
            if "-" in group:
                group = group.split("-")
            elif "'" in group:
                group = group.split("'")
            else:
                group = group.split(" ")

            # Join by an underscore to find file
            group = ("_").join(group)
            file = open("{0}/{1}.fa.aln".format(directory, group))

    for line in file:
        if line[0] == ">":
            ID = gff_dict[line.split("_")[0][1:]]
            fasta_dict.setdefault(ID, "")
        else:
            fasta_dict[ID] += line.rstrip()

    file.close()

    # Records list for Biopython write out
    records = []
    for i in fasta_dict.keys():
        cur_record = Bio.SeqRecord.SeqRecord(id = i, description="", seq=Bio.Seq.Seq(fasta_dict[i]))
        records.append(cur_record)

    # Write records to FASTA file
    Bio.SeqIO.write(records, "{0}/{1}.fa".format(final_out, group), "fasta")


def fastGEAR(final_out, directory):
    ############# Run fastGEAR on each individual gene #################
    print("Running fastGEAR on each renamed core gene alignment...")
    for alignment in os.scandir(final_out):
        # Make directory for gene and move to fastgear folder
        gene_dir = alignment.name[:-3]
        os.mkdir("{0}/{1}".format(directory, gene_dir))
        shutil.copy(alignment.path, "{0}/{1}/.".format(directory, gene_dir))

        # Run fastgear on core gene
        os.chdir("{0}/{1}/".format(directory, gene_dir))
        subprocess.run("fastGEAR ./{0} ./gene_fastgear /mnt/lustre/andam/shared/coopers_programs/fG_input_specs.txt".format(alignment.name), shell = True)
        os.chdir("../../")
        

def main(which_portion, portion_list, gff_dict):
    # Move core genes to new directory
    if which_portion == 0:
        label = "core"
    else:
        label = "accessory"
#    out = make_directory("{0}/sidekick_{1}_genome_alns".format(args.output, label))
#    file_finder(portion_list, args.alignments, out)

    # Rename gene alignments based on Gene ID
    print("Renaming ALN files based on protein IDs")
    final_out = make_directory("{0}/renamed_{1}_aln_files".format(args.output, label))
    for gene in portion_list:
        alignment_parser(gene, args.alignments, gff_dict, final_out)


    # <Optional> Run fastGEAR on genes
    if args.fastGEAR:
        directory = make_directory("{0}/{1}_gene_fastgears".format(args.output, label))
        fastGEAR(final_out, directory)



# Build pan-genome
genome_num = 0
gff_dict = {}
for file in os.scandir(args.gffs):
    gff_dict[gff_2_GeneID(file)] = file.name.split(".gff")[0]
    genome_num += 1

pan_genome = pan_genome_extraction(args.cluster_file, genome_num)

for num, portion in enumerate(pan_genome):
    if len(portion) > 0:
        main(num, portion, gff_dict)

