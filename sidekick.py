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
parser = argparse.ArgumentParser(description = "Use before HERO to extract individual alignment files for each core gene in pan-genome. Will also run fastGEAR individually on genes.", usage = "sidekick.py [options] alignments gffs")
parser.add_argument("alignments", help = "Filepath to individual gene alignments")
parser.add_argument("gffs", help = "Filepath to GFF files")
parser.add_argument("--cluster_file", metavar = '', default = "./clustered_proteins", help = "Filepath to the 'clustered_proteins' output of Roary [./clustered_proteins]")
parser.add_argument("--output", metavar = '', default = ".", help = "Filepath to output directory for core genome protein files. [Core_Genome_alignments]")
#parser.add_argument("--raxml", action = 'store_true', help = "Run RAxML on every core gene")
parser.add_argument("--fastGEAR", action = 'store_true', help = "Run fastGEAR on every core gene")
args = parser.parse_args()


def file_find(cluster_file, genome_num):
    # List of alignment filenames to pull
    core_cluster = []

    with open(args.cluster_file, "r") as cluster:
        for line in cluster:
            line = line.split(":")
        
            # Name of protein cluster
            group = line[0]
            if "/" in group:
                group = "_".join(group.split("/"))

            # Count number of occurances
            prots = line[1].split("\t")
            if len(prots) > math.floor(genome_num * 0.99):
                core_cluster.append(group)
    return core_cluster


def file_finder(core_cluster, alignments, output):
    # Copy alignment files from pan_genome_sequences folder to new folder
    for group in core_cluster:
        try:
            shutil.copy("{0}/{1}.fa.aln".format(alignments, group), output)
        except FileNotFoundError:

            # Split by the error causing character
            if "-" in group:
                group = group.split("-")
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

genome_num = 0
for file in os.scandir(args.gffs):
    genome_num += 1

core_cluster = file_find(args.cluster_file, genome_num)

print("{0} Groups found in core genome".format(len(core_cluster)))

out = make_directory("{0}/HiDe_Core_genome_alns".format(args.output))

file_finder(core_cluster, args.alignments, out)

########################## Parse GFF and ALN files  ###############################3
print("Renaming ALN files based on protein IDs")
def gff_dict_value(gff):
    # Read GFF file for ID
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


def alignment_parser(aln, gff_dict, final_out):
    # Read Alignment file (fasta) and sort into dictionary
    fasta_dict = {}
    with open("{0}/{1}".format(out, aln.name), "r") as file:
        for line in file:
            if line[0] == ">":
                ID = gff_dict[line.split("_")[0][1:]]
                fasta_dict.setdefault(ID, "")
            else:
                fasta_dict[ID] += line.rstrip()

    # Records dict for Biopython write out
    records = []
    for i in fasta_dict.keys():
        cur_record = Bio.SeqRecord.SeqRecord(id = i, description="", seq=Bio.Seq.Seq(fasta_dict[i]))
        records.append(cur_record)

    # Write records to FASTA file
    Bio.SeqIO.write(records, "{0}/{1}".format(final_out, aln.name[:-4]), "fasta")


gff_dict = {}

for file in os.scandir(args.gffs):
    gff_dict[gff_dict_value(file)] = file.name.split(".gff")[0]


final_out = make_directory("{0}/renamed_aln_files".format(args.output))

for file in os.scandir(out):
    alignment_parser(file, gff_dict, final_out)

def RAxML(final_out):
    ############# Run RAxML on all the renamed alignments ##############
    print("Running RAxML on each renamed alignment file.")
    for alignment in os.scandir(final_out):
        subprocess.run("srun raxmlHPC-PTHREADS -T 24 -s {0}/{1} -m GTRGAMMA -b 24961 -N 100 -p 14352 -n {2}".format(final_out, alignment.name, alignment.name[:-3]), shell = True)

    # Move files into appropriate folders
    try:
        os.mkdir("RAxML_info")
        os.mkdir("RAxML_trees")
        os.mkdir("RAxML_reduced_alns")
    except FileExistsError:
        print("RAxML Output Directories Already Made. This will override any trees already in the directory.")

    subprocess.run("mv RAxML_bootstrap.* RAxML_trees", shell = True)
    subprocess.run("mv RAxML_info.* RAxML_info", shell = True)
    subprocess.run("mv *.reduced RAxML_reduced_alns", shell = True)

    for file in os.scandir("RAxML_trees"):
        name = file.name.split(".")[1]
        os.rename(file.path, "{0}.newick".format(name))


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
        subprocess.run("fastGEAR ./{0} ./gene_fastgear ~/../shared/coopers_programs/fG_input_specs.txt".format(alignment.name), shell = True)
        os.chdir("../../")
        
#if args.raxml:
#    RAxML(final_out)

if args.fastGEAR:
    directory = make_directory("{0}/core_gene_fastgears".format(args.output))
    fastGEAR(final_out, directory)
