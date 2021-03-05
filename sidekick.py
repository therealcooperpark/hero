#! /usr/bin/env python3
'''
Convert Roary pan_genome_sequences FASTA headers into genome based headers
Optional ability to run fastGEAR on each new alignment
'''

import argparse
import multiprocessing
import os
import random
import shutil
import string
import subprocess
import sys
import time


def parse_args():
    parser = argparse.ArgumentParser(description = 'Use before HERO to convert Roary FASTA alignment headers to genome names',
                                     usage='sidekick.py [options] gff_table')
    parser.add_argument('gff_table', help='Tab-delimited file of GFF file location and associated genome for renaming')
    parser.add_argument('--alns', default='./pan_genome_sequences', help='Filepath to Roary pan_genome_sequences directory (requires -z argument) [./pan_genome_sequences]')
    parser.add_argument('--output', default='sidekick_genes', help='Output directory name [sidekick_genes]')
    parser.add_argument('--cpus', default=1, type=int, help='Number of CPUS to use [1]')
 
    fastgear_args = parser.add_argument_group('fastGEAR')
    fastgear_args.add_argument('--fastgear', help='Filepath to "run_fastGEAR.sh" script provided by fastGEAR. Must be used with --mcr')
    fastgear_args.add_argument('--mcr', help='Filepath to MCR executable. Must be used if using --fastgear')
    fastgear_args.add_argument('--fgout', default='fastgear_genes', help='Output directory name for fastgear runs [fastgear_genes]')
    fastgear_args.add_argument('--iters', default='15', help='Number of iterations [15]')
    fastgear_args.add_argument('--bounds', default='10', help='Upper bound for number of clusters [10]')
    fastgear_args.add_argument('--partition', default='-', help='File containing a partition for strains [NA]')
    fastgear_args.add_argument('--fg_output', default='1', help='1=reduced output, 0=complete output')
    return parser.parse_args()


def make_directory(name):
    '''
    Make output directory.
    If name already taken, make a random extension
    '''

    try:
        os.mkdir(name)
        out = name
    except FileExistsError:
        ext = ''.join(random.choice(string.ascii_uppercase + string.digits) for _ in range(6))
        out = name + '_' + ext
        os.mkdir(out)
        print('Output directory already found, writing to {0} instead'.format(out))
    return out


def parse_gff_table(gff_table):
    ''' 
    Parse each GFF file for associated CDS Protein code
    and assign code to the genome listed in user table
    '''

    cds_to_genome = {}

    with open(gff_table, 'r') as infile:
        for line in infile:
            line = line.strip().split()
            gff, genome = line[0], line[1]
            # Try to open the GFF
            try:
                gff_file = open(gff, 'r')
            except FileNotFoundError:
                print('GFF file not found, check location: {0}'.format(gff))
                sys.exit(1)

            # Get CDS value from the GFF
            # Assign CDS ID to dict with genome as key
            for gff_line in gff_file:
                if gff_line[0] == '#':
                    continue

                if gff_line.split()[2] == 'CDS':
                    try:
                        ID = gff_line.split('=')[1].split('_')[0]
                    except IndexError:
                        print('Naming error in {0}'.format(gff))
                        sys.exit(1)
                    cds_to_genome[ID] = genome
                    gff_file.close()
                    break
    return cds_to_genome


def rename_alignment(alignment_path):
    '''
    Create a copy of the fasta file with headers renamed
    using the cds_to_genome dictionary
    '''

    new_aln_name = alignment_path.split('/')[-1][:-4]
    with open(alignment_path, 'r') as old_aln:
        with open('{0}/{1}'.format(output, new_aln_name), 'w') as new_aln:
            for line in old_aln:

                # If header, rename and write to new
                if line[0] == '>':
                    try:
                        old_head = line.strip().split('_')[0][1:]
                        new_head = cds_dict[old_head]
                        new_aln.write('>{0}\n'.format(new_head))
                    except KeyError:
                        print('{0} has no genome equivalent in alignment {1}! Make sure all GFFs were provided'.format(old_head, alignment_path))
                        os._exit(1)
                # Else, just write the seq over
                else:
                    new_aln.write(line)


def make_specs_file(iters, bounds, partition, fg_output):
    '''
    Make fastgear specs file
    '''

    with open('fG_input_specs.txt', 'w') as output:
        output.write(iters + " # Number of iterations\n")
        output.write(bounds + ' # Upper bound for the number of clusters (possibly multiple values)\n')
        output.write('0 # Run clustering for all upper bounds (0=no/1=yes)\n')
        output.write(partition + ' # File containing a partition for strains\n')
        output.write(fg_output + ' # 1=produce reduced output, 0=produce complete output\n')
    return './fG_input_specs.txt'


def run_fastgear(fg_path, mcr_path, sidekick_dir, fg_dir, spec_file, cpus):
    '''
    Run fastgear on each renamed gene alignment
    '''

    # Get list of filepaths to give to parallel
    with open('renamed_gene_filepaths.txt', 'w') as infile:
        for gene in os.scandir(sidekick_dir):
            infile.write(gene.path + '\n')

    # Assume ~ 4 cpus per job. Also, just don't submit all jobs at once!
    j = cpus // 4 if cpus > 8 else 1

    # Build string to run everything and then do it!
    run_string = 'parallel -j {0} '.format(j)
    run_string += '"{0} {1}'.format(fg_path, mcr_path)
    run_string += ' ./{} ./' + fg_dir + '/{/.}/{/.} ' + spec_file + '"'
    run_string += ' < renamed_gene_filepaths.txt'
    print(run_string)
    subprocess.run(run_string, shell=True)


def make_hero_input(gene_output, fgout):
    '''
    Make the input file necessary for HERO to run
    '''

    with open('hero_input.txt', 'w') as outfile:
        for gene in os.scandir(gene_output):
            gene_name = gene.name.split('.')[0]
            fg_path = fgout + '/' + gene_name + '/'
            if os.path.isdir(fg_path):
                outfile.write('{0}\t{1}\t{2}\n'.format(gene_name,
                                                       gene.path,
                                                       fg_path))
            else:
                print('fastGEAR for gene {0} not found!'.format(gene.name))

def main():
    # Get arguments from command line
    args = parse_args()

    if args.fastgear:
        # Verify that filepath to fastgear exists
        if not os.path.exists(args.fastgear):
            print('Invalid filepath to run_fastgear.sh', flush = True)
            sys.exit(1)

        # Verify that MCR executable given
        if not args.mcr:
            print('MCR executable filepath required to run fastgear.', flush = True)
            sys.exit(1)

        # Verify that MCR executable path exists
        if not os.path.exists(args.mcr):
            print('Invalid filepath to MCR executable', flush = True)
            sys.exit(1)

    # Iterate over each gene alignment and build a renamed one
    args.output = make_directory(args.output)
    
    # Get dictionary of CDS IDs to genomes
    cds_to_genome = parse_gff_table(args.gff_table)

    # Set global variables for cpus
    global output
    output = args.output
    global cds_dict
    cds_dict = cds_to_genome

    # Get list of alignments to rename
    gene_alns = []
    for gene in os.scandir(args.alns):
        gene_alns.append(gene.path)
    
    # Iterate over each gene alignment and rename the headers
    # BONUS: We can use multithreading!
    print('Starting rename process!')
    pool = multiprocessing.Pool(processes=args.cpus)
    status = pool.map(rename_alignment, gene_alns)
    pool.close()
    pool.join()

    # Make fastGEAR specs file and run fastGEAR on each gene
    print('Starting batch fastGEARs!')
    if args.fastgear:
        args.fgout = make_directory(args.fgout)
        spec_file = make_specs_file(args.iters, args.bounds, args.partition, args.fg_output)
        run_fastgear(args.fastgear, args.mcr, output, args.fgout, spec_file, args.cpus)
        make_hero_input(output, args.fgout)
if __name__ == '__main__':
    main()
