#!/usr/bin/python
'''
Title: HCMV Mini-project for Comp 383 (main.py)
Author: Jasen M. Jackson, '19
Date: 2/20/19-

This scripts contains the main workflow for the Comp383 mini-projectself.
'''
from miniProjectDefs import *
import sys, getopt
import os
import argparse

## bowtie step: make sure reads map to the genome. remove those that don't
def bt(genome, genome_name, in_files, samples, log_message):
    print("\n## Downloading transcriptomes...")
    sra_downloader(in_files[0],in_files[1])
    print("\n## Creating a genome index from "+ genome)
    bowtie_build(genome, genome_name)
    print("\n## Mapping reads to genome index with Bowtie2...")
    log_message += bowtie_align(genome_name,in_files[0],samples[0]) + "\n"
    log_message += bowtie_align(genome_name,in_files[1],samples[1]) + "\n"
    print(log_message)
    return log_message

## tophat step: compare how expression is different between two samples (requires bt step)
def th(genome, genome_name, in_files, samples, log_message):
    print("\n## Mapping reads to genome index with TopHat...")
    tophat_align(genome_name,in_files[0],samples[0])
    tophat_align(genome_name, in_files[1], samples[1])
    print("\n## Computing change in gene expression between 2dpi & 6dpi..")
    make_merged_transcriptomes(genome, samples[0], samples[1])
    run_cuffdiff(genome, samples[0], samples[1])
    log_message += filter_genes("cds.diff", 12, "yes")
    print(log_message)
    return log_message

## spades step: which publically available strain is most similar to the patient's sample? (requires bt step)
def sp(genome, genome_name, in_files, samples):
    print("Assembling Bowtie2 reads in to: /hcmv_assembly/ ...")
    assembly_name = genome_name+"_assembly"
    fq_1_1 = in_files[0]+".mapped.1.fastq"; fq_1_2 = in_files[0]+".mapped.2.fastq"
    fq_2_1 = in_files[1]+".mapped.1.fastq"; fq_2_2 = in_files[1]+".mapped.2.fastq"
    spades_assemble(assembly_name, fq_1_1, fq_1_2, fq_2_1, fq_2_2)
    num_contigs, bp = filter_contigs(assembly_name, 1000)
    concatenate_fasta("contigs.1000.fasta", 50)
    log_message = ("There are "+str(num_contigs)+" contigs >1000 in the assembly.\n")
    log_message+= ("There are "+str(bp)+"bp in the assembly.\n")
    log_file.write(log_message); print(log_message)
    log_message = blastn_nr("contigs.1000.concatenated.fasta", 10)
    log_file.write(log_message)

if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('-m', '--mode', required = True, type=str)
    parser.add_argument('-g', '--genome', required = True, type=str)
    parser.add_argument('-i', '--input', required = True, type=str)
    parser.add_argument('-s', '--samples', type=str)
    args = parser.parse_args()
    samples_flag = (args.samples is not None)

    print("\n---User Preferences-----")
    print("Mode: "+args.mode) ## full, BT, BTandTH,BTandSP
    print("Genome (fasta file): "+args.genome)
    print("Input files: "+args.input)
    if samples_flag: print("Samples: "+args.samples)


    ## check genome for .fasta/.fa ending, get name of genome
    fasta_check = (args.genome[-6:] == ".fasta")
    fa_check = (args.genome[-3:] == ".fa")
    if not (fasta_check or fa_check):
        print("Error: Genome file must be in the FASTA format. (.fasta/.fa)")
        sys.exit()
    elif fasta_check: genome = args.genome; genome_name = args.genome[:-6]
    elif fa_check: genome = args.genome; genome_name = args.genome[:-3]

    ## split up input files/sample files and error check (only takes two arguments right now)
    in_files = args.input.split(",")
    if len(in_files) == 1: print("Error: Must provide at least two sequence IDs."); sys.exit()
    if len(in_files) > 2: print("Error: Must provide only two sequence IDs."); sys.exit()

    ## split up sample names into list. If none provided, use SRRxxxxx
    if samples_flag:
        samples = args.samples.split(",")
        if len(in_files) == 1: print("Error: Must provide at least two sample IDs."); sys.exit()
        if len(in_files) > 2: print("Error: Must provide only two sample IDs."); sys.exit()
    else: samples = [args.input[0], args.input[1]]

    ## setting up the "output" folder
    if not os.path.isdir("OptionB_Jasen_Jackson/"):
        print("\n## Made OptionB_Jasen_Jackson")
        os.system("mkdir OptionB_Jasen_Jackson")
    if os.path.isfile(genome):
        print("## Copying "+genome+" to OptionB_Jasen_Jackson/"+genome)
        os.system("cp "+genome+" OptionB_Jasen_Jackson/"+genome)
    print("## Moving to OptionB_Jasen_Jackson\n")
    prev_wd = os.getcwd()
    os.chdir("OptionB_Jasen_Jackson")

    ## get to work!
    log_message = ""
    log_message += bt(genome, genome_name, in_files, samples, log_message) # step 1 (bt = bowtie) // in all options
    if args.mode == "2" or args.mode == "full":
        log_message += th(genome, genome_name, in_files, samples, log_message) # step 2 (th = tophat)
    if args.mode == "3" or args.mode == "full":
        log_message += sp("HCMV.fa", "hcmv", in_files, samples, log_message) # step 3 (sp = spades)

    ## finish up
    print("## Writing to log file...")
    log_file = open("OptionB.log", "w+") ## start log file.
    log_file.write(log_message)
    print("## Finishing up...deleting unnecessary files")
    os.system("rm "+genome); print("## Removed copied genome folder")
    #os.system("rm *.sra")
    ## remove unncessary files
    print("")
    os.chdir(prev_wd)

    print("## All done! All generated data was stored in OptionB_Jasen_Jackson/")
