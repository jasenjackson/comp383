#!/usr/bin/python
'''
Title: HCMV Mini-project for Comp 383 (miniProjectDefs.py)
Author: Jasen M. Jackson, '19
Date: 2/20/19-

This scripts contains the function definitions for this project's python
wrapper.
'''
import os
from itertools import (takewhile,repeat)
from Bio import SeqIO
from Bio.Blast import NCBIWWW
import csv

def line_count(fi): ## TODO: better alternative?
    f = open(fi, 'rb')
    bufgen = takewhile(lambda x: x, (f.raw.read(1024*1024) for _ in repeat(None)))
    seq_count = int(sum( buf.count(b'\n') for buf in bufgen )/4)
    return seq_count

## Download and split list of SRA files
def sra_downloader(*inputs):
    from sys import platform
    if platform == "linux" or platform == "linux2": download_command = "wget -q"
    if platform == "darwin": download_command = "wget " # TODO: change to CURL.. get rid of WGET dependency
    for input in inputs:
        if os.path.isfile(input+"_1.fastq") and os.path.isfile(input+"_2.fastq"):
            print(input+" already downloaded. Skipping!")
        else:
            url = "\"ftp://ftp.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/SRR\""
            url += input[3:6]+'/'+input+'/'+input+'.sra'
            print("Downloading: "+input+" from url: "+url); os.system(download_command+url)
            print("Splitting "+input+".sra into 2 paired FASTQ files"); os.system("fastq-dump --split-3 "+input+".sra")

## Run Bowtie2 build function (i: .FASTA, o: bt2 index)
def bowtie_build(input,out):
    print("Building index for "+input+"...")
    bowtie_command = 'bowtie2'
    input_exists = os.path.isfile(input)
    if input_exists:
        bowtie_command += '-build --quiet -f '+input+' '+out
        print(bowtie_command)
    else:
        print("Error: could not locate "+input+"")
        return False
    os.system(bowtie_command)

## Run Bowtie2 alignment (idx: bt2 index, o: .SAM, i: paired FASTQs)
def bowtie_align(idx,input,out):
    r1, r2 = input+"_1.fastq", input+"_2.fastq"
    input_exists = os.path.isfile(r1) and os.path.isfile(r2)
    index_exists = os.path.isfile(idx+".1.bt2")
    ## Map input to index, if no mapped fastq file exists
    if not os.path.isfile(out+".mapped.1.fastq"):
        print("Aligning "+input+" to "+idx+" with Bowtie2")
        ## bowtie2 --al-conc -x <idx> -1 <m1> -2 <m2> -S
        if input_exists and index_exists:
            bowtie_command = "bowtie2 --al-conc "+out+".mapped.fastq -x "+idx+" -1 "+r1+" -2 "+r2+" -S "+out; print(bowtie_command)
        elif index_exists: print("Error: could not locate paired FASTQ files: "+r1+","+r2); return False
        elif input_exists: print("Error: could not locate index file: "+idx); return False
        else: print("Error: could not locate any files provided."); return False
        os.system(bowtie_command)
        print("Done!")
    else: print(input+" has already been mapped to "+idx+" with Bowtie2. Skipping!")
    ## Write to log file ..
    before = str(line_count(r1))
    after = str(line_count(out+".mapped.1.fastq"))
    log_message = out+" had "+before+" pairs before Bowtie2 filtering and "+after+" pairs after."
    print(log_message)
    return log_message

## Run TopHat alignment (idx: bt2 index, o: .SAM, i: paired FASTQs)
def tophat_align(idx,input,out):
    print("Aligning "+input+" to "+idx+" with TopHat")
    r1, r2 = input+"_1.fastq", input+"_2.fastq"
    input_exists = os.path.isfile(r1) and os.path.isfile(r2)
    index_exists = os.path.isfile(idx+".1.bt2")
    out_dir = out+"_tophat"
    ## Map input to index, if no mapped fastq file index_exists
    if not os.path.isfile(out_dir+"/accepted_hits.bam"):
        if input_exists and index_exists:
            tophat_command = "tophat -p 2 -o "+out_dir+" "+idx+" "+r1+" "+r2
        elif index_exists: print("Error: could not locate paired FASTQ files: "+r1+","+r2); return False
        elif input_exists: print("Error: could not locate index file: "+idx); return False
        else: print("Error: could not locate any files provided."); return False
        print(tophat_command)
        os.system(tophat_command)
    else: print(input+" has already been mapped to "+idx+" with TopHat. Skipping!")

## Run Cufflinks (accepted_hits.bam >> transcripts.gtf)
def run_cufflinks(input):
    print("Running Cufflinks on "+input+"...")
    output = input+"/transcripts.gtf" ## expected output
    if not os.path.isfile(output):
        os.system("cufflinks -p 4 -o " +input+" "+input+"/accepted_hits.bam")
        if os.path.isfile(output):
            print("Transcript assembly for "+input+"was a success. Results are located here: "+output)
            return output
        else:
            print("\tUh oh! Cufflinks failed. *shrug* ")
            return False
    else:
        print("\tCufflinks has already assembled transcripts for "+input+" in Skipping!")
        return output

## Run Cuffmerge (transcript.gtf >> merged.gtf)
def run_cuffmerge(genome):
    print("Merging transcripts...")
    output = "merged_asm/merged.gtf"
    assemblies_exists = os.path.isfile("assemblies.txt")
    merged_asm_exists = os.path.isfile(output)
    if not merged_asm_exists:
        if os.path.isfile("assemblies.txt"):
            cuffmerge_command = "cuffmerge -s "+genome+" -p 2 assemblies.txt"
            print(cuffmerge_command)
            os.system(cuffmerge_command)
        else: print("Could not locate assemblies.txt")
    else: print("\tTranscripts in assemblies.txt have already been merged. Skipping!")


## Prepare tophat files for cuffdiff
def make_merged_transcriptomes(genome,*inputs):
    if not os.path.isdir("merged_transcriptome"):
        assemblies = open("assemblies.txt", "w+")
        ## Assemble transcripts and collect in assemblies.txt ...
        for input in inputs:
            if os.path.isfile(input+"/accepted_hits.bam"):
                cl_output = run_cufflinks(input)
                assemblies.write(cl_output+"\n")
        ## Merge files together
        run_cuffmerge(genome)
    else: print("Error: Could not locate merged transcript file")

## Run cuffdiff (assumes pair-wise comparison)
def run_cuffdiff(genome, t1, t2):
    print("## Computing differences in transcript expression using Cuffdiff...")
    output = "cuffdiff_"+t1+"_"+t2
    if not os.path.isfile(output):
        cuffdiff_command = "cuffdiff -p 4 "+output+" -b "+genome+" -L "+t1+","+t2+" -u merged_asm/merged.gtf"+" ./"+t1+"/accepted_hits.bam"+",./"+t2+"/accepted_hits.bam"
        print(cuffdiff_command)
    else: print("Cuffdiff has already compares these transcriptomes: "+output)

def filter_genes(file, r, val):
    print("Compiling list of significant genes... ")
    with open(file,'r') as gene_file:
        log_message = ""
        reader = csv.reader(gene_file, delimiter='\t')
        for row in reader:
            if row[r] == val:
                log_message += " ".join(str(x) for x in row)+"\n"
    gene_file.close()
    return log_message

## Run spades (assumes "pair-wise" assembly)
def spades_assemble(out, *inputs):
    print("## Assembling "+str(inputs)+" using spades")
    in1_exists = os.path.isfile(inputs[0])
    in2_exists = os.path.isfile(inputs[2])
    output_exists = os.path.isdir(out)
    if not output_exists:
        if in1_exists and in2_exists:
            spades_command = "spades -k 55,77,99,127 -t 4 --only-assembler --pe1-1 "+inputs[0]+" --pe1-2 "+inputs[1]+" --pe2-1 "+inputs[2]+" --pe2-2 "+inputs[3]+" -o "+out+"/"
            print(spades_command)
            os.system(spades_command)
        else: print("Could not find input mapped fastqs")
    else: print("Files have already been assembled: "+out)

## Output filtered contigs with length > n (saves as contigs.filtered.n.fasta)
def filter_contigs(assembly, n):
    print("## Filtering contigs in "+assembly+" by size "+str(n))
    in_path = assembly+"/contigs.fasta"
    outfile = open("contigs."+str(n)+".fasta", "w+")
    num_contigs, size = 0,0
    for record in SeqIO.parse(in_path, "fasta"):
        s = record.id.find("length_")
        l = record.id[s+7:]; l = int(l[:l.find("_")])
        if l > n:
            outfile.write(">"+record.id+"\n"+str(record.seq)+"\n")
            num_contigs += 1
            size += l
    return num_contigs, size

## Concatenate FASTA into one string, adding k Ns inbetween
def concatenate_fasta(fasta, k):
    print("Concatenating "+fasta+" with "+str(k)+" 'N's in between")
    fasta_exists = os.path.isfile(fasta)
    if fasta_exists:
        output_name = fasta[:-5]+"concatenated.fasta"
        output_file = open(output_name, "w+")
        output_file.write(">"+output_name+" concatenated with "+str(k)+"stretches"+"\n")
        n_string = 'N'*k
        for record in SeqIO.parse(fasta, "fasta"):
            output_file.write(str(record.seq))
            output_file.write(n_string)
        ## remove last 50 Ns
        output_file.seek(0); output_file.readline()
        concatenated = output_file.readline(); output_file.seek(0)
        output_file.readline(); output_file.write(concatenated[:-50])

## Run BLASTN on FASTA file against NR database, report top n results
def blastn_nr(query, n):
    print("## Performing a BLAST search against the non-redundant nucleotide database using "+query+"as query")
    fasta_string = open(query).read()
    results_handle = NCBIWWW.qblast("blastn", "nr", fasta_string)
    blast_records = NCBIXML.read(result_handle)
    blast_records = list(blast_records)
    print("Printing top BLAST hits...")
    log_message = ""
    for i in range(0,n):
        log_message += "Sequence title: "+blast_records[i].title+"\n"
        log_message += "Alignment length: "+blast_records[i].length+"\n"
        log_message += "HSP identities"+blast_records[i].hsp.identities+"\n"
        log_message += "HSP query: "+blast_records[i].hsp.query+"\n"
        log_message += "HSP bits: "+blast_records[i].hsp.bits+"\n"
        log_message += "HSP Expect scores: "+blast_records[i].hsp.expect+"\n"
    print(log_message)
    return(log_message)
