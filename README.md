# Comp383 MiniProject
## Transcriptome-wide characterization of human cytomegalovirus in natural infection and experimental latency.

After infection, Human herpesvirus, also known as the Human cytomegalovirus (HCMV) remains latent within the body throughout life and can be reactivated at any time. This can be a life-threatening condition for immune-compromised patients. The recent paper of Cheng et al. 2017 (https://www.ncbi.nlm.nih.gov/pubmed/29158406) produced the transcriptomes of latent HCMV and as they enter and exit latency. 

The code in this repository functions as a pipeline to analyze the data from this paper. For an example, we use the transcriptomes (RNA-seq) of a single patient's viral load 2 days (2dpi) and 6 days post-infection (6dpi). The analysis pipeline can be broken down into three main steps:
1. Remove non-HCMV sequences by aligning the transcriptomes to the reference genome (HCMV.fa) and only keeping what aligns.
2. Analyze how expression levels of various genes differ between two provided samples.
3. Identify the specific strain, by assembling transcriptomes into one giant from both files and BLAST-aligning it to all nucleotides. 

## Data
A reference sequence of the HCMV genome is available for download on this repository (HCMV.fa). However, the transcriptome data sets are too large to store on github. Instead, the script handles SRA IDs as the file input for the transcripts and uses WGET to download them from the NCBI database. All other files will be generated locally on your computer, in the "OptionB_Jasen_Jackson" folder.

## Software Dependencies
This script depends on the following software tools: 
* Bowtie2
* Samtools
* TopHat
* Cufflinks 
* Spades
* Biopython
* wget for mac users ("brew install wget").
These must be installed to your path, for the script to work.

## Options
* *-m*, *--mode:* specifies which mode to run the program in ("full", "2" or "3")
  * *"full"* mode runs all three steps. 
  * *"2"* runs the differential expression analysis (1 & 2)
  * *"3"* runs the strain comparison (1 & 3)
*  -g,--genome: fasta file containg the HCMV reference genome (must be .fasta/.fa)
*  -i,--input: list of SRA *ids* for the transcriptome data
*  -s,--sample: (optional) a user-provided name for each sample. Must be in the same order as the samples!

## walkthroughs

To run everything from the beginning, with my original settings, just clone this repository and run the following command:

    python3 main.py -m full -g HCMVA.fasta --g HCMV.fa -i SRR5660044,SRR5660045 -s 2dpi,6dpi


  



