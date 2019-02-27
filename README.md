# comp383
Comp 383 Mini Project Code

The recent paper of Cheng et al. 2017 (https://www.ncbi.nlm.nih.gov/pubmed/29158406) produced the transcriptomes of latent HCMV and as they enter and exit latency. The code in this repository functions as a pipeline to analyze the data from this paper. 

# software dependencies
This script depends on the following software tools. These must be installed to your path, for the script to work.
* Bowtie2
* Samtools
* TopHat
* Cufflinks
* spades
* Biopython

# data 
A reference sequence of the HCMV genome is available on this repository (HCMV.fa). However, the transcriptome data sets are too large to store on github. Instead, this script takes as input 2 SRA IDs from the study (for comparisons) and downloads them from the NCBI database. All other files will be generated locally on your computer, in the "OptionB_Jasen_Jackson" folder.

# overview
This analysis can be broken down into three main steps:
1. align the transcriptomes to the reference genome (HCMV.fa), and remove anything that does not align.
2. analyze different expression levels between the provided transcriptomes using TopHat and CuffLinks.
3. assemble both transcriptomes using spades, and blast that against the non-repeat nucleotide.

# instructions
The main.py script takes 4 arguments:
*  '-m','--mode': specifies which mode to run the program in ("full", "2" or "3").
  * "full" mode runs all three steps. "2" runs steps 1 and 2. "3" runs steps 1 and 3. 
*  '-g','--genome': fasta file containg the HCMV reference genome (must be .fasta/.fa)
*  '-i','--input': list of SRA *ids* for the transcriptome data
*  '-s','--sample': (optional) a user-provided name for each sample. Must be in the same order as the samples!

# command line instruction


