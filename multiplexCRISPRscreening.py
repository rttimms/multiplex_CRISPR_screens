#!/usr/local/bin/python

import os
import pandas as pd
from subprocess import call

docstring= """
DESCRIPTION
    Python script to analyse Illumina sequence data from multiplex CRISPR screening

REQUIRES:
    Raw gzipped fastq paired-end Illumina data in rawData folder
    Bowtie2 indexes and library fasta files in bowtie2index folder
    Cutadapt installed and in PATH
    Bowtie2 installed and in PATH

USAGE:
    python3 multiplexCRISPRscreening.py
"""

# paths to library files
SUBSTRATE_LIBRARY_FASTA_FILENAME = 'bowtie2index/substrates.fasta'
SGRNA_LIBRARY_FASTA_FILENAME = 'bowtie2index/guides.fasta'
SUBSTRATE_BT2_INDEX = 'bowtie2index/substrates'
SGRNA_BT2_INDEX = 'bowtie2index/guides'

# sequences flanking the substrate and sgRNA for trimming
SUBSTRATE_F = 'CGAATTCCGCTGCGC' #sequence upstream of the GPS substrate
SUBSTRATE_R = 'CGCATCGCTCGAGCC' #sequence downstream of the GPS substrate
SGRNA_F = 'CAGCATAGCTCTTAAAC' #reverse complement of the start of the tracRNA
SGRNA_R = 'GGTGTTTCGTCCTTTC' #reverse complement of end of U6 promoter

# minimum bowtie2 mapping quality
MIN_MAPQ_SCORE = 20


######################################################################
def read_fasta_library(fasta_file):

    '''
    Reads in a fasta file and returns a list of the library members
    '''

    lib = []
    with open(fasta_file, 'r') as infile:
        for line in infile:
            if line.startswith('>'):
                id = line.strip().split('>')[1]
                lib.append(id)
    print(f'{len(lib)} lines read in from library fasta file')

    return lib
######################################################################


#######################################################################
def find_fastq_data():

    '''
    Examines ./rawData folder and returns a list of raw gzipped
    fastq files for analysis
    '''

    input_fastqs = []

    for filename in os.listdir('./rawData'):
        if filename.endswith('_R1.fastq.gz'):
            input_fastqs.append(filename)

    return input_fastqs
######################################################################


###########################################################################################################################################
def trim_with_cutadapt(input_R1, input_R2, output_R1, output_R2, \
                       substrate_flank5, substrate_flank3, sgrna_flank5, sgrna_flank3, \
                       substrate_min_len=0, substrate_max_len=1000, sgrna_min_len=0, sgrna_max_len=1000):

    '''
    Calls cutadapt to trim constant flanking regions
    '''

    send2cutadapt = f"cutadapt -j 0 -m {substrate_min_len}:{sgrna_min_len} -M {substrate_max_len}:{sgrna_max_len} --discard-untrimmed \
                     -g {substrate_flank5}...{substrate_flank3} -G {sgrna_flank5}...{sgrna_flank3} \
                     -o {output_R1} -p {output_R2} {input_R1} {input_R2}"
    call(send2cutadapt, shell=True)
###########################################################################################################################################


###############################################################################################
def map_with_bowtie2(input, output, index):

    '''
    Calls Bowtie 2 to map trimmed fastq reads to the reference index
    '''

    send2bowtie2 = f"bowtie2 --no-head --reorder -p 20 -t -x {index} -U {input} -S {output}"
    call(send2bowtie2, shell=True)
##############################################################################################


##############################################################
def simplify_sam(input, output, cols):

    '''
    Calls cut to extract specific columns from sam file
    '''

    send2cut = f"cut -f{cols} {input} > {output}"
    call(send2cut, shell=True)
##############################################################


#######################################################################################
def generate_substrate_guide_counts_dataframe(substrate, counts, guide_library):

    # create dataframe with all the possible guides
    dfp = pd.DataFrame(index=guide_library)

    # add the counts as a column and the substrate name
    dfp['Counts'] = counts
    dfp['Substrate'] = substrate

    # tidy the df
    dfp.index.name = 'Guide'
    dfp = dfp.reset_index()
    dfp = dfp[['Substrate','Guide','Counts']]
    dfp.fillna(0, inplace=True)

    return dfp
#######################################################################################




def main():

    # make a new directory to store the output files
    if not os.path.exists(f"./rawCounts"):
        os.makedirs(f"./rawCounts")

    # read in the library details from the fasta files
    lib1 = read_fasta_library(fasta_file=SUBSTRATE_LIBRARY_FASTA_FILENAME)
    lib2 = read_fasta_library(fasta_file=SGRNA_LIBRARY_FASTA_FILENAME)

    # populate a list with all the gzipped fastq file names to be analysed
    input_fastqs = find_fastq_data()

    # perform the analysis for each of the fastq files
    for filename in input_fastqs:

        # establish the sample name and update terminal
        sample_name = filename.split('_R1.fastq.gz')[0]
        print(f"Starting analysis of sample: {sample_name}")

        # trim the reads using Cutadapt
        print(f"Removing constant flanking sequences using Cutadapt from sample: {sample_name}")
        trim_with_cutadapt(input_R1 = f"rawData/{filename}",
                           input_R2 = f"rawData/{filename.replace('_R1.fastq.gz', '_R2.fastq.gz')}",
                           output_R1 = f"{sample_name}_R1_trimmed.fastq",
                           output_R2 = f"{sample_name}_R2_trimmed.fastq",
                           substrate_flank5 = SUBSTRATE_F,
                           substrate_flank3 = SUBSTRATE_R,
                           sgrna_flank5 = SGRNA_F,
                           sgrna_flank3 = SGRNA_R,
                           substrate_min_len = 72,
                           substrate_max_len = 72,
                           sgrna_min_len = 20,
                           sgrna_max_len = 21,
                           )

        # map to the reference index with Bowtie 2
        print(f"Using Bowtie2 to map reads from sample: {sample_name} to reference indexes: {SUBSTRATE_BT2_INDEX} and {SGRNA_BT2_INDEX}")
        map_with_bowtie2(input = f"{sample_name}_R1_trimmed.fastq",
                         output = f"{sample_name}_R1.sam",
                         index = SUBSTRATE_BT2_INDEX,
                         )

        map_with_bowtie2(input = f"{sample_name}_R2_trimmed.fastq",
                         output = f"{sample_name}_R2.sam",
                         index = SGRNA_BT2_INDEX,
                         )

        # prune the relevant columns from the two sam files
        print('Concatenating the two sam files...')
        simplify_sam(input=f"{sample_name}_R1.sam",
                     output=f"{sample_name}_R1_simplified.sam",
                     cols='3,5')

        simplify_sam(input=f"{sample_name}_R2.sam",
                     output=f"{sample_name}_R2_simplified.sam",
                     cols='3,5')

        # merge the two files in a pandas dataframe
        print('Generating merged counts dataframe...')
        df1 = pd.read_csv(f"{sample_name}_R1_simplified.sam", sep='\t', header=None, names=['RNAME_R1','MAPQ_R1'])
        df2 = pd.read_csv(f"{sample_name}_R2_simplified.sam", sep='\t', header=None, names=['RNAME_R2','MAPQ_R2'])
        df = pd.concat([df1,df2], axis=1)

        # filter for mapped reads of sufficient quality
        df = df[(df['MAPQ_R1'] >= MIN_MAPQ_SCORE) & (df['MAPQ_R2'] >= MIN_MAPQ_SCORE)].reset_index(drop=True)
        print(f"Number of mapped pairs: {len(df)}")

        # count occurences of each substrate-sgRNA combination
        dfmaster = pd.DataFrame()
        for i, data in df.groupby('RNAME_R1'):
            
            # convert into a dataframe in the format Substrate-Guide-Counts
            counts = data['RNAME_R2'].value_counts()
            dfc = generate_substrate_guide_counts_dataframe(substrate=i, counts=counts, guide_library=lib2)
            dfmaster = pd.concat([dfmaster,dfc], ignore_index=True)

        # output the master df to file
        dfmaster.to_csv(f"rawCounts/{sample_name}_rawCounts.csv", index=False)
        print(f"Finished analysis of sample: {sample_name}")




if __name__ == '__main__':
    main()

