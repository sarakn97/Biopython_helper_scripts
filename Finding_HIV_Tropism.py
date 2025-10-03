#!/usr/bin/env python3
"""
Title: Labeling HIV Sequence as X4 or R5-derived tropism, returns dictionary of annotations
Created By : Sara Nicholson
Date: February 24, 2025
Description: This program can be run on the command line by providing:
    [1] a fasta file with your sequences
    [2] path to output alignment file after clustal is run on fasta
    [3] path to fasta file after alignment file is converted to fasta
    [4] a length at which you would like to filter sequences at
"""
import sys
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import AlignIO
import subprocess

# First align sequences using subprocess module - subprocesses calls external packages
# (clustalw2 needs to be referenced or in path to be called here) - can use other MSAs
def clustal_align(infile, outfile, outfileFA):
    subprocess.run(["/home/sara/Downloads/clustalw-2.1-linux-x86_64-libcppstatic/clustalw2", "-infile=" + infile, "-outfile="+ outfile])
    # Convert Alignment file to Fasta
    outFA = AlignIO.convert(outfile, "clustal", outfileFA, "fasta")
    return outFA

# View Aligned Sequences
# alignment = AlignIO.read("/home/sara/GITHUB/Biopython_helper_scripts/aligned.fa", "fasta")
# print(alignment)
# for record in alignment:
#    print(f"Sequence ID: {record.id}, Sequence: {record.seq}")

# LOAD REFERENCES
refX4 = SeqRecord(
                Seq("TRPNNNTRKSIRIQRGPGRAFVTIGKI-GNMRQAHCNISRAKWNATLKQIASKLREQFGNNKTIIFKQSSGGDPEI"),
                id="NL43",
                name="REF_X4",
                description="sequence",
            )

refR5 = SeqRecord(
                Seq("TRPNNNTRKSINI--GPGRAFYTTGEIIGDIRQAHCNLSRAKWNDTLNKIVIKLREQFG-NKTIVFKHSSGGDPEI"),
                id="BAL",
                name="REF_R5",
                description="sequence",
            )

# prepare references so that each nucleotide is its own object in list
x4 = list(refX4.seq)
r5 = list(refR5.seq)

def list_sequences(recs, length):
    """ return list of sequences from SeqIO object of which are of certain length or greater """
    seq_rec = SeqIO.parse(recs, "fasta")
    fasta = []
    for record in seq_rec:
        if len(record.seq) >= length:  # remove sequences with <= # of nucleotides specified
            record_final = SeqRecord(
                record.seq,
                id=record.id,
                name="ASV",
                description="sequence",
            )
            fasta.append(record_final)  # add sequence to list
    return fasta


def list_list_sequences(fa):
    """ return list of list of sequences """
    listed_seqs = []
    for record in fa:
        record = list(record.seq)
        listed_seqs.append(record)

    return listed_seqs


def annotate_seqs(seqs):
    """ score alignments to references and determine if X4 or R5 sequence based on score"""
    scoresx4 = []
    scoresr5 = []

    for seq in seqs:
        score = 0  # set base score zero
        for i in range(len(seq)):  # iterate through each position in sequence
            if seq[i] != x4[i]:  # score based on if it matches the x4 position (higher score = less similarity)
                score += 1
        scoresx4.append(score)

    for seq in seqs:
        score = 0
        for i in range(len(seq)):
            if seq[i] != r5[i]:
                score += 1
        scoresr5.append(score)

    annotation = []  # open list
    for i in range(len(scoresx4)):  # iterate
        if scoresx4[i] > scoresr5[i]:
            trop = "R5"  # if the x4 score is greater than the r5 score then that sequence is designated R5
        elif scoresx4[i] == scoresr5[i]:
            trop = "NA"  # if the scores are equal then the sequence tropism is indistinguishable
        else:
            trop = "X4"  # otherwise, the sequence is designated X4
        annotation.append(trop)  # add annotation designation to list

    print(scoresx4)
    print(scoresr5)
    print(annotation)
    return annotation


def get_seq_ids(fa):
    seqs =[]
    for record in fa:
        record = record.id
        seqs.append(record)

    return seqs


if __name__ == "__main__":
    # Initializations
    records = str(sys.argv[1])
    outfile = str(sys.argv[2])
    outfileFA = str(sys.argv[3])
    leng = int(sys.argv[4])

    # align sequences
    clustal_align(records, outfile, outfileFA)
    # List sequences and list nucleotides within sequences
    fasta = list_sequences(outfileFA, leng)
    listed_seqs = list_list_sequences(fasta)
    # Adjust for any discrepancies between length of seqs and references by adding blanks to references
    max_len = max(len(sublist) for sublist in listed_seqs)
    while len(x4) < max_len:
        x4.append("-")
    while len(r5) < max_len:
        r5.append("-")
    # Compare Sequences by scoring alignments to reference and determining tropism
    annotation = annotate_seqs(listed_seqs)
    # Collect Sequence IDs to create dictionary
    seq_ids = get_seq_ids(fasta)

    # create dictionary using the list of sequence IDs and the designated annotations
    keys = seq_ids
    values = annotation

    # create dictionary using dict()
    seq_dict = dict(zip(keys, values))

    # print dictionary of annotations
    print(seq_dict)
