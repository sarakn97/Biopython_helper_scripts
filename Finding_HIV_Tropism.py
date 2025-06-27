#!/usr/bin/env python3
"""
Name: Sara K Nicholson
Title: Labeling HIV Sequence as X4 or R5 (X4 or R5-derived), returns dictionary of annotations
Date: February 24, 2025
Description:
"""

import sys
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import AlignIO
import subprocess

# First align sequences using subprocess module - subprocesses calls external packages
# (clustalw2 needs to be referenced or in path to be called here) - can use other MSAs
# subprocess.run(["/home/sara/Downloads/clustalw-2.1-linux-x86_64-libcppstatic/clustalw2", "-infile=/home/sara/GITHUB/Biopython_helper_scripts/sample_fasta.fa", "-outfile=/home/sara/GITHUB/Biopython_helper_scripts/aligned.aln"])

# Convert Alignment file to Fasta
# AlignIO.convert("/home/sara/GITHUB/Biopython_helper_scripts/aligned.aln", "clustal", "/home/sara/GITHUB/Biopython_helper_scripts/aligned.fa", "fasta")

# View Aligned Sequences
# alignment = AlignIO.read("/home/sara/GITHUB/Biopython_helper_scripts/sample_fasta.fa", "fasta")
# print(alignment)
# for record in alignment:
#     print(f"Sequence ID: {record.id}, Sequence: {record.seq}")

# second find divergence against a reference and have many if/then statements to conclude with label
refX4 = SeqIO.parse("/home/sara/GITHUB/Biopython_helper_scripts/NL4-3_ref.fa", "fasta")
refR5 = SeqIO.parse("/home/sara/GITHUB/Biopython_helper_scripts/BAL_ref.fa", "fasta")

# prepare references so that each nucleotide is its own object in list
for rec in refX4:
    seqX4 = rec.seq
x4 = list(seqX4)

for rec in refR5:
    seqR5 = rec.seq
r5 = list(seqR5)


def list_sequences(recs, length):
    """ return list of sequences from SeqIO object of which are of certain length or greater """
    seq_rec = SeqIO.parse(recs, "fasta")
    fasta = []
    for record in seq_rec:
        if len(record.seq) >= length:  # remove sequences with <= 60 nucleotides
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
    leng = int(sys.argv[2])

    fasta = list_sequences(records, leng)
    listed_seqs = list_list_sequences(fasta)
    annotation = annotate_seqs(listed_seqs)
    seq_ids = get_seq_ids(fasta)

    # create dictionary using the list of sequence IDs and the designated annotations
    keys = seq_ids
    values = annotation

    # create dictionary using dict()
    seq_dict = dict(zip(keys, values))

    # print dictionary of annotations
    print(seq_dict)
