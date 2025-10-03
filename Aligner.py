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

