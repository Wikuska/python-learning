from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqUtils import gc_fraction
import re

def records_num(filename):
    records = 0
    for seq_record in SeqIO.parse(filename, "fasta"):
        records += 1
    return records

def records_lenght(filename):
    for seq_record in SeqIO.parse(filename, "fasta"):
        print(f"\n{seq_record.description}")
        print(len(seq_record))
