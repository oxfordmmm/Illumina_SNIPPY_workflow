#!/usr/bin/env python3
import sys
from Bio import SeqIO

def _yieldPassSeqs(fasta, percent):
    for record in SeqIO.parse(fasta, "fasta"):
        seq = str(record.seq)
        pN=(seq.count('N')/len(seq))*100
        print(record.id, pN)
        if pN < float(percent):
            yield record

seqs=_yieldPassSeqs(sys.argv[1], float(sys.argv[2]))

with open(sys.argv[3], 'w') as f:
    SeqIO.write(seqs, f, "fasta")