"""
Pairwise Global Alignment
url: https://rosalind.info/problems/need/

Problem

An online interface to EMBOSS's Needle tool for aligning DNA and RNA strings can be found here.

Use:

The DNAfull scoring matrix; note that DNAfull uses IUPAC notation for ambiguous nucleotides.
Gap opening penalty of 10.
Gap extension penalty of 1.
For our purposes, the "pair" output format will work fine; this format shows the two strings aligned at the bottom of the output file beneath some statistics about the alignment.

Given: Two GenBank IDs.
Return: The maximum global alignment score between the DNA strings associated with these IDs.
"""
# found the easiest implementation ever on the internet
from Bio import Entrez, SeqIO, pairwise2

with open("./data/rosalind_need.txt", "r") as f:
    genbank_ids = ", ".join(f.readline().strip().split())
Entrez.email = "example@example.com"
handle = Entrez.efetch(db="nucleotide", id=[genbank_ids], rettype="fasta")
records = list(SeqIO.parse(handle, "fasta"))
print(pairwise2.align.globalms(records[0].seq, records[1].seq, 5, -4, -10, -1)[0][2])

# https://biopython.org/DIST/docs/api/Bio.pairwise2-module.html