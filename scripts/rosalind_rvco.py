"""
Complementing a Strand of DNA.
url: https://rosalind.info/problems/rvco/

Recall that in a DNA string s, 'A' and 'T' are complements of each other, as are 'C' and 'G'. Furthermore, the reverse complement of s is the string sc formed by reversing the symbols of s and then taking the complement of each symbol (e.g., the reverse complement of "GTCA" is "TGAC").

The Reverse Complement program from the SMS 2 package can be run online here.

Given: A collection of n (n≤10) DNA strings.
Return: The number of given strings that match their reverse complements.
"""
# completely disregarded the SMS package coz WHY
from Bio.Seq import Seq
from Bio import SeqIO
# from Bio.Alphabet import IUPAC
from typing import List
from tools import parse_fasta

def rc_match(seqs: List[str]) -> int:
    count = 0
    for seq in seqs:
        s = Seq(seq)
        rc_s = s.reverse_complement()
        if s == rc_s:
            count += 1
    return count

if __name__ == '__main__':
    with open('../data/rosalind_rvco.txt') as f:
        sequences = [str(record.seq) for record in SeqIO.parse(f, 'fasta')]
    counts = rc_match(sequences)
    print(counts)
