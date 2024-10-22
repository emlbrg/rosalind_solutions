"""
Complementing a Strand of DNA
url: https://rosalind.info/problems/revc/

In DNA strings, symbols 'A' and 'T' are complements of each other, as are 'C' and 'G'.

The reverse complement of a DNA string s is the string sc formed by reversing the symbols of s, then taking the
complement of each symbol (e.g., the reverse complement of "GTCA" is "TGAC").

Given: A DNA string s of length at most 1000 bp.
Return: The reverse complement sc of s.
"""

def reverse_complement(str):
    complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    reversed_s = s[::-1]
    reversed_complement_s = ''.join(complement[nucleotide] for nucleotide in reversed_s)
    
    return print(reversed_complement_s)


with open('../data/rosalind_revc.txt', 'r') as f:
    s = f.read()
# print(s)
reverse_complement(s)