"""
Counting DNA Nucleotides
url: http://rosalind.info/problems/dna/

A string is simply an ordered collection of symbols selected from some alphabet and formed into a word; the length
of a string is the number of symbols that it contains.

An example of a length 21 DNA string (whose alphabet contains the symbols 'A', 'C', 'G', and 'T') is "ATGCTTCAGAAAGGTCTTACG."

Given: A DNA string s of length at most 1000 nt.
Return: Four integers (separated by spaces) counting the respective number of times that the symbols 'A', 'C', 'G', and 'T' occur in s.
"""

def count_dna_nucleotides(sequence: str) -> None:
    count_A = 0
    count_C = 0
    count_G = 0
    count_T = 0
    
    for nucleotide in s:
        if nucleotide == 'A':
            count_A += 1
        elif nucleotide == 'C':
            count_C += 1
        elif nucleotide == 'G':
            count_G += 1
        elif nucleotide == 'T':
            count_T += 1
    print(count_A, count_C, count_G, count_T)

if __name__ == '__main__':
    with open('../data/rosalind_dna.txt', 'r') as f:
        s = f.read()

    count_dna_nucleotides(s)