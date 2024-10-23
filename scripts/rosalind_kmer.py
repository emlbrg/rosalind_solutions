"""
k-Mer Composition.
url: https://rosalind.info/problems/kmer/

For a fixed positive integer k, order all possible k-mers taken from an underlying alphabet lexicographically.
Then the k-mer composition of a string s can be represented by an array A for which A[m] denotes the number of times that the mth k-mer (with respect to the lexicographic order) appears in s.

Given: A DNA string s in FASTA format (having length at most 100 kbp).
Return: The 4-mer composition of s.
"""
from itertools import product
from collections import defaultdict
from tools import parse_fasta_str

def generate_kmers(k: int) -> list:
    """
    Generate all possible k-mers in lexicographic order from the DNA alphabet {A, C, G, T}.

    Args:
        k (int): The length of the k-mer.

    Returns:
        list: A list of all possible k-mers sorted lexicographically.
    """
    return [''.join(p) for p in product('ACGT', repeat=k)]

def count_kmers(sequence: str, k: int) -> list:
    """
    Count occurrences of each k-mer in the sequence.

    Args:
        sequence (str): The DNA sequence.
        k (int): The length of k-mers to count.

    Returns:
        list: A list of k-mer counts in lexicographic order.
    """
    kmer_dict = defaultdict(int)
    
    # sliding window of size k over the sequence
    for i in range(len(sequence) - k + 1):
        kmer = sequence[i:i+k]
        kmer_dict[kmer] += 1
    
    # Generate all k-mers in lexicographic order
    kmers = generate_kmers(k)
    
    return [kmer_dict[kmer] for kmer in kmers]  # counts for each k-mer in lexicographic order

if __name__ == '__main__':
    fasta_file_path = '../data/rosalind_kmer.txt'
    k = 4 
    sequence = parse_fasta_str(fasta_file_path)
    kmer_counts = count_kmers(sequence, k)
    print(' '.join(map(str, kmer_counts)))
