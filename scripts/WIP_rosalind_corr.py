"""
Error Correction in Reads.
url: https://rosalind.info/problems/corr/

As is the case with point mutations, the most common type of sequencing error occurs when a single nucleotide from a read is interpreted incorrectly.

Given: A collection of up to 1000 reads of equal length (at most 50 bp) in FASTA format. Some of these reads were generated with a single-nucleotide error. For each read s in the dataset, one of the following applies:

s was correctly sequenced and appears in the dataset at least twice (possibly as a reverse complement);
s is incorrect, it appears in the dataset exactly once, and its Hamming distance is 1 with respect to exactly one correct read in the dataset (or its reverse complement).
Return: A list of all corrections in the form "[old read]->[new read]". (Each correction must be a single symbol substitution, and you may return the corrections in any order.)
"""
from tools import parse_fasta
from rosalind_hamm import hamming_distance

from typing import Dict, List

def reverse_complement(s: str) -> str:
    """
    Compute the reverse complement of a given DNA string.
    
    Args:
        s (str): A DNA string composed of 'A', 'T', 'C', and 'G', with a length at most 1000 bp.

    Returns:
        str: The reverse complement of the input DNA string.
    """
    complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    reversed_s = s[::-1]
    return ''.join(complement[nucleotide] for nucleotide in reversed_s)

def find_corrections(fasta_file: str) -> List[str]:
    """
    Find corrections for sequencing errors in a FASTA file.

    Args:
        fasta_file (str): Path to the input FASTA file.

    Returns:
        List[str]: A list of corrections in the format "old_read->new_read".
    """
    sequences = parse_fasta(fasta_file)
    count = {}
    
    # Count occurrences of each sequence and its reverse complement
    for seq in sequences.values():
        count[seq] = count.get(seq, 0) + 1
        rev_seq = reverse_complement(seq)
        count[rev_seq] = count.get(rev_seq, 0) + 1

    corrections = set()  # Use a set to avoid duplicates
    # Check for sequences that appear exactly once
    for seq, num_occurrences in count.items():
        if num_occurrences == 1:
            # Check for Hamming distance with all other sequences
            for other_seq in sequences.values():
                if hamming_distance(seq, other_seq) == 1:
                    corrections.add(f"{seq}->{other_seq}")
                elif hamming_distance(reverse_complement(seq), other_seq) == 1:
                    corrections.add(f"{reverse_complement(seq)}->{other_seq}")

    return list(corrections)

# Example usage:
fasta_file = '../data/rosalind_corr.txt'
corrections = find_corrections(fasta_file)
print('\n'.join(corrections))

