"""
Finding Genes with ORFs
url: https://rosalind.info/problems/orfr/

Problem

An ORF begins with a start codon and ends either at a stop codon or at the end of the string. We will assume the standard genetic code for translating an RNA string into a protein string (i.e., see the standard RNA codon table).
ORF finder from the SMS 2 package can be run online here.

Given: A DNA string s of length at most 1 kbp.
Return: The longest protein string that can be translated from an ORF of s. If more than one protein string of maximum length exists, then you may output any solution.
"""
from Bio.Seq import Seq
from tools import parse_fasta

def find_orf_translate(sequences, table=1, min_pro_len=30):
    """
    Find long proteins from a dictionary of sequences and return detailed information.

    Args:
        sequences (dict): Dictionary of sequences where keys are IDs and values are sequences.
        table (int): Translation table to use (default is 1).
        min_pro_len (int): Minimum protein length to consider (default is 100).
    
    Returns:
        list: List of proteins with details on record ID, strand, frame, length, and sequence.
    """
    all_proteins = []

    for seq_id, seq in sequences.items():
        nuc = Seq(seq)  # Essentially, it's a string?

        for strand, nuc in [(+1, nuc), (-1, nuc.reverse_complement())]:  # Forward + reverse
            for frame in range(3):
                protein_sequence = nuc[frame:].translate(table)
                for pro in protein_sequence.split("*"):
                    start_index = pro.find("M")
                    if start_index != -1:  # There is at least one 'M'
                        trimmed_pro = pro[start_index:]
                        start = frame + start_index * 3
                        end = start + len(trimmed_pro) * 3
                        end += 3 # Stop codon 
                        # Check if the trimmed protein meets the min len
                        if len(trimmed_pro) >= min_pro_len:
                            all_proteins.append((seq_id, strand, frame, start, end, trimmed_pro))
    return all_proteins

def find_longest_orf(all_orfs):
    """
    Find the longest ORF for each sequence ID from a list of ORFs.

    Args:
        all_orfs (list): List of ORF details (ID, strand, frame, start, end, sequence).

    Returns:
        list: A list of the longest ORFs for each sequence ID.
    """
    longest_orfs = {}
    
    # Group ORFs by sequence ID and find the longest for each ID
    for orf in all_orfs:
        seq_id = orf[0]
        if seq_id not in longest_orfs or len(orf[-1]) > len(longest_orfs[seq_id][-1]):
            longest_orfs[seq_id] = orf  # Update with the longer ORF

    return list(longest_orfs.values())

if __name__ == '__main__':
    fasta_file = '../data/rosalind_orfr.txt'
    sequences = parse_fasta(fasta_file)
    # print(sequences)
    all_prots = find_orf_translate(sequences)
    # print(all_prots)
    longest_prot = find_longest_orf(all_prots)
    print(longest_prot)
