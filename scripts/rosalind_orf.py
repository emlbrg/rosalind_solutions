"""
Open Reading Frames.
url: https://rosalind.info/problems/orf/

Problem

Either strand of a DNA double helix can serve as the coding strand for RNA transcription. Hence, a given DNA string implies six total reading frames, or ways in which the same region of DNA can be translated into amino acids: three reading frames result from reading the string itself, whereas three more result from reading its reverse complement.

An open reading frame (ORF) is one which starts from the start codon and ends by stop codon, without any other stop codons in between. Thus, a candidate protein string is derived by translating an open reading frame into amino acids until a stop codon is reached.

Given: A DNA string s of length at most 1 kbp in FASTA format.
Return: Every distinct candidate protein string that can be translated from ORFs of s. Strings can be returned in any order.
"""
from Bio.Seq import Seq
from typing import List
from tools import parse_fasta

def find_orf_translate(sequences: dict, table: int=1, min_pro_len: int=1) -> List:
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

if __name__ == '__main__':
    fasta_file =  '../data/rosalind_orf.txt'
    dna_sequence = parse_fasta(fasta_file)
    # print(dna_sequence)

    proteins = find_orf_translate(dna_sequence)
    for protein in proteins:
        print(protein)
