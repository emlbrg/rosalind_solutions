"""
Edit Distance
url: http://rosalind.info/problems/edit/

Given: Two protein strings s and t in FASTA format (each of length at most 1000 aa).
Return: The edit distance dE(s,t).
"""
from Bio import SeqIO
from typing import List

def calculate_edit_distance(seq1: str, seq2: str) -> int:
    """
    Calculate the edit distance between two sequences.

    Args:
        seq1 (str): First protein sequence.
        seq2 (str): Second protein sequence.

    Returns:
        int: The edit distance between seq1 and seq2.
    """
    len1, len2 = len(seq1), len(seq2)

    if len1 == 0 or len2 == 0:  # if one sequence is empty
        return len1 + len2

    if len1 < len2:  # shorter sequence for columns -> less space uzage
        seq1, seq2 = seq2, seq1
        len1, len2 = len2, len1

    # initialize the two rows we need to keep track of
    previous_row = list(range(len2 + 1))
    current_row = [0] * (len2 + 1)

    ### Compute edit distance ###
    for i in range(1, len1 + 1):
        current_row[0] = i
        for j in range(1, len2 + 1):
            insert_cost = current_row[j - 1] + 1
            delete_cost = previous_row[j] + 1
            substitute_cost = previous_row[j - 1] + (seq1[i - 1] != seq2[j - 1])
            current_row[j] = min(insert_cost, delete_cost, substitute_cost)
        previous_row, current_row = current_row, previous_row

    return previous_row[len2]

def load_sequences(file_path: str) -> List:
    """
    Load sequences from a FASTA file. Used this becasue my usual parser expects a fasta format 

    Args:
        file_path (str): Path to the FASTA file.

    Returns:
        list: List of sequences extracted from the file.
    """
    sequences = []
    with open(file_path, "r") as fa:
        for seq_record in SeqIO.parse(fa, "fasta"):
            sequences.append(str(seq_record.seq))
    return sequences

if __name__ == "__main__":
    seqs = load_sequences("../data/rosalind_edit.txt")
    seq1, seq2 = seqs
    print(calculate_edit_distance(seq1, seq2))
