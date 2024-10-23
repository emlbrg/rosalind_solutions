"""
Creating a Distance Matrix 
url: https://rosalind.info/problems/pdst/

Problem

For two strings s1 and s2 of equal length, the p-distance between them, denoted dp(s1,s2), is the proportion of corresponding symbols that differ between s1 and s2.
For a general distance function d on n taxa s1,s2,…,sn (taxa are often represented by genetic strings), we may encode the distances between pairs of taxa via a distance matrix D in which Di,j=d(si,sj).

Given: A collection of n (n≤10) DNA strings s1,…,sn of equal length (at most 1 kbp). Strings are given in FASTA format.
Return: The matrix D corresponding to the p-distance dp on the given strings. As always, note that your answer is allowed an absolute error of 0.001.
"""
from Bio import SeqIO
from typing import List

def p_distance(s1: str, s2: str) -> float:
    """
    Calculate the p-distance between two DNA sequences.

    The p-distance is the proportion of differing symbols at corresponding positions 
    in the two sequences.

    Args:
        s1 (str): The first DNA sequence.
        s2 (str): The second DNA sequence.

    Returns:
        float: The p-distance between the two sequences, rounded to five decimal places.
    """
    differences = sum(1 for a, b in zip(s1, s2) if a != b)
    return round(differences / len(s1), 5)

def create_distance_matrix(sequences: List[str]) -> List[List[float]]:
    """
    Create a distance matrix for a collection of DNA sequences.

    Args:
        sequences (List[str]): A list of DNA sequences.

    Returns:
        List[List[float]]: A matrix where each element [i][j] is the p-distance between 
        sequences[i] and sequences[j].
    """
    n = len(sequences)
    matrix = [[p_distance(sequences[i], sequences[j]) for j in range(n)] for i in range(n)]
    return matrix

def read_fasta(file_path: str) -> List[str]:
    """
    Read sequences from a FASTA file.

    Args:
        file_path (str): Path to the FASTA file.

    Returns:
        List[str]: A list of DNA sequences extracted from the FASTA file.
    """
    with open(file_path, "r") as f:
        return [str(record.seq) for record in SeqIO.parse(f, "fasta")]

if __name__ == "__main__":
    file_path = "../data/rosalind_pdst.txt"
    sequences = read_fasta(file_path)

    distance_matrix = create_distance_matrix(sequences)
    for row in distance_matrix:
        print(" ".join(f"{distance:.5f}" for distance in row))
