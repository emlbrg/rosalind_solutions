"""
Consensus and Profile.
url: https://rosalind.info/problems/cons/

A matrix is a rectangular table of values divided into rows and columns. An m×n matrix has m rows and n columns. Given a matrix A, we write Ai,j to indicate the value found at the intersection of row i and column j.
Say that we have a collection of DNA strings, all having the same length n. Their profile matrix is a 4×n matrix P in which P1,j represents the number of times that 'A' occurs in the jth position of one of the strings, P2,j represents the number of times that C occurs in the jth position, and so on (see below).
A consensus string c is a string of length n formed from our collection by taking the most common symbol at each position; the jth symbol of c therefore corresponds to the symbol having the maximum value in the j-th column of the profile matrix. Of course, there may be more than one most common symbol, leading to multiple possible consensus strings.

Given: A collection of at most 10 DNA strings of equal length (at most 1 kbp) in FASTA format.
Return: A consensus string and profile matrix for the collection. (If several possible consensus strings exist, then you may return any one of them.)
"""
from tools import parse_fasta_list

def calculate_profile_and_consensus(sequences: str) -> tuple:
    """
    Calculate the profile matrix and consensus string from a list of DNA sequences.

    Args:
        sequences (list): List of DNA sequences.

    Returns:
        tuple: A tuple containing the consensus string and profile matrix.
    """
    n = len(sequences[0])  # Length of the sequences
    profile = {'A': [0] * n, 'C': [0] * n, 'G': [0] * n, 'T': [0] * n}

    # Count occurrences of each nucleotide at each position
    for seq in sequences:
        for i, nucleotide in enumerate(seq):
            profile[nucleotide][i] += 1

    # Generate the consensus string
    consensus = ''
    for j in range(n):
        consensus += max(profile, key=lambda x: profile[x][j])

    return consensus, profile


def main(fasta_file):
    """
    Main function to read a FASTA file and print the consensus string and profile matrix.

    Args:
        fasta_file (str): Path to the FASTA file.
    """
    sequences = parse_fasta_list(fasta_file)
    consensus, profile = calculate_profile_and_consensus(sequences)

    # Print the consensus string
    print(consensus)
    # Print the profile matrix
    for nucleotide in 'ACGT':
        print(f"{nucleotide}: {' '.join(map(str, profile[nucleotide]))}")


# Example usage:
# Replace 'your_fasta_file.fasta' with the path to your actual FASTA file
fasta_file_path = '../data/rosalind_cons.txt'
main(fasta_file_path)

