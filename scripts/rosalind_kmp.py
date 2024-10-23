"""
Speeding Up Motif Finding.
url: https://rosalind.info/problems/kmp/

A prefix of a length n string s is a substring s[1:j]; a suffix of s is a substring s[k:n].

The failure array of s is an array P of length n for which P[k] is the length of the longest substring s[j:k] that is equal to some prefix s[1:k−j+1], where j cannot equal 1 (otherwise, P[k] would always equal k). By convention, P[1]=0.

Given: A DNA string s (of length at most 100 kbp) in FASTA format.
Return: The failure array of s.
"""
# need to compute a failure array?!
from typing import List
from tools import parse_fasta_str

def failure_array(s: str) -> List[int]:
    """
    Compute the failure array (or prefix function) of a string.

    Args:
        s (str): Input DNA string.

    Returns:
        List[int]: The failure array of the string.
    """
    n = len(s)
    P = [0] * n
    j = 0  # Length of the previous longest prefix suffix

    for i in range(1, n):
        # Find the longest prefix which is also a suffix for s[0:i]
        while j > 0 and s[i] != s[j]:
            j = P[j - 1]

        if s[i] == s[j]:
            j += 1

        P[i] = j

    return P

if __name__ == '__main__':
    fasta_file = '../data/rosalind_kmp.txt'  # Adjust the path as needed
    dna_sequence = parse_fasta_str(fasta_file)
    result = failure_array(dna_sequence)
    
    # Print the failure array as space-separated values
    print(" ".join(map(str, result)))
