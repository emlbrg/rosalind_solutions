"""
Transcribing DNA into RNA
url: https://rosalind.info/problems/rna/

An RNA string is a string formed from the alphabet containing 'A', 'C', 'G', and 'U'.
Given a DNA string t corresponding to a coding strand, its transcribed RNA string u is formed by replacing all
occurrences of 'T' in t with 'U' in u.

Given: A DNA string t having length at most 1000 nt.
Return: The transcribed RNA string of t.
"""

def translate(sequence: str) -> None:
    """
    Convert a DNA sequence to its corresponding RNA sequence by replacing
    all occurrences of 'T' with 'U'.

    Args:
        sequence (str): A DNA sequence composed of 'A', 'T', 'C', and 'G'.

    Returns:
        str: The RNA sequence obtained from the input DNA sequence.
    """
    u = t.replace('T', 'U')
    print(u)

if __name__ == '__main__':
    with open('../data/rosalind_rna.txt', 'r') as f:
        t = f.read()

    translate(t)