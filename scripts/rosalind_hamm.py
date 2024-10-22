"""
Counting Point Mutations
url: https://rosalind.info/problems/hamm/

Problem
Given two strings s and t of equal length, the Hamming distance between s and t, denoted dH(s,t), is the number of corresponding symbols that differ in s and t. See Figure 2.

Given: Two DNA strings s and t of equal length (not exceeding 1 kbp).
Return: The Hamming distance dH(s,t).
"""
def hamming_distance(s: str, t: str) -> int:
    """
    Calculate the Hamming distance between two DNA strings of equal length.

    Args:
        s (str): The first DNA string.
        t (str): The second DNA string.

    Returns:
        int: The Hamming distance (number of differing characters).
    """
    if len(s) != len(t):  # 両方は同じ長さ
        raise ValueError("Strings must be of equal length")
    
    distance = sum(1 for a, b in zip(s, t) if a != b)
    return distance

if __name__ == '__main__':
    with open('../data/rosalind_hamm.txt', 'r') as f:
        lines = f.readlines()
        s = lines[0].strip()  # Read the first line and strip whitespace
        t = lines[1].strip()  # Read the second line and strip whitespace

    print(hamming_distance(s, t))