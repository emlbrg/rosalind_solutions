"""
Finding a Shared Spliced Motif
url: https://rosalind.info/problems/lcsq/

A string u is a common subsequence of strings s and t if the symbols of u appear in order as a subsequence of both s and t. For example, "ACTG" is a common subsequence of "AACCTTGG" and "ACACTGTGA".

Analogously to the definition of longest common substring, u is a longest common subsequence of s and t if there does not exist a longer common subsequence of the two strings. Continuing our above example, "ACCTTG" is a longest common subsequence of "AACCTTGG" and "ACACTGTGA", as is "AACTGG".

Given: Two DNA strings s and t (each having length at most 1 kbp) in FASTA format.
Return: A longest common subsequence of s and t. (If more than one solution exists, you may return any one.)
"""
from tools import parse_fasta_list

def longest_common_subsequence(s: str, t: str) -> str:
    """
    Compute the longest common subsequence of two strings.

    Args:
        s (str): First DNA string.
        t (str): Second DNA string.

    Returns:
        str: A longest common subsequence of s and t.
    """
    # Initialize the DP table
    m, n = len(s), len(t)
    dp = [[0] * (n + 1) for _ in range(m + 1)]

    # Fill the DP table
    for i in range(1, m + 1):
        for j in range(1, n + 1):
            if s[i - 1] == t[j - 1]:
                dp[i][j] = dp[i - 1][j - 1] + 1
            else:
                dp[i][j] = max(dp[i - 1][j], dp[i][j - 1])

    # Backtrack to find the actual LCS
    i, j = m, n
    lcs = []

    while i > 0 and j > 0:
        if s[i - 1] == t[j - 1]:
            lcs.append(s[i - 1])
            i -= 1
            j -= 1
        elif dp[i - 1][j] >= dp[i][j - 1]:
            i -= 1
        else:
            j -= 1

    # The LCS will be in reverse order, so reverse it before returning
    return ''.join(reversed(lcs))

if __name__ == '__main__':
    # Input file path
    fasta_file = '../data/rosalind_lcsq.txt'  # Adjust the path to your input file

    # Parse the input sequences from the FASTA file
    sequences = parse_fasta_list(fasta_file)

    # We expect two sequences from the file
    s = sequences[0]
    t = sequences[1]

    # Find and print the longest common subsequence
    result = longest_common_subsequence(s, t)
    print(result)
