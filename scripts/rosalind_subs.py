"""
Finding a Motif in DNA.
url: https://rosalind.info/problems/subs/

Given two strings s and t, t is a substring of s if t is contained as a contiguous collection of symbols in s (as a result, t must be no longer than s).

The position of a symbol in a string is the total number of symbols found to its left, including itself (e.g., the positions of all occurrences of 'U' in "AUGCUUCAGAAAGGUCUUACG" are 2, 5, 6, 15, 17, and 18). The symbol at position i of s is denoted by s[i].

A substring of s can be represented as s[j:k], where j and k represent the starting and ending positions of the substring in s; for example, if s = "AUGCUUCAGAAAGGUCUUACG", then s[2:5] = "UGCU".

The location of a substring s[j:k] is its beginning position j; note that t will have multiple locations in s if it occurs more than once as a substring of s (see the Sample below).

Given: Two DNA strings s and t (each of length at most 1 kbp).
Return: All locations of t as a substring of s.
"""
def find_substring_locations(s: str, t: str) -> list:
    """
    Find all starting positions of substring t in string s.

    Args:
        s (str): The main string to search in.
        t (str): The substring to find.

    Returns:
        list: A list of starting positions (1-indexed) where t occurs in s.
    """
    positions = []
    t_length = len(t)
    
    for i in range(len(s) - t_length + 1):  # Iterate through possible starting points
        if s[i:i + t_length] == t:  # Check for substring match
            positions.append(i + 1)  # Convert to 1-indexed position
            
    return positions

if __name__ == '__main__':
    with open('../data/rosalind_subs.txt') as f:
        lines = f.readlines()
        s = lines[0].strip()
        t = lines[1].strip()

    positions = find_substring_locations(s, t)
    print(" ".join(map(str, positions)))
