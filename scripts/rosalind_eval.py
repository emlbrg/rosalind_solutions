"""
Expected Number of Restriction Sites
url: http://rosalind.info/problems/eval/

Problem

Say that you place a number of bets on your favorite sports teams. If their chances of winning are 0.3, 0.8, and 0.6, then you should expect on average to win 0.3 + 0.8 + 0.6 = 1.7 of your bets (of course, you can never win exactly 1.7!)
More generally, if we have a collection of events A1,A2,…,An, then the expected number of events occurring is Pr(A1)+Pr(A2)+⋯+Pr(An) (consult the note following the problem for a precise explanation of this fact). In this problem, we extend the idea of finding an expected number of events to finding the expected number of times that a given string occurs as a substring of a random string.

Given: A positive integer n (n≤1,000,000), a DNA string s of even length at most 10, and an array A of length at most 20, containing numbers between 0 and 1.
Return: An array B having the same length as A in which B[i] represents the expected number of times that s will appear as a substring of a random DNA string t of length n, where t is formed with GC-content A[i] (see “Introduction to Random Strings”).
"""
def expected_restriction_sites(file_path: str) -> None:
    """
    Calculate the expected number of times that a DNA string appears
    as a substring in a random DNA string of a given length with specific GC content.

    Args:
        file_path (str): Path to the input file containing values for n, s, and the GC content array A.

    Returns:
        None: Prints the expected number of restriction sites for each GC content in A.
    """
    with open(file_path, "r") as f:
        n = int(f.readline().strip())  # len of DNA
        s = f.readline().strip()  # DNA substring
        A = list(map(float, f.readline().strip().split()))  # array of GC contents

    at_count = s.count('A') + s.count('T')  # count occurencies
    gc_count = s.count('G') + s.count('C')
    
    # expected number of restriction sites for each GC content
    expected_sites = [
        (pow((1 - gc_content) / 2, at_count) * pow(gc_content / 2, gc_count) * (n - len(s) + 1))
        for gc_content in A
    ]
    
    print(" ".join(f"{site:.6f}" for site in expected_sites))


if __name__ == '__main__':
    file_path = "../data/rosalind_eval.txt"
    expected_restriction_sites(file_path)
