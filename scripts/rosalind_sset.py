"""
Counting Subsets
url: http://rosalind.info/problems/sset/

Problem

A set is the mathematical term for a loose collection of objects, called elements. Examples of sets include {the moon, the sun, Wilford Brimley} and ℝ, the set containing all real numbers. We even have the empty set, represented by ∅ or {}, which contains no elements at all. Two sets are equal when they contain the same elements. In other words, in contrast to permutations, the ordering of the elements of a set is unimportant (e.g., {the moon, the sun, Wilford Brimley} is equivalent to {Wilford Brimley, the moon, the sun}). Sets are not allowed to contain duplicate elements, so that {Wilford Brimley, the sun, the sun} is not a set. We have already used sets of 2 elements to represent edges from a graph.
A set A is a subset of B if every element of A is also an element of B, and we write A⊆B. For example, {the sun, the moon}⊆{the sun, the moon, Wilford Brimley}, and ∅ is a subset of every set (including itself!).
As illustrated in the biological introduction, we can use subsets to represent the collection of taxa possessing a character. However, the number of applications is endless; for example, an event in probability can now be defined as a subset of the set containing all possible outcomes.
Our first question is to count the total number of possible subsets of a given set.

Given: A positive integer n (n≤1000).
Return: The total number of subsets of {1,2,…,n} modulo 1,000,000.
"""
def total_subsets(n: int) -> int:
    """
    Calculate the total number of subsets of the set {1, 2, ..., n}
    and return the result modulo 1,000,000.

    Args:
        n (int): The number of elements in the set.
    
    Returns:
        int: The total number of subsets modulo 1,000,000.
    """
    return pow(2, n, 1000000)  # Efficient power mod calculation

if __name__ == '__main__':
    with open('../data/rosalind_sset.txt', 'r') as f:
        n = int(f.read().strip())
    
    print(total_subsets(n))