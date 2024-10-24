"""
Introduction to Alternative Splicing
url: http://rosalind.info/problems/aspc/

Problem

In “Counting Subsets”, we saw that the total number of subsets of a set S containing n elements is equal to 2n.
However, if we intend to count the total number of subsets of S having a fixed size k, then we use the combination statistic C(n,k), also written (nk).

Given: Positive integers n and m with 0≤m≤n≤2000.
Return: The sum of combinations C(n,k) for all k satisfying m≤k≤n, modulo 1,000,000. In shorthand, ∑nk=m(nk).
"""
import math

with open("../data/rosalind_aspc.txt", "r") as f:
    n, m = map(int, f.readline().strip().split(" "))
# print("Positive integers n and  m are: {}, {}".format(n,m))

# The sum of combinations C(n,k) for all k satisfying m≤k≤n, modulo 1,000,000.
c = 0
for k in range(m, n+1):
    c += math.factorial(n) // (math.factorial(k) * math.factorial(n-k))

print(c%1000000)
# print("The sum of combinations: {}".format(c))
# print("The sum of combinations (modulo 1,000,000): {}".format(c%1000000))