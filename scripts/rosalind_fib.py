"""
Rabbits and Recurrence Relations
url: https://rosalind.info/problems/fib/

A sequence is an ordered collection of objects (usually numbers), which are allowed to repeat. Sequences can be finite or infinite. Two examples are the finite sequence (π,−2‾√,0,π) and the infinite sequence of odd numbers (1,3,5,7,9,…). We use the notation an to represent the n-th term of a sequence.

A recurrence relation is a way of defining the terms of a sequence with respect to the values of previous terms. In the case of Fibonacci's rabbits from the introduction, any given month will contain the rabbits that were alive the previous month, plus any new offspring. A key observation is that the number of offspring in any month is equal to the number of rabbits that were alive two months prior. As a result, if Fn represents the number of rabbit pairs alive after the n-th month, then we obtain the Fibonacci sequence having terms Fn that are defined by the recurrence relation Fn=Fn−1+Fn−2 (with F1=F2=1 to initiate the sequence). Although the sequence bears Fibonacci's name, it was known to Indian mathematicians over two millennia ago.

When finding the n-th term of a sequence defined by a recurrence relation, we can simply use the recurrence relation to generate terms for progressively larger values of n. This problem introduces us to the computational technique of dynamic programming, which successively builds up solutions by using the answers to smaller cases.

Given: Positive integers n≤40 and k≤5.
Return: The total number of rabbit pairs that will be present after n months, if we begin with 1 pair and in each generation, every pair of reproduction-age rabbits produces a litter of k rabbit pairs (instead of only 1 pair).
"""
from typing import List

def rabbit_pairs(n: int, k: int) -> List[int]:
    """
    Calculate the number of rabbit pairs after n months given that each pair produces k new pairs each month.

    Args:
        n (int): Number of months.
        k (int): Number of new pairs produced by each pair each month.

    Returns:
        List[int]: Total number of rabbit pairs after n months.
    """
    if n == 1:
        return 1
    elif n == 2:
        return 1

    dp = [0] * (n + 1)
    dp[1] = 1
    dp[2] = 1
    
    for i in range(3, n + 1):  # fill the DP table
        dp[i] = dp[i - 1] + k * dp[i - 2]
    
    return dp[n]

def read_input(file_path: str) -> tuple:
    """
    Reads the input from a text file and returns n and k.

    Args:
        file_path (str): Path to the input text file.

    Returns:
        tuple: A tuple containing two integers (n, k) from the file.
    """
    with open(file_path, 'r') as file:
        # Read the first line of the file and split by space
        data = file.readline().strip().split()
        n = int(data[0])
        k = int(data[1])
        return n, k

if __name__ == '__main__':
    input_file = '../data/rosalind_fib.txt'
    n, k = read_input(input_file)
    result = rabbit_pairs(n, k)

    print(result)