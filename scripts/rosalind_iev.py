"""
Calculating Expected Offspring
url: https://rosalind.info/problems/iev/

For a random variable X taking integer values between 1 and n, the expected value of X is E(X)=∑nk=1k×Pr(X=k). The expected value offers us a way of taking the long-term average of a random variable over a large number of trials.

As a motivating example, let X be the number on a six-sided die. Over a large number of rolls, we should expect to obtain an average of 3.5 on the die (even though it's not possible to roll a 3.5). The formula for expected value confirms that E(X)=∑6k=1k×Pr(X=k)=3.5.

More generally, a random variable for which every one of a number of equally spaced outcomes has the same probability is called a uniform random variable (in the die example, this "equal spacing" is equal to 1). We can generalize our die example to find that if X is a uniform random variable with minimum possible value a and maximum possible value b, then E(X)=a+b2. You may also wish to verify that for the dice example, if Y is the random variable associated with the outcome of a second die roll, then E(X+Y)=7.

Given: Six nonnegative integers, each of which does not exceed 20,000. The integers correspond to the number of couples in a population possessing each genotype pairing for a given factor. In order, the six given integers represent the number of couples having the following genotypes:

AA-AA
AA-Aa
AA-aa
Aa-Aa
Aa-aa
aa-aa
Return: The expected number of offspring displaying the dominant phenotype in the next generation, under the assumption that every couple has exactly two offspring.
"""
def expected_dominant_offspring(couples: list) -> float:
    """
    Calculate the expected number of offspring displaying the dominant phenotype.

    Args:
        couples (list): A list of six integers representing the number of couples 
                        for each genotype pairing.

    Returns:
        float: The expected number of dominant offspring.
    """
    probabilities = [1, 1, 1, 0.75, 0.5, 0]  # Probabilities of dominant phenotype for each genotype pairing
    
    # expected number of dominant offspring
    expected_value = 2 * sum(couples[i] * probabilities[i] for i in range(6))
    
    return expected_value

if __name__ == '__main__':
    with open('../data/rosalind_iev.txt') as f:
        couples = list(map(int, f.read().strip().split()))
# couples = [19679, 17901, 18859, 19845, 19780, 16248]
        result = expected_dominant_offspring(couples)
        print(result)
