"""
Protein Translation.
url: https://rosalind.info/problems/ptra/

Problem

The 20 commonly occurring amino acids are abbreviated by using 20 letters from the English alphabet (all letters except for B, J, O, U, X, and Z). Protein strings are constructed from these 20 symbols. The RNA codon table shows the encoding from each RNA codon to the amino acid alphabet.

The Translate tool from the SMS 2 package can be found here in the SMS 2 package

A detailed list of genetic code variants (codon tables) along with indexes representing these codes (1 = standard genetic code, etc.) can be obtained here.

For now, when translating DNA and RNA strings, we will start with the first letter of the string and ignore stop codons.

Given: A DNA string s of length at most 10 kbp, and a protein string translated by s.
Return: The index of the genetic code variant that was used for translation. (If multiple solutions exist, you may return any one.)
"""
from Bio.Seq import translate # type: ignore

def genetic_code(dna: str, protein: str) -> list[int]:
    """
    Identifies the genetic code variant used for translating a given DNA string into a protein string.

    Args:
        dna (str): A DNA string of length at most 10 kbp.
        protein (str): A protein string that is translated from the DNA string.

    Returns:
        list[int]: A list of genetic code variant indices that could have been used for the translation.
                   If multiple solutions exist, any one of them can be returned.
    
    This function compares the translated DNA string with the given protein string using several genetic code
    tables (identified by their indices). If the translation matches the protein string, the index of the 
    genetic code variant is stored in the result list.
    """
    res = []
    table_ids = [1, 2, 3, 4, 5, 6, 9, 10, 11, 12, 13, 14, 15]
    for t in table_ids:
        if translate(dna, table=t, to_stop=True) == protein:
            res.append(t)
    return res

if __name__ == '__main__':
    with open('../data/rosalind_ptra.txt', 'r') as f:
        dna = f.readline().strip()
        protein = f.readline().strip()
    res = genetic_code(dna, protein)
    print(res)