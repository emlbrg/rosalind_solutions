"""
Translating RNA into Protein.
url: https://rosalind.info/problems/prot/

The 20 commonly occurring amino acids are abbreviated by using 20 letters from the English alphabet (all letters except for B, J, O, U, X, and Z). Protein strings are constructed from these 20 symbols. Henceforth, the term genetic string will incorporate protein strings along with DNA strings and RNA strings.

The RNA codon table dictates the details regarding the encoding of specific codons into the amino acid alphabet.

Given: An RNA string s corresponding to a strand of mRNA (of length at most 10 kbp).
Return: The protein string encoded by s.
"""
from constants import RNA_CODON_TABLE

def translate_rna_to_protein(rna: str) -> str:
    """
    Translates an RNA string into a protein string based on the provided RNA codon table.

    Args:
        rna (str): RNA string to be translated.

    Returns:
        str: Protein string corresponding to the RNA string.
    """
    protein = []
    
    # Iterate through the RNA string in chunks of 3 (codons)
    for i in range(0, len(rna), 3):
        codon = rna[i:i+3]
        amino_acid = RNA_CODON_TABLE.get(codon, '*')
        
        if amino_acid == '*':  # Stop translation if a stop codon is encountered
            break
        protein.append(amino_acid)
    
    return ''.join(protein)

if __name__ == '__main__':
    with open('../data/rosalind_prot.txt', 'r') as f:
        
        print(translate_rna_to_protein(f.readline()))