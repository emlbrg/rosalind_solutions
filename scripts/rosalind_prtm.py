"""
Calculating Protein Mass.
url: https://rosalind.info/problems/prtm/

Problem

In a weighted alphabet, every symbol is assigned a positive real number called a weight. A string formed from a weighted alphabet is called a weighted string, and its weight is equal to the sum of the weights of its symbols.
The standard weight assigned to each member of the 20-symbol amino acid alphabet is the monoisotopic mass of the corresponding amino acid.

Given: A protein string P of length at most 1000 aa.
Return: The total weight of P. Consult the monoisotopic mass table.
"""
from constants import monoisotopic_mass_table

def calculate_protein_mass(protein_string: str) -> float:
    """Calculate the total weight of the protein string based on the monoisotopic mass table.
   
    Args:
        protein_string (str): the protein sequence
    
    Returns:
        total_mass (float): a float representing the calculated mass
    """
    total_mass = 0.0
    for aa in protein_string:
        total_mass += monoisotopic_mass_table[aa]
    return total_mass

if __name__ == '__main__':
    with open('rosalind_prtm.txt', 'r') as f:
        prt_list = f.readlines()
        # prt_string = ''.join(prt_list)
        prt_string = ''.join(prt_list).replace('\n', '').replace(' ', '')
        print(prt_string)
        # print(type(prt_string))
        total_weight = calculate_protein_mass(prt_string)
        print(f'{total_weight:.3f}')
