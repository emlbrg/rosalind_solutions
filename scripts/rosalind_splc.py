"""
RNA Splicing.
url: https://rosalind.info/problems/splc/

Problem

After identifying the exons and introns of an RNA string, we only need to delete the introns and concatenate the exons to form a new string ready for translation.

Given: A DNA string s (of length at most 1 kbp) and a collection of substrings of s acting as introns. All strings are given in FASTA format.
Return: A protein string resulting from transcribing and translating the exons of s. (Note: Only one solution will exist for the dataset provided.)
"""
from Bio.Seq import Seq
from Bio import SeqIO

def parse_fasta(file_path: str) -> tuple:
    """
    Parse the input FASTA file to extract the main DNA string and the introns.

    Args:
        file_path (str): Path to the FASTA file.

    Returns:
        s, introns (tuple): A tuple containing the main DNA string (s) and a list of introns.
    """
    records = list(SeqIO.parse(file_path, "fasta"))
    s = str(records[0].seq)  # Main DNA sequence
    introns = [str(record.seq) for record in records[1:]]  # List of introns
    return s, introns

def remove_introns(dna_seq: str, introns: list) -> str:
    """
    Remove the introns from the main DNA sequence.

    Args:
        dna_seq (str): The main DNA sequence.
        introns (list): List of introns to be removed from the DNA sequence.

    Returns:
        dna_seq (str): The exon sequence after removing introns.
    """
    for intron in introns:
        dna_seq = dna_seq.replace(intron, '')
    return dna_seq

def translate_dna_to_protein(dna_seq: str) -> str:
    """
    Transcribe the exon DNA sequence into mRNA and translate it into a protein.

    Args:
        dna_seq (str): The exon DNA sequence.

    Returns:
        protein (str): The translated protein string.
    """
    exon_seq = Seq(dna_seq)
    # Transcribe to mRNA + translate to protein
    protein = exon_seq.translate(to_stop=True)  # Translate
    return str(protein)

def main(file_path):
    """
    Main function to process the input FASTA file and output the protein string.

    Args:
        file_path (str): Path to the FASTA file.
    """
    s, introns = parse_fasta(file_path)
    exon_seq = remove_introns(s, introns)
    protein = translate_dna_to_protein(exon_seq)
    print(protein)

if __name__ == '__main__':
    file_path = '../data/rosalind_splc.txt'
    main(file_path)
