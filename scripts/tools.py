from typing import Dict

def parse_fasta(fasta_file_path: str) -> Dict:
    """
    Parse the sequence from a FASTA formatted file.

    Args:
        fasta_file (str): path to the inport fasta file

    Returns:
        dict: a dictionary where seq_id the the key and seq the value
    """
    sequences = {}
    current_sequence = None
    with open(fasta_file_path, 'r') as file:
        for line in file:
            line = line.strip()
            if line.startswith('>'):
                current_sequence = line[1:]
                sequences[current_sequence] = ''
            else:
                sequences[current_sequence] += line
    return sequences

def parse_fasta_str(fasta_file_path: str) -> str:
    """
    Parse the sequence from a FASTA formatted file.
    
    Args:
        fasta_file (str): path to the inport fasta file

    Returns:
        str: The RNA sequence.
    """
    with open(fasta_file_path, 'r') as file:
        lines = file.readlines()
        # Remove the first line (header) and concatenate the rest
        rna_sequence = ''.join(line.strip() for line in lines if not line.startswith('>'))
    return rna_sequence