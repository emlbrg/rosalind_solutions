from typing import Dict, List

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

def parse_fasta_list(fasta_file_path: str) -> List:
    """
    Parse a FASTA file and return a list of sequences.

    Args:
        file_path (str): Path to the FASTA file.

    Returns:
        list: A list of DNA sequences.
    """
    with open(fasta_file_path, 'r') as file:
        sequences = []
        current_seq = ''
        
        for line in file:
            line = line.strip()
            if line.startswith('>'):  # Skip the header line
                if current_seq:  # Save the previous sequence
                    sequences.append(current_seq)
                    current_seq = ''
            else:
                current_seq += line  # Append the line to the current sequence

        if current_seq:  # Append the last sequence
            sequences.append(current_seq)
    
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
        sequence = ''.join(line.strip() for line in lines if not line.startswith('>'))
    return sequence