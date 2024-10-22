"""
Computing GC Content.
url: https://rosalind.info/problems/gc/

The GC-content of a DNA string is given by the percentage of symbols in the string that are 'C' or 'G'. For example, the GC-content of "AGCTATAG" is 37.5%. Note that the reverse complement of any DNA string has the same GC-content.

DNA strings must be labeled when they are consolidated into a database. A commonly used method of string labeling is called FASTA format. In this format, the string is introduced by a line that begins with '>', followed by some labeling information. Subsequent lines contain the string itself; the first line to begin with '>' indicates the label of the next string.

In Rosalind's implementation, a string in FASTA format will be labeled by the ID "Rosalind_xxxx", where "xxxx" denotes a four-digit code between 0000 and 9999.

Given: At most 10 DNA strings in FASTA format (of length at most 1 kbp each).
Return: The ID of the string having the highest GC-content, followed by the GC-content of that string. Rosalind allows for a default error of 0.001 in all decimal answers unless otherwise stated; please see the note on absolute error below.
"""
from typing import Optional, Any
from tools import parse_fasta

# read_fasta('/home/molcure/storageA/emilia/perso/rosalind/fasta.fa')

def calculate_gc_content(sequence: str) -> float:
    """
    Calculate the GC content of a given DNA sequence.

    Args:
        sequence (str): The DNA sequence to calculate the GC content for.

    Returns:
        float: The GC content as a percentage of 'G' and 'C' in the sequence.
    """
    g_count = sequence.count('G')
    c_count = sequence.count('C')
    total_length = len(sequence)
    
    if total_length == 0:
        return 0.0  #  empty sequenceの場合は
    
    gc_content = (g_count + c_count) / total_length * 100
    return gc_content

def find_highest_gc(fasta_file_path: str) -> tuple[Optional[Any], float]:
    """
    Find the sequence with the highest GC content from a FASTA file.

    Args:
        fasta_file_path (str): The file path to the FASTA file.

    Returns:
        tuple[str, float]: A tuple containing the sequence ID with the highest GC content and the GC content percentage.
    """
    fasta_dict = parse_fasta(fasta_file_path)
    
    highest_gc_id = None
    highest_gc_content = 0.0
    
    for seq_id, sequence in fasta_dict.items():
        gc_content = calculate_gc_content(sequence)
        if gc_content > highest_gc_content:
            highest_gc_content = gc_content
            highest_gc_id = seq_id
    
    return highest_gc_id, highest_gc_content

if __name__ == '__main__':
    fasta_file_path = '../data/rosalind_gc.txt'
    id_with_highest_gc, gc_content = find_highest_gc(fasta_file_path)
    print(f"{id_with_highest_gc}\n{gc_content:.6f}")