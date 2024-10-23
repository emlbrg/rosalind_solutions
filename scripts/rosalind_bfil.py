"""
Base Filtration by Quality
url: http://rosalind.info/problems/bfil/

Bad quality bases can be easily trimmed out using certain threshold (defined by quality plot similar to what we did in “Base Quality Distribution”) There is a lot of trimming tools, you can try one of following:
- FASTQ Quality Trimmer tool on the Galaxy. It uses a "sliding window" approach so for a simple trimming of the ends you should set window size 1.
- Trimmomatic. It is a command-line java-based tool, detail description and download link can be found here. For a simple trimming from both ends you should specify parameters LEADING and TRAILING.

Given: FASTQ file, quality cut-off value q, Phred33 quality score assumed.
Return: FASTQ file trimmed from the both ends (removed leading and trailing bases with quality lower than q)
"""
from Bio import SeqIO
from typing import List, Tuple

def trim_bases(phred: List[int], threshold: int) -> Tuple[int, int]:
    """
    Find the start and end indices for trimming based on the quality threshold.

    Args:
        phred (List[int]): List of Phred quality scores for the sequence.
        threshold (int): The quality threshold below which bases will be trimmed.

    Returns:
        Tuple[int, int]: A tuple containing the start and end indices for trimming.
    """
    start = next((i for i, q in enumerate(phred) if q >= threshold), len(phred))
    end = next((i for i in range(len(phred) - 1, -1, -1) if phred[i] >= threshold), -1)
    return start, end + 1  # end is exclusive

def base_filtration_quality(file_path: str) -> None:
    """
    Trim FASTQ entries based on quality scores, removing low-quality bases from both ends.

    Args:
        file_path (str): Path to the input FASTQ file containing sequences.

    Returns:
        None: The function prints trimmed FASTQ entries to the standard output.
    """
    with open(file_path, "r") as f:
        threshold = int(f.readline().strip())
        for record in SeqIO.parse(f, "fastq"):
            phred = record.letter_annotations["phred_quality"]
            start, end = trim_bases(phred, threshold)
            if start < end:  # Ensure there's something to trim
                trimmed_record = record[start:end]
                print(trimmed_record.format('fastq'))

if __name__ == "__main__":
    base_filtration_quality("../data/rosalind_bfil.txt")
