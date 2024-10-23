"""
Base Quality Distribution.
url: https://rosalind.info/problems/bphr/

Quality of the bases can vary depends on position in read due to nature of the sequencing procedure. One can check this quality distribution using "Per Base Sequence Quality" module of the FastQC program.
Average accepted quality values is a 10 for the lower quartile and 25 for median. If the values falls below this limit, then the module returns a warning.
Note that for the reads >50bp long FastQC will group the bases. To show data for every base in the read use "--nogroup" option.

Given: FASTQ file, quality threshold q
Return: Number of positions where mean base quality falls below given threshold
"""
from Bio import SeqIO
from typing import List, Any

def load_quality_data(file_path: str) -> tuple[int, list[Any]]:
    """
    Load quality scores from a FASTQ file.

    Args:
        file_path (str): Path to the input FASTQ file.

    Returns:
        List[List[int]]: A list of quality scores for each record in the FASTQ file.
    """
    qualities = []
    with open(file_path, "r") as f:
        # Read the threshold value
        threshold = int(f.readline().strip())
        for record in SeqIO.parse(f, "fastq"):
            quality = record.letter_annotations["phred_quality"]
            qualities.append(quality)
    return threshold, qualities

def count_below_threshold(qualities: List[List[int]], threshold: int) -> int:
    """
    Count the number of positions with average quality below a given threshold.

    Args:
        qualities (List[List[int]]): A list of quality scores for each record.
        threshold (int): The quality threshold value.

    Returns:
        int: The count of positions below the threshold.
    """
    count = 0
    num_records = len(qualities)
    num_positions = len(qualities[0])

    for i in range(num_positions):
        average_quality = sum(q[i] for q in qualities) / num_records
        if average_quality < threshold:
            count += 1

    return count

def bphr(file_path: str) -> int:
    """
    Main function to calculate the number of positions with average quality below a threshold.

    Args:
        file_path (str): Path to the input FASTQ file.

    Returns:
        int: The count of positions below the quality threshold.
    """
    threshold, qualities = load_quality_data(file_path)
    return count_below_threshold(qualities, threshold)

if __name__ == "__main__":
    data = "../data/rosalind_bphr.txt"
    count = bphr(data)
    print(count)
