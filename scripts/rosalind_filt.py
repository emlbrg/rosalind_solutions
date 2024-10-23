"""
Read Filtration by Quality.
url: https://rosalind.info/problems/filt/

Problem

Poor-quality reads can be filtered out using the FASTQ Quality Filter tool from the FASTX toolkit. A command-line version of FASTX can be downloaded for Linux or MacOS from its website. An online interface for the FASTQ Quality Filter is also available here within the Galaxy web platform.

Given: A quality threshold value q, percentage of bases p, and set of FASTQ entries.
Return: Number of reads in filtered FASTQ entries
"""
from Bio import SeqIO

def passes_quality_threshold(record, threshold: int, percentage: int) -> bool:
    """
    Check if a sequence passes the quality threshold.

    Args:
        record: A SeqIO record containing sequence and quality data.
        threshold (int): Minimum acceptable PHRED quality score.
        percentage (int): Required percentage of bases meeting the threshold.

    Returns:
        bool: True if the sequence passes, False otherwise.
    """
    phred_scores = record.letter_annotations["phred_quality"]
    passing_bases = sum(phred >= threshold for phred in phred_scores)
    return (passing_bases / len(phred_scores)) * 100 >= percentage

def quality_filtration(data: str) -> int:
    """
    Filters reads based on quality threshold and percentage.

    Args:
        data (str): Path to the input file containing FASTQ entries.

    Returns:
        int: Number of reads passing the quality filtration.
    """
    with open(data, "r") as f:
        threshold, percentage = map(int, f.readline().strip().split())
        # Use generator expression to count passing reads
        return sum(
            passes_quality_threshold(record, threshold, percentage)
            for record in SeqIO.parse(f, 'fastq')
        )

if __name__ == '__main__':
    data = '../data/rosalind_filt.txt'
    count = quality_filtration(data)
    print(count)
