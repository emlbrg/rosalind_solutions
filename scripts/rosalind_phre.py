"""
Read Quality Distribution.
url: https://rosalind.info/problems/phre/

A version of FastQC can be downloaded here and run locally on any operating system with a suitable Java Runtime Environment (JRE) installed.
An online version of FastQC is also available here in the "Andromeda" Galaxy instance.

Given: A quality threshold, along with FASTQ entries for multiple reads.
Return: The number of reads whose average quality is below the threshold.
"""
from Bio import SeqIO

def count_low_quality_reads(fastq_file):
    """
    Count the number of reads whose average quality score is below the given threshold.
    
    Args:
        fastq_file (str): Path to the input FASTQ file where the first line is the threshold.
        
    Returns:
        int: The number of reads with average quality below the threshold.
    """
    low_quality_count = 0
    
    with open(fastq_file, 'r') as handle:
        threshold = int(handle.readline().strip())  # first line == thread threshold
        for record in SeqIO.parse(handle, "fastq"):  # the rest is a nomral fastq WTF why did they choose violence
            avg_quality = sum(record.letter_annotations["phred_quality"]) / len(record)
            # This is quite self explanatory
            if avg_quality < threshold:
                low_quality_count += 1
    
    return low_quality_count

if __name__ == "__main__":
    fastq_file = "../data/rosalind_phre.txt"
    
    result = count_low_quality_reads(fastq_file)
    print(result)
