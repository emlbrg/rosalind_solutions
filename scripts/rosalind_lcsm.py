"""
Finding a Shared Motif.
url: https://rosalind.info/problems/lcsm/

Problem

A common substring of a collection of strings is a substring of every member of the collection. We say that a common substring is a longest common substring if there does not exist a longer common substring. For example, "CG" is a common substring of "ACGTACGT" and "AACCGTATA", but it is not as long as possible; in this case, "CGTA" is a longest common substring of "ACGTACGT" and "AACCGTATA".
Note that the longest common substring is not necessarily unique; for a simple example, "AA" and "CC" are both longest common substrings of "AACC" and "CCAA".

Given: A collection of k (k≤100) DNA strings of length at most 1 kbp each in FASTA format.
Return: A longest common substring of the collection. (If multiple solutions exist, you may return any single solution.)
"""
from Bio import SeqIO
    
def shortest_seq(seq):
    min_len = 10000
    shortest_seq = ''
    for i in seq.keys():
        if len(seq[i]) < min_len:
            min_len = len(seq[i])
            shortest_seq = seq[i]
    return shortest_seq

def shared_motif(seq):
    s_seq = shortest_seq(seq)
    motif = set()
    for i in range(len(s_seq)):
        for j in range(i+1,len(s_seq)+1):
            motif.add(s_seq[i:j])
    for s in seq.values():
        update_motif = list(motif)
        for m in update_motif:
            if m not in s:
                motif.remove(m)
    n = 0
    longest_motif = ''
    for i in motif:
        if len(i) > n:
            longest_motif = i
            n = len(i)
    return longest_motif

if __name__ == '__main__':
    seq_name, seq_string = [], []
    with open ('../data/rosalind_lcsm.txt','r') as fa:
        for seq_record  in SeqIO.parse(fa,'fasta'):
            seq_name.append(str(seq_record.name))
            seq_string.append(str(seq_record.seq))
    seq = {seq_name[i]:seq_string[i] for i in range(len(seq_name))}
    print(shared_motif(seq))