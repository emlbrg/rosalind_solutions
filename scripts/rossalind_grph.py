"""
Overlap Graphs.
url: https://rosalind.info/problems/grph/

A graph whose nodes have all been labeled can be represented by an adjacency list, in which each row of the list contains the two node labels corresponding to a unique edge.

A directed graph (or digraph) is a graph containing directed edges, each of which has an orientation. That is, a directed edge is represented by an arrow instead of a line segment; the starting and ending nodes of an edge form its tail and head, respectively. The directed edge with tail v and head w is represented by (v,w) (but not by (w,v)). A directed loop is a directed edge of the form (v,v).

For a collection of strings and a positive integer k, the overlap graph for the strings is a directed graph Ok in which each string is represented by a node, and string s is connected to string t with a directed edge when there is a length k suffix of s that matches a length k prefix of t, as long as s≠t; we demand s≠t to prevent directed loops in the overlap graph (although directed cycles may be present).

Given: A collection of DNA strings in FASTA format having total length at most 10 kbp.
Return: The adjacency list corresponding to O3. You may return edges in any order.
"""
from tools import parse_fasta

def overlap_graph(sequences, k):
    """
    Generate the overlap graph for the given sequences based on the specified overlap length k.

    Args:
        sequences (dict): A dictionary of sequences.
        k (int): The length of the suffix and prefix to compare.

    Returns:
        list: A list of directed edges representing the overlap graph.
    """
    edges = []
    
    # Convert keys to a list for iteration
    headers = list(sequences.keys())
    
    # Compare each pair of sequences
    for i in range(len(headers)):
        for j in range(len(headers)):
            if i != j:  # Prevent self-loops
                suffix = sequences[headers[i]][-k:]  # k-length suffix of s
                prefix = sequences[headers[j]][:k]  # k-length prefix of t
                if suffix == prefix:
                    edges.append((headers[i], headers[j]))

    return edges

def main(fasta_file, k):
    """
    Main function to read a FASTA file and generate the overlap graph.

    Args:
        fasta_file (str): Path to the FASTA file.
        k (int): Length of the overlap to consider.
    """
    sequences = parse_fasta(fasta_file)
    edges = overlap_graph(sequences, k)

    for edge in edges:
        print(f"{edge[0]} {edge[1]}")

if __name__ == '__main__':
     fasta_file_path = '../data/rosalind_grph.txt'
     overlap_length = 3
     main(fasta_file_path, overlap_length)
