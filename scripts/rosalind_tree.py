"""
Completing a Tree.
url: https://rosalind.info/problems/tree/

An undirected graph is connected if there is a path connecting any two nodes. A tree is a connected (undirected) graph containing no cycles; this definition forces the tree to have a branching structure organized around a central core of nodes, just like its living counterpart.
We have already grown familiar with trees in “Mendel's First Law”, where we introduced the probability tree diagram to visualize the outcomes of a random variable.
In the creation of a phylogeny, taxa are encoded by the tree's leaves, or nodes having degree 1. A node of a tree having degree larger than 1 is called an internal node.

Given: A positive integer n (n≤1000n≤1000) and an adjacency list corresponding to a graph on n nodes that contains no cycles.
Return: The minimum number of edges that can be added to the graph to produce a tree.
"""
def count_edges_to_form_tree(file_path: str) -> int:
    """
    Count the minimum number of edges required to convert a graph into a tree.

    Args:
        file_path (str): Path to the input file containing the adjacency list.

    Returns:
        int: Minimum number of edges needed to form a tree.
    """
    with open(file_path, 'r') as f:
        lines = f.readlines()
    
    n = int(lines[0].strip())  # Number of nodes
    # Number of edges is equal to the total lines - 1 (since first line is node count)
    number_of_edges = len(lines) - 1

    required_edges = n - 1  # A tree with n nodes must have exactly n - 1 edges
    
    edges_to_add = required_edges - number_of_edges  # Edges needed to be added?
    
    return edges_to_add

if __name__ == '__main__':
    file_path = '../data/rosalind_tree.txt'
    edges_needed = count_edges_to_form_tree(file_path)
    print(edges_needed)
