"""
Distances in Trees
url: http://rosalind.info/problems/nwck/

Problem

Newick format is a way of representing trees even more concisely than using an adjacency list, especially when dealing with trees whose internal nodes have not been labeled.
First, consider the case of a rooted tree T. A collection of leaves v1,v2,…,vn of T are neighbors if they are all adjacent to some internal node u. Newick format for T is obtained by iterating the following key step: delete all the edges {vi,u} from T and label u with (v1,v2,…,vn)u. This process is repeated all the way to the root, at which point a semicolon signals the end of the tree.
A number of variations of Newick format exist. First, if a node is not labeled in T, then we simply leave blank the space occupied by the node. In the key step, we can write (v1,v2,…,vn) in place of (v1,v2,…,vn)u if the vi are labeled; if none of the nodes are labeled, we can write (,,…,).
A second variation of Newick format occurs when T is unrooted, in which case we simply select any internal node to serve as the root of T. A particularly peculiar case of Newick format arises when we choose a leaf to serve as the root.

Given: A collection of n trees (n≤40) in Newick format, with each tree containing at most 200 nodes; each tree Tk is followed by a pair of nodes xk and yk in Tk.
Return: A collection of n positive integers, for which the kth integer represents the distance between xk and yk in Tk.
"""

def dis_tree(T: str, x: str, y: str) -> int:
    """
    Calculate the distance between two nodes in a tree represented in Newick format.

    Args:
        T (str): Newick format string representing the tree.
        x (str): Name of the first node.
        y (str): Name of the second node.

    Returns:
        int: The distance between the two nodes.
    """
    x_index = T.find(x)
    y_index = T.find(y)
    t = [i for i in T[min(x_index, y_index):max(x_index, y_index)] if i in [')','(',',']]
    bracket = ''
    for i in t:
        bracket += i
    while '(,)' in bracket:
        bracket = bracket.replace('(,)','')
    if bracket.count('(') == len(bracket):
        return len(bracket)
    elif bracket.count(')') == len(bracket):
        return len(bracket)
    elif bracket.count(',') == len(bracket):
        return 2
    else:
        return bracket.count(')') + bracket.count('(') + 2

if __name__ == '__main__':
    with open('../data/rosalind_nwck.txt', 'r') as f:
        tree = [line.strip().replace(';','') for line in f.readlines() if len(line.strip()) > 0]
    for i in range(0, len(tree), 2):
        T = tree[i]
        x, y = tree[i+1].split(' ')
        print(dis_tree(T, x, y), end=" ")
