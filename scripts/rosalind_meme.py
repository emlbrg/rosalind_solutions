"""
New Motif Discovery
url: https://rosalind.info/problems/meme/

Problem

The novel-motif finding tool MEME can be found here.

Given: A set of protein strings in FASTA format that share some motif with minimum length 20.
Return: Regular expression for the best-scoring motif.
"""

# meme ../data/rosalind_meme.txt -protein -oc ./rosalinf_meme -nostatus -time 14400 -mod zoops -nmotifs 1 -minw 6 -maxw 50 -objfun classic -markov_order 0
# Results: SATWDGLKSNYNEYS[NVT][FHC]GIDK[IG]DKWMYAPAAYRMY