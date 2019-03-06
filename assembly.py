#!/usr/bin/env python

from Bio import SeqIO
import argparse
import networkx as nx
from itertools import product

library = "You are using %s networkx version. This program depends on version 2.2." % nx.__version__
parser = argparse.ArgumentParser(description='Assembly contigs based on reads from input file.', epilog=library)
parser.add_argument('infile', metavar=('input'), help='input fasta file with reads')
parser.add_argument('outfile', metavar=('output'), help='output fasta file with contigs')
args = parser.parse_args()

class DeBruijnGraph():
      def __init__(self, seq_iter, k):
         self.sequences = seq_iter
         self.k = k
         self.G = nx.DiGraph()

         for seq in self.sequences:
            for kmer in self.extract_kmers(seq):
               if kmer.l_node not in self.G.nodes():
                  self.G.add_node(kmer.l_node)
               if kmer.r_node not in self.G.nodes():
                  self.G.add_node(kmer.r_node)
               if (kmer.l_node, kmer.r_node) in self.G.edges():
                  self.G.edges[kmer.l_node, kmer.r_node]['weight'] += 1 
               else:
                  self.G.add_edge(kmer.l_node, kmer.r_node, weight = 1)

         self.subgraphs = self.extract_subgraphs()

      def extract_kmers(self, sequence):
         return [Kmer(sequence[i:i + self.k]) for i in range(len(sequence) - self.k + 1)]
   
      def extract_subgraphs(self):
         return [self.G.subgraph(c).copy() for c in nx.weakly_connected_components(self.G)]   

      def extract_contigs(self):
         contigs = []
         for sub in self.subgraphs:
            starts = get_start_nodes(sub)
            ends = get_end_nodes(sub)
            pairs = product(starts, ends)
            paths = []
            for pair in pairs:
               try:
                  path = nx.dijkstra_path(sub, *pair, weight = 'weight')
                  paths.append(''.join(short_seq[-1] for short_seq in path))
               except Exception as e:
                  pass
            if paths:
               contigs.append(max(paths, key=len))
         return contigs

class Kmer():
   def __init__(self, kmer_sequence):
      self.kmer_sequence = kmer_sequence
      self.k = len(self.kmer_sequence)
      self.l_node = self.kmer_sequence[:-1]
      self.r_node = self.kmer_sequence[1:]


def get_start_nodes(g):
   return [n for n in g if g.in_degree(n) == 0]

def get_end_nodes(g):
   return [n for n in g if g.out_degree(n) == 0]

with open(args.infile, "rU") as handle:
   records =  SeqIO.parse(handle, "fasta")
   sequences = [str(record.seq) for record in records]
   graph = DeBruijnGraph(sequences, k = 31) 
   contigs = graph.extract_contigs()
with open(args.outfile, 'w') as f:
    for item in contigs:
        f.write(">contig\n")
        f.write("%s\n" % item)

