#!/Users/hekate/anaconda3/bin/python3.6

from Bio import SeqIO
from collections import Counter, defaultdict

#reference_seq = "/Users/hekate/2018Z/TSG/ex2/assignment2/reference/reference_changed.fasta"
reads_file = "/Users/hekate/2018Z/TSG/ex2/assignment2/reads/reads0.fasta"
#reads_file = reference_seq

class DeBruijnGraph:
   ''' De Bruijn directed multigraph built from a collection of
      strings. User supplies strings and k-mer length k.  Nodes
      are k-1-mers.  An Edge corresponds to the k-mer that joins
      a left k-1-mer to a right k-1-mer. '''
 
   @staticmethod
   def chop(st, k):
      ''' Chop string into k-mers of given length '''
      for i in range(len(st)-(k-1)):
         yield (st[i:i+k], st[i:i+k-1], st[i+1:i+k])
   
   class Node:
      ''' Node representing a k-1 mer.  Keep track of # of
         incoming/outgoing edges so it's easy to check for
         balanced, semi-balanced. '''
      
      def __init__(self, km1mer):
         self.km1mer = km1mer
         self.nin = 0
         self.nout = 0
      
      def isSemiBalanced(self):
         return abs(self.nin - self.nout) == 1
      
      def isBalanced(self):
         return self.nin == self.nout
      
      def __hash__(self):
         return hash(self.km1mer)
      
      def __str__(self):
         return self.km1mer
   
   def __init__(self, strIter, k, circularize=False):
      ''' Build de Bruijn multigraph given string iterator and k-mer
         length k '''
      self.G = {}    # multimap from nodes to neighbors
      self.nodes = {} # maps k-1-mers to Node objects
      ###
      self.weights = defaultdict(int)
      for st in strIter:
         if circularize:
            st += st[:k-1]
         for kmer, km1L, km1R in self.chop(st, k):
            nodeL, nodeR = None, None
            if km1L in self.nodes:
               nodeL = self.nodes[km1L]
            else:
               nodeL = self.nodes[km1L] = self.Node(km1L)
            if km1R in self.nodes:
               nodeR = self.nodes[km1R]
            else:
               nodeR = self.nodes[km1R] = self.Node(km1R)
            ###
            #nodeL.nout += 1 
            #nodeR.nin += 1
            # self.G.setdefault(nodeL, []).append(nodeR)
            self.G.setdefault(nodeL, [])
            if nodeR not in self.G[nodeL]:
               self.G[nodeL].append(nodeR)
               nodeL.nout += 1
               nodeR.nin += 1
               self.weights[nodeL] = 1
            else:
               self.weights[nodeL] += 1
      # Iterate over nodes; tally # balanced, semi-balanced, neither
      self.nsemi, self.nbal, self.nneither = 0, 0, 0
      # Keep track of head and tail nodes in the case of a graph with
      # Eularian walk (not cycle)
      self.head, self.tail = None, None
      for node in iter(self.nodes.values()):
         if node.isBalanced():
            self.nbal += 1
         elif node.isSemiBalanced():
            if node.nin == node.nout + 1:
               self.tail = node
            if node.nin == node.nout - 1:
               self.head = node
            self.nsemi += 1
         else:
            self.nneither += 1
   
   def nnodes(self):
      ''' Return # nodes '''
      return len(self.nodes)
   
   def nedges(self):
      ''' Return # edges '''
      return len(self.G)
   
   def hasEulerianWalk(self):
      ''' Return true iff graph has Eulerian walk. '''
      return self.nneither == 0 and self.nsemi == 2
   
   def hasEulerianCycle(self):
      ''' Return true iff graph has Eulerian cycle. '''
      return self.nneither == 0 and self.nsemi == 0
   
   def isEulerian(self):
      ''' Return true iff graph has Eulerian walk or cycle '''
      # technically, if it has an Eulerian walk
      return self.hasEulerianWalk() or self.hasEulerianCycle()
   
   def eulerianWalkOrCycle(self):
      ''' Find and return sequence of nodes (represented by
         their k-1-mer labels) corresponding to Eulerian walk
         or cycle '''
      assert self.isEulerian()
      g = self.G
      if self.hasEulerianWalk():
         g = g.copy()
         g.setdefault(self.tail, []).append(self.head)
      # graph g has an Eulerian cycle
      tour = []
      src = next(iter(g.keys())) # pick arbitrary starting node
      
      def __visit(n):
         while len(g[n]) > 0:
            dst = g[n].pop()
            __visit(dst)
         tour.append(n)
      __visit(src)
      tour = tour[::-1][:-1] # reverse and then take all but last node
         
      if self.hasEulerianWalk():
         # Adjust node list so that it starts at head and ends at tail
         sti = tour.index(self.head)
         tour = tour[sti:] + tour[:sti]
      
      # Return node list
      return list(map(str, tour))

   def condense(self):
      g = self.G.copy()
      for nodeL, nodesR in iter(self.G.items()):
         if nodeL.nout == 1:
            nodeR = nodesR[0]
            nodeL.km1mer += nodeR.km1mer[-1]
            nodeL.nout = nodeR.nout
            self.nodes[nodeL.km1mer] = nodeL
            if nodeR in self.G:
               g[nodeL] = self.G[nodeR]
               if nodeR in g:
                  g.pop(nodeR)
               if nodeR.km1mer in self.nodes:
                  self.nodes.pop(nodeR.km1mer)
      return g

def scout(start, graph, path = []):
   path += [start]
   next_node = graph[start][0]
   if next_node.nout == 1:
      new_path = scout(next_node, graph, path)
      if new_path: return new_path
   else:
      new_path += [next_node]
   return new_path

with open(reads_file, "rU") as handle:
   records =  SeqIO.parse(handle, "fasta")
   graph = DeBruijnGraph([str(record.seq) for record in records], k = 30)
   print(graph.nnodes())
   for i in range(100):
      #singles = [k for k,v in graph.G.items() if k.nout == 1]
      #print(len(singles)) 
      new_g = graph.condense()
      graph.G = new_g
   print(graph.nnodes())
   for nodeL, nodeR in iter(graph.G.items()):
      print(nodeL, nodeR)
   #start_nodes = [node for node in iter(graph.G.keys()) if node.nin > 1 or node.nout > 1]
   '''
   for nodeL, Rnodes in iter(graph.G.items()):
      if nodeL.nout == 1:
         nodeR = Rnodes[0]
         nodeL.km1mer += Rnodes[0].km1mer[-1]
         graph.G[nodeL] = graph.G[nodeR]
         del(graph.G[nodeR])
         #del(graph.nodes[nodeR.km1mer] jeśli nie ma już wchodzących i wychodzących to nara
      #for nodeR in Rnodes:
      #   if nodeR.nout == 0:
      #      print(nodeR)    
   #   print("Node: ", node)
   #   print("Wieght: ", graph.weights[node])
   #   print("in/out: ", node.nin, node.nout)
   #for Lnode, Rnodes in iter(graph.G.items()):
   #   print(Lnode, Rnodes[0])
   '''

   print(graph.isEulerian(), graph.hasEulerianWalk(), graph.hasEulerianCycle())
   #walk = graph.eulerianWalkOrCycle()
   #print(walk)
