##
# Minimum cost flow problem instance
#
from node import Node
from arc import Arc

class Graph:
  def __init__(self, filename):
    self.id = filename
    self.nodes = {}
    self.sources = {}
    self.sinks = {}
    self.arcs = {}
    self.readDimacs(filename)

  def readDimacs(self, filename):
    """
    Reads instance data (m = arc cardinality, n = node cardinality, h=capacity, b=supply(+)/demand(-), c=cost) from a text file on DIMACS format
    """
    reader = open(filename,'r')
    string = reader.readline().strip()
    while string[0] != 'p':
      string = reader.readline().strip()
    self.n = int(string.split()[2].strip())
    self.m = int(string.split()[3].strip())
    for i in range(self.n):
      self.nodes[i] = Node(self, i, 0)
    k=0
    while k<self.m:
      string = reader.readline().strip()
      if string[0] == 'n':
        i = int(string.split()[1].strip())-1
        b = int(string.split()[2].strip())
        node = self.nodes[i]
        node.b = b
        if b>0:
          self.sources[i] = node
        elif b<0:
          self.sinks[i] = node
      elif string[0] == 'a':
        i = int(string.split()[1].strip())-1
        j = int(string.split()[2].strip())-1
        h = int(string.split()[4].strip())
        c = int(string.split()[5].strip())
        self.arcs[i,j] = Arc(self, self.nodes[i], self.nodes[j], h, c)
        self.nodes[i].addOut(self.arcs[i,j])
        self.nodes[j].addIn(self.arcs[i,j])
        k += 1
    reader.close()

  def writeDimacs(self, filename):
    """
    Writes instance data to a file on DIMACS format
    """
    writer = open(filename,'w')
    writer.write('c Instance %s\n' % self.id)
    writer.write('p min %d %d\n' % (self.n, self.m))
    for i in self.nodes:
      writer.write('n %d %d' % (i+1, self.nodes[i].b))
    for (i,j) in self.arcs:
      writer.write('a %d %d 0 %d %d' % (i+1, j+1, self.arcs[i,j].h, self.arcs[i,j].c))
    writer.write('c End of file\n')
    writer.close()

