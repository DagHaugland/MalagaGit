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
    print('Reading %d nodes and %d arcs.' % (self.n,self.m))
    for i in range(self.n):
      self.nodes[i] = Node(self, i, 0)
    k=0
    while k<self.m:
      string = reader.readline().strip()
      if string[0] == 'n':
        i = int(string.split()[1].strip())-1
        b = int(string.split()[2].strip())
        node = self.nodes[i]
        nmax = max(self.nodes)
        if b>0:
          source = Node(self, nmax+1, b)
          self.nodes[nmax+1] = source
          self.sources[nmax+1] = source
          arc = Arc(self, source, node, b, 0)
          self.arcs[nmax+1,i] = arc
          source.addOut(arc)
          node.addIn(arc)
          k += 1
          self.n += 1
          self.m += 1
        elif b<0:
          sink = Node(self, nmax+1, b)
          self.nodes[nmax+1] = sink
          self.sinks[nmax+1] = sink
          arc = Arc(self, node, sink, -b, 0)
          self.arcs[i,nmax+1] = arc
          sink.addIn(arc)
          node.addOut(arc)
          k += 1
          self.n += 1
          self.m += 1
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
    print('Created graph with %d nodes, %d arcs, %d sources, and %d sinks.'
          % (self.n,self.m,len(self.sources),len(self.sinks)))
    self.connectSources()
    self.connectSinks()

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

  def connectSources(self): 
    nconnections = 0
    for s in self.sources.values():
      nconnections += self.connectSource(s)
    print('Average %.1f source connections per node' % (nconnections / len(self.nodes)))

  def connectSinks(self): 
    nconnections = 0
    for t in self.sinks.values():
      nconnections += self.connectSink(t)
    print('Average %.1f sink connections per node' % (nconnections / len(self.nodes)))

  def connectSource(self, s): 
    s.sources[s.id] = s
    q = [s]
    head = 0
    while head < len(q):
      node = q[head]
      for neigh in node.outarcs:
        if s.id not in neigh.sources:
          neigh.sources[s.id] = s
          q.append(neigh)
      head += 1
    return len(q)

  def realSource(self, s):
    if s.id in self.sources:
      for neigh in s.outarcs:
        return neigh
    return None

  def connectSink(self, t): 
    t.sinks[t.id] = t
    q = [t]
    head = 0
    while head < len(q):
      node = q[head]
      for neigh in node.inarcs:
        if t.id not in neigh.sinks:
          neigh.sinks[t.id] = t
          q.append(neigh)
      head += 1
    return len(q)

  def realSink(self, t):
    if t.id in self.sinks:
      for neigh in t.inarcs:
        return neigh
    return None

