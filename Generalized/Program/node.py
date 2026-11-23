class Node:
  def __init__(self, g, id, bL=0.0, b=None):
    self.g = g
    self.id = id
    self.b = b
    self.bL = bL
    self.inarcs = {}
    self.outarcs = {}
    self.sources = {}
    self.pools = {}
    self.sinks = {}

  def addIn(self, arc):
    self.inarcs[arc.fromnode] = arc

  def addOut(self, arc):
    self.outarcs[arc.tonode] = arc

  def isSource(self):
    return self.id in self.g.sources

  def isPool(self):
    return self.id in self.g.pools

  def isSink(self):
    return self.id in self.g.sinks

class Source(Node):
  def __init__(self, g, id, bL=0.0, b=None, q=None):
    super().__init__(g, id, bL, b)
    self.sources = {id:self}
    self.q = list()
    if q is not None:
      self.q = list(q)

class Pool(Node):
  def __init__(self, g, id, bL=0.0, b=None):
    super().__init__(g, id, bL, b)
    self.ylb = {}
    self.yub = {}

class Sink(Node):
  def __init__(self, g, id, bL=0.0, b=None, qL=None, q=None):
    super().__init__(g, id, bL, b)
    self.sinks = {id:self}
    self.q = list()
    if q is not None:
      self.q = list(q)
    self.qL = list()
    if qL is not None:
      self.qL = list(qL)
    self.ylb = {}
    self.yub = {}

