class Node:
  def __init__(self, g, id, b):
    self.g = g
    self.id = id
    self.b = b
    self.inarcs = {}
    self.outarcs = {}
    self.sources = {}
    self.sinks = {}
    self.yOpt = {}

  def addIn(self, arc):
    self.inarcs[arc.fromnode] = arc

  def addOut(self, arc):
    self.outarcs[arc.tonode] = arc

  def isSource(self):
    return self in self.g.sources

  def isPool(self):
    return self in self.g.pools

  def isSink(self):
    return self in self.g.sinks
