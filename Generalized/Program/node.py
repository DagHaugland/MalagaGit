class Node:
  def __init__(self, g, id, b):
    self.g = g
    self.id = id
    self.b = b
    self.inarcs = {}
    self.outarcs = {}

  def addIn(self, arc):
    self.inarcs[arc.fromnode] = arc

  def addOut(self, arc):
    self.outarcs[arc.tonode] = arc

  def isSource(self):
    return self.b > 0

  def isSink(self):
    return self.b < 0
