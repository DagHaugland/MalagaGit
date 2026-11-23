class Arc:
  def __init__(self, g, fromnode, tonode, c=None, h=None):
    self.g = g
    self.fromnode = fromnode
    self.tonode = tonode
    self.c = c
    self.h = h

  def updateCap(self):
    if self.h is None:
      self.h = min(self.fromnode.b, self.tonode.b)
    else:
      self.h = min(self.h, self.fromnode.b, self.tonode.b)
