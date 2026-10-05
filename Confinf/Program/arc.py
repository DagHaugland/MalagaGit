class Arc:
  def __init__(self, g, i, j, h, c):
    self.g = g
    self.fromnode = i
    self.tonode = j
    self.h = h
    self.c = c

