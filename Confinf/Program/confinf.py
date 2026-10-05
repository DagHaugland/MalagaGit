##
# Formulations of the confluent-influent minimum cost flow problem
#
import sys
import gurobipy
from gurobipy import GRB
from gurobipy import LinExpr
from node import Node
from arc import Arc
from graph import Graph
status = {2:'OPT', 3:'INF', 6:'CUT', 9:'TIME'}
formulation = {False:'STRONG', True:'WEAK'}

class CIMCF:
  def __init__(self, filename):
    self.g = Graph(filename)

  def makeModel(self):
    # Create a minimization model:
    self.model = gurobipy.Model('cimcf')
    self.model.modelSense = GRB.MINIMIZE
    # self.model.Params.outputFlag = 0
    # Create decision variables:
    self.makeVars()
    # Create constraints:
    self.makeConstrs()

  def makeVars(self):
    # Flow variables:
    self.x = {}
    self.y = {}
    for ((i,j),arc) in self.g.arcs.items():
      self.x[arc] = self.model.addVar(lb=0.0, ub=arc.h, obj=arc.c, name = 'x_%d_%d' % (i,j))
    # Binary selection variables:
    for (i,node) in self.g.nodes.items():
      for (s,source) in self.g.sources.items():
        self.y[node,source] = self.model.addVar(vtype=GRB.BINARY, obj=0, name = 'y_%d^%d' % (i,s))
        if s==i:
          self.y[node,source].setAttr('lb',1)
      for (t,sink) in self.g.sinks.items():
        self.y[node,sink] = self.model.addVar(vtype=GRB.BINARY, obj=0, name = 'y_%d^%d' % (i,t))
        if t==i:
          self.y[node,sink].setAttr('lb',1)

  def consFlow(self):
    # Conservation of flow:
    self.cons = {}
    for (i,node) in self.g.nodes.items():
      lhs = LinExpr()
      for (j,arc) in node.outarcs.items():
        lhs += self.x[arc]
      for (j,arc) in node.inarcs.items():
        lhs -= self.x[arc]
      if node.b:
        self.cons[i] = self.model.addConstr(lhs<=node.b, name = 'cons%d' % i)
      else:
        self.cons[i] = self.model.addConstr(lhs==0, name = 'cons%d' % i)

  def convexity(self):
    # Associate a source or a sink
    self.conv = {}
    for (i,node) in self.g.nodes.items():
      lhs = LinExpr()
      for (s,source) in self.g.sources.items():
        lhs += self.y[node,source]
      for (t,sink) in self.g.sinks.items():
        lhs += self.y[node,sink]
      self.conv[i] = self.model.addConstr(lhs==1, name = 'conv%d' % i)

  def close(self):
    # Close arc if their end nodes are not assigned to compatible sources/sinks:
    self.close = {}
    for ((i,j),arc) in self.g.arcs.items():
      for (s,source) in self.g.sources.items():
        ysum = LinExpr()
        ysum += self.y[arc.fromnode,source]
        for (ss,ext) in self.g.sources.items():
          if ss!=s:
            ysum += self.y[arc.tonode,ext]
        for (t,sink) in self.g.sinks.items():
          ysum += self.y[arc.tonode,sink]
        self.close[arc,source] = self.model.addConstr(self.x[arc] <= arc.h*ysum, name = 'close%d_%d^%d' % (i,j,s))
      for (t,sink) in self.g.sinks.items():
        ysum = LinExpr()
        ysum += self.y[arc.tonode,sink]
        for (tt,ext) in self.g.sinks.items():
          if tt!=t:
            ysum += self.y[arc.fromnode,ext]
        for (s,source) in self.g.sources.items():
          ysum += self.y[arc.fromnode,source]
        self.close[arc,sink] = self.model.addConstr(self.x[arc] <= arc.h*ysum, name = 'close%d_%d^%d' % (i,j,t))

  def makeConstrs(self):
    # Create all constraints:
    self.consFlow()
    self.convexity()
    self.close()

  def getSol(self):
    # Retrieve the optimal solution from the solver
    for (arc,x) in self.x.items():
      arc.flow = x.getAttr('x')
    for ((node,st),y) in self.y.items():
      if y.getAttr('x') > 0.5:
        node.st = st

  def solve(self):
    self.makeModel()
    self.model.optimize()
    self.getSol()
    print('Optimal objective function value: %d found in \t%.1f seconds.\n' % (self.model.getObjective().getValue(), self.model.Runtime))

  def stat(self):
    # Statistics on the values of y and the cut
    self.flow = [0,0,0]
    self.single = [0,0]
    for node in self.g.nodes.values():
      if node.isSource():
        for arc in node.outarcs.values():
          self.flow[0] += arc.flow
      elif node.isSink():
        for arc in node.inarcs.values():
          self.flow[1] += arc.flow
      else:
        self.single[node.st.isSink()] += 1
    self.cut = {}
    for arc in self.g.arcs.values():
      (sti,stj) = (arc.fromnode.st, arc.tonode.st)
      if sti.isSource() and stj.isSink():
        self.cut[arc] = arc.flow
        self.flow[2] += arc.flow
  
if __name__ == '__main__':
  if len(sys.argv) > 1:
    cimcf = CIMCF('../Data/' + sys.argv[1])
    cimcf.solve()
    cimcf.stat()
    print('%d nodes, %d sources, %d sinks, %d arcs' % (cimcf.g.n, len(cimcf.g.sources), len(cimcf.g.sinks), cimcf.g.m))
    print('%d source associations, %d sink associations, %d cut arcs' % (cimcf.single[0], cimcf.single[0], len(cimcf.cut))) 
    print('flow check: ', cimcf.flow)
