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
    self.model = None

  def makeModel(self,relax=False):
    # Create a minimization model:
    self.model = gurobipy.Model('cimcf')
    self.model.modelSense = GRB.MINIMIZE
    # self.model.Params.outputFlag = 0
    # Create decision variables:
    self.makeVars(relax)
    # Create constraints:
    self.makeConstrs()

  def makeVars(self,relax=False):
    # Flow variables:
    self.x = {}
    self.y = {}
    for ((i,j),arc) in self.g.arcs.items():
      self.x[arc] = self.model.addVar(lb=0.0, ub=arc.h, obj=arc.c, name = 'x_%d_%d' % (i,j))
    # Binary selection variables:
    for (i,node) in self.g.nodes.items():
      for (s,source) in node.sources.items():
        if relax:
          self.y[node,source] = self.model.addVar(lb=0, ub=1, obj=0, name = 'y_%d^%d' % (i,s))
        else:
          self.y[node,source] = self.model.addVar(vtype=GRB.BINARY, obj=0, name = 'y_%d^%d' % (i,s))
        if s==i:
          self.y[node,source].setAttr('lb',1)
      for (t,sink) in node.sinks.items():
        if relax:
          self.y[node,sink] = self.model.addVar(lb=0, ub=1, obj=0, name = 'y_%d^%d' % (i,t))
        else:
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
      for (s,source) in node.sources.items():
        lhs += self.y[node,source]
      for (t,sink) in node.sinks.items():
        lhs += self.y[node,sink]
      self.conv[i] = self.model.addConstr(lhs==1, name = 'conv%d' % i)

  def close(self):
    # Close arc if their end nodes are not assigned to compatible sources/sinks:
    self.closeArc = {}
    for arc in self.g.arcs.values():
      (fromnode,tonode) = (arc.fromnode,arc.tonode)
      for (s,source) in tonode.sources.items():
        ysum = LinExpr()
        if s in fromnode.sources:
          ysum += self.y[fromnode,source]
        for (ss,ext) in tonode.sources.items():
          if ss!=s:
            ysum += self.y[tonode,ext]
        for sink in tonode.sinks.values():
          ysum += self.y[tonode,sink]
        self.closeArc[arc,source] = self.model.addConstr(self.x[arc] <= arc.h*ysum, name = 'close%d_%d^%d'
          % (fromnode.id,tonode.id,s))
      for (t,sink) in self.g.sinks.items():
        ysum = LinExpr()
        if t in tonode.sinks:
          ysum += self.y[tonode,sink]
        for (tt,ext) in fromnode.sinks.items():
          if tt!=t:
            ysum += self.y[fromnode,ext]
        for source in fromnode.sources.values():
          ysum += self.y[fromnode,source]
        self.closeArc[arc,sink] = self.model.addConstr(self.x[arc] <= arc.h*ysum, name = 'close%d_%d^%d'
          % (fromnode.id,tonode.id,t))

  def makeConstrs(self):
    # Create all constraints:
    self.consFlow()
    self.convexity()
    self.close()

  def nullObj(self):
    for arc in self.g.arcs.values():
      self.x[arc].setAttr('obj',0)

  def restoreObj(self):
    for arc in self.g.arcs.values():
      self.x[arc].setAttr('obj',arc.c)

  def maxFlow(self):
    done = False
    it = 0
    if self.model is None:
      self.makeModel(True)
    self.nullObj()
    while not (done or it):
      it += 1
      done = True
      aCount=0
      for ((i,j),arc) in sorted(self.g.arcs.items()):
        self.model.Params.outputFlag = 0
        aCount += 1
        self.x[arc].setAttr('obj',-1)
        self.model.optimize()
        flow = self.x[arc].getAttr('x')
        if arc.h - flow > 0.001:
          if arc.h - flow > 0.1*arc.h:
            print('max flow iteration %d,%d: (%d,%d) %.1f -> %.1f'
                  % (it, aCount, arc.fromnode.id, arc.tonode.id, arc.h, flow))
          arc.h = flow
          self.makeModel(True)
          self.nullObj()
          done = False
        else:
          self.x[arc].setAttr('obj',0)
      print('max flow iteration %d done' % it)
    self.restoreObj()
    print('max flow done')

  def getSol(self):
    # Retrieve the optimal solution from the solver
    for (arc,x) in self.x.items():
      arc.flow = x.getAttr('x')
    for ((node,st),y) in self.y.items():
      node.yOpt[st] = y.getAttr('x')
      if node.yOpt[st] > 0.5:
        node.st = st

  def solve(self,relax=False):
    # self.maxFlow()
    self.makeModel(relax)
    # self.addCuts()
    self.model.optimize()
    self.getSol()
    if relax:
      print('RELAXED optimal objective function value: %d found in \t%.1f seconds.\n' % (self.model.getObjective().getValue(), self.model.Runtime))
    else:
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
    self.fcut = {}
    for arc in self.g.arcs.values():
      (sti,stj) = (arc.fromnode.st, arc.tonode.st)
      if sti.isSource() and stj.isSink():
        self.cut[arc] = arc.flow
        if arc.flow > 0.01:
          self.fcut[arc] = arc.flow
        self.flow[2] += arc.flow

  def checkArcs(self):
    print('Arc checks (sources):')
    bad = 0
    for ((i,j),arc) in sorted(self.g.arcs.items()):
      found = 0
      if arc.flow > 0.01:
        (fromnode, tonode) = (arc.fromnode, arc.tonode)
        for source in tonode.sources.values():
          fromy = 0.0
          if source.id in fromnode.sources:
            fromy = fromnode.yOpt[source]
          toy = tonode.yOpt[source]
          if toy - fromy > 0.001:
            found = 1
            print('Arc check: %d %d %d: %.1f %.1f %.2f %.2f' % (fromnode.id, tonode.id, self.g.realSource(source).id, arc.flow, arc.h, fromy, toy))
      bad += found
    print('Found %d bad arcs.' % bad)
    print('Arc checks (sinks):')
    bad = 0
    for ((i,j),arc) in sorted(self.g.arcs.items()):
      found = 0
      if arc.flow > 0.001:
        (fromnode, tonode) = (arc.fromnode, arc.tonode)
        for sink in tonode.sinks.values():
          toy = 0.0
          if sink.id in tonode.sinks:
            toy = tonode.yOpt[sink]
          fromy = fromnode.yOpt[sink]
          if fromy - toy > 0.01:
            found = 1
            print('Arc check: %d %d %d: %.1f %.1f %.2f %.2f' % (fromnode.id, tonode.id, self.g.realSink(sink).id, arc.flow, arc.h, fromy, toy))
      bad += found
    print('Found %d bad arcs.' % bad)
  
  def collectArcs(self):
    collection = {}
    for ((i,j),arc) in sorted(self.g.arcs.items()):
      if arc.flow > 0.001:
        (fromnode, tonode) = (arc.fromnode, arc.tonode)
        for source in tonode.sources.values():
          fromy = 0.0
          if source.id in fromnode.sources:
            fromy = fromnode.yOpt[source]
          toy = tonode.yOpt[source]
          if toy - fromy > 0.001:
            (stlist,ysum) = collection.get(arc,([],1.0))
            stlist.append(source)
            ysum -= (toy-fromy)
            collection[arc] = (stlist,ysum)
    for ((i,j),arc) in sorted(self.g.arcs.items()):
      if arc.flow > 0.001:
        (fromnode, tonode) = (arc.fromnode, arc.tonode)
        for sink in tonode.sinks.values():
          toy = 0.0
          if sink.id in tonode.sinks:
            toy = tonode.yOpt[sink]
          fromy = fromnode.yOpt[sink]
          if fromy - toy > 0.001:
            (stlist,ysum) = collection.get(arc,([],1.0))
            stlist.append(sink)
            ysum -= (fromy-toy)
            collection[arc] = (stlist,ysum)
    return collection

  def makevi(self):
    collection = self.collectArcs()
    for (arc,(stlist,ysum)) in collection.items():
      print(arc.fromnode.id, arc.tonode.id, arc.flow)
      (fromnode, tonode) = (arc.fromnode, arc.tonode)
      if arc.flow > arc.h*ysum + 0.001:
        yvarsum = LinExpr(1.0)
        for u in stlist:
          if u.id in self.g.sources:
            if u in fromnode.sources.values():
              yvarsum += self.y[fromnode,u]
            yvarsum -= self.y[tonode,u]
          elif u.id in self.g.sinks:
            if u in tonode.sinks.values():
              yvarsum += self.y[tonode,u]
            yvarsum -= self.y[fromnode,u]
          else:
            print('Illegal source/sink', u.id)
        self.vi.append(self.model.addConstr(self.x[arc] <= arc.h*yvarsum, name = 'vi%d' % len(self.vi)))

  def addCuts(self):
    self.vi = []
    if self.model is None:
      self.makeModel(True)
    else:
      for y in self.y.values():
        y.setAttr('vtype', 'C')
    done = False
    it = 0
    # self.maxFlow()
    flag = self.model.Params.outputFlag
    self.model.Params.outputFlag = 0
    while not done and it<999:
      it += 1
      self.model.optimize()
      self.getSol()
      if it==1:
        self.checkArcs()
      oldlen = len(self.vi)
      self.makevi()
      print('Iteration: %d. Lower bound: %.1f. Cuts: %d' % (it, self.model.getObjective().getValue(), len(self.vi)))
      done = len(self.vi) == oldlen
    self.checkArcs()
    for y in self.y.values():
      y.setAttr('vtype', 'B')
    self.model.Params.outputFlag = flag
  
if __name__ == '__main__':
  if len(sys.argv) > 2 and sys.argv[2] == 'cut':
    cimcf = CIMCF('../Data/' + sys.argv[1])
    cimcf.addCuts()
  elif len(sys.argv) == 2 or len(sys.argv) > 2 and sys.argv[2] != 'cut':
    cimcf = CIMCF('../Data/' + sys.argv[1])
    relax = False
    cimcf.solve(relax)
    if relax:
      cimcf.checkArcs()
    else:
      cimcf.stat()
      print('%d nodes, %d sources, %d sinks, %d arcs' % (cimcf.g.n, len(cimcf.g.sources), len(cimcf.g.sinks), cimcf.g.m))
      print('%d source associations, %d sink associations, %d cut arcs, %d flow cut arcs.' % (cimcf.single[0], cimcf.single[0], len(cimcf.cut), len(cimcf.fcut))) 
      print('flow check: ', cimcf.flow)
  else:  
    print('Usage: python3 confinf.py filename ["cut"]')
