##
# Alfaki's model for the generalized pooling problem
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

class GPP:
  def __init__(self, filename, formulation='s'):
    self.g = Graph(filename)
    self.instance = filename
    while '/' in self.instance:
      pos = self.instance.index('/')
      self.instance = self.instance[pos+1:]
    if '.' in self.instance:
      pos = self.instance.index('.')
      self.instance = self.instance[:pos]
    self.formulation = formulation.lower()
    self.model = None

  def makeModel(self,relax=False):
    if 's' not in self.formulation and 't' not in self.formulation:
      print('Invalid formulation: %s. Valid formulations: "s", "t", and "st"' % self.formulation)
      return
    # Create a minimization model:
    self.model = gurobipy.Model('gpp'+self.instance)
    self.model.modelSense = GRB.MINIMIZE
    # self.model.Params.outputFlag = 0
    self.model.Params.timeLimit = 300
    # Create decision variables:
    self.makeVars(relax)
    # Create constraints:
    self.makeConstrs()

  def makeVars(self,relax=False):
    # Flow variables:
    self.x = {}
    self.f = {}
    for ((i,j),arc) in self.g.arcs.items():
      self.f[arc] = self.model.addVar(lb=0.0, ub=arc.h, obj=arc.c, name = 'f_%d_%d' % (i,j))
      if 's' in self.formulation:
        for source in arc.fromnode.sources.values():
          self.x[arc,source] = self.model.addVar(lb=0.0, ub=arc.h, obj=0.0, name = 'x_%d_%d^%d' % (i,j,source.id))
      if 't' in self.formulation:
        for sink in arc.tonode.sinks.values():
          self.x[arc,sink] = self.model.addVar(lb=0.0, ub=arc.h, obj=0.0, name = 'x_%d_%d^%d' % (i,j,sink.id))
    # Proportion variables:
    self.y = {}
    for pool in self.g.pools.values():
      if 's' in self.formulation:
        for source in pool.sources.values():
          self.y[pool,source] = self.model.addVar(lb=0, ub=1, obj=0, name = 'y_%d^%d' % (pool.id,source.id))
      if 't' in self.formulation:
        for sink in pool.sinks.values():
          self.y[pool,sink] = self.model.addVar(lb=0, ub=1, obj=0, name = 'y_%d^%d' % (pool.id,sink.id))

  def capacity(self):
    # Capacity constraints:
    self.cap = {}
    self.capL = {}
    for node in self.g.nodes.values():
      lhs = LinExpr()
      if node.isSink():
        for arc in node.inarcs.values():
          lhs += self.f[arc]
      else:
        for arc in node.outarcs.values():
          lhs += self.f[arc]
      self.cap[node] = self.model.addConstr(lhs<=node.b, name = 'cap%d' % node.id)
      if node.bL > 0.0005:
        self.capL[node] = self.model.addConstr(lhs>=node.bL, name = 'capL%d' % node.id)

  def consFlow(self):
    # Conservation of flow:
    self.cons = {}
    for node in self.g.nodes.values():
      if 's' in self.formulation:
        for source in node.sources.values():
          lhs = LinExpr()
          for arc in node.outarcs.values():
            lhs -= self.x[arc,source]
          for (fromnode,arc) in node.inarcs.items():
            if source.id in fromnode.sources:
              lhs += self.x[arc,source]
          if node==source:
            self.cons[node] = self.model.addConstr(lhs<=0, name = 'cons%d^%d' % (node.id,source.id))
          elif node.isSink():
            self.cons[node] = self.model.addConstr(lhs>=0, name = 'cons%d^%d' % (node.id,source.id))
          else:
            self.cons[node] = self.model.addConstr(lhs==0, name = 'cons%d^%d' % (node.id,source.id))
      if 't' in self.formulation:
        for sink in node.sinks.values():
          lhs = LinExpr()
          for (tonode,arc) in node.outarcs.items():
            if sink.id in tonode.sinks:
              lhs -= self.x[arc,sink]
          for arc in node.inarcs.values():
            lhs += self.x[arc,sink]
          if node.isSource():
            self.cons[node] = self.model.addConstr(lhs<=0, name = 'cons%d^%d' % (node.id,sink.id))
          elif node==sink:
            self.cons[node] = self.model.addConstr(lhs>=0, name = 'cons%d^%d' % (node.id,sink.id))
          else:
            self.cons[node] = self.model.addConstr(lhs==0, name = 'cons%d^%d' % (node.id,sink.id))

  def concentration(self):
    # Upper concentration bounds at sinks:
    self.conc = {}
    self.concL = {}
    for sink in self.g.sinks.values():
      for k in range(self.g.nK):
        if 's' in self.formulation:
          lhs = LinExpr()
          lhsL = LinExpr()
          for arc in sink.inarcs.values():
            for source in arc.fromnode.sources.values():
              lhs += (source.q[k]-sink.q[k]) * self.x[arc,source]
              if len(sink.qL):
                lhsL += (source.q[k]-sink.qL[k]) * self.x[arc,source]
          self.conc[sink,k] = self.model.addConstr(lhs<=0, name = 'concT%d^%d' % (sink.id,k))
          if len(sink.qL):
            self.concL[sink,k] = self.model.addConstr(lhsL>=0, name = 'concTL%d^%d' % (sink.id,k))
        if 't' in self.formulation:
          lhs = LinExpr()
          lhsL = LinExpr()
          for source in sink.sources.values():
            for (tonode,arc) in source.outarcs.items():
              if sink in tonode.sinks.values():
                lhs += (source.q[k]-sink.q[k]) * self.x[arc,sink]
                if len(sink.qL):
                  lhsL += (source.q[k]-sink.qL[k]) * self.x[arc,sink]
          self.conc[sink,k] = self.model.addConstr(lhs<=0, name = 'concS%d^%d' % (sink.id,k))
          if len(sink.qL):
            self.concL[sink,k] = self.model.addConstr(lhsL>=0, name = 'concSL%d^%d' % (sink.id,k))

  def convexity(self):
    # Proportions add up to 1.0
    self.conv = {}
    for pool in self.g.pools.values():
      if 's' in self.formulation:
        lhs = LinExpr()
        for source in pool.sources.values():
          lhs += self.y[pool,source]
        self.conv[pool] = self.model.addConstr(lhs==1, name = 'convS%d' % pool.id)
      if 't' in self.formulation:
        lhs = LinExpr()
        for sink in pool.sinks.values():
          lhs += self.y[pool,sink]
        self.conv[pool] = self.model.addConstr(lhs==1, name = 'convT%d' % pool.id)

  def accumulation(self):
    # Arc flow is an accumulation of commodity flow
    self.accu = {}
    for arc in self.g.arcs.values():
      (fromnode,tonode) = (arc.fromnode,arc.tonode)
      if 's' in self.formulation:
        lhs = LinExpr()
        for source in fromnode.sources.values():
          lhs += self.x[arc,source]
        self.accu[arc] = self.model.addConstr(lhs==self.f[arc], name = 'accuS%d_%d' % (fromnode.id,tonode.id))
      if 't' in self.formulation:
        lhs = LinExpr()
        for sink in tonode.sinks.values():
          lhs += self.x[arc,sink]
        self.accu[arc] = self.model.addConstr(lhs==self.f[arc], name = 'accuT%d_%d' % (fromnode.id,tonode.id))

  def validIneq(self):
    # Redundant valid inequalities
    self.vi = {}
    for pool in self.g.pools.values():
      if 's' in self.formulation:
        for source in pool.sources.values():
          lhs = LinExpr()
          for arc in pool.outarcs.values():
            lhs += self.x[arc,source]
          self.vi[pool,source] = self.model.addConstr(lhs<=pool.b*self.y[pool,source], name = 'vi%d^%d' % (pool.id,source.id))
      if 't' in self.formulation:
        for sink in pool.sinks.values():
          lhs = LinExpr()
          for arc in pool.inarcs.values():
            lhs += self.x[arc,sink]
          self.vi[pool,sink] = self.model.addConstr(lhs<=pool.b*self.y[pool,sink], name = 'vi%d^%d' % (pool.id,sink.id))

  def bilinear(self):
    # Connection between arc flow and commodity flow
    self.bl = {}
    for pool in self.g.pools.values():
      if 's' in self.formulation:
        for source in pool.sources.values():
          for (tonode,arc) in pool.outarcs.items():
            self.bl[arc,source] = self.model.addConstr(self.x[arc,source]==self.y[pool,source]*self.f[arc], name = 'bl%d_%d^%d'
                                % (pool.id,tonode.id,source.id))
      if 't' in self.formulation:
        for sink in pool.sinks.values():
          for (fromnode,arc) in pool.inarcs.items():
            self.bl[arc,sink] = self.model.addConstr(self.x[arc,sink]==self.y[pool,sink]*self.f[arc], name = 'bl%d_%d^%d'
                                % (pool.id,fromnode.id,sink.id))

  def fracbound(self):
    # Grothey: Bounds on (source,pool) and (source,sink) fractions
    self.ylbound = {}
    self.yubound = {}
    for ((i,j),arc) in self.g.arcs.items():
      (fromnode, tonode) = (arc.fromnode, arc.tonode)
      if fromnode.isSource():
        (ylb, yub) = (tonode.ylb.get(fromnode,0.0), tonode.yub.get(fromnode,1.0))
        rhs = LinExpr()
        for a in tonode.inarcs.values():
          rhs += self.f[a]
        if ylb >= 0.001:
          self.ylbound[arc] = self.model.addConstr(self.f[arc]>=ylb*rhs, name = 'grotheyl%d^%d' % (i,j))
        if yub <= 0.999:
          self.yubound[arc] = self.model.addConstr(self.f[arc]<=yub*rhs, name = 'grotheyu%d^%d' % (i,j))

  def makeConstrs(self):
    # Create all constraints:
    self.capacity()
    self.consFlow()
    self.concentration()
    self.convexity()
    self.accumulation()
    self.validIneq()
    self.bilinear()
    if self.g.grothey:
      self.fracbound()

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
                  % (it, aCount, i, j, arc.h, flow))
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
    for (arc,f) in self.f.items():
      arc.flow = f.getAttr('x')
      arc.compflow = {}
      """
      if arc.flow > 0.05:
        print('Flow: %d %d %.1f' % (arc.fromnode.id, arc.tonode.id, arc.flow))
      """
    for ((arc,source),x) in self.x.items():
      arc.compflow[source] = x.getAttr('x')
    for node in self.g.nodes.values():
      node.yOpt = {}
    for ((node,source),y) in self.y.items():
      node.yOpt[source] = y.getAttr('x')

  def zeroSol(self):
    # puts the solution to zero
    for arc in self.g.arcs.values():
      arc.flow = 0.0
      arc.compflow = {}
      for source in arc.fromnode.sources.values():
        arc.compflow[source] = 0.0
      for sink in arc.tonode.sinks.values():
        arc.compflow[sink] = 0.0
    for node in self.g.nodes.values():
      node.yOpt = {}
      for source in node.sources.values():
        node.yOpt[source] = 0.0
      for sink in node.sinks.values():
        node.yOpt[sink] = 0.0

  def checkSol(self):
    # Check the computed solution for correctness
    incount = {}
    inflow = {}
    outflow = {}
    qual = {}
    (queue,h) = ([],0)
    for node in self.g.nodes.values():
      incount[node] = 0
      inflow[node] = 0.0
      outflow[node] = 0.0
      if node.isSource():
        qual[node] = node.q
        queue.append(node)
      else:
        qual[node] = [0.0]*self.g.nK
    cost = 0.0
    for arc in self.g.arcs.values():
      cost += arc.c * arc.flow
      if arc.flow > 0.001:
        if arc.flow > arc.h + 0.001:
          print('Failed capacity check at arc %d: %.2f > %.2f.' % (arc.fromnode.id, arc.tonodeid, arc.flow,arc.h))
        incount[arc.tonode] += 1
        inflow[arc.tonode] += arc.flow
        outflow[arc.fromnode] += arc.flow
    while h < len(queue):
      node = queue[h]
      h += 1
      for (tonode,arc) in node.outarcs.items():
        if arc.flow > 0.001:
          q = qual[tonode]
          for k in range(self.g.nK):
            q[k] += arc.flow * qual[node][k] / inflow[tonode]
          incount[arc.tonode] -= 1
          if incount[arc.tonode] == 0:
            queue.append(tonode)
    ok = True
    for source in self.g.sources.values():
      if outflow[source] > source.b + 0.001:
        print('Failed capacity check at source %d: %.2f > %.2f.' % (source.id, outflow[source], source.b))
        ok = False
      if outflow[source] < source.bL - 0.001:
        print('Failed capacity check at source %d: %.2f < %.2f.' % (source.id, outflow[source], source.bL))
        ok = False
    for pool in self.g.pools.values():
      if abs(inflow[pool]-outflow[pool]) > 0.001:
        print('Failed flow conservation check at pool %d: %.2f != %.2f.' % (pool.id, inflow[pool], outflow[pool]))
        ok = False
      if inflow[pool] > pool.b + 0.001:
        print('Failed capacity check at pool %d: %.2f > %.2f.' % (pool.id, inflow[pool], pool.b))
        ok = False
      if inflow[pool] < pool.bL - 0.001:
        print('Failed capacity check at pool %d: %.2f < %.2f.' % (pool.id, inflow[pool], pool.bL))
        ok = False
    for sink in self.g.sinks.values():
      if inflow[sink] > sink.b + 0.001:
        print('Failed capacity check at sink %d: %.2f > %.2f.' % (sink.id, inflow[sink], sink.b))
        ok = False
      if inflow[sink] < sink.bL - 0.001:
        print('Failed capacity check at sink %d: %.2f < %.2f.' % (sink.id, inflow[sink], sink.bL))
        ok = False
      for k in range(self.g.nK):
        if len(sink.q) and qual[sink][k]*inflow[sink] > sink.q[k]*inflow[sink] + 0.01:
          print('Failed quality check of component %d at sink %d: %.6f > %.6f (inflow: %.4f).' % (k, sink.id, qual[sink][k], sink.q[k], inflow[sink]))
          ok = False
        if len(sink.qL) and qual[sink][k]*inflow[sink] < sink.qL[k]*inflow[sink] - 0.01:
          print('Failed quality check of component %d at sink %d: %.6f < %.6f (inflow: %.4f).' % (k, sink.id, qual[sink][k], sink.qL[k], inflow[sink]))
          ok = False
    return (ok,cost)

  def solve(self,relax=False):
    self.makeModel(relax)
    if self.model is None:
      return
    """
    pre = self.model.presolve()
    pre.write('presolve.lp')
    for var in pre.getVars():
      if var.VType != GRB.CONTINUOUS:
        print('Integer variable: ', var.varName)
    print('We generated %d integer variables.' % self.model.NumIntVars)
    """
    self.model.optimize()
    if relax:
      print('gpp %s: RELAXED optimal objective function value: %.1f found in \t%.1f seconds.\n'
         % (self.instance, self.model.getObjective().getValue(), self.model.Runtime))
      self.getSol()
    elif self.model.SolCount:
      if status[self.model.Status] == 'OPT':
        print('gpp %s: Optimal objective function value: %.1f found in \t%.1f seconds.\n'
           % (self.instance, self.model.getObjective().getValue(), self.model.Runtime))
      elif status[self.model.Status] == 'TIME':
        print('gpp %s: Best objective function value: %.1f found in \t%.1f seconds.\n'
           % (self.instance, self.model.getObjective().getValue(), self.model.Runtime))
      self.getSol()
    else:
      self.zeroSol()
      print('gpp %s: No solution found in %.1f seconds.\n' % (self.instance, self.model.Runtime))

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
  """
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
  """
  if len(sys.argv) > 1:
    if len(sys.argv) > 2:
      gpp = GPP(sys.argv[1],sys.argv[2])
    else:
      gpp = GPP(sys.argv[1])
    gpp.solve()
    if gpp.model is not None:
      (ok,cost) = gpp.checkSol()
      if ok:
        print('Solution with profit %.2f ok.' % -cost)
      else:
        print('Solution with profit %.2f not ok.' % -cost)
    """
    if sys.argv[1][-4:] == '.dat':
      newfile = sys.argv[1][:-4] + '.txt'
      print(newfile)
      gpp.g.write(newfile)
    elif sys.argv[1][-5:] == '.json':
      newfile = sys.argv[1][:-5] + '.txt'
      print(newfile)
      gpp.g.write(newfile)
    elif sys.argv[1][-2:] == '.p':
      newfile = sys.argv[1][:-2] + '.txt'
      print(newfile)
      gpp.g.write(newfile)
    """
  else:  
    print('Usage: python3 gpp.py filename')
