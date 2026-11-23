##
# Minimum cost flow problem instance
#
from node import Node, Source, Pool, Sink
from arc import Arc
from json import load
import pickle

class Graph:
  def __init__(self, filename):
    self.id = filename
    self.nodes = {}
    self.sources = {}
    self.pools = {}
    self.sinks = {}
    self.arcs = {}
    self.nK = 0
    self.grothey = 'Grothey' in filename
    if filename[-4:] == '.txt':
      self.read(filename)
    elif 'Schweiger' in filename:
      self.readSchweiger(filename)
    elif 'Alfaki' in filename or 'Gestsson' in filename:
      self.readAlfaki(filename)
    elif 'Ceccon' in filename:
      self.readCeccon(filename)
    elif 'Dey' in filename:
      self.readDey(filename)
    elif 'Grothey' in filename:
      self.readGrothey(filename)
    else:
      self.read(filename)
    for arc in self.arcs.values():
      arc.updateCap()

  def read(self, filename):
    """
    Reads instance data from a text file with the simple format
    """
    reader = open(filename,'r')
    reader.readline()
    string = reader.readline().strip()
    data = string.split()
    self.nS = int(data[0].strip())
    self.nP = int(data[1].strip())
    self.nT = int(data[2].strip())
    self.n = self.nS + self.nP + self.nT
    self.nA = int(data[3].strip())
    self.nK = int(data[4].strip())
    for i in range(self.nS):
      self.sources[i] = Source(self, i)
    for i in range(self.nS,self.nS+self.nP):
      self.pools[i] = Pool(self, i)
    for i in range(self.nS+self.nP, self.n):
      self.sinks[i] = Sink(self, i)
    self.nodes = self.sources | self.pools | self.sinks
    # print('Reading %d arcs.' % self.nA)
    reader.readline()
    for k in range(self.nA):
      string = reader.readline().strip()
      data = string.split()
      i = int(data[0].strip())
      j = int(data[1].strip())
      c = float(data[2].strip())
      h = None
      if len(data) > 3:
        h = float(data[3].strip())
        if len(data) > 4:
          (ylb,yub) = (float(data[4].strip()), 1.0)
          if len(data) > 5:
            yub = float(data[5].strip())
          self.nodes[j].ylb[self.nodes[i]] = ylb
          self.nodes[j].yub[self.nodes[i]] = yub
      arc = Arc(self, self.nodes[i], self.nodes[j], c, h)
      self.arcs[i,j] = arc
      self.nodes[i].addOut(arc)
      self.nodes[j].addIn(arc)
    reader.readline()
    for i in range(self.n):
      string = reader.readline().strip()
      data = string.split()
      i = int(data[0].strip())
      b = float(data[1].strip())
      self.nodes[i].b = b
      if len(data) > 2:
        bL = float(data[2].strip())
        self.nodes[i].bL = bL
    reader.readline()
    for i in range(self.nS):
      string = reader.readline().strip()
      data = string.split()
      i = int(data[0].strip())
      for k in range(self.nK):
        q = float(data[k+1].strip())
        self.sources[i].q.append(q)
    for i in range(self.nS+self.nP, self.n):
      string = reader.readline().strip()
      data = string.split()
      i = int(data[0].strip())
      for k in range(self.nK):
        q = float(data[k+1].strip())
        self.sinks[i].q.append(q)
      if len(data) > 1+self.nK:
        for k in range(self.nK):
          qL = float(data[k+1+self.nK].strip())
          self.sinks[i].qL.append(qL)
    reader.close()
    # print('Created graph with %d sources, %d pools, %d sinks, %d arcs, and %d flow components.'
    #       % (self.nS, self.nP, self.nT, self.nA, self.nK))
    self.connectSources()
    self.connectSinks()

  def readArcs(self, reader):
    arcs = []
    string = reader.readline().strip()
    while '.' not in string:
      string = reader.readline().strip()
    while '.' in string:
      string.strip(',')
      string.strip('/;')
      newarc = string.split('.')
      arcs.append((newarc[0].strip(' ,'), newarc[1].strip(' ,')))
      string = reader.readline().strip()
    return (arcs,'A')

  def readSet(self, reader):
    string = reader.readline().strip()
    while len(string)==0 or string[0]=='*':
      string = reader.readline().strip()
    while 'set ' not in string or '/' not in string:
      string = reader.readline().strip()
    pos = string.find('set ')
    string = string[pos+len('set '):].strip()
    t = string[0]
    if t=='A':
      return self.readArcs(reader)
    pos = string.find('/')
    string = string[pos+1:]
    pos = string.find('/')
    elements = []
    for entry in string[:pos].split(','):
      elements.append(entry.strip(' ,'))
    return (elements,t)

  def readParam(self, reader):
    string = reader.readline().strip()
    while len(string)==0 or string[0]=='*':
      string = reader.readline().strip()
    while 'parameter ' not in string or '/' not in string:
      string = reader.readline().strip()
    pos = string.find('parameter ')
    string = string[pos+len('parameter '):].strip()
    pos = string.find('(')
    param = string[:pos]
    values = {}
    pos = string.find('/')
    string = string[pos+1:].strip()
    if len(string) == 0:
      string = reader.readline().strip()
    done = False
    while not done:
      if '/' in string:
        done = True
        pos = string.find('/')
        string = string[:pos].strip()
      if len(string):
        entries = string.split()
        values[entries[0].strip()] = entries[1].strip()
        string = reader.readline().strip()
    return (values, param)

  def readSchweiger(self, filename):
    """
    Reads instance data from a text file with Schweiger's gms-format
    """
    sets = {'V', 'K', 'I', 'L', 'J', 'A'}
    lists = {}
    reader = open(filename,'r')
    for i in range(len(sets)):
      (elements, t) = self.readSet(reader)
      lists[t] = elements
    self.nS = len(lists['I'])
    self.nP = len(lists['L'])
    self.nT = len(lists['J'])
    self.n = self.nS + self.nP + self.nT
    self.nA = len(lists['A'])
    self.nK = len(lists['K'])
    # print(self.nS, self.nP, self.nT, self.n, self.nA, self.nK) 
    nodelist = lists['I'] + lists['L'] + lists['J']
    for i in range(self.nS):
      self.sources[i] = Source(self, i)
      self.sources[i].q = [None]*self.nK
    for i in range(self.nS,self.nS+self.nP):
      self.pools[i] = Pool(self, i)
    for i in range(self.nS+self.nP, self.n):
      self.sinks[i] = Sink(self, i)
      self.sinks[i].q = [None]*self.nK
    self.nodes = self.sources | self.pools | self.sinks
    for k in range(self.nA):
      (iname,jname) = (lists['A'][k][0].strip(' ,'), lists['A'][k][1].strip(','))
      (i,j) = (nodelist.index(iname), nodelist.index(jname))
      arc = Arc(self, self.nodes[i], self.nodes[j],None)
      self.arcs[i,j] = arc
    params = {'cost', 'C', 'lambda', 'overbeta'}
    tables = {}
    for i in range(len(params)):
      (values, param) = self.readParam(reader)
      tables[param] = values
    data = tables['cost']
    for (a,c) in data.items():
      names = a.strip().split('.')
      (i,j) = (nodelist.index(names[0].strip(' ,')), nodelist.index(names[1].strip(' ,')))
      self.arcs[i,j].c = float(c)
    data = tables['C']
    for (i,b) in data.items():
      iname = i.strip()
      i = nodelist.index(iname)
      self.nodes[i].b = float(b)
    data = tables['lambda']
    for (ki,q) in data.items():
      (kname, iname) = ki.strip().split('.')
      (k,i) = (lists['K'].index(kname), nodelist.index(iname))
      self.sources[i].q[k] = float(q)
    data = tables['overbeta']
    for (ki,q) in data.items():
      (kname, iname) = ki.strip().split('.')
      (k,i) = (lists['K'].index(kname), nodelist.index(iname))
      self.sinks[i].q[k] = float(q)
    reader.close()
    # print('Created graph with %d sources, %d pools, %d sinks, %d arcs, and %d flow components.'
    #       % (self.nS, self.nP, self.nT, self.nA, self.nK))
    self.connectSources()
    self.connectSinks()

  def getNumbers(self, string):
    entries = string.split()
    numbers = []
    for entry in entries:
      numbers.append(float(entry.strip(' ;')))
    return numbers

  def readTable(self, reader):
    string = reader.readline().strip()
    while 'table ' not in string or '(' not in string or ')' not in string:
      string = reader.readline().strip()
    pos = string.find('table ') + len('table ')
    t = string[pos]
    string = reader.readline().strip()
    cols = []
    entries = string.split()
    for entry in entries:
      cols.append(int(entry.strip())-1)
    values = []
    rows = []
    done = False
    while not done:
      string = reader.readline().strip()
      done = ';' in string
      string = string.strip(' ;')
      if len(string):
        row = self.getNumbers(string)
        rows.append(int(row[0])-1)
        values.append(row[1:])
      else:
        done = True
    return (rows, cols, values, t)

  def getRange(self, string):
    if '*' in string:
      entries = string.split('*')
      return int(entries[1].strip()) - int(entries[0].strip()) + 1
    entries = string.split()
    (low,up) = (int(entries[0].strip()), int(entries[0].strip()))
    for i in range(1, len(entries)):
      (low, up) = (min(low, int(entries[i].strip())), max(up, int(entries[i].strip())))
    return up - low + 1

    
  def getCost(self, row, col, table):
    rows = table[0]
    cols = table[1]
    values = table[2]
    i = rows.index(row)
    j = cols.index(col)
    return values[i][j]

  def readAlfaki(self, filename):
    """
    Reads instance data from a text file with Alfaki's gms-format
    """
    sets = {'i', 's', 't', 'k'}
    lists = {}
    reader = open(filename,'r')
    for i in range(len(sets)):
      (elements, t) = self.readSet(reader)
      lists[t] = elements
    self.nS = self.getRange(lists['s'][0])
    self.nT = self.getRange(lists['t'][0])
    self.n = self.getRange(lists['i'][0])
    self.nP = self.n - self.nS - self.nT
    self.nK = self.getRange(lists['k'][0])
    # print(self.nS, self.nP, self.nT, self.n, self.nK) 
    tableNames = {'c', 'a', 'q'}
    tables = {}
    for i in range(len(tableNames)):
      (rows, cols, values, t) = self.readTable(reader)
      tables[t] = (rows, cols, values)
    paramNames = {'bl(i)', 'bu(i)'}
    for i in range(len(paramNames)):
      (values, param) = self.readParam(reader)
      tables[param] = values
    data = tables['bu']
    b = {}
    for (i,value) in data.items():
      b[int(i.strip())-1] = float(value.strip())
    q = {}
    (rows,cols,values) = tables['q']
    for i in range(len(rows)):
      q[rows[i]] = values[i]
    for i in range(self.nS):
      self.sources[i] = Source(self, i, 0.0, b[i], None, q[i])
    for i in range(self.nS,self.nS+self.nP):
      self.pools[i] = Pool(self, i, 0.0, b[i])
    for i in range(self.nS+self.nP, self.n):
      self.sinks[i] = Sink(self, i, 0.0, b[i], None, q[i])
    self.nodes = self.sources | self.pools | self.sinks
    (rows, cols, values) = tables['a']
    for i in range(len(values)):
      for j in range(len(values[i])):
        if values[i][j] > 0.5:
          c = self.getCost(rows[i], cols[j], tables['c'])
          arc = Arc(self, self.nodes[rows[i]], self.nodes[cols[j]], c)
          self.arcs[rows[i],cols[j]] = arc
          self.nodes[rows[i]].addOut(arc)
          self.nodes[cols[j]].addIn(arc)
    self.nA = len(self.arcs)
    # print(self.nS, self.nP, self.nT, self.n, self.nA, self.nK) 
    reader.close()
    # print('Created graph with %d sources, %d pools, %d sinks, %d arcs, and %d flow components.'
    #       % (self.nS, self.nP, self.nT, self.nA, self.nK))
    self.connectSources()
    self.connectSinks()

  def readCeccon(self, filename):
    """
    Reads instance data from a text file with Ceccon's json-format
    """
    file = open(filename, 'r')
    data = load(file)
    file.close()
    self.nS = len(data['components'])
    self.nP = len(data['pool_size'])
    self.nT = len(data['products'])
    self.n = self.nS + self.nP + self.nT
    self.nK = len(data['components'][0]['quality'])
    ids = {}
    for s in range(self.nS):
      ids[data['components'][s]['name']] = s
      b = data['components'][s]['upper']
      q = list(data['components'][s]['quality'].values())
      self.sources[s] = Source(self, s, 0.0, b, q)
    p = self.nS
    for (name,b) in data['pool_size'].items():
      ids[name] = p
      self.pools[p] = Pool(self, p, 0.0, b)
      p += 1
    for t in range(self.nT):
      sinkId = self.nS+self.nP+t
      ids[data['products'][t]['name']] = sinkId
      b = data['products'][t]['upper']
      q = list(data['products'][t]['quality_upper'].values())
      qL = None
      if data['products'][t]['quality_lower'] is not None:
        qL = list(data['products'][t]['quality_lower'].values())
      self.sinks[sinkId] = Sink(self, sinkId, 0.0, b, qL, q)
    self.nodes = self.sources | self.pools | self.sinks
    for value in data['component_to_pool_fraction']:
      (i,j,c) = (ids[value['component']], ids[value['pool']], value.get('cost',0))
      c += data['components'][i].get('price',0)
      arc = Arc(self, self.nodes[i], self.nodes[j], c)
      self.arcs[i,j] = arc
      self.nodes[i].addOut(arc)
      self.nodes[j].addIn(arc)
    for value in data['component_to_product_bound']:
      (i,j,c,h) = (ids[value['component']], ids[value['product']], value.get('cost',0), value['bound'])
      c += data['components'][i].get('price',0)
      c -= data['products'][j-self.nS-self.nP].get('price',0)
      arc = Arc(self, self.nodes[i], self.nodes[j], c, h)
      self.arcs[i,j] = arc
      self.nodes[i].addOut(arc)
      self.nodes[j].addIn(arc)
      if h < min(self.nodes[i].b, self.nodes[j].b):
        print('Warning: Small arc capacity on file', filename, i, j, self.nodes[i].b, self.nodes[j].b, h)
    for value in data['pool_to_product_bound']:
      (i,j,c,h) = (ids[value['pool']], ids[value['product']], value.get('cost',0), value['bound'])
      c -= data['products'][j-self.nS-self.nP].get('price',0)
      arc = Arc(self, self.nodes[i], self.nodes[j], c, h)
      self.arcs[i,j] = arc
      self.nodes[i].addOut(arc)
      self.nodes[j].addIn(arc)
      if h < min(self.nodes[i].b, self.nodes[j].b):
        print('Warning: Small arc capacity on file', filename, i, j, self.nodes[i].b, self.nodes[j].b, h)
    self.nA = len(self.arcs)
    self.connectSources()
    self.connectSinks()

  def readDey(self, filename):
    """
    Reads instance data from a text file with Dey's pickle-format
    """
    data = pickle.load(open(filename, 'rb'))
    (self.nS, self.nP, self.nT, self.nK) = (len(data['I']), len(data['L']), len(data['J']), len(data['K']))
    matrix = data['adj_matrix']
    for s in range(self.nS):
      bL = data['node_lb'][s]
      b = data['node_ub'][s]
      q = [None]*self.nK
      for k in range(self.nK):
        q[k] = data['speclevel'][s,k]
      self.sources[s] = Source(self, s, bL, b, q)
    for p in range(self.nP):
      bL = data['node_lb'][p]
      b = data['node_ub'][p]
      self.pools[self.nS+p] = Pool(self, self.nS+p, bL, b)
    for t in range(self.nT):
      bL = data['node_lb'][t]
      b = data['node_ub'][t]
      q = [None]*self.nK
      for k in range(self.nK):
        q[k] = data['maxspec'][self.nS+self.nP+t,k]
      self.sinks[self.nS+self.nP+t] = Sink(self, self.nS+self.nP+t, bL, b, None, q)
    self.nodes = self.sources | self.pools | self.sinks
    for (i,j) in matrix:
      if matrix[i,j]:
        arc = Arc(self, self.nodes[i], self.nodes[j], data['arc_cost'][i,j])
        self.arcs[i,j] = arc
        self.nodes[i].addOut(arc)
        self.nodes[j].addIn(arc)
    self.nA = len(self.arcs)
    self.connectSources()
    self.connectSinks()

  def skip(reader, term):
    data = reader.readline()
    while term not in data:
      data = reader.readline()

  def readFraction(self, reader, torange, lower=False):
    """
    Reads bounds on the entering flow fractions
    """
    for j in torange:
      tonode = self.nodes[j]
      ybound = tonode.yub
      if lower:
        ybound = tonode.ylb
      data = reader.readline().split()
      for i in range(self.nS+self.nP):
        fromnode = self.nodes[i]
        ybound[fromnode] = float(data[i+1])
        if ybound[fromnode] > 0.001 and (i,j) not in self.arcs:
          if fromnode.isPool() and tonode.isPool():
            print('Paper contradiction: pool-to-pool arc (%d,%d).' % (i,j))
          arc = Arc(self, fromnode, tonode, 0.0)
          self.arcs[i,j] = arc
          fromnode.addOut(arc)
          tonode.addIn(arc)
        if ybound[fromnode] > 0.001 and (ybound[fromnode] < 0.999 or lower) and i >= self.nS:
          print('Warning: Unexpected fraction bound (%d,%d): %.4f' % (i,j,ybound[fromnode]))

  def readGrothey(self, filename):
    """
    Reads instance data from a text file with Grothey's AMPL-format
    """
    reader = open(filename,'r')
    self.nK = len(reader.readline().split())-3
    self.nS = len(reader.readline().split())-3
    self.nP = len(reader.readline().split())-3
    self.nT = len(reader.readline().split())-3
    self.n = self.nS + self.nP + self.nT
    for i in range(self.nS):
      self.sources[i] = Source(self, i)
    for i in range(self.nS,self.nS+self.nP):
      self.pools[i] = Pool(self, i)
    for i in range(self.nS+self.nP, self.n):
      self.sinks[i] = Sink(self, i)
    self.nodes = self.sources | self.pools | self.sinks
    Graph.skip(reader, 'param')
    for k in range(self.nK):
      data = reader.readline().split()
      for s in range(self.nS):
        q = float(data[s+1])
        self.sources[s].q.append(q)
    Graph.skip(reader, 'param')
    for t in range(self.nT):
      data = reader.readline().split()
      for k in range(self.nK):
        q = float(data[k+1])
        self.sinks[self.nS+self.nP+t].q.append(q)
    Graph.skip(reader, 'param')
    for t in range(self.nT):
      data = reader.readline().split()
      for k in range(self.nK):
        q = float(data[k+1])
        self.sinks[self.nS+self.nP+t].qL.append(q)
    Graph.skip(reader, 'param')
    bTot = 0
    for t in range(self.nT):
      data = reader.readline().split()
      b = float(data[1])
      self.sinks[self.nS+self.nP+t].b = b
      self.sinks[self.nS+self.nP+t].bL = b
      bTot += b
    for i in range(self.nS+self.nP):
      self.nodes[i].b = bTot
    Graph.skip(reader, 'param')
    cost = {}
    for s in range(self.nS):
      data = reader.readline().split()
      cost[s] = float(data[1])
    Graph.skip(reader, 'param')
    self.readFraction(reader, range(self.nS, self.nS+self.nP), True)
    Graph.skip(reader, 'param')
    self.readFraction(reader, range(self.nS, self.nS+self.nP), False)
    Graph.skip(reader, 'param')
    self.readFraction(reader, range(self.nS+self.nP, self.nS+self.nP+self.nT), True)
    Graph.skip(reader, 'param')
    self.readFraction(reader, range(self.nS+self.nP, self.nS+self.nP+self.nT), False)
    for ((i,j), arc) in self.arcs.items():
      if i < self.nS:
        arc.c = cost[i]
    self.nA = len(self.arcs)
    self.connectSources()
    self.connectSinks()

  def write(self, filename):
    """
    Writes instance data to a file on compact format
    """
    writer = open(filename,'w')
    writer.write('Number of sources, pools, sinks, arcs, and flow components (|S|, |P|, |T|, |A|, |K|):\n')
    writer.write('%d %d %d %d %d\n' % (self.nS, self.nP, self.nT, self.nA, self.nK))
    writer.write('Arcs (i, j, c_{ij}[, h_{ij}]):\n')
    for ((i,j),arc) in self.arcs.items():
      if arc.fromnode.isSource() and arc.fromnode in arc.tonode.ylb | arc.tonode.yub:
        ylb = arc.tonode.ylb.get(arc.fromnode, 0.0)
        yub = arc.tonode.yub.get(arc.fromnode, 1.0)
        writer.write('%d %d %.6f %.6f %.6f %.6f\n' % (i, j, arc.c, arc.h, ylb, yub))
      elif arc.h < min(arc.fromnode.b, arc.tonode.b):
        writer.write('%d %d %.6f %.6f\n' % (i, j, arc.c, arc.h))
      else:
        writer.write('%d %d %.6f\n' % (i, j, arc.c))
    writer.write('Node capacities (i, b_i):\n')
    for node in self.nodes.values():
      if node.bL > 0.0005:
        writer.write('%d %.6f %.6f\n' % (node.id, node.b, node.bL))
      else:
        writer.write('%d %.6f\n' % (node.id, node.b))
    writer.write('Source concentrations and sink concentration upper (+lower) bounds (i=s/t, [q_i^k: k=0,...,|K|-1]):\n')
    for source in self.sources.values():
      string = '%d' % source.id
      for k in range(self.nK):
        string += ' %.6f' % source.q[k]
      writer.write(string+'\n')
    for sink in self.sinks.values():
      string = '%d' % sink.id
      for k in range(self.nK):
        string += ' %.6f' % sink.q[k]
      if len(sink.qL):
        for k in range(self.nK):
          string += ' %.6f' % sink.qL[k]
      writer.write(string+'\n')
    writer.close()

  def connectSources(self): 
    nconnections = 0
    for s in self.sources.values():
      nconnections += self.connectSource(s)
    # print('Average %.1f source connections per node' % (nconnections / len(self.nodes)))

  def connectSinks(self): 
    nconnections = 0
    for t in self.sinks.values():
      nconnections += self.connectSink(t)
    # print('Average %.1f sink connections per node' % (nconnections / len(self.nodes)))

  def connectSource(self, s): 
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

  def connectSink(self, t): 
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

