##
# Reads data file names, writes a jobfile running all experiments
#

import sys

def makejob(datafilenames, jobfilename):
  """
  Reads instance data from a text file with the simple format
  """
  reader = open(datafilenames,'r')
  filelist = reader.readlines()
  reader.close()
  writer = open(jobfilename,'w')
  for string in filelist:
    data = string.split(' ')
    file = data[-1].strip()
    writer.write('python3 gpp.py %s | grep ": Optimal"\n' % file)
    writer.write('python3 wcigpp.py %s | grep ": Optimal"\n' % file)

if __name__ == '__main__':
  if len(sys.argv) > 2:
    makejob(sys.argv[1],sys.argv[2])
  else:
    print('Usage: makejob datafilenames jobfilename')
