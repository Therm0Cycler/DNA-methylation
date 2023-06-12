import seaborn as sns
import numpy as np
import matplotlib.pyplot as plt
import sys

fileOut=''
if len(sys.argv) == 4: fileOut=sys.argv[3]
elif len(sys.argv) != 3: sys.exit("* usage: plotState.py filename n_run [filename.pdf]")
fileName = sys.argv[1]
plotRun = int(sys.argv[2])

with open(fileName) as f:
    lines = f.read().splitlines()


nRuns = int(lines[-1].split()[0]) + 1
print("nRuns =",nRuns)

dyn = []
times = []

for l in lines:
    iRun = int(l.split()[0])    
    time = float(l.split()[1])
    state = l.split()[2]

    if iRun==plotRun:
      dyn.append( list(state)  )
      times.append(time)
                            
dyn = [list(map(int, i)) for i in dyn]

s = sns.heatmap( dyn, cmap="vlag", cbar=False  )
s.set_xlabel('locus')
s.set_ylabel(["Time: {}".format(times[-1])])
s.set_title(fileName)
if fileOut=='':
   plt.show()
else:
   plt.savefig(fileOut)
