from matplotlib.backends.backend_pdf import PdfPages
import seaborn as sns
import math
import numpy as np
import matplotlib.pyplot as plt
import sys

fileOut=''
LogData='LogData.dat'
MethLog='MethLog.dat'

if len(sys.argv) != 5: sys.exit("* usage: Analiser.py filename n_run filename.pdf filename_conts")
fileName = sys.argv[1]
plotRun = int(sys.argv[2])
fileOut=sys.argv[3]
file_conts = sys.argv[4]

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
rows = len(dyn)
nsomes = len(dyn[0])

#DataArr
DataArr = []


#Compute mean Methylation

Meth = []

for i in range(rows):
    temp = 0 
    for j in range(nsomes):
        temp += dyn[i][j]

    Meth.append((temp / nsomes))

Meth_tot = 0

for i in range(rows):
    Meth_tot += (Meth[i] / rows)

print("Mean Metylation: ", Meth_tot)

DataArr.append(Meth_tot)


#Compute persistences

pers = []
CEmpty = 0

for i in range(nsomes):
    temp = 0
    starter = 0
    prst = []
    for j in range(rows):
    	
        if j < (rows-1):
            if dyn[j][i]== 0 and not prst: 
                starter +=1
            if dyn[j][i] == 1 and dyn[j+1][i] == 1:
                temp += 1
	
            if dyn[j][i] == 1 and dyn[j+1][i] == 0:
                temp += 1
                prst.append(temp/(rows-starter))
                temp = 0
            else:
                None
        if j == (rows -1) and dyn[j][i] == 1:
            temp += 1
            if not prst: prst.append(1)
            prst.append(temp/(rows-starter))
            temp = 0
            break
        else:
            None

    if not prst: CEmpty += 1
    pers.append(prst)

mean_persistence = []

for i in range(nsomes):
    tmp = 0
    for j in range(len(pers[i])):
        tmp += pers[i][j]/len(pers[i])

    mean_persistence.append(tmp)

mean_perstot = 0.

for i in range(len(mean_persistence)):

    mean_perstot += mean_persistence[i]/(nsomes - CEmpty)

print("Mean Persistence: ", mean_perstot)

DataArr.append(mean_perstot)

#Compute Null-persistences

Nullpers = []
CNullEmpty = 0

for i in range(nsomes):
    temp = 0
    starter = 0
    Nprst = []
    for j in range(rows):
        if j < (rows-1):
            if dyn[j][i]== 1 and not Nprst: 
                starter +=1
            if dyn[j][i] == 0 and dyn[j+1][i] == 0:
                temp += 1
	
            if dyn[j][i] == 0 and dyn[j+1][i] == 1:
                temp += 1
                Nprst.append(temp/(rows-starter))
                temp = 0
            
            else:
                None
	
        if j == (rows -1) and dyn[j][i] == 0:
            temp += 1
            if not Nprst: Nprst.append(1)
            Nprst.append(temp/(rows-starter))
            temp = 0
            break
        else:
            None

    if not Nprst: CNullEmpty += 1
    Nullpers.append(Nprst)

mean_null_persistence = []

for i in range(nsomes):
    tmp = 0
    for j in range(len(Nullpers[i])):
        tmp += Nullpers[i][j]/len(Nullpers[i])

    mean_null_persistence.append(tmp)

mean_null_perstot = 0.

for i in range(len(mean_null_persistence)):

    mean_null_perstot += mean_null_persistence[i]/(nsomes - CNullEmpty)

print("Mean Unmeth Persistence: ", mean_null_perstot)

DataArr.append(mean_null_perstot)

#Compute blocks dimension

blocks = []

for i in range(rows):
    temp = 0
    blct = []
    for j in range(nsomes):
        if j < (nsomes-1):
            if dyn[i][j] == 1 and dyn[i][j]*dyn[i][j+1] == 1:
                temp += 1
	
            if dyn[i][j] == 1 and dyn[i][j]*dyn[i][j+1] == 0:
                temp += 1
                blct.append(temp)
                temp = 0
            else:
                None

        if j == (nsomes -1) and dyn[i][j] == 1:
            temp += 1
            blct.append(temp)
            temp = 0
            break
        else:
            None

    blocks.append(blct)

mean_dim = []

for i in range(rows):
    tmp = 0
    for j in range(len(blocks[i])):
        tmp += blocks[i][j]/len(blocks[i])

    mean_dim.append(tmp)
    
#Contacter

with open(file_conts) as fc:
    lines_c = fc.read().splitlines()

length_c = len(lines_c)
Contacts = []
times_c = []
for l in lines_c:

    Contacts.append(int(l.split()[2]))
    times_c.append(float(l.split()[1]))
    
Mean_c = 0
temp_c = 0 
for i in range(length_c):
    
    temp_c += Contacts[i]

Mean_c = temp_c / length_c


print("Mean number of contacts: ", Mean_c )


#PrintLog
s0 = fileName.split("(")
s1 = s0[1].split(")")
sval = s1[0].split("-")

with open("LogData.dat", "a") as text_file:
    text_file.write(str(sval[0]) + " " +str(sval[1]) + " " + str(sval[2]) + " " + str(sval[3]) + " " + str(sval[4]) + " " + str(sval[5]) + " " + str(sval[6]) + " "  +   str(Meth_tot) + " " +  str(mean_perstot) + " " +  str(mean_null_perstot) + "\n")

with open("MethLog.dat", "a") as text_file:
    text_file.write(str(Meth) + "\n")

#Plotting
import matplotlib.pyplot as plt

#plt.figure(1)	
#s = sns.heatmap( dyn, cmap="vlag", cbar=False  )
#s.set_xlabel('locus')
#s.set_ylabel('time')
#s.set_title(fileName)

plt.figure(2)
plt.plot(times,Meth)
plt.ylabel('Methylation rate')
plt.legend(["Mean methylation: {}".format(Meth_tot)])
plt.xlabel('time')
#plt.savefig('Methylation.pdf')

plt.figure(3)
plt.plot(times,mean_dim)
plt.xlabel('time')
plt.ylabel('Mean blocks dimension')
#plt.savefig('Blocks.pdf')	

plt.figure(4)
plt.plot(mean_persistence)
plt.ylabel('Mean Persistence')
plt.xlabel('Nucleosome')
plt.legend(["Mean persistence: {}".format(mean_perstot)])
#plt.savefig('Persistence.pdf')

plt.figure(5)
plt.plot(mean_null_persistence)
plt.ylabel('Mean Unmeth Persistence')
plt.xlabel('Nucleosome')
plt.legend(["Mean unmeth persistence: {}".format(mean_null_perstot)])


plt.figure(6)
bins = np.linspace(0, 20, 50)
plt.hist(blocks[math.floor(rows/4)], bins, alpha=0.5, label='1/4 time')
plt.hist(blocks[math.floor(rows/2)], bins, alpha=0.5, label='1/2 time')
plt.hist(blocks[math.floor(rows/1.33)], bins, alpha=0.5, label='3/4 time')
plt.hist(blocks[rows-1], bins, alpha=0.5, label='end time')
plt.legend(loc='upper right')
#plt.savefig('Blocks-distribution.pdf')

plt.figure(7)
plt.plot(times_c,Contacts)
plt.ylabel('Number of contacts')
plt.legend(["Mean number of contacts: {}".format(Mean_c)])
plt.xlabel('time')
plt.tight_layout()

def save_image(filenom):
    

    p = PdfPages(filenom)
      
    fig_nums = plt.get_fignums()  
    figs = [plt.figure(n) for n in fig_nums]

    for fig in figs: 

        fig.savefig(p, format='pdf') 

    p.close()  
   
  
save_image(fileOut)


