import subprocess
import matplotlib.pyplot as plt

# Creating an empty dictionary
myDict = {}
words = ['ST', 'CG', 'LB']

str = "head -n -3  | awk 'NR >= 10 {print $4}'"
for i in words:
    word = 'iterate_' + i + '.dat'
    process = str[:11] + word + str[11:]
    proc = subprocess.Popen(process, shell=True,
                            stdin=subprocess.PIPE,
                            stdout=subprocess.PIPE,
                            stderr=subprocess.PIPE)
    stdout_value = proc.communicate()[0]
    lst = [item for item in stdout_value.decode().split('\n')]
    del lst[-1]
    myDict[i] = [float(i) for i in lst]

plt.figure()
plt.plot(myDict['ST'], fillstyle='none', linestyle='-', marker='x',
         color='k', markersize=5, linewidth=1, label='Steepest descent')
plt.plot(myDict['CG'], fillstyle='none', linestyle='-', marker='^',
         color='k', markersize=5, linewidth=1, label='Nonlinear conjugate gradient')
plt.plot(myDict['LB'], fillstyle='none', linestyle='-', marker='s',
         color='k', markersize=5, linewidth=1, label='l-BFGS')
plt.xlabel('Iterations')
plt.ylabel('Relative objective function (log scale)')
plt.yscale('log')
plt.grid(linestyle=':', linewidth=0.5)
plt.legend(loc='upper right', frameon=True, prop={'size': 10})
plt.savefig('convergence_curves.svg')
#
nwDict = {}
str = "head -n -3  | awk 'NR >= 10 {print $(NF)}'"
for i in words:
    word = 'iterate_' + i + '.dat'
    process = str[:11] + word + str[11:]
    proc = subprocess.Popen(process, shell=True,
                            stdin=subprocess.PIPE,
                            stdout=subprocess.PIPE,
                            stderr=subprocess.PIPE)
    stdout_value = proc.communicate()[0]
    lst = [item for item in stdout_value.decode().split('\n')]
    del lst[-1]
    nwDict[i] = [float(i) for i in lst]

plt.figure()
plt.plot(nwDict['ST'], myDict['ST'], fillstyle='none', linestyle='-', marker='x',
         color='k', markersize=5, linewidth=1, label='Steepest descent')
plt.plot(nwDict['CG'], myDict['CG'], fillstyle='none', linestyle='-', marker='^',
         color='k', markersize=5, linewidth=1, label='Nonlinear conjugate gradient')
plt.plot(nwDict['LB'], myDict['LB'], fillstyle='none', linestyle='-', marker='s',
         color='k', markersize=5, linewidth=1, label='l-BFGS')
plt.xlabel('Computed gradients')
plt.ylabel('Relative objective function (log scale)')
plt.yscale('log')
plt.grid(linestyle=':', linewidth=0.5)
plt.legend(loc='upper right', frameon=True, prop={'size': 10})
plt.savefig('computationalcost_curves.svg')
