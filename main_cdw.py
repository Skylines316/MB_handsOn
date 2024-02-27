from scipy.integrate import dblquad
import numpy as np
from math import sqrt

import matplotlib.pyplot as plt

def Delta(d, g):
  return dblquad(lambda x, y: np.power(2*np.pi,-2)*d*g/np.sqrt(d**2*g**2/4 + 4*(np.cos((x+y)/sqrt(2))+np.cos((x-y)/sqrt(2)))**2), -np.pi/sqrt(2), np.pi/sqrt(2), lambda x: -np.pi/sqrt(2), lambda x: np.pi/sqrt(2))

def recursive(d, g):
  d1 = Delta(d, g)[0]
  dif = d1 - d 
  if abs(dif) <= 1e-2:
    return d1
  return recursive(d1, g)

g = np.linspace(0.1, 5, 50)
d = []
with open('data_cdw.dat', 'w') as file:  
  for i in g:
    d.append(recursive(i+1, i))
    file.write(str(i)+","+str(d[-1])+"\n")
    print(i, d[-1],"and", Delta(d[-1],i))

plt.plot(g,d)
plt.show()

############################################################

i3, data3 = np.loadtxt("data_cdw.dat", delimiter=',', unpack=True)
i2, data2 = np.loadtxt("data2.dat", delimiter=',', unpack=True)

plt.plot(i3,data3, label='CDW')
plt.plot(i2,2*data2, label='SC')
plt.legend()
plt.grid(True)
plt.title("Non interacting case: comparison between SC and CDW")
plt.xlabel('U/t')
plt.ylabel(r'$\Delta$')
plt.show()