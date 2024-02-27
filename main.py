from scipy.integrate import dblquad
import numpy as np

import matplotlib.pyplot as plt

def Delta(d, g):
  return dblquad(lambda x, y: np.power(2*np.pi,-2)*0.5*d*g/np.sqrt(d**2*g**2 + 4*(np.cos(x)+np.cos(y))**2), -np.pi, np.pi, lambda x: -np.pi, lambda x: np.pi)

def recursive(d, g):
  d1 = Delta(d, g)[0]
  dif = d1 - d 
  if abs(dif) <= 1e-2:
    return d1
  return recursive(d1, g)

g = np.linspace(0.1, 5, 50)
d = []
with open('data2.dat', 'w') as file:  
  for i in g:
    d.append(recursive(i+1, i))
    file.write(str(i)+","+str(d[-1])+"\n")
    print(i, d[-1],"and", Delta(d[-1],i))

plt.plot(g,d)
plt.show()
