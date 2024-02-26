from scipy.integrate import dblquad
import numpy as np

import matplotlib.pyplot as plt

def Delta(d, u, j, a):
  return dblquad(lambda x, y: np.power(2*np.pi,-2)*d*(u-j)/np.sqrt(4*d**2*(u-j)**2 + a**2*u**2 + 16*(np.cos(x)+np.cos(y))**2), -np.pi, np.pi, lambda x: -np.pi, lambda x: np.pi)

def alpha(d, u, j, a):
  return dblquad(lambda x, y: np.power(2*np.pi,-2)*d*u/np.sqrt(4*d**2*(u-j)**2 + a**2*u**2 + 16*(np.cos(x)+np.cos(y))**2), -np.pi, np.pi, lambda x: -np.pi, lambda x: np.pi)

def recursive(d, u, j, a):
  d1 = Delta(d, u, j, a)[0]
  a1 = alpha(d, u, j, a)[0]
  difd = d1 - d 
  difa = a1 - a
  if abs(difd) <= 1e-2 and abs(difa) <= 1e-2:
    return d1, a1
  return recursive(d1, u, j, a1)

u_arr = np.linspace(0.01, 15, 60)
d = []
c = 1
with open('data3PnJ0p5.dat', 'w') as file: 
  for i in u_arr:
    d, a = recursive(i+1, i, 0.1, i+1)
    file.write(str(i)+","+str(a)+","+str(d)+"\n")
    print(c,"/",len(u_arr), d-Delta(d, i, 0.1, a)[0], a-alpha(d, i, 0.1, a)[0])
    c+=1

# plt.plot(u_arr,d)
# plt.show()
