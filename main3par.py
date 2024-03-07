from scipy.integrate import dblquad
import numpy as np

import matplotlib.pyplot as plt

import random as rand

def IntDelta(d, u, j, a, start, N, var):
  """ d: BCS parameter, u = U/t, j = J/t, a = CDW parameter"""
  acc = 0
  
  x1 = 0
  y1 = 0
    
  delta = lambda x,y: (2/(2*np.pi)**2)*d*(u-j)/np.sqrt(16*(np.cos(x) + np.cos(y))**2 + 4*(d**2)*(u-j)**2 + (a**2)*u**2)
  for i in range(N):
    rand1 = rand.random()
    rand2 = rand.random()
    rand3 = rand.random()
    x_try = ((x1 + var*(1 - 2*rand1))- np.pi)%(2*np.pi) + np.pi
    y_try = ((y1 + var*(1 - 2*rand2))- np.pi)%(2*np.pi) + np.pi
    z1 = delta(x_try,y_try)/delta(x1,y1)
    valInt = 0
    if(rand3 < z1):
      x1 = x_try
      y1 = y_try
      valInt += delta(x1,y1)
      acc +=1
    else:
      valInt += delta(x1,y1) 
  valInt = valInt/N
  return valInt, acc/N
    


def IntAlpha(d, u, j, a, start, N, var):
  """ d: BCS parameter, u = U/t, j = J/t, a = CDW parameter"""
  acc = 0
  x1 = 0
  y1 = 0
  alpha =   lambda x,y: (2/(2*np.pi)**2)*a*(u)/np.sqrt(16*(np.cos(x) + np.cos(y))**2 + 4*(d**2)*(u-j)**2 + (a**2)*u**2)
  for i in range(N):
    rand1 = rand.random()
    rand2 = rand.random()
    rand3 = rand.random()
    x_try = ((x1 + var*(1 - 2*rand1))- np.pi)%(2*np.pi) + np.pi
    y_try = ((y1 + var*(1 - 2*rand2))- np.pi)%(2*np.pi) + np.pi
    z1 = alpha(x_try,y_try)/alpha(x1,y1)
    valInt = 0
    if(rand3 < z1):
      x1 = x_try
      y1 = y_try
      valInt += alpha(x1,y1)
      acc +=1
    else:
      valInt += alpha(x1,y1) 
  valInt = valInt/N
  return valInt, acc/N       


def Delta(d, u, j, a):
  return dblquad(lambda x, y: 2*np.power(2*np.pi,-2)*d*(u-j)/np.sqrt(4*d**2*(u-j)**2 + a**2*u**2 + 16*(np.cos(x)+np.cos(y))**2), -np.pi, np.pi, lambda x: -np.pi, lambda x: np.pi)

def Delta_zero(d, u, j):
  return dblquad(lambda x, y: 1*np.power(2*np.pi,-2)*d*np.abs(u-j)/np.sqrt(4*d**2*(u-j)**2 + 16*(np.cos(x)+np.cos(y))**2), -np.pi, np.pi, lambda x: -np.pi, lambda x: np.pi)

def alpha_zero(u, a):
  return dblquad(lambda x, y: 1*np.power(2*np.pi,-2)*a*u/np.sqrt(a**2*u**2 + 16*(np.cos(x)+np.cos(y))**2), -np.pi, np.pi, lambda x: -np.pi, lambda x: np.pi)

def alpha(d, u, j, a):
  return dblquad(lambda x, y: 2*np.power(2*np.pi,-2)*a*u/np.sqrt(4*d**2*(u-j)**2 + a**2*u**2 + 16*(np.cos(x)+np.cos(y))**2), -np.pi, np.pi, lambda x: -np.pi, lambda x: np.pi)

def recursive(d, u, j, a):
  d1 = Delta_zero(d, u, j)[0]
  a1 = alpha_zero(u, a)[0]
  if a1<0:
    a1 = 0
  if d1<0:
    d1 = 0
  difd = d1 - d
  difa = a1 - a
  # print(d1[1])
  # print(difd, d1[0])
  if abs(difd) <= 1e-2 and abs(difa) <= 1e-2:
    en_Delta = -dblquad(lambda x, y: np.power(2*np.pi,-2)*np.sqrt(4*d1**2*(u-j)**2 + 16*(np.cos(x)+np.cos(y))**2), -np.pi, np.pi, lambda x: -np.pi, lambda x: np.pi)[0] - 2*(j-u)*d1**2
    
    en_alpha = -dblquad(lambda x, y: np.power(2*np.pi,-2)*np.sqrt( a1**2*u**2 + 16*(np.cos(x)+np.cos(y))**2), -np.pi, np.pi, lambda x: -np.pi, lambda x: np.pi)[0] + 0.5*u*a1**2

    en_base = -dblquad(lambda x, y: np.power(2*np.pi,-2)*np.sqrt(16*(np.cos(x)+np.cos(y))**2), -np.pi, np.pi, lambda x: -np.pi, lambda x: np.pi)[0]

    if en_Delta > en_alpha:
      if en_alpha > en_base:
        return 0, 0, en_Delta, en_alpha, en_base
      else:
        return 0, a1, en_Delta, en_alpha, en_base
    else:
      if en_Delta > en_base:
        return 0, 0, en_Delta, en_alpha, en_base
      else:
        return d1, 0, en_Delta, en_alpha, en_base     
  return recursive(d1, u, j, a1)

u_arr = np.linspace(0.01, 10, 400)
d = []
c = 1
J = 9
filename = 'dataJ='+str(J)+'_U=0.01-10_400_w_energies.dat'
with open(filename, 'w') as file: 
  for i in u_arr:
    d, a , eD, eA, eB= recursive(i+1, i, J, i+1)
    file.write(str(i)+","+str(a)+","+str(d)+","+str(eA)+","+str(eD)+","+str(eB)+"\n")
    print(c,"/",len(u_arr), round(a-alpha_zero(i, a)[0],4), round(d-Delta_zero(d, i, J)[0],4), "U="+str(i), "delta="+str(round(d,4)), "alpha="+str(round(a,4)))
    # print(d-IntDelta(d, i, 0.1, a, rand.random(), 300, 0.5)[0], a-IntAlpha(d, i, 0.1, a, rand.random(), 300, 0.5)[0] )
    c+=1

# plt.plot(u_arr,d)
# plt.show()
