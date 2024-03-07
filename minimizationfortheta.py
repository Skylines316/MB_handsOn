import numpy as np
import matplotlib.pyplot as plt
import math
from mpl_toolkits.mplot3d import Axes3D
from scipy.optimize import minimize

def xi(kx, ky):
    return -2 * (np.cos(kx) + np.cos(ky))

def energy(params, *args):
    a = params
    U, J, kx, ky, Nk = args
    return 2 / Nk**2 * np.sum(- np.sqrt(np.power(xi(kx, ky), 2) + (U - J)**2*U**2/(np.cos(a)**2*(U  - J)**2 + U**2*np.sin(a)**2) )) +2 * (U*(U-J)*(U - J*np.cos(a)**2))/((U-J)**2*np.cos(a)**2 + U**2*np.sin(a)**2)

def energy1(params, *args):
    a = params
    U, J, kx, ky, Nk = args
    return 2 / Nk**2 * np.sum(- np.sqrt(np.power(xi(kx, ky), 2) + (J-U)**2*np.sin(a)**2 + U**2*np.cos(a)**2 )) +2*((U-J)*(np.sin(a))**2 +U*np.cos(a)**2)
Nk = 500
kx, ky = np.meshgrid(np.linspace(-math.pi, math.pi, Nk), np.linspace(-math.pi, math.pi, Nk))
J_value = 5.
U_value = 4
theta=np.linspace(0,np.pi/2,600)
en=np.zeros_like(theta)
for i in range(len(en)): 
    en[i]=energy1(theta[i], U_value, J_value, kx, ky, Nk)
plt.figure()
plt.plot(theta, en)
plt.xlabel(r'$\theta$')
plt.ylabel('Energy')
plt.title(f'Boundary energy, J={J_value}, U={U_value} ')
#plt.savefig(f"Plots/Minimisation/angle_J={J_value}_U={U_value}.svg")
plt.show()


Nk = 500
kx, ky = np.meshgrid(np.linspace(-math.pi, math.pi, Nk), np.linspace(-math.pi, math.pi, Nk))
theta=np.linspace(0,np.pi/2,600)
en=np.zeros_like(theta)
for i in range(len(en)): 
    en[i]=energy(theta[i], U_value, J_value, kx, ky, Nk)
plt.figure()
plt.plot(theta, en)
plt.xlabel(r'$\theta$')
plt.ylabel('Energy')
plt.title(f'Boundary energy, J={J_value}, U={U_value} ')
#plt.savefig(f"Plots/Minimisation/angle_J={J_value}_U={U_value}.svg")
plt.show()

guess = 1.172
bounds = [(0, np.pi/2)]
Nk_value = 700
kx, ky = np.meshgrid(np.linspace(-np.pi, np.pi, Nk_value), np.linspace(-np.pi, np.pi, Nk_value))
U_list = np.linspace(0.1, J_value, 100)
a_min = np.zeros(len(U_list))

for i in range(len(U_list)):
    result = minimize(energy, guess, args=(U_list[i], J_value, kx, ky, Nk_value), bounds=bounds, options={'xtol':1e-6, 'ftol':1e-6,})
    a_min[i] = result.x
    

np.savetxt(f'Data/minimisation/angle_dataJ={J_value}_U={U_list[0]}-{U_list[-1]}_{len(U_list)}.dat', np.c_[U_list,a_min], delimiter=',') 

plt.plot(U_list,a_min, label='minimum angle')
plt.legend()
plt.grid(True)
plt.title(f'minimum angle on the boundaries, U < J={J_value} ' )
plt.xlabel('U/t')
plt.ylabel(r'$\alpha$')
plt.title(f'J={J_value}')
plt.savefig(f"Plots/Minimisation/angle_J={J_value}_U={U_list[0]}-{U_list[-1]}_{len(U_list)}.svg")
plt.show()

en=np.zeros_like(U_list)
for i in range(len(en)): 
    en[i]=energy(a_min[i], U_list[i], J_value, kx, ky, Nk)
plt.plot(U_list,en)
plt.legend()
plt.grid(True)
plt.title(f'minimum energy on the boundaries, U < J={J_value} ' )
plt.xlabel('U/t')
plt.ylabel(r'Energy')
#plt.savefig(f"Plots/Minimisation/Boundary_energy_angle_J={J_value}_U={U_list[0]}-{U_list[-1]}_{len(U_list)}.svg")
plt.show()