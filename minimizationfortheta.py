import numpy as np
import matplotlib.pyplot as plt
import math
from mpl_toolkits.mplot3d import Axes3D
from scipy.optimize import minimize_scalar

def xi(kx, ky):
    return -2 * (np.cos(kx) + np.cos(ky))

#def energy(params, *args):
    a = params
    U, J, kx, ky, Nk = args
    return 2 / Nk**2 * np.sum(- np.sqrt(np.power(xi(kx, ky), 2) + (U - J)**2*U**2/(np.cos(a)**2*(U  - J)**2 + U**2*np.sin(a)**2) )) +2 * (U*(U-J)*(U - J*np.cos(a)**2))/((U-J)**2*np.cos(a)**2 + U**2*np.sin(a)**2)
def energy(params, *args):
    Alpha = params
    U, kx, ky, Nk = args
    return 2 / Nk**2 * np.sum(- np.sqrt(np.power(xi(kx, ky), 2)  + U**2 * np.power(Alpha, 2))) + 2 * ( U * ( np.power(Alpha, 2)))

def energy1(params, *args):
    a = params
    U, J, kx, ky, Nk = args
    return 2 / Nk**2 * np.sum(- np.sqrt(np.power(xi(kx, ky), 2) + (U - J)**2 * np.power(np.sin(a), 2) + U**2 * np.power(np.cos(a), 2))) + 2 * ((U - J) * np.power(np.sin(a), 2) + U * (np.power(np.cos(a), 2)))



Nk = 500
kx, ky = np.meshgrid(np.linspace(-math.pi, math.pi, Nk), np.linspace(-math.pi, math.pi, Nk))
J_value = 5.
U_value = 4.9
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

'''
J_value = 5
Nk = 500
kx, ky = np.meshgrid(np.linspace(-math.pi, math.pi, Nk), np.linspace(-math.pi, math.pi, Nk))
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
'''
############################################
J_list = [0.0,1.0,2.0,3.0,4.0,5.0,6.0,7.0]
guess = 1
bounds1 = (0, np.pi/2)
bounds2 = (0,1)
Nk_value = 700
kx, ky = np.meshgrid(np.linspace(-np.pi, np.pi, Nk_value), np.linspace(-np.pi, np.pi, Nk_value))
U_list = np.linspace(0.1, 10, 100)
en = np.zeros((len(J_list),len(U_list)))
a_min = np.zeros(len(U_list))
for i in range(len(J_list)):
    for j in range(len(U_list)):
        if(U_list[j]<J_list[i]):
            result = minimize_scalar(energy1, guess, args=(U_list[j], J_list[i], kx, ky, Nk_value), bounds=bounds1, options={'xtol':1e-6, 'ftol':1e-6,})
            a_min[j] = result.x
            en[i,j] = energy1(a_min[j],U_list[j], J_list[i], kx, ky, Nk_value)
        else:
            result = minimize_scalar(energy, guess, args=(U_list[j], kx, ky, Nk_value), bounds=bounds2, options={'xtol':1e-4, 'ftol':1e-4,})
            a_min[j] = result.x
            en[i,j] = energy(a_min[j],U_list[j], kx, ky, Nk_value)
    plt.plot(U_list[:],en[i,:],label=f'J={J_list[i]}',linestyle='-', marker='o',markersize=3)
    np.savetxt(f'Data/minimisation/datascalarJ={J_list[i]}_U={U_list[0]}-{U_list[-1]}_{len(U_list)}.dat', np.c_[U_list,a_min], delimiter=',') 

plt.legend()
plt.grid(True)
plt.xlabel('U')
plt.ylabel('Energy')
plt.title('Ground state energy')
plt.savefig('Plots/Minimisation/Energyscalar_vs_U.svg')
plt.show()

'''
plt.plot(U_list,a_min, label='minimum angle')
plt.legend()
plt.grid(True)
plt.title(f'minimum angle on the boundaries, U < J={J_value} ' )
plt.xlabel('U/t')
plt.ylabel(r'$\alpha$')
plt.title(f'J={J_value}')
plt.savefig(f"Plots/Minimisation/angle_J={J_value}_U={U_list[0]}-{U_list[-1]}_{len(U_list)}.svg")
plt.show()
'''
'''
en=np.zeros_like(U_list)
for i in range(len(en)): 
    en[i]=energy1(a_min[i], U_list[i], J_value, kx, ky, Nk)
plt.plot(U_list,en)
plt.legend()
plt.grid(True)
plt.title(f'minimum energy on the boundaries, U < J={J_value} ' )
plt.xlabel('U/t')
plt.ylabel(r'Energy')
plt.savefig(f"Plots/Minimisation/Boundary_energy_angle_J={J_value}_U={U_list[0]}-{U_list[-1]}_{len(U_list)}.svg")
plt.show()
'''
J_list = [0.0,1.0,2.0,3.0,4.0,5.0,6.0,7.0]
fig, axs = plt.subplots(2, len(J_list)//2, figsize=(15, 8))

for i in range(len(J_list)):
    U, a = np.loadtxt(f"Data/minimisation/datascalarJ={J_list[i]}_U=0.1-10.0_100.dat", delimiter=',', unpack=True)
    alpha = np.zeros(len(U))
    delta = np.zeros(len(U))
    for j in range(len(U)):
        if(U[j]<J_list[i]):
            alpha[j] = np.sin(a[j])
            delta[j] = np.cos(a[j])
        else:
            alpha[j] = 0
            delta[j] = a[j]
    axs[i // (len(J_list) // 2), i % (len(J_list) // 2)].plot(U, alpha, label=r'$\Delta_{SC}$', linestyle='-', marker='o', markersize=3)
    axs[i // (len(J_list) // 2), i % (len(J_list) // 2)].plot(U, delta, label=r'$\Delta_{CDW}$', linestyle='-', marker='o', markersize=3)
    axs[i // (len(J_list) // 2), i % (len(J_list) // 2)].set_title(f'J={J_list[i]}')

for ax in axs.flat:
    ax.set_xlabel(r'U')
    ax.set_ylabel(r'$\Delta$')
    ax.grid(True)
    ax.legend()

fig.tight_layout()
fig.show()
fig.savefig('Plots/Minimisation/Deltascalar_vs_U.svg')

