# import numpy as np
# import matplotlib.pyplot as plt
# from scipy.optimize import minimize
# import math
# from mpl_toolkits.mplot3d import Axes3D

# def xi(kx, ky):
#     return -2*(np.cos(kx)+np.cos(ky))

# Nk = 100
# kx, ky = np.meshgrid(np.linspace(-math.pi, math.pi, Nk), np.linspace(-math.pi, math.pi, Nk))
# dk = kx[0,1]-kx[0,0]

# def energy(Delta, Alpha, U, J):
#     return -2/(Nk**2)*np.sum(np.sqrt(np.power(xi(kx,ky),2) + (U-J)**2*np.power(Delta,2) + 0.25*U**2*np.power(Alpha,2)))*dk**2 + 2*Nk**2*((U-J)*np.power(Delta,2)+U*(1+0.25*np.power(Alpha,2)))

# Delta, Alpha = np.meshgrid(np.linspace(0, 1, 20), np.linspace(0, 1, 20))
# z = energy(Delta, Alpha, 0.1, 0)
# z.shape
# Z = z.reshape(Delta.shape)
# fig = plt.figure(figsize=(8, 6))
# ax = fig.add_subplot(111, projection="3d")
# ax.plot_surface(Delta, Alpha, Z, cmap='viridis')
# plt.show()

import numpy as np
import matplotlib.pyplot as plt
import math
from mpl_toolkits.mplot3d import Axes3D
from scipy.optimize import minimize

def xi(kx, ky):
    return -2 * (np.cos(kx) + np.cos(ky))

def energy(params, *args):
    Delta, Alpha = params
    U, J, kx, ky, Nk = args
    return 2 / Nk**2 * np.sum(- np.sqrt(np.power(xi(kx, ky), 2) + (U - J)**2 * np.power(Delta, 2) + U**2 * np.power(Alpha, 2))) + 2 * ((U - J) * np.power(Delta, 2) + U * ( np.power(Alpha, 2)))

def energy_circle(params, *args):
    Rho, Theta = params
    U, J, kx, ky, Nk = args
    return 2 / Nk**2 * np.sum(- np.sqrt(np.power(xi(kx, ky), 2) + np.power(Rho, 2))) + 2 /((U - J)*U) * np.power(Rho, 2) * (U - J * np.power(np.cos(Theta), 2))

def con(t):
    return - t[0]**2 - t[1]**2 + 1

Nk = 500
kx, ky = np.meshgrid(np.linspace(-math.pi, math.pi, Nk), np.linspace(-math.pi, math.pi, Nk))

Delta, Alpha = np.meshgrid(np.linspace(-1, 1, 50), np.linspace(-1, 1, 50))
U_value = 4.9
J_value = 5.
# Calculate the energy for each point in the mesh
z = np.zeros_like(Delta)
for i in range(Delta.shape[0]):
    for j in range(Delta.shape[1]):
        z[i, j] = energy([Delta[i, j], Alpha[i, j]], U_value, J_value, kx, ky, Nk)

fig = plt.figure(figsize=(8, 6))
ax = fig.add_subplot(111, projection="3d")
ax.plot_surface(Delta, Alpha, z, cmap='viridis')
ax.set_xlabel(r'$\Delta_{SC}$')
ax.set_ylabel(r'$\Delta_{CDW}$')
ax.set_zlabel(r'Energy')
ax.set_title(f'Energy surface for U={U_value} and J={J_value}')
plt.savefig(f"Plots/Surfaces/J={J_value}_U={U_value}.svg")
plt.show()

#####################################################

bounds = [(0, 1), (0, 1)]
cons = {'type':'ineq', 'fun': con}
J_value = 0.
Nk_value = 600
kx, ky = np.meshgrid(np.linspace(-np.pi, np.pi, Nk_value), np.linspace(-np.pi, np.pi, Nk_value))
U_list = np.linspace(0.1, 10, 40)
Delta_min = np.zeros(len(U_list))
alpha_min = np.zeros(len(U_list))
guess = np.zeros((2,len(U_list)))

for i in range(len(guess[0,:])):
    if (U_list[i]<J_value):
        guess[:,i]=[0.9,0.1]
    if (U_list[i]>=J_value):
        guess[:,i]=[0.01,0.01] # [0.,0.2]

for i in range(len(U_list)):
    # guess = [0.9, max((U_list[i]-J_value)/U_list[i],0)]
    result = minimize(energy, guess[:,i], args=(U_list[i], J_value, kx, ky, Nk_value), method='SLSQP', bounds=bounds, constraints=cons, options={'xtol':1e-6})#, 'ftol':1e-6,})
    Delta_min[i] = result.x[0]
    alpha_min[i] = result.x[1]

# np.savetxt(f'Data/minimisation/datanewJ={J_value}_U={U_list[0]}-{U_list[-1]}_{len(U_list)}.dat', np.c_[U_list,Delta_min,alpha_min], delimiter=',') 

plt.figure()
plt.plot(U_list, Delta_min, label='SC', marker = 'o',markersize=4)
plt.plot(U_list, alpha_min, label='CDW', marker = 'o',markersize=4)
plt.legend()
plt.grid(True)
plt.xlabel('U')
plt.ylabel(r'$\Delta$')
plt.title(f'J={J_value}')
# plt.savefig(f"Plots/Minimisation/J={J_value}_U={U_list[0]}-{U_list[-1]}_{len(U_list)}.svg")
plt.show()

J_value = 2.
Nk_value = 500
kx, ky = np.meshgrid(np.linspace(-np.pi, np.pi, Nk_value), np.linspace(-np.pi, np.pi, Nk_value))
U_list = np.linspace(0.1, 10, 40)
Rho_min = np.zeros(len(U_list))
Theta_min = np.zeros(len(U_list))

for i in range(len(U_list)):
    bounds = [(0, np.sqrt(np.power(abs(U_list[i]-J_value),2)+np.power(U_list[i]/2,2))), (0, 0.5*np.pi)]
    guess = [np.sqrt(np.power(abs(U_list[i]-J_value),2)+np.power(U_list[i]/2,2)), 0.25*np.pi]
    print(guess)
    result = minimize(energy_circle, guess, args=(U_list[i], J_value, kx, ky, Nk_value), bounds=bounds, options={'xtol':1e-6, 'ftol':1e-6,})
    Rho_min[i] = result.x[0]
    Theta_min[i] = result.x[1]

##############################################

Nk_value = 1000
kx, ky = np.meshgrid(np.linspace(-np.pi, np.pi, Nk_value), np.linspace(-np.pi, np.pi, Nk_value))
J_list = [0.0,1.0,2.0,3.0,4.0,5.0,6.0,7.0]
data = []
for i in range(len(J_list)):
    data.append(np.loadtxt(f'Data/minimisation/dataJ={J_list[i]}_U=0.1-10.0_40.dat', delimiter=',',unpack=True))
data = np.array(data)

en = np.zeros((len(J_list),len(data[0,0,:])))
plt.figure(figsize=(6,4))
for i in range(len(J_list)):
    for j in range(len(en[0,:])):
        en[i,j] = energy([data[i,1,j], data[i,2,j]], data[i,0,j], J_list[i], kx, ky, Nk_value)
    plt.plot(data[i,0,:],en[i,:],label=f'J={J_list[i]}',linestyle='-', marker='o',markersize=3)

# u5, d5, a5 = np.loadtxt(f'Data/minimisation/dataJ=5.0_U=0.1-10.0_40.dat', delimiter=',',unpack=True)
# en5=np.zeros(len(u5))
# for j in range(len(en5)):
#     en5[j] = energy([d5[j], a5[j]], u5[j], 5.0, kx, ky, Nk_value)
# plt.plot(u5,en5,label=f'new',linestyle='-', marker='o',markersize=3)

plt.legend()
plt.grid(True)
plt.xlabel('U')
plt.ylabel('Energy')
plt.title('Ground state energy')
# plt.savefig('Plots/Minimisation/EnergyOld_vs_U.svg')
plt.show()

##############################################

J_list = [0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0]

fig, axs = plt.subplots(2, len(J_list)//2, figsize=(15, 8))

for i in range(len(J_list)):
    g, delta, alpha = np.loadtxt(f"Data/minimisation/dataJ={J_list[i]}_U=0.1-10.0_40.dat", delimiter=',', unpack=True)
    axs[i // (len(J_list) // 2), i % (len(J_list) // 2)].plot(g, delta, label=r'$\Delta_{SC}$', linestyle='-', marker='o', markersize=3)
    axs[i // (len(J_list) // 2), i % (len(J_list) // 2)].plot(g, alpha, label=r'$\Delta_{CDW}$', linestyle='-', marker='o', markersize=3)
    axs[i // (len(J_list) // 2), i % (len(J_list) // 2)].set_title(f'J={J_list[i]}')

for ax in axs.flat:
    ax.set_xlabel(r'U')
    ax.set_ylabel(r'$\Delta$')
    ax.grid(True)
    ax.legend()

plt.tight_layout()
plt.savefig('Plots/Minimisation/DeltaOld_vs_U.svg')
plt.show()
################################################################################
