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
    return 2 / Nk**2 * np.sum(- np.sqrt(np.power(xi(kx, ky), 2) + (U - J)**2 * np.power(Delta, 2) + 0.25 * U**2 * np.power(Alpha, 2))) + 2 * ((U - J) * np.power(Delta, 2) + U * (1 + 0.25 * np.power(Alpha, 2)))

Nk = 500
kx, ky = np.meshgrid(np.linspace(-math.pi, math.pi, Nk), np.linspace(-math.pi, math.pi, Nk))
dk = kx[0, 1] - kx[0, 0]

Delta, Alpha = np.meshgrid(np.linspace(0, 1, 20), np.linspace(0, 1, 20))
U_value = 1
J_value = 0.1
# Calculate the energy for each point in the mesh
z = np.zeros_like(Delta)
for i in range(Delta.shape[0]):
    for j in range(Delta.shape[1]):
        z[i, j] = energy([Delta[i, j], Alpha[i, j]], U_value, J_value, kx, ky, Nk)

fig = plt.figure(figsize=(8, 6))
ax = fig.add_subplot(111, projection="3d")
ax.plot_surface(Delta, Alpha, z, cmap='viridis')
ax.set_xlabel(r'$\Delta$')
ax.set_ylabel(r'$\alpha$')
ax.set_zlabel(r'Energy')
plt.show()

guess = [0.1, 0.1]
bounds = [(0, 1), (0, 1)]
J_value = 0.05
Nk_value = 500
kx, ky = np.meshgrid(np.linspace(-np.pi, np.pi, Nk_value), np.linspace(-np.pi, np.pi, Nk_value))
U_list = np.linspace(1, 20, 40)
Delta_min = np.zeros(len(U_list))
alpha_min = np.zeros(len(U_list))

for i in range(len(U_list)):
    result = minimize(energy, guess, args=(U_list[i], J_value, kx, ky, Nk_value), bounds=bounds, options={'xtol':1e-6, 'ftol':1e-6,})
    Delta_min[i] = result.x[0]
    alpha_min[i] = result.x[1]

plt.plot(U_list,Delta_min, label='SC')
plt.plot(U_list,alpha_min, label='CDW')
plt.legend()
plt.grid(True)
plt.xlabel('U/t')
plt.ylabel(r'$\Delta$')
plt.title(f'J={J_value}')
#plt.savefig(f"Plots/Minimisation/J={J_value}_U=0.1-20_20.svg")
plt.show()