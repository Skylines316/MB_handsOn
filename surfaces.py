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
    return 2 / Nk**2 * np.sum(- np.sqrt(np.power(xi(kx, ky), 2) + (U - J)**2 * np.power(Delta, 2) + U**2 * np.power(Alpha, 2))) + 2 * ((U - J) * np.power(Delta, 2) + U * (np.power(Alpha, 2)))

Nk = 500
kx, ky = np.meshgrid(np.linspace(-math.pi, math.pi, Nk), np.linspace(-math.pi, math.pi, Nk))
dk = kx[0, 1] - kx[0, 0]

Delta, Alpha = np.meshgrid(np.linspace(-1, 1, 30), np.linspace(-1, 1, 30))
U_value = 5.5
J_value = 5
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
plt.savefig("Plots/Surfaces/U=5,5_J=5.svg")
plt.show()