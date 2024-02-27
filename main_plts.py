import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

######## Fermion Dispersion ########

kx = np.linspace(-np.pi, np.pi, 100)
ky = np.linspace(-np.pi, np.pi, 100)
def xi_k(kx,ky):
    return -4*(np.cos(kx)+np.cos(ky))

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
X, Y = np.meshgrid(kx, ky)
zs = np.array(xi_k(np.ravel(X), np.ravel(Y)))
Z = zs.reshape(X.shape)

ax.plot_surface(X, Y, Z)

ax.set_title(r'Single orbital free electrons dispersion')
ax.set_xlabel(r'$k_x$')
ax.set_ylabel(r'$k_y$')
ax.set_zlabel(r'$\xi_k$')

plt.savefig("Plots/FreeBand_Equal_t.svg")
plt.show()

######## J=0, SC vs. CDW ########

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

######## J, Energy ########

g1, e1 = np.loadtxt("Data/energy/dataJ=0.1_U=0.01-18_100_alpha=0_energy.dat", delimiter=',', unpack=True)
g2, e2 = np.loadtxt("Data/energy/dataJ=0.1_U=0.01-18_100_delta=0_energy.dat", delimiter=',', unpack=True)
g3, e3 = np.loadtxt("Data/energy/dataJ=0.1_U=0.01-18_100_energy.dat", delimiter=',', unpack=True)

plt.scatter(g1,e1, label='Only SC', marker='.',s=8)
plt.scatter(g2,e2, label='Only CDW', marker='*',s=8)
plt.scatter(g3,e3, label='SC + CDW', marker='v',s=8)

plt.legend()
plt.grid(True)
plt.title("Energy for J=0.1")
plt.xlabel('U/t')
plt.ylabel(r'$\Delta$')

plt.savefig("Plots/Energy_J=0.1.svg")
plt.show()

######## J, SC ########

g, alpha, delta = np.loadtxt("Data/J=0.1/dataJ=0.1_U=0.01-18_100.dat", delimiter=',', unpack=True)

plt.plot(g, np.power(0.5*alpha,2)+np.power((1-0.1/g)*delta,2), label=r'$\Delta_{CDW}^2 + (1-J/U)\Delta_{SC}^2$')
plt.plot(g, np.power(0.5*alpha,2)+np.power(delta,2), label=r'$\Delta_{CDW}^2 + \Delta_{SC}^2$')
plt.legend()
plt.grid(True)
plt.title("J=0.1")
plt.xlabel('U/t')
plt.ylabel(r'$\Delta$')
plt.savefig('Plots/Delta_vs_U_J=0.1.svg')
plt.show()

plt.scatter(0.5*alpha, delta, c=g, cmap='viridis', alpha=0.7)
cbar = plt.colorbar()
cbar.set_label('U/t')
plt.grid(True)
plt.title("J=0.1")
plt.xlabel(r'$\Delta_{CDW}$')
plt.ylabel(r'$\Delta_{SC}$')
plt.savefig('Plots/Delta_vs_Delta_J=0.1.svg')
plt.show()

#######################

g1, alpha1, delta1 = np.loadtxt("Data/J=0.1/dataJ=0.1_U=0.01-18_100_delta=0.dat", delimiter=',', unpack=True)
g2, alpha2, delta2 = np.loadtxt("Data/J=0.1/dataJ=0.1_U=0.01-18_100_alpha=0.dat", delimiter=',', unpack=True)

plt.plot(g1,0.5*alpha1, label='CDW')
plt.plot(g2,delta2, label='SC')
plt.legend()
plt.grid(True)
plt.title("J=0.1")
plt.xlabel('U/t')
plt.ylabel(r'$\Delta$')
plt.show()