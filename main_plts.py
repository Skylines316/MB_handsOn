import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

######## Fermion Dispersion ########

# kx = np.linspace(-np.pi, np.pi, 100)
# ky = np.linspace(-np.pi, np.pi, 100)
# def xi_k(kx,ky):
#     return -4*(np.cos(kx)+np.cos(ky))

# fig = plt.figure()
# ax = fig.add_subplot(111, projection='3d')
# X, Y = np.meshgrid(kx, ky)
# zs = np.array(xi_k(np.ravel(X), np.ravel(Y)))
# Z = zs.reshape(X.shape)

# ax.plot_surface(X, Y, Z)

# ax.set_title(r'Single orbital free electrons dispersion')
# ax.set_xlabel(r'$k_x$')
# ax.set_ylabel(r'$k_y$')
# ax.set_zlabel(r'$\xi_k$')

# plt.savefig("Plots/FreeBand_Equal_t.svg")
# plt.show()

######## J=0, SC vs. CDW ########

i3, data3 = np.loadtxt("Data/J=0/data_cdw.dat", delimiter=',', unpack=True)
i2, data2 = np.loadtxt("Data/J=0/data2.dat", delimiter=',', unpack=True)
g1, alpha1, delta1 = np.loadtxt("Data/J=0/dataJ=0_U=0.01-20_200_alpha=0.dat", delimiter=',', unpack=True)
g2, alpha2, delta2 = np.loadtxt("Data/J=0/dataJ=0_U=0.01-20_200_delta=0.dat", delimiter=',', unpack=True)
g3, alpha3, delta3 = np.loadtxt("Data/J=0/dataJ=0_U=0.01-20_200.dat", delimiter=',', unpack=True)

plt.plot(i3,np.sqrt(2)*data3, label='CDW')
plt.plot(i2,2*np.sqrt(2)*data2, label='SC')
plt.plot(g1,delta1, label='SC_W')
plt.plot(g2,0.5*alpha2, label='CDW_W')
plt.legend()
plt.grid(True)
plt.title("Non interacting case: comparison between SC and CDW")
plt.xlabel('U/t')
plt.ylabel(r'$\Delta$')
plt.savefig('SC_vs_CDW_J=0.svg')
plt.show()

######## J, Energy ########

J_values = [0.1,0.3,0.5,0.7,1]
fig, axs = plt.subplots(1, len(J_values), figsize=(15, 5))

for i in range(len(J_values)):
    g1, e1 = np.loadtxt(f"Data/energy/J={J_values[i]}/dataJ={J_values[i]}_U=0.01-20_200_alpha=0_energy.dat", delimiter=',', unpack=True)
    g2, e2 = np.loadtxt(f"Data/energy/J={J_values[i]}/dataJ={J_values[i]}_U=0.01-20_200_delta=0_energy.dat", delimiter=',', unpack=True)
    g3, e3 = np.loadtxt(f"Data/energy/J={J_values[i]}/dataJ={J_values[i]}_U=0.01-20_200_energy.dat", delimiter=',', unpack=True)

    axs[i].plot(g1,e1, label='Only SC', linestyle='-', marker='o', markersize=2)
    axs[i].plot(g2,e2, label='Only CDW', linestyle='-', marker='v', markersize=2)
    axs[i].plot(g3,e3, label='SC + CDW', linestyle='-', marker='*', markersize=2)
    axs[i].set_title(f'Energy, J={J_values[i]}')

for ax in axs:
    ax.set_xlabel(r'U/t')
    ax.set_ylabel(r'Energy')
    ax.grid(True)
    ax.legend()

plt.tight_layout()
#plt.savefig('Plots/Delta_vs_Delta.svg')
plt.show()


# g1, e1 = np.loadtxt("Data/energy/dataJ=0.1_U=0.01-18_100_alpha=0_energy.dat", delimiter=',', unpack=True)
# g2, e2 = np.loadtxt("Data/energy/dataJ=0.1_U=0.01-18_100_delta=0_energy.dat", delimiter=',', unpack=True)
# g3, e3 = np.loadtxt("Data/energy/dataJ=0.1_U=0.01-18_100_energy.dat", delimiter=',', unpack=True)

# plt.scatter(g1,e1, label='Only SC', marker='.',s=8)
# plt.scatter(g2,e2, label='Only CDW', marker='*',s=8)
# plt.scatter(g3,e3, label='SC + CDW', marker='v',s=8)

# plt.legend()
# plt.grid(True)
# plt.title("Energy for J=0.1")
# plt.xlabel('U/t')
# plt.ylabel(r'$\Delta$')

# plt.savefig("Plots/Energy_J=0.1.svg")
# plt.show()

######## J, SC ########

J_values = [0,0.1,0.15,0.2,0.3,0.5,0.7,1]

fig, (ax1, ax2) = plt.subplots(1, 2)

# g, alpha, delta = np.loadtxt("Data/J=0.1/dataJ=0.1_U=0.01-20_200.dat", delimiter=',', unpack=True)
# ax1.plot(g, 0.5*alpha, label=r'J=0.1')
# ax2.plot(g, delta, label=r'J=0.1')

for i in range(np.size(J_values)):
    g, alpha, delta = np.loadtxt("Data/J="+str(J_values[i])+"/dataJ="+str(J_values[i])+"_U=0.01-20_200.dat", delimiter=',', unpack=True)
    #plt.plot(g, np.power(0.5*alpha,2)+np.power(delta,2), label=r'J='+str(J_values[i]))
    ax1.plot(g, 0.5*alpha, label=r'J='+str(J_values[i]))
    ax2.plot(g, delta, label=r'J='+str(J_values[i]))

#plt.plot(g, np.power(0.5*alpha,2)+np.power((1-0.1/g)*delta,2), label=r'$\Delta_{CDW}^2 + (1-J/U)\Delta_{SC}^2$')
ax1.legend()
ax2.legend()
ax1.grid(True)
ax2.grid(True)
ax1.set_title(r'$\Delta_{CDW}$')
ax2.set_title(r"$\Delta_{SC}$")
ax1.set_xlabel('U/t')
ax2.set_xlabel('U/t')
ax1.set_ylabel(r'$\Delta$')
#plt.savefig('Plots/Delta_vs_U.svg')
plt.show()

fig, axs = plt.subplots(1, 2, figsize=(15, 5))

# g, alpha, delta = np.loadtxt("Data/J=0.1/dataJ=0.1_U=0.01-18_100.dat", delimiter=',', unpack=True)
# axs[0].plot(g, np.power(0.5*alpha,2)+np.power(delta,2), label=r'J=0.1')
# axs[1].plot(g, np.power(0.5*alpha,2)+np.power((1-0.1/g)*delta,2), label=r'J=0.1')

for i in range(np.size(J_values)):
    g, alpha, delta = np.loadtxt("Data/J="+str(J_values[i])+"/dataJ="+str(J_values[i])+"_U=0.01-20_200.dat", delimiter=',', unpack=True)
    axs[0].plot(g, np.power(0.5*alpha,2)+np.power(delta,2), label=r'J='+str(J_values[i]))
    axs[1].plot(g, np.power(0.5*alpha,2)+np.power((1-J_values[i]/g)*delta,2), label=r'J='+str(J_values[i]))

for ax in axs:
    ax.legend()
    ax.set_xlabel(r'U/t')
    ax.set_ylabel(r'$\Delta$')
    ax.grid(True)

axs[0].set_title(r'$\Delta_{CDW}^2 + \Delta_{SC}^2$')
axs[1].set_title(r'$\Delta_{CDW}^2 + (1-J/U)^2\Delta_{SC}^2$')
plt.savefig('Plots/Delta_Sum_vs_U.svg')
plt.show()

num_subplots = len(J_values)
fig, axs = plt.subplots(1, num_subplots, figsize=(15, 5))

# g, alpha, delta = np.loadtxt("Data/J=0.1/dataJ=0.1_U=0.01-18_100.dat", delimiter=',', unpack=True)
# scatter = axs[0].scatter(0.5 * alpha, delta, c=g, cmap='viridis', alpha=0.7, label='J=0.1')
# axs[0].set_title('J=0.1')

for i in range(len(J_values)):
    g, alpha, delta = np.loadtxt(f"Data/J={J_values[i]}/dataJ={J_values[i]}_U=0.01-20_200.dat", delimiter=',', unpack=True)
    scatter = axs[i].scatter(0.5 * alpha, delta, c=g, cmap='viridis', alpha=0.7, label=f'J={J_values[i]}')
    axs[i].set_title(f'J={J_values[i]}')
    #cbar = fig.colorbar(scatter, ax=axs[i+1])
    #cbar.set_label('U/t')

cbar = fig.colorbar(scatter, ax=axs[-1])
cbar.set_label('U/t')

for ax in axs:
    ax.set_xlabel(r'$\Delta_{CDW}$')
    ax.set_ylabel(r'$\Delta_{SC}$')
    ax.grid(True)

plt.tight_layout()
plt.savefig('Plots/Delta_vs_Delta.svg')
plt.show()

#######################

g1, alpha1, delta1 = np.loadtxt("Data/J=0.1/dataJ=0.1_U=0.01-20_200.dat", delimiter=',', unpack=True)
g2, alpha2, delta2 = np.loadtxt("Data/J=0.1/dataJ=0.1_U=0.01-10_200.dat", delimiter=',', unpack=True)

plt.plot(g1[0:100],delta1[0:100], label='old')
plt.plot(g2,delta2, label='new')
plt.legend()
plt.grid(True)
plt.title("J=0.1")
plt.xlabel('U/t')
plt.ylabel(r'$\Delta$')
plt.show()