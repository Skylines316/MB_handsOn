import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D


J_values = [0,0.5,0.7,1]


for i in range(len(J_values)):
    g1, e1 = np.loadtxt(f"Data/energy/J={J_values[i]}/dataJ={J_values[i]}_U=0.01-20_200_alpha=0_energy.dat", delimiter=',', unpack=True)
    g2, e2 = np.loadtxt(f"Data/energy/J={J_values[i]}/dataJ={J_values[i]}_U=0.01-20_200_delta=0_energy.dat", delimiter=',', unpack=True)
   # g3, e3 = np.loadtxt(f"Data/energy/J={J_values[i]}/dataJ={J_values[i]}_U=0.01-20_200_energy.dat", delimiter=',', unpack=True)

    plt.plot(g1[g1>J_values[i]],e1[g1>J_values[i]], label=f"J={J_values[i]}", linestyle='-', marker='o', markersize=2)
   # axs[i].plot(g2,e2, label='Only CDW', linestyle='-', marker='v', markersize=2)
    #axs[i].plot(g3,e3, label='SC + CDW', linestyle='-', marker='*', markersize=2)
    #axs[i].set_title(f'Energy, J={J_values[i]}')
    plt.title('Ground state energy')
    plt.xlabel(r'U/t')
    plt.ylabel(r'Energy')
    plt.grid(True)
    plt.xlim(0,10)
    plt.legend()

plt.tight_layout()
plt.savefig('Plots/Energyminimum_alpha=0.svg')
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