import numpy as np
import matplotlib.pyplot as plt

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

######## J ########

