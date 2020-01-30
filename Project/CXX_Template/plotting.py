import numpy as np
import matplotlib.pyplot as plt
import matplotlib.pyplot as plt
plt.rcParams.update({'font.size': 22})
from PyConstants import Constants

data = np.loadtxt("../data/cosmology.txt")
print(data.shape)

x = data[:,0]
eta_of_x = data[:,1]#/Constants.Mpc
Hp_of_x = data[:,2]
OmegaB = data[:,4]
OmegaCDM = data[:,5]
OmegaLambda = data[:,6]
OmegaR = data[:,7]


# plt.figure(figsize=(20,14))
# plt.plot(x, OmegaB, lw=4, label="OmegaB")
# plt.plot(x, OmegaCDM, lw=4, label="OmegaCDM")
# plt.plot(x, OmegaR, lw=4, label="OmegaR")
# plt.plot(x, OmegaLambda, lw=4, label="OmegaLambda")
# plt.plot(x, OmegaB + OmegaCDM + OmegaR + OmegaLambda, lw=4, label="Total")
# plt.legend()
# plt.show()

plt.figure(figsize=(20,14))
plt.plot(x, eta_of_x, lw=4)
plt.xlabel("x")
plt.ylabel("eta(x)")
plt.show()

plt.figure(figsize=(20,14))
plt.plot(x, Hp_of_x, lw=4)
plt.xlabel("x")
plt.ylabel("Hp(x)")
plt.show()