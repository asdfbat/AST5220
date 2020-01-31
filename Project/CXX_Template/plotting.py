import numpy as np
import matplotlib.pyplot as plt
import matplotlib.pyplot as plt
plt.rcParams.update({'font.size': 22})
from PyConstants import Constants

data = np.loadtxt("../data/cosmology.txt")
print(data.shape)

x = data[:,0]
eta_of_x = data[:,1]/Constants.Mpc/1000  # Eta(x) in Gpc
Hp_of_x = data[:,2]/1000*Constants.Mpc
H_of_x = Hp_of_x/np.exp(x)
OmegaB = data[:,4]
OmegaCDM = data[:,5]
OmegaLambda = data[:,6]
OmegaR = data[:,7]

a = np.exp(x)
z = 1/a - 1

print(z[-2])
print(H_of_x[-1])

fig, ax = plt.subplots(figsize=(20,14))
ax.plot(x, OmegaB, lw=4, label="OmegaB")
ax.plot(x, OmegaCDM, lw=4, label="OmegaCDM")
ax.plot(x, OmegaR, lw=4, label="OmegaR")
ax.plot(x, OmegaLambda, lw=4, label="OmegaLambda")
ax.plot(x, OmegaB + OmegaCDM + OmegaR + OmegaLambda, lw=4, label="Total")
ax.fill_between(x, 0, 1, OmegaR>(OmegaCDM+OmegaB), alpha=0.2)
ax.fill_between(x, 0, 1, (OmegaR<(OmegaCDM+OmegaB) * ((OmegaCDM+OmegaB) > OmegaLambda)), alpha=0.2)
ax.fill_between(x, 0, 1, (OmegaCDM+OmegaB) < OmegaLambda, alpha=0.2)
plt.legend()
plt.show()

fig, ax = plt.subplots(2, 2, figsize=(24, 18))
ax[0,0].plot(x, eta_of_x, lw=4)
ax[0,0].set_xlabel("x")
ax[0,0].set_ylabel("eta(x) [Gpc]")

ax[0,1].plot(x, Hp_of_x, lw=4)
ax[0,1].set_xlabel("x")
ax[0,1].set_ylabel("Hp(x)")

ax[1,0].semilogy(x, H_of_x, lw=4)
ax[1,0].set_xlabel("x")
ax[1,0].set_ylabel("H(x)")

ax[1,1].loglog(z, H_of_x, lw=4)
ax[1,1].set_xlim(1.5*z[0], 0.01)
ax[1,1].set_xlabel("z")
ax[1,1].set_ylabel("H(z)")
plt.show()