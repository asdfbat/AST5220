import numpy as np
import matplotlib.pyplot as plt
import matplotlib.pyplot as plt
plt.rcParams.update({'font.size': 22})
from PyConstants import Constants


data = np.loadtxt("../data/cosmology.txt")
print(data.shape)
N = data.shape[0]

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

Eta_rad_dom = 3*Constants.c/H_of_x/a

print(z[-2])
print(H_of_x[-1])


fig, ax = plt.subplots(figsize=(12, 8))
ax.plot(x, OmegaB, lw=4, label=r"$\Omega_{B}$", color="green")
ax.plot(x, OmegaCDM, lw=4, label=r"$\Omega_{CDM}$", color="blue")
ax.plot(x, OmegaCDM+OmegaB, lw=4, ls=":", label=r"$\Omega_{CDM+B}$", color="purple")
ax.plot(x, OmegaR, lw=4, label=r"$\Omega_{R}$", color="red")
ax.plot(x, OmegaLambda, lw=4, label=r"$\Omega_\Lambda$", color="black")
ax.plot(x, OmegaB + OmegaCDM + OmegaR + OmegaLambda, lw=4, label=r"$\sum \Omega$", color="orange", ls="--")
ax.fill_between(x, 0, 1, OmegaR>(OmegaCDM+OmegaB), alpha=0.2, color="red")
ax.fill_between(x, 0, 1, (OmegaR<(OmegaCDM+OmegaB) * ((OmegaCDM+OmegaB) > OmegaLambda)), alpha=0.2, color="blue")
ax.fill_between(x, 0, 1, (OmegaCDM+OmegaB) < OmegaLambda, alpha=0.2, color="black")
ax.set_xlabel("x")
ax.set_ylabel(r"$\Omega_i(x)$")
plt.legend()
plt.tight_layout()
plt.savefig("../figs/Omegas.pdf", bbox_inches="tight")
plt.clf()

rad_mat_equal = np.argmin(np.abs(OmegaR-(OmegaCDM+OmegaB)))
DE_mat_equal = np.argmin(np.abs(OmegaLambda[N//2:]-(OmegaCDM[N//2:]+OmegaB[N//2:]))) + N//2
print(x[rad_mat_equal], x[DE_mat_equal])
print(z[rad_mat_equal], z[DE_mat_equal])

fig, ax = plt.subplots(2, 2, figsize=(16, 12))
ax[0,0].semilogy(x, eta_of_x, lw=4)
#ax[0,0].plot(x, Eta_rad_dom)
ax[0,0].axvline(x=x[rad_mat_equal], ls="--", lw=2, c="r", label="rad/mat equality")
ax[0,0].axvline(x=x[DE_mat_equal], ls="--", lw=2, c="k", label="mat/DE equality")
ax[0,0].set_xlabel("x")
ax[0,0].set_ylabel(r"$\eta(x)$")
ax[0,0].set_title("$\eta \ - \ [Gpc]$")
#ax[0,0].set_ylim(0, 12.5)
ax[0,0].legend()

ax[0,1].semilogy(x, Hp_of_x, lw=4)
ax[0,1].axvline(x=x[rad_mat_equal], ls="--", lw=2, c="r")
ax[0,1].axvline(x=x[DE_mat_equal], ls="--", lw=2, c="k")
ax[0,1].set_xlabel("x")
ax[0,1].set_ylabel(r"$\mathcal{H}(x)$")
ax[0,1].set_title("$\mathcal{H} \ - \ [km/s/Mpc]$")

ax[1,0].semilogy(x, H_of_x, lw=4)
ax[1,0].axvline(x=x[rad_mat_equal], ls="--", lw=2, c="r")
ax[1,0].axvline(x=x[DE_mat_equal], ls="--", lw=2, c="k")
ax[1,0].set_xlabel("x")
ax[1,0].set_ylabel("H(x)")
ax[1,0].set_ylabel(r"$H(x)$")
ax[1,0].set_title("$H \ - \ [km/s/Mpc]$")

ax[1,1].loglog(z, H_of_x, lw=4)
ax[1,1].axvline(x=z[rad_mat_equal], ls="--", lw=2, c="r")
ax[1,1].axvline(x=z[DE_mat_equal], ls="--", lw=2, c="k")
ax[1,1].set_xlim(1.5*z[0], 0.01)
ax[1,1].set_xlabel("z")
ax[1,1].set_ylabel(r"$H(z)$")
ax[1,1].set_title("$H \ - \ [km/s/Mpc]$")

plt.tight_layout()
plt.savefig("../figs/Eta.pdf", bbox_inches="tight")
plt.clf()