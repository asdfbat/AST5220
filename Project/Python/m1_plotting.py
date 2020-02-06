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
Hp_of_x_SI = data[:,2]
Hp_of_x = Hp_of_x_SI/1000*Constants.Mpc
H_of_x_SI = Hp_of_x_SI/np.exp(x)
H_of_x = Hp_of_x/np.exp(x)
OmegaB = data[:,4]
OmegaCDM = data[:,5]
OmegaLambda = data[:,6]
OmegaR = data[:,7]

a = np.exp(x)
z = 1/a - 1



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
fig.legend()
fig.tight_layout()
fig.savefig("../m1_figs/Omegas.pdf", bbox_inches="tight")
plt.close(fig)

rad_mat_equal = np.argmin(np.abs(OmegaR-(OmegaCDM+OmegaB)))
DE_mat_equal = np.argmin(np.abs(OmegaLambda[N//2:]-(OmegaCDM[N//2:]+OmegaB[N//2:]))) + N//2
print(x[rad_mat_equal], x[DE_mat_equal])
print(z[rad_mat_equal], z[DE_mat_equal])
print(a[rad_mat_equal], a[DE_mat_equal])

rad_max_idx = 0
mat_max_idx = np.argmax(OmegaCDM + OmegaB)
DE_max_idx = N-1


a_star = a[mat_max_idx]
a_star_idx = mat_max_idx
a_lambda = a[DE_max_idx]
a_lambda_idx = DE_max_idx

Eta_rad_dom = Constants.c/H_of_x_SI/a/Constants.Mpc/1000  # Analytical Eta in Gpc.
Eta_rad_dom_a_star = Constants.c/H_of_x_SI[a_star_idx]/a_star/Constants.Mpc/1000

Eta_mat_dom = eta_of_x[mat_max_idx] + 2*Constants.c*(1/(a*H_of_x_SI) - a_star**0.5/(a**1.5*H_of_x_SI))/Constants.Mpc/1000
Eta_mat_dom_a_lambda = Eta_rad_dom_a_star + 2*Constants.c*(1/(a_lambda*H_of_x_SI) - a_star**0.5/(a_lambda**1.5*H_of_x_SI))/Constants.Mpc/1000

Eta_DE_dom = eta_of_x[DE_max_idx] + Constants.c/(H_of_x_SI)*(1/a_lambda - 1/a)/Constants.Mpc/1000

fig, ax = plt.subplots(figsize=(12, 8))
ax.semilogy(x, eta_of_x)
ax.semilogy(x, Eta_rad_dom)
ax.semilogy(x, Eta_mat_dom)
ax.semilogy(x[N//2:], Eta_DE_dom[N//2:])
ax.axvline(x=x[rad_mat_equal], ls="--")
ax.axvline(x=x[DE_mat_equal], ls="--")
plt.show()


fig, ax = plt.subplots(2, 2, figsize=(16, 12))
ax[0,0].semilogy(x, eta_of_x, lw=4)
#ax[0,0].plot(x, Eta_rad_dom)
ax[0,0].axvline(x=x[rad_mat_equal], ls="--", lw=2, c="r", label="rad/mat equality")
ax[0,0].axvline(x=x[DE_mat_equal], ls="--", lw=2, c="k", label="mat/DE equality")
ax[0,0].set_xlabel("x")
ax[0,0].set_ylabel(r"$\eta(x)$")
ax[0,0].set_title(r"$\eta \ - \ [Gpc]$")
#ax[0,0].set_ylim(0, 12.5)
ax[0,0].legend()

ax[0,1].semilogy(x, Hp_of_x, lw=4)
ax[0,1].axvline(x=x[rad_mat_equal], ls="--", lw=2, c="r")
ax[0,1].axvline(x=x[DE_mat_equal], ls="--", lw=2, c="k")
ax[0,1].set_xlabel("x")
ax[0,1].set_ylabel(r"$\mathcal{H}(x)$")
ax[0,1].set_title(r"$\mathcal{H} \ - \ [km/s/Mpc]$")

ax[1,0].semilogy(x, H_of_x, lw=4)
ax[1,0].axvline(x=x[rad_mat_equal], ls="--", lw=2, c="r")
ax[1,0].axvline(x=x[DE_mat_equal], ls="--", lw=2, c="k")
ax[1,0].set_xlabel("x")
ax[1,0].set_ylabel("H(x)")
ax[1,0].set_ylabel(r"$H(x)$")
ax[1,0].set_title(r"$H \ - \ [km/s/Mpc]$")

ax[1,1].loglog(z, H_of_x, lw=4)
ax[1,1].axvline(x=z[rad_mat_equal], ls="--", lw=2, c="r")
ax[1,1].axvline(x=z[DE_mat_equal], ls="--", lw=2, c="k")
ax[1,1].set_xlim(1.5*z[0], 0.01)
ax[1,1].set_xlabel("z")
ax[1,1].set_ylabel(r"$H(z)$")
ax[1,1].set_title(r"$H \ - \ [km/s/Mpc]$")

fig.tight_layout()
fig.savefig("../m1_figs/Eta.pdf", bbox_inches="tight")
plt.close(fig)