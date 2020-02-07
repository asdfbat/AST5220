import numpy as np
import matplotlib.pyplot as plt
import matplotlib.pyplot as plt
plt.rcParams.update({'font.size':18})
from PyConstants import Constants


data = np.loadtxt("../data/cosmology.txt")
N = data.shape[0]

x = data[:,0]
eta_of_x = data[:,1]/Constants.Mpc/Constants.km     # Eta(x) in Gpc
Hp_of_x_SI = data[:,2]                              # aH(x) in SI
Hp_of_x = Hp_of_x_SI/Constants.km*Constants.Mpc     # aH(x) in km/s/Mpc
H_of_x_SI = Hp_of_x_SI/np.exp(x)                    # H(x) in SI  
H_of_x = Hp_of_x/np.exp(x)                          # H(x) in km/s/Mpc
OmegaB = data[:,4]
OmegaCDM = data[:,5]
OmegaLambda = data[:,6]
OmegaR = data[:,7]

a = np.exp(x)
z = 1/a - 1
today_idx = np.argmin(np.abs(x))

H_of_x_SI_rad_dom = Constants.H0_over_h*0.7*np.sqrt(a**(-4)*OmegaR[today_idx])
H_of_x_rad_dom = H_of_x_SI_rad_dom/Constants.km*Constants.Mpc
H_of_x_SI_mat_dom = Constants.H0_over_h*0.7*np.sqrt(a**(-3)*(OmegaCDM[today_idx] + OmegaB[today_idx]))
H_of_x_mat_dom = H_of_x_SI_mat_dom/Constants.km*Constants.Mpc
H_of_x_SI_DE_dom = Constants.H0_over_h*0.7*np.sqrt(OmegaLambda[today_idx])
H_of_x_DE_dom = H_of_x_SI_DE_dom/Constants.km*Constants.Mpc

# --- Finding the equality points --- #
rad_mat_equal = np.argmin(np.abs(OmegaR[N//4:3*N//4]-(OmegaCDM[N//4:3*N//4]+OmegaB[N//4:3*N//4]))) + N//4
DE_mat_equal = np.argmin(np.abs(OmegaLambda[N//2:]-(OmegaCDM[N//2:]+OmegaB[N//2:]))) + N//2
print("   Rad/Mat equal   Mat/DE equal")
print("x =  %.3e      %.3e" % (x[rad_mat_equal], x[DE_mat_equal]))
print("a =  %.3e      %.3e" % (a[rad_mat_equal], a[DE_mat_equal]))
print("z =  %.3e      %.3e" % (z[rad_mat_equal], z[DE_mat_equal]))

# --- Finding the maxima points --- #
rad_max_idx = 0
mat_max_idx = np.argmax(OmegaCDM + OmegaB)
DE_max_idx = N-1


# --- Omega(x) plots --- #
fig, ax = plt.subplots(figsize=(10, 6))
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
fig.legend(loc=(0.1, 0.4))
fig.tight_layout()
fig.savefig("../m1_figs/Omegas.pdf", bbox_inches="tight")
plt.close(fig)


# Anchor points for the analytical solutions of Eta(a).
a_star = a[mat_max_idx]
a_lambda = a[DE_max_idx]

# --- Analytical solutions of Eta, with and without H_i = H replacement. --- #
# Analytic Eta(a) in the radiation dominated regime, units of Gpc.
Eta_rad_dom_H = Constants.c/H_of_x_SI/a/Constants.Mpc/1000
Eta_rad_dom_Hi = Constants.c/H_of_x_SI_rad_dom/a/Constants.Mpc/1000

# Analytical Eta(a) in the matter dominated regime, units of Gpc.
Eta_mat_dom_H = eta_of_x[mat_max_idx] + 2*Constants.c*(1/(a*H_of_x_SI) - a_star**0.5/(a**1.5*H_of_x_SI))/Constants.Mpc/1000
Eta_mat_dom_Hi = eta_of_x[mat_max_idx] + 2*Constants.c*(1/(a*H_of_x_SI_mat_dom) - a_star**0.5/(a**1.5*H_of_x_SI_mat_dom))/Constants.Mpc/1000

# Analytical Eta(a) in the DE dominated regime, units of Gpc.
Eta_DE_dom_H = eta_of_x[DE_max_idx] + Constants.c/(H_of_x_SI)*(1/a_lambda - 1/a)/Constants.Mpc/1000
Eta_DE_dom_Hi = eta_of_x[DE_max_idx] + Constants.c/(H_of_x_SI_DE_dom)*(1/a_lambda - 1/a)/Constants.Mpc/1000


# --- Plotting Eta --- #
fig, ax = plt.subplots(1, 2, figsize=(13.5, 6), sharey=True)
ax[0].semilogy(x, eta_of_x, lw=4, c="orange", label=r"$\eta(x)$")
ax[0].semilogy(x, Eta_rad_dom_Hi, lw=4, ls=":", c="r", label=r"$\eta_r(x)$")
ax[0].semilogy(x, Eta_mat_dom_Hi, lw=3, ls="--", c="b", label=r"$\eta_m(x)$")
ax[0].semilogy(x, Eta_DE_dom_Hi, lw=3, ls="-.", c="k", label=r"$\eta_\Lambda(x)$")
ax[1].semilogy(x, eta_of_x, lw=4, c="orange", label=r"$\eta(x)$")
ax[1].semilogy(x, Eta_rad_dom_H, lw=4, ls=":", c="r", label=r"$\eta_r(x)$")
ax[1].semilogy(x, Eta_mat_dom_H, lw=3, ls="--", c="b", label=r"$\eta_m(x)$")
ax[1].semilogy(x, Eta_DE_dom_H, lw=3, ls="-.", c="k", label=r"$\eta_\Lambda(x)$")
for i in range(2):
    ax[i].fill_between(x, 0, 100, OmegaR>(OmegaCDM+OmegaB), alpha=0.2, color="red")
    ax[i].fill_between(x, 0, 100, (OmegaR<(OmegaCDM+OmegaB) * ((OmegaCDM+OmegaB) > OmegaLambda)), alpha=0.2, color="blue")
    ax[i].fill_between(x, 0, 100, (OmegaCDM+OmegaB) < OmegaLambda, alpha=0.2, color="black")
    ax[i].set_xlabel("x")
    ax[i].set_xlim(-15, 4)
ax[0].set_ylabel(r"$\eta(x) \ - \ [Gpc]$")
ax[0].set_title(r"$\eta(x)$")
ax[1].set_title(r"$\eta(x)\ - \ H_i = H$")
plt.ylim(1e-4, 1e2)
plt.legend()
fig.tight_layout()
fig.savefig("../m1_figs/Eta.pdf", bbox_inches="tight")
plt.close(fig)


# --- Chained Analytical expressions of Eta --- #
Eta_rad_dom2_H = Constants.c/H_of_x_SI/a/Constants.Mpc/1000
Eta_mat_dom2_H = Eta_rad_dom2_H[rad_mat_equal] + 2*Constants.c*(1/(a*H_of_x_SI) - a[rad_mat_equal]**0.5/(a**1.5*H_of_x_SI))/Constants.Mpc/1000
Eta_DE_dom2_H = Eta_mat_dom2_H[DE_mat_equal] + Constants.c/(H_of_x_SI)*(1/a[DE_mat_equal] - 1/a)/Constants.Mpc/1000

Eta_rad_dom2_Hi = Constants.c/H_of_x_SI_rad_dom/a/Constants.Mpc/1000
Eta_mat_dom2_Hi = Eta_rad_dom2_Hi[rad_mat_equal] + 2*Constants.c*(1/(a*H_of_x_SI_mat_dom) - a[rad_mat_equal]**0.5/(a**1.5*H_of_x_SI_mat_dom))/Constants.Mpc/1000
Eta_DE_dom2_Hi = Eta_mat_dom2_Hi[DE_mat_equal] + Constants.c/(H_of_x_SI_DE_dom)*(1/a[DE_mat_equal] - 1/a)/Constants.Mpc/1000


fig, ax = plt.subplots(1, 2, figsize=(13.5, 6), sharey=True)
ax[0].semilogy(x, eta_of_x, lw=3, c="orange", label=r"$\eta(x)$")
ax[0].semilogy(x[:rad_mat_equal], Eta_rad_dom2_Hi[:rad_mat_equal], lw=4, ls=":", c="r", label=r"$\eta_r(x)$")
ax[0].semilogy(x[rad_mat_equal:DE_mat_equal], Eta_mat_dom2_Hi[rad_mat_equal:DE_mat_equal], lw=3, ls="--", c="b", label=r"$\eta_m(x)$")
ax[0].semilogy(x[DE_mat_equal:], Eta_DE_dom2_Hi[DE_mat_equal:], lw=3, ls="-.", c="purple", label=r"$\eta_\Lambda(x)$")
ax[i].semilogy(x, eta_of_x, lw=3, c="orange", label=r"$\eta(x)$")
ax[i].semilogy(x[:rad_mat_equal], Eta_rad_dom2_H[:rad_mat_equal], lw=4, ls=":", c="r", label=r"$\eta_r(x)$")
ax[i].semilogy(x[rad_mat_equal:DE_mat_equal], Eta_mat_dom2_H[rad_mat_equal:DE_mat_equal], lw=3, ls="--", c="b", label=r"$\eta_m(x)$")
ax[i].semilogy(x[DE_mat_equal:], Eta_DE_dom2_H[DE_mat_equal:], lw=3, ls="-.", c="purple", label=r"$\eta_\Lambda(x)$")
for i in range(2):
    ax[i].fill_between(x, 0, 100, OmegaR>(OmegaCDM+OmegaB), alpha=0.2, color="red")
    ax[i].fill_between(x, 0, 100, (OmegaR<(OmegaCDM+OmegaB) * ((OmegaCDM+OmegaB) > OmegaLambda)), alpha=0.2, color="blue")
    ax[i].fill_between(x, 0, 100, (OmegaCDM+OmegaB) < OmegaLambda, alpha=0.2, color="black")
    ax[i].set_xlabel("x")
    ax[i].set_xlim(-15, 4)
ax[0].set_ylabel(r"$\eta(x)\ - \ [Gpc]$")
ax[0].set_title(r"$\eta(x)$")
ax[1].set_title(r"$\eta(x)\ - \ H_i = H$")
fig.tight_layout()
plt.legend()
fig.savefig("../m1_figs/Eta2.pdf", bbox_inches="tight")
plt.close(fig)

fig, ax = plt.subplots(2, 2, figsize=(13.5, 10))
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
ax[1,1].set_xlim(1.5*z[0], 0.001)
ax[1,1].set_xlabel("z")
ax[1,1].set_ylabel(r"$H(z)$")
ax[1,1].set_title(r"$H \ - \ [km/s/Mpc]$")

fig.tight_layout()
fig.savefig("../m1_figs/H.pdf", bbox_inches="tight")
plt.close(fig)