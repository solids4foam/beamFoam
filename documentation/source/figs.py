import numpy as np, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from model import loop_matrix, mono_matrix, rho
plt.rcParams.update({"font.size": 10, "font.family": "serif"})

Om = np.logspace(-3, np.log10(3.4), 300)
fig, ax = plt.subplots(figsize=(6.4, 4.4))
for K, w, c, ls, lab in [(1,0.5,"b","-","Loop, 1 corrector (relaxation 0.5)"),
                         (1,1.0,"c","--","Loop, 1 corrector (no relaxation)"),
                         (2,0.5,"g","-","Loop, 2 correctors"),
                         (4,0.5,"orange","-","Loop, 4 correctors"),
                         (8,0.5,"m","-","Loop, 8 correctors")]:
    g = [rho(loop_matrix(o,K,w))**(2*np.pi/o) - 1 for o in Om]
    ax.loglog(Om, np.maximum(g, 1e-12), color=c, ls=ls, label=lab)
ax.text(1.2e-3, 2.5e-7, "Monolithic: growth exactly 0 for every $\\Omega$ (not visible on the log axis)", color="r", fontsize=8)
for lo, hi, lab, y in [(0.0016,0.019,"base",3e-1),(0.0053,0.062,"lightBody",1e-4),(0.017,0.2,"stiffLine",3e-2),(0.05,0.62,"stiffLight",3e-3)]:
    ax.annotate("", xy=(lo,y), xytext=(hi,y), arrowprops=dict(arrowstyle="<->", color="k", lw=0.8))
    ax.text(np.sqrt(lo*hi), y*1.4, lab, ha="center", fontsize=8)
ax.set_ylim(1e-7, 20)
ax.set_xlabel(r"$\Omega = \omega_n \Delta t$  (line frequency $\times$ time step)")
ax.set_ylabel("Amplitude growth per period $-1$")
ax.grid(which="both", ls=":", lw=0.4)
ax.legend(fontsize=7.5, loc="upper center", bbox_to_anchor=(0.5, -0.17), ncol=3, frameon=False)
plt.tight_layout(); plt.savefig("growth.png", dpi=200)

fig, ax = plt.subplots(figsize=(6.4, 3.0))
k = np.arange(1, 9)
loop = [4.99,4.98,3.00,3.00,2.97,2.87,2.83,2.36]
mono = [3.99,2.00,2.00,2.00,2.00,2.00,2.00,2.00]
ax.bar(k-0.2, loop, 0.4, color="b", label="Loop, 8 correctors (2402 in total)")
ax.bar(k+0.2, mono, 0.4, color="r", label="Monolithic, 8 correctors (1601 in total)")
ax.set_xlabel("PIMPLE outer corrector within the time step")
ax.set_ylabel("Mean beam Newton iterations")
ax.set_xticks(k); ax.set_ylim(0, 6); ax.grid(axis="y", ls=":", lw=0.4); ax.legend(fontsize=8)
plt.tight_layout(); plt.savefig("newton.png", dpi=200)
