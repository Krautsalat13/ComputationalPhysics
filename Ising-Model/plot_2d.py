"""Plot the 2D Ising observables against the exact Onsager theory.

Loads the Monte-Carlo results produced by ``ising_2d.py`` from ``Data/``
(internal energy U, specific heat C and absolute magnetisation |M| per spin, for
N = 10, 50, 100, free and periodic boundaries, N_s = 1e3 / 1e4 sweeps) and plots
each observable against temperature, overlaying the exact infinite-lattice
Onsager result. Figures go to ``Plots/``.

The Onsager curves are for the infinite lattice, so some disagreement near the
critical temperature T_c is expected and physical: a finite lattice rounds off
the transition (|M| leaks above T_c, the specific-heat peak is finite and
slightly shifted) instead of being sharp. The agreement well away from T_c is
what shows the simulation is correct.

Run from the project root:  ``python plot_2d.py``
"""
import matplotlib
matplotlib.use("Agg")   # non-interactive backend
import matplotlib.pyplot as plt
import numpy as np
from scipy.special import ellipk

# plain mathtext labels, no LaTeX needed
plt.rcParams["font.family"] = "serif"
plt.rcParams["font.size"] = "18"


def Color(color_name):
    """RWTH Aachen corporate colour palette, returned as a 0-1 RGB tuple."""
    colors = {
        "RWTH": (0, 84, 159),
        "Schwarz": (0, 0, 0),
        "Petrol": (0, 97, 101),
        "Türkis": (0, 152, 161),
        "Grün": (87, 171, 39),
        "Maigrün": (189, 205, 0),
        "Gelb": (255, 237, 0),
        "Orange": (246, 168, 0),
        "Rot": (204, 7, 30),
        "Magenta": (227, 0, 102),
        "Bordeaux": (161, 16, 53),
        "Violett": (97, 33, 88),
        "Lila": (122, 111, 172),
        "RWTHlight": (142, 186, 229),
        "Bordeauxlight": (205, 139, 135),
    }
    if color_name not in colors:
        raise ValueError(f"Color '{color_name}' not found.")
    r, g, b = colors[color_name]
    return (r / 255, g / 255, b / 255)


# colours for the simulation curves, RWTH and Bordeaux first
SIM_COLORS = [Color(c) for c in ("RWTH", "Bordeaux", "Petrol", "Orange")]

Tc = 2 / np.log(1 + np.sqrt(2))     # critical temperature ~ 2.2692
T_GRID = np.linspace(0.2, 4.0, 40)  # must match T_GRID in ising_2d.py
T_fine = np.linspace(0.2, 4.0, 400) # smooth grid for the theory lines


def M_theory(T):
    """Onsager spontaneous magnetisation per spin (0 above T_c)."""
    return np.where(T < Tc, np.abs(1 - np.sinh(2 / T) ** -4) ** (1 / 8), 0.0)


def U_theory(T):
    """Onsager internal energy per spin."""
    b = 1 / T
    k1 = 2 * np.sinh(2 * b) / np.cosh(2 * b) ** 2
    K = ellipk(k1 ** 2)                 # scipy uses the parameter m = k^2
    return -1 / np.tanh(2 * b) * (1 + (2 / np.pi) * (2 * np.tanh(2 * b) ** 2 - 1) * K)


def C_theory(T):
    """Onsager specific heat per spin, C = dU/dT (diverges at T_c).

    Taken as the exact derivative of the Onsager energy above, which is the
    thermodynamic definition of the specific heat.
    """
    h = 1e-4
    return (U_theory(T + h) - U_theory(T - h)) / (2 * h)


# the four simulation settings shown on each plot
CURVES = [("per", 1000, "periodic, $N_s=10^3$"),
          ("per", 10000, "periodic, $N_s=10^4$"),
          ("free", 1000, "free, $N_s=10^3$"),
          ("free", 10000, "free, $N_s=10^4$")]


def load(obs, N, Ns, bc):
    return np.load(f"Data/2D_{obs}_{N}N_{Ns}NS_{bc}.npy")


def make_plot(obs, N, ylabel, title, theory_fn, clip_to_data=False):
    plt.figure(figsize=(10, 6))
    plt.plot(T_fine, theory_fn(T_fine), "k--", lw=2, label="Onsager theory")
    if obs == "C":
        plt.axvline(Tc, color="grey", ls=":", lw=1.5, label="$T_c$")
    ymax = 0.0
    for i, (bc, Ns, lab) in enumerate(CURVES):
        data = load(obs, N, Ns, bc)
        ymax = max(ymax, data.max())
        ls = "-" if bc == "per" else "--"
        plt.plot(T_GRID, data, color=SIM_COLORS[i], ls=ls, marker="o", ms=4,
                 lw=1.3, label=lab)
    if clip_to_data:
        # the theory specific heat diverges at T_c; cap the axis so the finite
        # simulation peaks are still readable (the theory line runs off the top).
        plt.ylim(0, ymax * 1.4)
    plt.xlabel("temperature $T$")
    plt.ylabel(ylabel)
    plt.title(title)
    plt.grid(ls="--", alpha=0.5)
    plt.legend(fontsize=12)
    plt.savefig(f"Plots/2D_{obs}_{N}.png", dpi=150, bbox_inches="tight")
    plt.close()


for N in [10, 50, 100]:
    make_plot("U", N, r"$U/N^2$", f"Internal energy per spin, N={N}", U_theory)
    make_plot("C", N, r"$C/N^2$", f"Specific heat per spin, N={N}", C_theory,
              clip_to_data=True)
    make_plot("M", N, r"$|M|/N^2$", f"Absolute magnetisation per spin, N={N}",
              M_theory)
