#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Plot the 1D Ising-model internal energy and specific heat.

Reads the Monte-Carlo averages tabulated in ``Data/1D/U{10,100,1000}.csv``
(produced by ``1D/1D.py``) and plots internal energy U/N and specific heat
C/N against temperature for chains of N = 10, 100, 1000 spins, overlaying the
exact 1D results U/N = -(N-1)/N * tanh(1/T) and the corresponding C/N. Each
figure is written as PNG to ``Plots/``.

Run from the project root:  ``python 1D/plot.py``
"""

import numpy as np
import matplotlib
matplotlib.use("Agg")   # non-interactive backend
import matplotlib.pyplot as plt
from cycler import cycler

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


# simulation curves cycle through the RWTH palette (RWTH and Bordeaux first);
# the theory curves are drawn in black below.
plt.rcParams["axes.prop_cycle"] = cycler(
    color=[Color(c) for c in ("RWTH", "Bordeaux", "Petrol")])


T = np.linspace(0.2,4.2,21)
N = np.array([10,100,1000])


a, b, c, d= np.loadtxt('Data/1D/U10.csv', delimiter=';')

E10_1000 = np.array([a[::-1],b[::-1]])
E10_10000 = np.array([c[::-1],d[::-1]])


a, b, c, d = np.loadtxt('Data/1D/U100.csv', delimiter=';')

E100_1000 = np.array([a[::-1],b[::-1]])
E100_10000 = np.array([c[::-1],d[::-1]])


a, b, c, d= np.loadtxt('Data/1D/U1000.csv', delimiter=';')

E1000_1000 = np.array([a[::-1],b[::-1]])
E1000_10000 = np.array([c[::-1],d[::-1]])



#10
plt.plot(T, E10_1000[0]/10, label = r"10 Spins $N_{S} = 1000$")
plt.plot(T, E10_10000[0]/10, label = r"10 Spins $N_{S} = 10000$")

plt.plot(T, -(N[0]-1)/N[0] *np.tanh(1/T), color="black", ls="--", label="Theory")

plt.title("Average Energy")
plt.xlabel("Temperature")
plt.ylabel(r"U/N")
#plt.legend(bbox_to_anchor=(1.04, 0.5), loc="center left", borderaxespad=0)
plt.legend()
plt.grid()
plt.tight_layout()
plt.savefig("Plots/1D_U10.png", dpi=150, bbox_inches='tight')

plt.close()


plt.plot(T, 1/T**2*(E10_1000[1]- E10_1000[0]**2)/10, label = r"10 Spins $N_{S} = 1000$")
plt.plot(T, 1/T**2*(E10_10000[1]- E10_10000[0]**2)/10, label = r"10 Spins $N_{S} = 10000$")

plt.plot(T, (N[0]-1)/N[0] *1/(T*np.cosh(1/T))**2, color="black", ls="--", label="Theory")

plt.title("Specific Heat")
plt.xlabel("Temperature")
plt.ylabel(r"C/N")
#plt.legend(bbox_to_anchor=(1.04, 0.5), loc="center left", borderaxespad=0)
plt.legend()
plt.grid()
plt.tight_layout()
plt.savefig("Plots/1D_C10.png", dpi=150, bbox_inches='tight')

plt.close()

#100
plt.plot(T, E100_1000[0]/100, label = r"100 Spins $N_{S} = 1000$")
plt.plot(T, E100_10000[0]/100, label = r"100 Spins $N_{S} = 10000$")

plt.plot(T, -(N[1]-1)/N[1] *np.tanh(1/T), color="black", ls="--", label="Theory")

plt.title("Average Energy")
plt.xlabel("Temperature")
plt.ylabel(r"U/N")
#plt.legend(bbox_to_anchor=(1.04, 0.5), loc="center left", borderaxespad=0)
plt.legend()
plt.grid()
plt.tight_layout()
plt.savefig("Plots/1D_U100.png", dpi=150, bbox_inches='tight')

plt.close()
plt.plot(T, 1/T**2*(E100_1000[1]- E100_1000[0]**2)/100, label = r"100 Spins $N_{S} = 1000$")
plt.plot(T, 1/T**2*(E100_10000[1]- E100_10000[0]**2)/100, label = r"100 Spins $N_{S} = 10000$")

plt.plot(T, (N[1]-1)/N[1] *1/(T*np.cosh(1/T))**2, color="black", ls="--", label="Theory")

plt.title("Specific Heat")
plt.xlabel("Temperature")
plt.ylabel(r"C/N")
#plt.legend(bbox_to_anchor=(1.04, 0.5), loc="center left", borderaxespad=0)
plt.legend()
plt.grid()
plt.tight_layout()
plt.savefig("Plots/1D_C100.png", dpi=150, bbox_inches='tight')

plt.close()


#1000
plt.plot(T, E1000_1000[0]/1000, label = r"1000 Spins $N_{S} = 1000$")
plt.plot(T, E1000_10000[0]/1000, label = r"1000 Spins $N_{S} = 10000$")

plt.plot(T, -(N[2]-1)/N[2] *np.tanh(1/T), color="black", ls="--", label="Theory")

plt.title("Average Energy")
plt.xlabel("Temperature")
plt.ylabel(r"U/N")
#plt.legend(bbox_to_anchor=(1.04, 0.5), loc="center left", borderaxespad=0)
plt.legend()
plt.grid()
plt.tight_layout()
plt.savefig("Plots/1D_U1000.png", dpi=150, bbox_inches='tight')

plt.close()

plt.plot(T, 1/T**2*(E1000_1000[1]- E1000_1000[0]**2)/1000, label = r"10 Spins $N_{S} = 1000$")
plt.plot(T, 1/T**2*(E1000_10000[1]- E1000_10000[0]**2)/1000, label = r"10 Spins $N_{S} = 10000$")

plt.plot(T, (N[2]-1)/N[2] *1/(T*np.cosh(1/T))**2, color="black", ls="--", label="Theory")

plt.title("Specific Heat")
plt.xlabel("Temperature")
plt.ylabel(r"C/N")
#plt.legend(bbox_to_anchor=(1.04, 0.5), loc="center left", borderaxespad=0)
plt.legend()
plt.grid()
plt.tight_layout()
plt.savefig("Plots/1D_C1000.png", dpi=150, bbox_inches='tight')

plt.close()

