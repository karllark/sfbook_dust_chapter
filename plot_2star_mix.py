import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl

import astropy.units as u

from dust_extinction.parameter_averages import G16

# setup the plots
fontsize = 16
font = {"size": fontsize}

mpl.rc("font", **font)

mpl.rc("lines", linewidth=2)
mpl.rc("axes", linewidth=2)
mpl.rc("xtick.major", width=2)
mpl.rc("ytick.major", width=2)
mpl.rc("xtick.minor", width=2)
mpl.rc("ytick.minor", width=2)

fig, ax = plt.subplots()

x = np.arange(0.3, 10.0, 0.1) / u.micron

# input extinction curve
ext_model = G16()
ax.plot(1.0 / x, ext_model(x), "k--", label="input extinction")

# stellar intrinsic fluxes
f1 = np.full(len(x), 1.0)
f2 = np.full(len(x), 1.0)

f12 = f1 + f2

av1 = [1.0, 0.5, 0.2, 0.0]
av2 = [2.0, 2.0, 2.0, 2.0]

# get the index for the "V" band
sindxs = np.argsort(np.abs(x - (1.0 / 0.55) / u.micron))

for k in range(len(av1)):

    # attenuated
    rf1 = f1 * ext_model.extinguish(x, Av=av1[k])
    rf2 = f2 * ext_model.extinguish(x, Av=av2[k])

    rf12 = rf1 + rf2

    # attenuation curve
    att12 = -2.5 * np.log10(rf12 / f12)

    label_str = "$A(V)_1 = {:3.1f}$, $A(V)_2 = {:3.1f}$, $Att(V) = {:3.2f}$"
    ax.plot(
        1.0 / x,
        att12 / att12[sindxs[0]],
        label=label_str.format(av1[k], av2[k], att12[sindxs[0]]),
    )

ax.set_xlabel(r"$\lambda$ [$\mu m$]")
ax.set_ylabel(r"$Att(\lambda)/Att(V)$")

ax.set_xscale("log")
ax.set_xlim(1.0 / 10.0, 1.0 / 0.3)

ax.legend(loc="best", fontsize=14)

ax.tick_params(which="major", length=7)
ax.tick_params(which="minor", length=4)

plt.tight_layout()

fig.savefig("2star_mix.pdf")
