import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from ResearchTools import SpecTools

filenames = np.loadtxt("ALLSPEC.txt", dtype="str")

spec = fits.open(filenames[0])[1]


fig = SpecTools.plot_SDSSspec(spec.data['Lam'], spec.data['Lum'], figsize=(8, 3), specType="", title="", xmin=3800, xmax=10000, spec_box_size=10, kind="object", major_tick_space=1000, minor_tick_space=100)

ax = plt.gca()
ax.lines[0].set_linewidth(1.0)
ax.set_ylabel("Luminosity [erg s$^{-1}$ \AA$^{-1}$]")

plt.savefig(f"PLOTS/{filenames[0].split('/')[1].replace('fits', 'pdf')}", dpi=600)
plt.clf()
plt.close(fig)




for filename in filenames:
    spec = fits.open(filename)[1]


    fig = SpecTools.plot_SDSSspec(spec.data['Lam'], spec.data['Lum'], figsize=(8, 3), specType="", title="", xmin=3800, xmax=10000, spec_box_size=10, kind="object", major_tick_space=1000, minor_tick_space=100)
    
    ax = plt.gca()
    ax.lines[0].set_linewidth(1.0)
    ax.set_ylabel("Luminosity [erg s$^{-1}$ \AA$^{-1}$]")

    plt.savefig(f"PLOTS/{filename.split('/')[1].replace('fits', 'pdf')}", dpi=600)
    plt.clf()
    plt.close(fig)