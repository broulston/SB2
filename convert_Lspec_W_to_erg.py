import numpy as np
from astropy.table import Table
import astropy.units as u
from tqdm import tqdm

indvid_dir = "SB2_IndividualSpec/"
invid_spec_list = np.loadtxt(f"{indvid_dir}ALLSPEC.txt", dtype="U")

for ii, filename in tqdm(enumerate(invid_spec_list)):
    spec = Table.read(f"{indvid_dir}{filename}")

    spec['Lum'] *= 1e7
    spec['LumErr'] *= 1e7

    spec.write(f"{indvid_dir}{filename}", format='fits', overwrite=True)

SB2_lists = ['list_of_SB2_04-30-2020.txt', 'list_of_SB2_05-08-2020.txt']

for SB2list in SB2_lists:
    spec_list = np.loadtxt(SB2list, dtype="U")

    for ii, filename in enumerate(tqdm(spec_list)):
        spec = Table.read(f"{filename}")

        spec['Lum'] *= 1e7
        spec['LumErr'] *= 1e7

        spec.write(f"{filename}", format='fits', overwrite=True)