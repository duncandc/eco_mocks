"""
define galaxy properties and build PDF for ECO survey
"""

from __future__ import (division, print_function, absolute_import,
                        unicode_literals)
import numpy as np
import os
#from gaussian_kde.gaussian_kde import gaussian_kde
from astropy.io import ascii
from scipy.stats import norm
from astropy.table import Table

# load ECO data
data_path = os.path.dirname(os.path.realpath(__file__))+'/data/'
data = ascii.read(data_path + "eco_data.csv")

# trim data
mask = (data['LOGMSTAR'] > 0.0)
data = data[mask]

# add random noise to stellar mass 
N = len(data)
f = 1.0+(np.random.random(N)*2.0-1.0)/10.0
f = np.random.normal(1, 0.1, N)

# process values
h = 0.7
mstar = f*(10.0**data['LOGMSTAR'])*h**2  # stellar mass
mgas = (10.0**data['LOGMGAS'])*h**2  # gas mass
mbary = mgas + mstar  # baryonic mass
fsmgr = np.log10(data['MEANFSMGR'])  # fractional stellar mass growth rate
abs_rmag = data['ABSRMAG']  # absolute r-band magnitude
b_to_a = data['B_A']  # axis ratios
dgr = data['DGR']  # (g-r) color gradient
morph_type = data['MORPHEL']  # quantitative morphological type
u_minus_r = data['MODELU_RCORR']  # (u-r) color
g_minus_r = data['MODELG_RCORR']  # (g-r) color

# add some derived quantities
abs_umag = u_minus_r + abs_rmag  # absolute u-band magnitude
abs_gmag = g_minus_r + abs_rmag  # absolute g-band magnitude
fgas = mgas/(mbary)  # gas fraction

# create b_to_a variable with better properties
bijected_b_to_a = np.log(b_to_a/(1.0-b_to_a))
# replace b_to_a==0 with a small value
mask = (b_to_a == 0.0)
min_b_to_a = np.min(b_to_a[~mask])
bijected_b_to_a[mask] = np.log(min_b_to_a/(1.0-min_b_to_a))

# create fgas variable with better properties
bijected_fgas = np.log(fgas/(1.0-fgas))

# convert morphological type to numerical value
# L=0, E=1, NA=-1
num_morph_type = np.zeros(len(data))
mask = (morph_type == 'E')
num_morph_type[mask] = 1
mask = (morph_type == 'NA')
num_morph_type[mask] = -1

# trim dgr dist
mask = np.array((dgr > -1) & (dgr < 1))
params = norm.fit(dgr[mask])
model = norm(loc=params[0], scale=params[1])
dgr[~mask] = model.rvs(size=np.sum(~mask))

eco_values = np.vstack([mstar, mgas, mbary,
                       abs_rmag, abs_umag, abs_gmag,
                       u_minus_r, g_minus_r, fsmgr, fgas, bijected_fgas,
                       b_to_a, bijected_b_to_a,
                       num_morph_type, dgr]).T
eco_names = ['stellar_mass', 'gas_mass', 'baryonic_mass',
             'abs_rmag', 'abs_umag', 'abs_gmag',
             'u_minus_r', 'g_minus_r', 'fsmgr', 'fgas', 'bijected_fgas',
             'b_to_a', 'bijected_b_to_a',
             'num_morph_type', 'dgr']

eco_table = Table(eco_values, names=eco_names)

"""
# build KDE model for galaxy property PDF
eco_kernel = gaussian_kde(eco_values)

# build KDE 2D kernels
values = np.vstack([mstar, fsmgr])
mstar_fsmgr_kernel = gaussian_kde(values)


def rvs_fsmgr_mstar(x):
    x = np.log10(x)
    result = mstar_fsmgr_kernel.conditional_sample(x, [True, False])
    return result.flatten()
"""


