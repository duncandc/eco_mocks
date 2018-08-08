"""
Define Fiducial Models for ECO Mocks
"""

import os
import numpy as np

# set default models
from halotools.empirical_models import SubhaloModelFactory
from eco_mocks.sham_model import DeconvolveSHAM, HaloProps
from eco_mocks.galaxy_abundance_functions import Eckert_2016_phi

Lbox = np.array([130.0,130.0,130.0])  # Vishnu box size

####################################################################################

phi = Eckert_2016_phi(sample='RESOLVE-B SMF double')
sm_model =  DeconvolveSHAM(stellar_mass_function = phi, scatter=0.15, prim_haloprop='halo_vpeak', Lbox=Lbox)
additional_halo_properties = HaloProps()
model_1 = SubhaloModelFactory(stellar_mass = sm_model, haloprops = additional_halo_properties)

####################################################################################

phi = Eckert_2016_phi(sample='RESOLVE-B BMF double')
sm_model =  DeconvolveSHAM(stellar_mass_function = phi, scatter=0.15, prim_haloprop='halo_vpeak', Lbox=Lbox)
additional_halo_properties = HaloProps()
model_2 = SubhaloModelFactory(baryonic_mass = sm_model, haloprops = additional_halo_properties)

####################################################################################