"""
Define Fiducial Models for ECO Mocks
"""

import os
import numpy as np

# set default models
from halotools.empirical_models import SubhaloModelFactory
from eco_mocks.sham_model import DeconvolveSHAM, HaloProps, CAMGalProp
from eco_mocks.galaxy_abundance_functions import Eckert_2016_phi
from eco_mocks.galaxy_secondary_functions import prim_prop_nn

# some parameters to set for all models
Lbox = np.array([130.0,130.0,130.0])  # Vishnu box size
mstar_bins = 10**np.arange(5,12.5,0.1)  # stellar mass bins for CAM
mbary_bins = 10**np.arange(5,12.5,0.1)  # baryonic mass bins for CAM

####################################################################################
#### stellar mass based models
####################################################################################

#define galaxy selection function
def galaxy_selection_func(table, min_mass=10**8.0, max_mass=np.inf, prim_gal_prop='stellar_mass'):
    """
    define which galaxies should appear in the mock galaxy table

    Parameters
    ----------
    table : astropy.table object
        table containing mock galaxies
    """

    mask = (table[prim_gal_prop] >= min_mass) & (table[prim_gal_prop] < max_mass)
    return mask


####################################################################################

phi = Eckert_2016_phi(sample='RESOLVE-B SMF double')
sm_model =  DeconvolveSHAM(stellar_mass_function = phi, scatter=0.15,
                           prim_galprop='stellar_mass', prim_haloprop='halo_vpeak', Lbox=Lbox)
additional_halo_properties = HaloProps()
color_dist = prim_prop_nn(prim_prop='stellar_mass', sec_prop='u_minus_r')
color_model = CAMGalProp('stellar_mass', mstar_bins, rho=-1.0,
                         secondary_haloprop='halo_half_mass_scale',
                         secondary_galprop='u_minus_r',
                         conditional_rvs=color_dist.rvs)
model_1a = SubhaloModelFactory(stellar_mass = sm_model,
                               haloprops = additional_halo_properties,
                               galaxy_color = color_model,
                               galaxy_selection_func = galaxy_selection_func,
                               model_feature_calling_sequence = ('stellar_mass', 'galaxy_color', 'haloprops'))

####################################################################################

phi = Eckert_2016_phi(sample='RESOLVE-B SMF double')
sm_model =  DeconvolveSHAM(stellar_mass_function = phi, scatter=0.15,
                           prim_galprop='stellar_mass', prim_haloprop='halo_vpeak', Lbox=Lbox)
additional_halo_properties = HaloProps()
fsmgr_dist = prim_prop_nn(prim_prop='stellar_mass', sec_prop='fsmgr')
fsmgr_model = CAMGalProp('stellar_mass', mstar_bins, rho=1.0,
                         secondary_haloprop='halo_half_mass_scale',
                         secondary_galprop='fsmgr',
                         conditional_rvs=fsmgr_dist.rvs)
model_1b = SubhaloModelFactory(stellar_mass = sm_model,
                               haloprops = additional_halo_properties,
                               galaxy_fsmgr = fsmgr_model,
                               galaxy_selection_func = galaxy_selection_func,
                               model_feature_calling_sequence = ('stellar_mass', 'galaxy_fsmgr', 'haloprops'))


####################################################################################

phi = Eckert_2016_phi(sample='RESOLVE-B SMF double')
sm_model =  DeconvolveSHAM(stellar_mass_function = phi, scatter=0.15,
                           prim_galprop='stellar_mass', prim_haloprop='halo_vpeak', Lbox=Lbox)
additional_halo_properties = HaloProps()
fgas_dist = prim_prop_nn(prim_prop='stellar_mass', sec_prop='fgas')
fgas_model = CAMGalProp('stellar_mass', mstar_bins, rho=1.0,
                         secondary_haloprop='halo_half_mass_scale',
                         secondary_galprop='fgas',
                         conditional_rvs=fgas_dist.rvs)
model_1c = SubhaloModelFactory(stellar_mass = sm_model,
                               haloprops = additional_halo_properties,
                               galaxy_fgas = fgas_model,
                               galaxy_selection_func = galaxy_selection_func,
                               model_feature_calling_sequence = ('stellar_mass', 'galaxy_fgas', 'haloprops'))

####################################################################################

phi = Eckert_2016_phi(sample='RESOLVE-B SMF double')
sm_model =  DeconvolveSHAM(stellar_mass_function = phi, scatter=0.15,
                           prim_galprop='stellar_mass', prim_haloprop='halo_vpeak', Lbox=Lbox)
additional_halo_properties = HaloProps()
color_dist = prim_prop_nn(prim_prop='stellar_mass', sec_prop='u_minus_r')
color_model = CAMGalProp('stellar_mass', mstar_bins, rho=-1.0,
                         secondary_haloprop='halo_acc_rate_100myr',
                         secondary_galprop='u_minus_r',
                         conditional_rvs=color_dist.rvs)
model_1d = SubhaloModelFactory(stellar_mass = sm_model,
                               haloprops = additional_halo_properties,
                               galaxy_color = color_model,
                               galaxy_selection_func = galaxy_selection_func,
                               model_feature_calling_sequence = ('stellar_mass', 'galaxy_color', 'haloprops'))

####################################################################################

phi = Eckert_2016_phi(sample='RESOLVE-B SMF double')
sm_model =  DeconvolveSHAM(stellar_mass_function = phi, scatter=0.15,
                           prim_galprop='stellar_mass', prim_haloprop='halo_vpeak', Lbox=Lbox)
additional_halo_properties = HaloProps()
fsmgr_dist = prim_prop_nn(prim_prop='stellar_mass', sec_prop='fsmgr')
fsmgr_model = CAMGalProp('stellar_mass', mstar_bins, rho=1.0,
                         secondary_haloprop='halo_acc_rate_100myr',
                         secondary_galprop='fsmgr',
                         conditional_rvs=fsmgr_dist.rvs)
model_1e = SubhaloModelFactory(stellar_mass = sm_model,
                               haloprops = additional_halo_properties,
                               galaxy_fsmgr = fsmgr_model,
                               galaxy_selection_func = galaxy_selection_func,
                               model_feature_calling_sequence = ('stellar_mass', 'galaxy_fsmgr', 'haloprops'))


####################################################################################

phi = Eckert_2016_phi(sample='RESOLVE-B SMF double')
sm_model =  DeconvolveSHAM(stellar_mass_function = phi, scatter=0.15,
                           prim_galprop='stellar_mass', prim_haloprop='halo_vpeak', Lbox=Lbox)
additional_halo_properties = HaloProps()
fgas_dist = prim_prop_nn(prim_prop='stellar_mass', sec_prop='fgas')
fgas_model = CAMGalProp('stellar_mass', mstar_bins, rho=1.0,
                         secondary_haloprop='halo_acc_rate_100myr',
                         secondary_galprop='fgas',
                         conditional_rvs=fgas_dist.rvs)
model_1f = SubhaloModelFactory(stellar_mass = sm_model,
                               haloprops = additional_halo_properties,
                               galaxy_fgas = fgas_model,
                               galaxy_selection_func = galaxy_selection_func,
                               model_feature_calling_sequence = ('stellar_mass', 'galaxy_fgas', 'haloprops'))



####################################################################################
#### baryonic mass based models
####################################################################################

#define galaxy selection function
def galaxy_selection_func(table, min_mass=10**8.5, max_mass=np.inf, prim_gal_prop='baryonic_mass'):
    """
    define which galaxies should appear in the mock galaxy table

    Parameters
    ----------
    table : astropy.table object
        table containing mock galaxies
    """

    mask = (table[prim_gal_prop] >= min_mass) & (table[prim_gal_prop] < max_mass)
    return mask

####################################################################################

phi = Eckert_2016_phi(sample='RESOLVE-B BMF double')
sm_model =  DeconvolveSHAM(stellar_mass_function = phi, scatter=0.15,
                           prim_galprop='baryonic_mass', prim_haloprop='halo_vpeak', Lbox=Lbox)
additional_halo_properties = HaloProps()
color_dist = prim_prop_nn(prim_prop='baryonic_mass', sec_prop='u_minus_r')
color_model = CAMGalProp('baryonic_mass', mbary_bins, rho=-1,
                         secondary_haloprop='halo_half_mass_scale',
                         secondary_galprop='u_minus_r',
                         conditional_rvs=color_dist.rvs)
model_2a = SubhaloModelFactory(baryonic_mass = sm_model,
                               haloprops = additional_halo_properties,
                               galaxy_color=color_model,
                               galaxy_selection_func = galaxy_selection_func,
                               model_feature_calling_sequence = ('baryonic_mass', 'galaxy_color', 'haloprops'))

####################################################################################

phi = Eckert_2016_phi(sample='RESOLVE-B BMF double')
sm_model =  DeconvolveSHAM(stellar_mass_function = phi, scatter=0.15,
                           prim_galprop='baryonic_mass', prim_haloprop='halo_vpeak', Lbox=Lbox)
additional_halo_properties = HaloProps()
fsmgr_dist = prim_prop_nn(prim_prop='baryonic_mass', sec_prop='fsmgr')
fsmgr_model = CAMGalProp('baryonic_mass', mstar_bins, rho=1.0,
                         secondary_haloprop='halo_half_mass_scale',
                         secondary_galprop='fsmgr',
                         conditional_rvs=fsmgr_dist.rvs)
model_2b = SubhaloModelFactory(baryonic_mass = sm_model,
                               haloprops = additional_halo_properties,
                               galaxy_fsmgr = fsmgr_model,
                               galaxy_selection_func = galaxy_selection_func,
                               model_feature_calling_sequence = ('baryonic_mass', 'galaxy_fsmgr', 'haloprops'))


####################################################################################

phi = Eckert_2016_phi(sample='RESOLVE-B BMF double')
sm_model =  DeconvolveSHAM(stellar_mass_function = phi, scatter=0.15,
                           prim_galprop='baryonic_mass', prim_haloprop='halo_vpeak', Lbox=Lbox)
additional_halo_properties = HaloProps()
fgas_dist = prim_prop_nn(prim_prop='baryonic_mass', sec_prop='fgas')
fgas_model = CAMGalProp('baryonic_mass', mstar_bins, rho=1.0,
                         secondary_haloprop='halo_half_mass_scale',
                         secondary_galprop='fgas',
                         conditional_rvs=fgas_dist.rvs)
model_2c = SubhaloModelFactory(baryonic_mass = sm_model,
                               haloprops = additional_halo_properties,
                               galaxy_fgas = fgas_model,
                               galaxy_selection_func = galaxy_selection_func,
                               model_feature_calling_sequence = ('baryonic_mass', 'galaxy_fgas', 'haloprops'))

####################################################################################

phi = Eckert_2016_phi(sample='RESOLVE-B BMF double')
sm_model =  DeconvolveSHAM(stellar_mass_function = phi, scatter=0.15,
                           prim_galprop='baryonic_mass', prim_haloprop='halo_vpeak', Lbox=Lbox)
additional_halo_properties = HaloProps()
color_dist = prim_prop_nn(prim_prop='baryonic_mass', sec_prop='u_minus_r')
color_model = CAMGalProp('baryonic_mass', mbary_bins, rho=-1,
                         secondary_haloprop='halo_acc_rate_100myr',
                         secondary_galprop='u_minus_r',
                         conditional_rvs=color_dist.rvs)
model_2d = SubhaloModelFactory(baryonic_mass = sm_model,
                               haloprops = additional_halo_properties,
                               galaxy_color = color_model,
                               galaxy_selection_func = galaxy_selection_func,
                               model_feature_calling_sequence = ('baryonic_mass', 'galaxy_color', 'haloprops'))

####################################################################################

phi = Eckert_2016_phi(sample='RESOLVE-B BMF double')
sm_model =  DeconvolveSHAM(stellar_mass_function = phi, scatter=0.15,
                           prim_galprop='baryonic_mass', prim_haloprop='halo_vpeak', Lbox=Lbox)
additional_halo_properties = HaloProps()
fsmgr_dist = prim_prop_nn(prim_prop='baryonic_mass', sec_prop='fsmgr')
fsmgr_model = CAMGalProp('baryonic_mass', mstar_bins, rho=1.0,
                         secondary_haloprop='halo_acc_rate_100myr',
                         secondary_galprop='fsmgr',
                         conditional_rvs=fsmgr_dist.rvs)
model_2e = SubhaloModelFactory(baryonic_mass = sm_model,
                               haloprops = additional_halo_properties,
                               galaxy_fsmgr=fsmgr_model,
                               galaxy_selection_func = galaxy_selection_func,
                               model_feature_calling_sequence = ('baryonic_mass', 'galaxy_fsmgr', 'haloprops'))


####################################################################################

phi = Eckert_2016_phi(sample='RESOLVE-B BMF double')
sm_model =  DeconvolveSHAM(stellar_mass_function = phi, scatter=0.15,
                           prim_galprop='baryonic_mass', prim_haloprop='halo_vpeak', Lbox=Lbox)
additional_halo_properties = HaloProps()
fgas_dist = prim_prop_nn(prim_prop='baryonic_mass', sec_prop='fgas')
fgas_model = CAMGalProp('baryonic_mass', mstar_bins, rho=1.0,
                         secondary_haloprop='halo_acc_rate_100myr',
                         secondary_galprop='fgas',
                         conditional_rvs=fgas_dist.rvs)
model_2f = SubhaloModelFactory(baryonic_mass = sm_model,
                               haloprops = additional_halo_properties,
                               galaxy_fgas=fgas_model,
                               galaxy_selection_func = galaxy_selection_func,
                               model_feature_calling_sequence = ('baryonic_mass', 'galaxy_fgas', 'haloprops'))
