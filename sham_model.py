"""
Module containing classes used in SHAM models
"""

from __future__ import (division, print_function, absolute_import)
import numpy as np
from astropy.table import Table
from halotools.empirical_models.smhm_models.smhm_helpers import safely_retrieve_redshift
from halotools.empirical_models import model_helpers as model_helpers
from halotools.empirical_models.component_model_templates import PrimGalpropModel
from AbundanceMatching import AbundanceFunction, calc_number_densities

__author__=('Duncan Campbell')
__all__=('DeconvolveSHAM', 'CAMGalProp')


class DeconvolveSHAM(object):
    """
    class to assign stellar mass, or some other primary galaxy property, 
    based on the de-convolution SHAM method, based on Yau-Yun Mao's
    implementation of Peter Behroozi's deconvolution package.
    """
    def __init__(self,
                 stellar_mass_function,
                 prim_haloprop = 'halo_mpeak',
                 prim_galprop = 'stellar_mass',
                 scatter = 0.0,
                 gal_log_prop = True,
                 gal_min_mass = 5.0,
                 gal_max_mass = 12.5,
                 gal_nbins = 100,
                 Lbox = 130.0,
                 **kwargs):
        """
        Parameters
        ----------
        stellar_mass_function : callable object

        prim_haloprop : string, optional
            key indicating halo property used to assign the primary galaxy property

        prim_galprop : string, optional
            key indicating primary galaxy property
        
        gal_log_prop : boolean, optional
            boolean indicating whther the scatter in the primary galaxy property
            is modelled as a log-normal.

        gal_min_mass, gal_max_mass, gal_nbins: float, float, int
            parameters used to bin and span galaxy abundance function. 
            If gal_log_prop is True, these should also be logged, 
            and the spacing will be even in log.

        redshift : float
             redshit for which the model is applicable
        """

        if 'redshift' in list(kwargs.keys()):
            self.redshift = kwargs['redshift']
        else:
            self.redshift=0.0
        
        self.prim_haloprop = prim_haloprop
        self.prim_galprop = prim_galprop
        self._galprop_dtypes_to_allocate = np.dtype([(str(self.prim_galprop), 'f4')])
        self._mock_generation_calling_sequence = ['inherit_halocat_properties', 'assign_stellar_mass']
        self.list_of_haloprops_needed = [prim_haloprop]

        # set default box size.
        self._Lbox = Lbox

        # update Lbox if a halo catalog object is passed.
        self._additional_kwargs_dict = dict(inherit_halocat_properties=['Lbox'])
        
        # store model parametrs
        self.gal_log_prop = gal_log_prop
        self.param_dict = ({'scatter': scatter})
        self.gal_min_mass = gal_min_mass
        self.gal_max_mass = gal_max_mass
        self.gal_nbins = gal_nbins
        
        # deconvolve abundance function
        self._stellar_mass_function = stellar_mass_function
        self.deconvolve_scatter(self.param_dict['scatter'])

    def inherit_halocat_properties(self, seed=None, **kwargs):
        """
        inherit the box size during mock population
        """
        try:
            Lbox = kwargs['Lbox']
            self._Lbox = Lbox
        except KeyError:
            print("Lbox not found in halocat.")
    
    def deconvolve_scatter(self, scatter):
        """
        """
        # tabulate stellar mass function
        if self.gal_log_prop == True:
            msample = np.logspace(self.gal_min_mass, self.gal_max_mass, self.gal_nbins)
            nsample = self._stellar_mass_function(msample)
            self.af = AbundanceFunction(np.log10(msample), nsample, faint_end_first=True)
        else:
            msample = np.linspace(self.gal_min_mass, self.gal_max_mass, self.gal_nbins)
            nsample = self._stellar_mass_function(msample)
            self.af = AbundanceFunction(msample, nsample, faint_end_first=True)
        
        remainder = self.af.deconvolute(scatter, 100)

    def assign_stellar_mass(self,  **kwargs):
        """
        assign stellar mass
        """

        if 'table' in kwargs.keys():
            table = kwargs['table']
            try:
                Lbox = kwargs['Lbox']
            except KeyError:
                Lbox = self._Lbox
        else:
            try:
                Lbox = kwargs['Lbox']
            except KeyError:
                Lbox = self._Lbox

        Lbox = np.atleast_1d(Lbox)
        if len(Lbox)==3:
            Lbox = (np.prod(Lbox))**(1.0/3.0)
        else:
            Lbox = Lbox[0]

        nd_halos = calc_number_densities(table[self.prim_haloprop], Lbox)
        mstar = self.af.match(nd_halos, self.param_dict['scatter'])

        if self.gal_log_prop==True:
            mstar = 10.0**mstar
            table[self.prim_galprop] = mstar
        else:
            table[self.prim_galprop] = mstar

        if 'table' in kwargs.keys():
            return table
        else:
            return mstar


class CAMGalProp(object):
    """
    A class for assigning a secondary galaxy property
    based on a secondary halo property using the conditional
    abundance matching (CAM) techinque.
    """
    def __init__(self,
                 conditonal_rvs,
                 prim_galprop,
                 prim_galprop_bins,
                 secondary_haloprop = 'halo_halfmass_scale',
                 secondary_galprop = 'ssfr',
                 rho = 1.0,
                 **kwargs):
        """
        Parameters
        ----------
        """

        if 'redshift' in list(kwargs.keys()):
            self.redshift = kwargs['redshift']
        else:
            self.redshift=0.0
        
        self.secondary_haloprop = secondary_haloprop
        self.secondary_galprop = secondary_galprop
        self._galprop_dtypes_to_allocate = np.dtype([(str(self.secondary_galprop), 'f4')])
        self._mock_generation_calling_sequence = ['assign_secondary_galprop']
        self.list_of_haloprops_needed = [prim_haloprop, secondary_haloprop]

        # store model parametrs
        self.param_dict = ({'rho': rho})
        self.prim_galprop_bins = prim_galprop_bins
        
        # deconvolve abundance function
        self._stellar_mass_function = stellar_mass_function
        self.deconvolve_scatter(self.param_dict['scatter'])

    def assign_secondary_galprop(self):
        """
        Assign galaxy properties using CAM technique
        """





