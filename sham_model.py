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
from halotools.utils import fuzzy_digitize
from scipy import interpolate

__author__ = ('Duncan Campbell')
__all__ = ('DeconvolveSHAM', 'CAMGalProp')


class DeconvolveSHAM(object):
    """
    class to assign stellar mass, or some other primary galaxy property,
    based on the de-convolution SHAM method, based on Yau-Yun Mao's
    implementation of Peter Behroozi's deconvolution package.
    """
    def __init__(self,
                 stellar_mass_function,
                 prim_haloprop='halo_mpeak',
                 prim_galprop='stellar_mass',
                 scatter=0.0,
                 gal_log_prop=True,
                 gal_min_mass=2.0,
                 gal_max_mass=12.5,
                 gal_nbins=200,
                 **kwargs):
        """
        Parameters
        ----------
        stellar_mass_function : callable object

        prim_haloprop : string, optional
            key indicating halo property used to assign the primary
            galaxy property

        prim_galprop : string, optional
            key indicating primary galaxy property

        gal_log_prop : boolean, optional
            boolean indicating whther the scatter in the primary
            galaxy property is modelled as a log-normal.

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
            self.redshift = 0.0

        self.prim_haloprop = prim_haloprop
        self.prim_galprop = prim_galprop
        self._galprop_dtypes_to_allocate = np.dtype([(str(self.prim_galprop), 'f4')])
        self._mock_generation_calling_sequence = ['inherit_halocat_properties', 'assign_stellar_mass']
        self.list_of_haloprops_needed = [prim_haloprop]

        # set default box size.
        if 'Lbox' in kwargs.keys():
            self._Lbox = kwargs['Lbox']
        else:
            self._Lbox = np.inf
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
        Inherit the box size during mock population
        """
        try:
            Lbox = kwargs['Lbox']
            self._Lbox = Lbox
        except KeyError:
            print("Error automatically detecting Lbox.")

    def deconvolve_scatter(self, scatter):
        """
        """
        # tabulate stellar mass function
        if self.gal_log_prop is True:
            msample = np.logspace(self.gal_min_mass, self.gal_max_mass, self.gal_nbins)
            nsample = self._stellar_mass_function(msample)
            self.af = AbundanceFunction(np.log10(msample), nsample, faint_end_first=True)
        else:
            msample = np.linspace(self.gal_min_mass, self.gal_max_mass, self.gal_nbins)
            nsample = self._stellar_mass_function(msample)
            self.af = AbundanceFunction(msample, nsample, faint_end_first=True)

        remainder = self.af.deconvolute(scatter, 100)
        return remainder

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
        if len(Lbox) == 3:
            Lbox = (np.prod(Lbox))**(1.0/3.0)
        else:
            Lbox = Lbox[0]

        nd_halos = calc_number_densities(table[self.prim_haloprop], Lbox)
        mstar = self.af.match(nd_halos, self.param_dict['scatter'])

        if self.gal_log_prop is True:
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
                 prim_galprop,
                 prim_galprop_bins,
                 secondary_galprop='ssfr',
                 secondary_haloprop='halo_half_mass_scale',
                 rho=1.0,
                 conditional_rvs=None,
                 **kwargs):
        """
        Parameters
        ----------

        rho : float
             correlation strength
        """

        if 'redshift' in list(kwargs.keys()):
            self.redshift = kwargs['redshift']
        else:
            self.redshift = 0.0

        self.prim_galprop = prim_galprop
        self.secondary_haloprop = secondary_haloprop
        self.secondary_galprop = secondary_galprop
        self._galprop_dtypes_to_allocate = np.dtype([(str(self.secondary_galprop), 'f4')])
        self._mock_generation_calling_sequence = ['assign_secondary_galprop']
        self.list_of_haloprops_needed = [secondary_haloprop]

        # store model parametrs
        self.param_dict = ({'rho': rho})
        self.prim_galprop_bins = prim_galprop_bins

        # set rvs function
        if conditional_rvs is None:
            self.conditional_rvs = self.default_conditional_rvs
        else:
            self.conditional_rvs = conditional_rvs

    def assign_secondary_galprop(self, **kwargs):
        """
        Assign galaxy properties using CAM technique
        """

        table = kwargs['table']

        binning_inds = fuzzy_digitize(table[self.prim_galprop],
                                      self.prim_galprop_bins)
        binning_inds = np.digitize(table[self.prim_galprop],
                                   bins=self.prim_galprop_bins)

        # loop through bins
        result = np.zeros(len(table))
        inds = np.arange(0, len(table)).astype('int')
        for i in range(0, len(self.prim_galprop_bins)):
            # mask of galaxies in the bin
            mask = (binning_inds == i)

            if (np.sum(mask) == 0.0):
                continue

            # sort galaxies by the value of the secondary halo property
            sort_inds = np.argsort(table[self.secondary_haloprop][mask])

            # apply index shuffling for required correlation strength
            sigma = rho_to_sigma_function(np.fabs(self.param_dict['rho']))
            sort_inds = scatter_ranks(sort_inds, sigma)
            # flip the direction of correlation if rho<0.0
            if self.param_dict['rho'] < 0:
                sort_inds = sort_inds[::-1]

            # draw from secondary galaxy property distribution
            y = np.sort(self.conditional_rvs(table[self.prim_galprop][mask]))
            inds_in_bin = inds[mask]
            result[inds_in_bin[sort_inds]] = y

        table[str(self.secondary_galprop)] = result

    def default_conditional_rvs(self, x):
        """
        """
        return np.random.random(size=len(x))


class HaloProps(object):
    """
    class to carry over halo properties to galaxy mock
    """
    def __init__(self,
                 haloprop_keys=['halo_mpeak', 'halo_vpeak', 'halo_half_mass_scale'],
                 **kwargs):
        """
        Parameters
        ----------
        haloprop_keys : list
        """

        self._mock_generation_calling_sequence = []
        self._galprop_dtypes_to_allocate = np.dtype([])
        self.list_of_haloprops_needed = haloprop_keys


def scatter_ranks(arr, sigma):
    """
    Scatter the index of values in an array.

    Parameters
    ----------
    arr : array_like
        array of values to scatter

    sigma : array_like
        scatter relative to len(arr)

    Returns
    -------
    scatter_array : numpy.array
        array with same values as ``arr``, but the locations of those values
        have been scatter.
    """

    sigma = np.atleast_1d(sigma)
    if len(sigma) == 1:
        sigma = np.repeat(sigma, len(arr))
    elif len(sigma) != len(arr):
        raise ValueError("sigma must be same length as ``arr``.")

    # get array of indicies before scattering
    N = len(arr)
    inds = np.arange(0, N)

    mask = (sigma > 1000.0)
    sigma[mask] = 1000.0

    # get array to scatter positions
    mask = (sigma > 0.0)
    dn = np.zeros(N)
    dn[mask] = np.random.normal(0.0, sigma[mask]*N)

    # get array of new indicies
    new_inds = inds + dn
    new_inds = np.argsort(new_inds, kind='mergesort')

    return arr[new_inds]


def sigma_to_rho_function(x, x1=0.64, x2=0.30, alpha=-1.07, beta=-1.97):
    """
    """

    return (2.0-1.0*np.exp(-1.0*(x/x1)**alpha)-1.0/(1.0+(x/x2)**beta))/2.0


def rho_to_sigma_function(rho):
    """
    """

    x = np.logspace(-2, 1, 30)
    y = sigma_to_rho_function(x)
    f_yx = interpolate.interp1d(y, x, fill_value="extrapolate")

    sigma = f_yx(rho)

    mask = (sigma < 0.0)
    sigma[mask] = 0.0

    return sigma








