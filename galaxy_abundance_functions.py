# -*- coding: utf-8 -*-

"""
callable stellar mass functions from the literature
"""

from __future__ import (division, print_function, absolute_import, unicode_literals)
import numpy as np

from astropy.table import Table
from astropy.modeling.models import custom_model

__all__ = ['Eckert_2016_phi']


class Eckert_2016_phi(object):
    """
    stellar mass function from Eckert et al. 2016, arXiv:1604.03957
    """
    
    def __init__(self, sample='RESOLVE-B SMF single', **kwargs):
        """
        Paramaters
        ==========
        sample : string
            string indicting which fit to use
        """
        
        self.littleh = 0.7
        
        # parameters from table 3
        self.set_params(**kwargs)
        
        #define components of double Schechter function
        s1 = Log_Schechter(phi0=self.phi1, x0=self.x1, alpha=self.alpha1)
        s2 = Log_Schechter(phi0=self.phi2, x0=self.x2, alpha=self.alpha2)
        
        #create piecewise model
        self.s = s1 + s2

    def set_params(self, **kwargs):
        """
        Set the parameters for model.
        The values are taken from table 3 in Eckert et al. 2016, arXiv:1604.03957
        """

        if 'sample' not in kwargs.keys():
            sample = 'RESOLVE-B SMF single'
        else:
            sample = kwargs['sample']

        if sample == 'RESOLVE-B SMF single':
            self.params = {'M1': 11.25,
                           'phi1': 4.47*10**3,
                           'alpha1': -1.28,
                           'M2': 11.25,
                           'phi2': 0.0,
                           'alpha2': 0.0,
                           }
        elif sample == 'RESOLVE-B SMF double':
            self.params = {'M1': 10.87,
                           'phi1': 9.00*10**3,
                           'alpha1': -0.52,
                           'M2': 10.87,
                           'phi2': 3.25,
                           'alpha2':-1.38,
                           }
        elif sample == 'ECO SMF single':
            self.params = {'M1': 10.92,
                           'phi1': 5.95*10**3,
                           'alpha1': -1.19,
                           'M2': 10.87,
                           'phi2': 0.0,
                           'alpha2':0.0,
                           }
        elif sample == 'ECO SMF double':
            self.params = {'M1': 10.87,
                           'phi1': 3.44*10**3,
                           'alpha1': -0.91,
                           'M2': 10.87,
                           'phi2': 3.62,
                           'alpha2':-1.26,
                           }
        elif sample == 'RESOLVE-B BMF single':
            self.params = {'M1': 10.11,
                           'phi1': 6.93*10**3,
                           'alpha1': -1.30,
                           'M2': 10.11,
                           'phi2': 0.0,
                           'alpha2':0.0,
                           }
        elif sample == 'RESOLVE-B BMF double':
            self.params = {'M1': 11.11,
                           'phi1': 6.93*10**3,
                           'alpha1': -1.30,
                           'M2': 11.11,
                           'phi2': 0.0,
                           'alpha2':0.0,
                           }
        elif sample == 'ECO BMF single':
            self.params = {'M1': 10.98,
                           'phi1': 2.74*10**3,
                           'alpha1': -0.52,
                           'M2': 10.87,
                           'phi2': 3.25,
                           'alpha2':-1.38,
                           }
        elif sample == 'ECO BMF double':
            self.params = {'M1': 10.87,
                           'phi1': 9.00*10**3,
                           'alpha1': -0.52,
                           'M2': 10.87,
                           'phi2': 3.25,
                           'alpha2':-1.38,
                           }
        else:
            msg = ('Unknown sample passed!')
            raise ValueError(msg)

        self.phi1 = self.params['phi1']
        self.x1 = self.params['M1']
        self.alpha1 = self.params['alpha1']
        self.phi2 = self.params['phi2']
        self.x2 = self.params['M2']
        self.alpha2 = self.params['alpha2']
    
    def __call__(self, mstar):
        """
        stellar mass function
        
        Parameters
        ----------
        mstar : array_like
            stellar mass in units Msol/h^2
        
        Returns
        -------
        phi : nunpy.array
            number density in units h^3 Mpc^-3 dex^-1
        """
        
        #convert from h=1 to h=0.7
        mstar = mstar / self.littleh**2
        
        #take log of stellar masses
        mstar = np.log10(mstar)
        
        #convert from h=0.7 to h=1.0
        return self.s(mstar) / self.littleh**3


@custom_model
def Log_Schechter(x, phi0=0.001, x0=10.5, alpha=-1.0):
    """
    log schecter x function
    """
    x = np.asarray(x)
    x = x.astype(float)
    norm = np.log(10.0)*phi0
    val = norm*(10.0**((x-x0)*(1.0+alpha)))*np.exp(-10.0**(x-x0))
    return val