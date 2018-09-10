"""
objects to model secondary galaxy properties
"""

from __future__ import (division, print_function, absolute_import,
                        unicode_literals)
import numpy as np
from astropy.table import Table
from sklearn import mixture
from scipy.stats import norm

from eco_mocks.utils import nearest_nieghbors_1d
from eco_mocks.eco_galaxy_properties import eco_table as default_data


class prim_prop_nn(object):
    """
    primary galaxy property nearest neighbor model for secondary property.
    """

    def __init__(self, prim_prop='stellar_mass', sec_prop='u_minus_r', reference_table=None):
        """
        Parameters
        ----------
        prim_prop : string
            primary galaxy property

        sec_prop : string
            secondary galaxy property

        reference_table : astropy.table object
            reference table to search for nearest neighbors.  This table must have keys for
            `prim_prop` and `sec_prop`.
        """

        # set the reference table to search for nearest neighbors and
        # assiging the asscoaited vales for the secondary property
        if reference_table is None:
            self.data = default_data
        else:
            self.data = reference_table

        # make sure the reference table has the requested keys
        if prim_prop not in self.data.keys():
            msg = ('primary galaxy property not found in reference table.')
            raise ValueError(msg)
        else:
            self.prim_prop_key = prim_prop

        if sec_prop not in self.data.keys():
            msg = ('secondary galaxy property not found in reference table.')
            raise ValueError(msg)
        else:
            self.sec_prop_key = sec_prop

    def rvs(self, x, seed=None, return_index=False):
        """
        Return values of secondary galaxy property
        given values of the primary galaxy property.

        Parameters
        ----------
        x : array_like
            array of primary galaxy properties

        Returns
        -------
        y : numpy.array
            array of secondary galaxy properties

        idx : numpy.array
            if return_index is True, return the index from the reference catalog
        """

        idx = nearest_nieghbors_1d(self.data[self.prim_prop_key], x, seed=seed)

        if return_index:
            return self.data[self.sec_prop_key][idx], idx
        else:
            return self.data[self.sec_prop_key][idx]


class color_model_2(object):
    """
    a double gaussian mixture model for galaxy color
    """

    def __init__(self, data=None):
        """
        """

        if data is not None:
            self.mstar = data[0, :]
            self.color = data[1, :]
        else:
            self.mstar = default_data['stellar_mass']
            self.color = default_data['u_minus_r']

        self.set_params(None)

    def set_params(self, params):
        """
        """

        names = ['B', 'Q', 'nu',  # red fraction params
                 'm1', 'm1',      # blue mean params
                 'm2', 'm2',      # red mean params
                 's1', 's2']      # scatters

        # set default params
        self.params = {'B': 0.6248,
                       'Q': 35.00,
                       'nu': 0.089,
                       'm1': 0.267,
                       'b1': -1.459,
                       'm2': 0.399,
                       'b2': -2.306,
                       's1': 0.2,
                       's2': 0.15
                       }

        # replace paramaters that are passed
        if params is not None:
            for key in params.keys():
                if key in names:
                    self.params[key] = params[key]

    def f_blue(self, mstar):
        """
        blue fraction of galaxies
        """

        return 1.0 - self.f_red(mstar)

    def f_red(self, mstar):
        """
        red fraction of galaxies
        """

        x = np.log10(mstar)

        B = self.params['B']
        Q = self.params['Q']
        nu = self.params['nu']

        return gen_logistic(x, B, Q, nu)

    def blue_mean(self, mstar):
        """
        """
        x = np.log10(mstar)
        m = self.params['m1']
        b = self.params['b1']
        return linear(x, m, b)

    def red_mean(self, mstar):
        """
        """

        x = np.log10(mstar)
        m = self.params['m2']
        b = self.params['b2']
        return linear(x, m, b)

    def blue_scatter(self, mstar):
        """
        width of blue disttribution
        """
        x = np.log10(mstar)
        return self.params['s1'] + 0.0*x

    def red_scatter(self, mstar):
        """
        width of red disttribution
        """
        x = np.log10(mstar)
        return self.params['s2'] + 0.0*x

    def pdf(self, mstar, color):
        """
        the probability of color given stellar mass
        """

        p1 = self.f_blue(mstar)*norm.pdf(color,
                                         loc=self.blue_mean(mstar),
                                         scale=self.blue_scatter(mstar))
        p2 = self.f_red(mstar)*norm.pdf(color,
                                        loc=self.red_mean(mstar),
                                        scale=self.red_scatter(mstar))

        return p1+p2

    def lnlike(self, params=None):
        """
        sum of the log liklihoods
        """

        if params is not None:
            self.set_params(params)

        return np.sum(np.log(self.pdf(self.mstar, self.color)))


class color_model_3(object):
    """
    a triple gaussian mixture model for galaxy color
    """

    def __init__(self):
        """
        """
        self.set_params(None)

    def set_params(self, params):
        """
        """

        names = ['B', 'Q', 'nu',    # red fraction params
                 'mf', 'bf',        # blue fraction params
                 'm1', 'm1',        # blue mean params
                 'm2', 'm2',        # green mean params
                 'm3', 'm3',        # red mean params
                 's1', 's2', 's3']  # scatters

        # set default params
        self.params = {'B': 0.6248,
                       'Q': 35.00,
                       'nu': 0.089,
                       'mf': -0.33,  #-0.330,
                       'bf': 3.0, # 2.8248,
                       'm1': 0.267,
                       'b1': -1.459,
                       'm2': 0.399,
                       'b2': -2.306,
                       'm3': 0.233,
                       'b3': -0.210,
                       's1': 0.141,
                       's2': 0.157,
                       's3': 0.147
                       }

        # replace paramaters that are passed
        if params is not None:
            for key in params.keys():
                if key in names:
                    self.params[key] = params[key]

    def f_blue(self, mstar):
        """
        blue fraction of galaxies
        """

        x = np.log10(mstar)

        m = self.params['mf']
        B = self.params['bf']

        num = log_linear(x, m, B)
        den = 1.0 + log_linear(x, m, B)

        return num/den

    def f_green(self, mstar):
        """
        green fraction of galaxies
        """
        return 1.0 - self.f_blue(mstar) - self.f_red(mstar)

    def f_red(self, mstar):
        """
        red fraction of galaxies
        """

        x = np.log10(mstar)

        B = self.params['B']
        Q = self.params['Q']
        nu = self.params['nu']

        return gen_logistic(x, B, Q, nu)

    def blue_mean(self, mstar):
        """
        """
        x = np.log10(mstar)
        m = self.params['m1']
        b = self.params['b1']
        return linear(x, m, b)

    def green_mean(self, mstar):
        """
        """
        x = np.log10(mstar)
        m = self.params['m2']
        b = self.params['b2']
        return linear(x, m, b)

    def red_mean(self, mstar):
        """
        """

        x = np.log10(mstar)
        m = self.params['m3']
        b = self.params['b3']
        return linear(x, m, b)

    def blue_scatter(self, mstar):
        """
        width of blue disttribution
        """
        x = np.log10(mstar)
        return self.params['s1'] + 0.0*x

    def green_scatter(self, mstar):
        """
        width of green disttribution
        """
        x = np.log10(mstar)
        return self.params['s2'] + 0.0*x

    def red_scatter(self, mstar):
        """
        width of red disttribution
        """
        x = np.log10(mstar)
        return self.params['s3'] + 0.0*x


def gen_logistic(x, B, Q, nu):
    """
    generalized logistic function model
    """
    A = 0.0
    K = 1.0
    C = 1.0
    return A+(K-A)/(C+Q*np.exp(-B*x))**(1.0/nu)


def log_linear(x, m, b):
    """
    log-linear model
    """
    return 10**(m*x + b)


def linear(x, m, b):
    """
    linear model
    """
    return m*x + b
