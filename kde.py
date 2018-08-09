"""
custom implementation of kde for conditional sampling
of multidimensional distributions
"""

from __future__ import division, print_function, absolute_import
import numpy as np
from scipy.stats import gaussian_kde as scipy_gaussian_kde
# from scipy.stats import multivariate_normal
from numpy.random import multivariate_normal
from sklearn.neighbors import BallTree


class gaussian_kde(scipy_gaussian_kde):
    """
    a gaussian kde class with an additonal method to allow
    for approximate conditional sampling of the pdf based on
    the scipy.stats.gaussian_kde class
    """

    def __init__(self, dataset, bw_method=None):
        """
        """
        super(gaussian_kde, self).__init__(dataset, bw_method=bw_method)

    def conditional_sample(self, x, c):
        """
        perform an approximate conditonal sampling of the pdf of the
        N dimenstional of pdf.

        Parameters
        ----------
        x : array_like
            array of dependent variables of dinemsion q

        c : array_like
            boolean array defining conditonal dimensions with
            q elements set to True

        Returns
        -------
        y : numpy.array
            array of sampled values of dimantion N-q
        """

        x = np.atleast_1d(x)
        c = np.atleast_1d(c)

        if x.ndim == 1:
            x = x[..., np.newaxis]

        # find nearest neighbor
        V = np.diagonal(self.covariance)
        d = self.dataset.T
        tree = BallTree(d[:, c], metric='seuclidean', V=V[c])

        cov = np.atleast_2d(self.covariance[c, c])
        d = np.sum(c)
        size = len(x)
        x = x + multivariate_normal(np.zeros((d,), float),
                                    cov, size=size)

        indices = tree.query(x, k=1, return_distance=False)

        cov = np.atleast_2d(self.covariance[~c, ~c])

        size = len(x)
        d = self.d - np.sum(c)
        norm = multivariate_normal(np.zeros((d,), float),
                                   cov, size=size)
        means = self.dataset[~c, indices]
        print(np.shape(means), np.shape(norm))
        return means + norm

