"""
custom implementation of kde for conditional sampling of multidimensional distributions
"""

from __future__ import division, print_function, absolute_import
import numpy as np
from scipy.stats import scipy_gaussian_kde
from sklearn.neighbors import NearestNeighbors


class gaussian_kde(scipy_gaussian_kde):
    """
    a gaussian kde class with an additonal method to allow
    for approximate conditional sampling of the pdf based on
    the scipy.stats.gaussian_kde class
    """

    def __init__(self, dataset, bw_method=None):
        """
        """
        super(Y, self).__init__(dataset, bw_method=bw_method)

    def conditional_sample(self, x, c):
        """
        perform an approximate conditonal sampling of the pdf 

        Parameters
        ----------
        x : array_like
            array of dependent variables

        c : array_like
            boolean array defining conditonal dimensions

        Returns
        -------
        y : numpy.array
            array of sampled values
        """
        
        # define distance metric using kernel
        inv_var = 1.0/np.diagonal(self.covariance)[c]
        metric_dict = {'w': inv_var}
        
        nn = NearestNeighbors(metric='seuclidean', n_neighbors=1, metric_params=metric_dict)
        nn = nn.fit(self.dataset[:,c]) 
        
        # find nearest neighbor
        inds = nn.kneighbors(x, return_distance=False)

        # define covaraince matrix
        cov = self.covariance[~c,~c]

        norms = transpose(multivariate_normal(zeros((self.d[inds],), float),
                          cov, size=size))
        means = self.dataset[inds]

        return means + norms

