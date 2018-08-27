"""
A set of utilites for the ECO/RESOLVE surveys
"""

from __future__ import (division, print_function, absolute_import,
                        unicode_literals)
import numpy as np
from halotools.utils import sample_spherical_surface

__author__ = ['Duncan Campbell']
__all__ = ['inside_survey_area', 'mc_survey_area']


def inside_survey_area(ra, dec, survey='resolve_a'):
    r"""
    Determine if a set of coordinates are within the RESOLVE/ECO
    survey boundaries.

    Parameters
    ==========
    ra : array_like
        right acension of galaxies in degrees (-180,180)

    dec : array_like
        declination of galaxies in degrees (-90,90)

    survey : string

    Returns
    =======
    mask : numpy.array
        boolean array
    """

    ra = np.atleast_1d(ra)
    dec = np.atleast_1d(dec)

    # survey boundaries taken from Eckert + 2015, ArXiv:1507.08669
    if survey == 'resolve_b':
        params = {'min_ra': 3.0*(360.0/24.0),
                  'max_ra': 22.0*(360.0/24.0),
                  'min_dec': -1.25,
                  'max_dec': 1.25,
                  'wrap_ra': True,
                  'wrap_dec': False
                  }
    elif survey == 'resolve_a':
        params = {'min_ra': 236.25,
                  'max_ra': 131.25,
                  'min_dec': 0.0,
                  'max_dec': 5.0,
                  'wrap_ra': False,
                  'wrap_dec': False
                  }
    elif survey == 'eco_a':  # by eye looking at fig. 1, Moffett + 2015
        params = {'min_ra': 8.5*(360.0/24.0),
                  'max_ra': 16.0*(360.0/24.0),
                  'min_dec': 0.0,
                  'max_dec': 50.0,
                  'wrap_ra': False,
                  'wrap_dec': False
                  }
    else:
        msg = ('Survey string not recognized.')
        raise ValueError(msg)

    # note that I define ra to be bounded between [-180, 180]
    params['min_ra'] = params['min_ra']-180.0
    params['max_ra'] = params['max_ra']-180.0

    # make ra-dec cut
    radec_mask = ra_dec_box_mask(ra, dec, params)

    return radec_mask


def ra_dec_box_mask(ra, dec, params):
    r"""
    Determine if a set of coordinates are within a box on the sky
    defined by pairs of ra and dec boundaries.

    Parameters
    ==========
    ra : array_like
        right acension of galaxies in degrees (-180,180)

    dec : array_like
        declination of galaxies in degrees (-90,90)

    params : dictionary
        dictionary indicating ra-dec boundaries of box

    Returns
    =======
    mask : numpy.array
        boolean array
    """

    if params['wrap_ra']:
        ra_mask = (ra > params['max_ra']) | (ra < params['min_ra'])
    else:
        ra_mask = (ra < params['max_ra']) & (ra > params['min_ra'])

    dec_mask = (dec > params['min_dec']) & (dec < params['max_dec'])

    mask = ra_mask & dec_mask

    return mask


def mc_survey_area(N_points, survey='resolve_a', seed=None):
    r"""
    Estimate the area of the survey.

    Paramaters
    ==========
    N_points : float

    survey : string

    seed : integer

    Returns
    =======
    area : float
        an estimate of the area of the survey in steradians
    """

    ran_coords = np.asarray(sample_spherical_surface(N_points, seed))
    ran_ra = ran_coords[:, 0]
    ran_dec = ran_coords[:, 1]

    mask = inside_survey_area(ran_ra, ran_dec, survey=survey)

    area_estimate = 1.0 * (np.sum(mask)/N_points) * (4.0 * np.pi)

    return area_estimate
