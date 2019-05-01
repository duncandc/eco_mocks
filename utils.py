"""
A set of utilites for the ECO/RESOLVE surveys
"""

from __future__ import (division, print_function, absolute_import,
                        unicode_literals)
import numpy as np
import sys
from halotools.mock_observables.mock_survey import ra_dec_z
from halotools.utils import rotation_matrices_from_angles, normalized_vectors
from eco_mocks.survey_utils import inside_survey_area
from astropy.table import Table


__author__ = ['Duncan Campbell']
__all__ = ['random_survey']


def nearest_nieghbors_1d(arr1, arr2, seed=None):
    """
    Find the indices of the nearest neighbors in `arr1` for each element in `arr2`.

    Parameters
    ----------
    arr1 : array_like

    arr2 : array_like

    seed : int, optional
        random seed

    Returns
    -------
    idx : numpy.array
        an array of length len(arr2) with indices into arr1

    Notes
    -----
    If there are repeat values in `arr1`, the elements to match to
    are chosen randomly amongst the repeat values.  Setting the
    `seed` argument will make the results for such cases reproducible.
    """

    np.random.seed(seed)

    arr1 = np.atleast_1d(arr1)
    arr2 = np.atleast_1d(arr2)

    n = len(arr1)

    sort_inds = np.argsort(arr1)

    # find repeat elements in arr1
    unq_x, n_x = np.unique(arr1[sort_inds], return_counts=True)
    n_x = np.repeat(n_x,n_x)  # number of times each element repeats

    # find the nearest ***larger*** value in arr1
    matching_inds_a = np.searchsorted(arr1[sort_inds], arr2)
    mask = (matching_inds_a>=n)
    matching_inds_a[mask] = n-1

    # find the nearest ***smaller*** value in arr1
    matching_inds_b = np.searchsorted(-1.0*arr1[sort_inds][::-1], -1.0*arr2)
    matching_inds_b = n-1-matching_inds_b
    mask = (matching_inds_b<0)
    matching_inds_b[mask] = 0

    # choose the match wehich results in the smallest difference
    da = np.fabs(arr1[matching_inds_a] - arr2)
    db = np.fabs(arr1[matching_inds_b] - arr2)
    matching_inds = np.where(da<=db, matching_inds_a, matching_inds_b)
    #matching_inds = matching_inds_a

    # account for repeats in arr1
    ran_int = np.round(np.random.random(len(matching_inds))*(n_x[matching_inds]-1))
    ran_int = ran_int.astype(int)
    matching_inds = matching_inds+ran_int

    return sort_inds[matching_inds]


def pad_box(coords, Lbox=130):
    """
    Tile a PBC simulation cube into 3x3x3 cube.
    Note the memory requirements to store this will by 27 times
    larger than that required to store the original data!
    """

    N = len(coords)
    new_coords = np.zeros((N*(3**3), 3))

    l = 0
    for i in range(0, 3):
        for j in range(0, 3):
            for k in range(0, 3):
                dx = Lbox*i
                dy = Lbox*j
                dz = Lbox*k
                new_coords[l*N:(l+1)*N, 0] = coords[:, 0]+dx
                new_coords[l*N:(l+1)*N, 1] = coords[:, 1]+dy
                new_coords[l*N:(l+1)*N, 2] = coords[:, 2]+dz
                l = l+1

    return new_coords, np.tile(np.arange(0, N), 27)


def random_rotation(coords, vels):
    """
    apply a random rotation to a set of coordinates and velocities
    """

    ran_angle = np.random.random(size=1)*(np.pi)
    ran_direction = normalized_vectors(np.random.random((3,)))*2.0 - 1.0
    ran_rot = rotation_matrices_from_angles(ran_angle, ran_direction)

    new_coords = rotate_vector_collection(ran_rot, coords)

    new_vels = rotate_vector_collection(ran_rot, vels)

    return new_coords, new_vels


def rotate_vector_collection(rotation_matrices, vectors, optimize=False):
    r""" Given a collection of rotation matrices and a collection of 3d vectors,
    apply each matrix to rotate the corresponding vector.

    Parameters
    ----------
    rotation_matrices : ndarray
        Numpy array of shape (npts, 3, 3) storing a collection of rotation matrices.
        If an array of shape (3, 3) is passed, all the vectors
        are rotated using the same rotation matrix.

    vectors : ndarray
        Numpy array of shape (npts, 3) storing a collection of 3d vectors

    Returns
    -------
    rotated_vectors : ndarray
        Numpy array of shape (npts, 3) storing a collection of 3d vectors

    Examples
    --------
    In this example, we'll randomly generate two sets of unit-vectors, `v0` and `v1`.
    We'll use the `rotation_matrices_from_vectors` function to generate the
    rotation matrices that rotate each `v0` into the corresponding `v1`.
    Then we'll use the `rotate_vector_collection` function to apply each
    rotation, and verify that we recover each of the `v1`.

    >>> npts = int(1e4)
    >>> v0 = normalized_vectors(np.random.random((npts, 3)))
    >>> v1 = normalized_vectors(np.random.random((npts, 3)))
    >>> rotation_matrices = rotation_matrices_from_vectors(v0, v1)
    >>> v2 = rotate_vector_collection(rotation_matrices, v0)
    >>> assert np.allclose(v1, v2)
    """

    # apply same rotation matrix to all vectors
    if np.shape(rotation_matrices) == (3, 3):
        return np.dot(rotation_matrices, vectors.T).T
    # rotate each vector by associated rotation matrix
    else:
        try:
            return np.einsum('ijk,ik->ij', rotation_matrices, vectors, optimize=optimize)
        except TypeError:
            return np.einsum('ijk,ik->ij', rotation_matrices, vectors)


def random_survey(mock, Lbox=130.0, survey='eco_a', seed=None):
    """
    create a random mock survey within a periodic box
    """

    np.random.seed(seed=seed)

    # put mock coordinates in halotools format
    coords = np.vstack((mock['x'], mock['y'], mock['z'])).T
    vels = np.vstack((mock['vx'], mock['vy'], mock['vz'])).T

    # pad the simulation
    new_coords, inds = pad_box(coords)
    new_vels = np.tile(vels, (27, 1))

    # choose a new coordionate center within the central box
    new_coords[:, 0] = new_coords[:, 0] - (1.5*Lbox) + ((np.random.random(size=1) *2-1) * (0.5*Lbox))
    new_coords[:, 1] = new_coords[:, 1] - (1.5*Lbox) + ((np.random.random(size=1) *2-1) * (0.5*Lbox))
    new_coords[:, 2] = new_coords[:, 2] - (1.5*Lbox) + ((np.random.random(size=1) *2-1) * (0.5*Lbox))

    # randomly rotate coordinate system
    new_coords, new_vels = random_rotation(new_coords, new_vels)

    # put into ra, dec, redshift space
    ra, dec, z = ra_dec_z(new_coords, new_vels)
    ra = np.degrees(ra)
    dec = np.degrees(dec)

    # keep only galaxies within the footprint
    mask = inside_survey_area(ra, dec, z, survey=survey)
    new_coords = new_coords[mask]
    new_vels = new_vels[mask]
    inds = inds[mask]
    ra = ra[mask]
    dec = dec[mask]
    z = z[mask]

    # create new table to store mock
    new_mock = Table([new_coords[:, 0],
                      new_coords[:, 1],
                      new_coords[:, 2],
                      new_vels[:, 0],
                      new_vels[:, 1],
                      new_vels[:, 2]],
                     names=['x', 'y', 'z',
                            'vx', 'vy', 'vz'])

    # add the remaining galaxy properties in the mock
    cols = list(mock.dtype.names)
    cols.remove('x')
    cols.remove('y')
    cols.remove('z')
    cols.remove('vx')
    cols.remove('vy')
    cols.remove('vz')

    for name in cols:
        new_mock[name] = mock[name][inds]

    # add ra, dec, and redshift
    new_mock['ra'] = ra
    new_mock['dec'] = dec
    new_mock['redshift'] = z

    # remove excess
    r = np.sqrt(new_mock['x']**2+new_mock['y']**2+new_mock['z']**2)
    mask = (r < Lbox)

    return new_mock[mask]


def fractional_polar_axes(f, thlim=(0, 180), rlim=(0, 1), step=(30, 0.01),
                          thlabel=r'$\alpha$', rlabel='z', ticklabels=True):
    import matplotlib.pyplot as plt

    from matplotlib.transforms import Affine2D
    from matplotlib.projections import PolarAxes
    from mpl_toolkits.axisartist import angle_helper
    from mpl_toolkits.axisartist.grid_finder import MaxNLocator
    from mpl_toolkits.axisartist.floating_axes import GridHelperCurveLinear, FloatingSubplot
    from mpl_toolkits.axisartist.grid_finder import (FixedLocator, MaxNLocator,
                                                 DictFormatter)

    """Return polar axes that adhere to desired theta (in deg) and r limits. steps for theta
    and r are really just hints for the locators. Using negative values for rlim causes
    problems for GridHelperCurveLinear for some reason"""
    th0, th1 = thlim # deg
    r0, r1 = rlim
    thstep, rstep = step

    # scale degrees to radians:
    tr_scale = Affine2D().scale(np.pi/180., 1.)
    tr = tr_scale + PolarAxes.PolarTransform()
    theta_grid_locator = angle_helper.LocatorDMS((th1-th0) // thstep)
    r_grid_locator = MaxNLocator((r1-r0) // rstep)
    theta_tick_formatter = angle_helper.FormatterDMS()
    theta_tick_formatter = None
    theta_ticks = [(0, r"$90^{\circ}$"),
                   (30, r"$120^{\circ}$"),
                   (60, r"$150^{\circ}$"),
                   (90, r"$180^{\circ}$"),
                   (120, r"$210^{\circ}$"),
                   (150, r"$270^{\circ}$"),
                   (180, r"$0^{\circ}$")]
    theta_tick_formatter = DictFormatter(dict(theta_ticks))
    grid_helper = GridHelperCurveLinear(tr,
                                        extremes=(th0, th1, r0, r1),
                                        grid_locator1=theta_grid_locator,
                                        grid_locator2=r_grid_locator,
                                        tick_formatter1=theta_tick_formatter,
                                        tick_formatter2=None)

    a = FloatingSubplot(f, 111, grid_helper=grid_helper)
    f.add_subplot(a)

    # adjust x axis (theta):
    a.axis["bottom"].set_visible(False)
    a.axis["top"].set_axis_direction("bottom") # tick direction
    a.axis["top"].toggle(ticklabels=ticklabels, label=bool(thlabel))
    a.axis["top"].major_ticklabels.set_axis_direction("top")
    a.axis["top"].label.set_axis_direction("top")

    # adjust y axis (r):
    a.axis["left"].set_axis_direction("bottom") # tick direction
    a.axis["right"].set_axis_direction("top") # tick direction
    a.axis["left"].toggle(ticklabels=ticklabels, label=bool(rlabel))

    # add labels:
    a.axis["top"].label.set_text(thlabel)
    a.axis["left"].label.set_text(rlabel)

    # create a parasite axes whose transData is theta, r:
    auxa = a.get_aux_axes(tr)
    # make aux_ax to have a clip path as in a?:
    auxa.patch = a.patch
    # this has a side effect that the patch is drawn twice, and possibly over some other
    # artists. So, we decrease the zorder a bit to prevent this:
    a.patch.zorder = -2

    # add sector lines for both dimensions:
    thticks = grid_helper.grid_info['lon_info'][0]
    rticks = grid_helper.grid_info['lat_info'][0]
    for th in thticks[1:-1]: # all but the first and last
        auxa.plot([th, th], [r0, r1], '--', c='grey', zorder=-1)
    for ri, r in enumerate(rticks):
        # plot first r line as axes border in solid black only if it isn't at r=0
        if ri == 0 and r != 0:
            ls, lw, color = 'solid', 2, 'black'
        else:
            ls, lw, color = 'dashed', 1, 'grey'
        # From http://stackoverflow.com/a/19828753/2020363
        auxa.add_artist(plt.Circle([0, 0], radius=r, ls=ls, lw=lw, color=color, fill=False,
                        transform=auxa.transData._b, zorder=-1))

    return auxa


