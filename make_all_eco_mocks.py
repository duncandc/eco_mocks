"""
create and save all fiducial mocks
"""
import os
from fiducial_models import model_1a, model_1b, model_1c, model_1d, model_1e, model_1f
from fiducial_models import model_2a, model_2b, model_2c, model_2d, model_2e, model_2f
import numpy as np
from astropy.io import ascii
import sys

cwd = os.getcwd()
save_path = cwd + '/mocks/'

def main():

    # load halo catalog
    import warnings
    warnings.simplefilter(action='ignore', category=FutureWarning)

    from halotools import sim_manager
    simname = 'vishnu_130'
    halocat = sim_manager.CachedHaloCatalog(simname = simname, redshift=0.0, dz_tol = 0.1,
                                           version_name='custom', halo_finder='Rockstar')

    # populate mock catalogs
    model_1a.populate_mock(halocat)
    mock_1 = model_1a.mock.galaxy_table
    mock_1.write(save_path+'mock_1a.hdf5', format='hdf5', path='data', overwrite=True)

    model_1b.populate_mock(halocat)
    mock_1 = model_1b.mock.galaxy_table
    mock_1.write(save_path+'mock_1b.hdf5', format='hdf5', path='data', overwrite=True)

    model_1c.populate_mock(halocat)
    mock_1 = model_1c.mock.galaxy_table
    mock_1.write(save_path+'mock_1c.hdf5', format='hdf5', path='data', overwrite=True)

    model_1d.populate_mock(halocat)
    mock_1 = model_1d.mock.galaxy_table
    mock_1.write(save_path+'mock_1d.hdf5', format='hdf5', path='data', overwrite=True)

    model_1e.populate_mock(halocat)
    mock_1 = model_1e.mock.galaxy_table
    mock_1.write(save_path+'mock_1e.hdf5', format='hdf5', path='data', overwrite=True)

    model_1f.populate_mock(halocat)
    mock_1 = model_1f.mock.galaxy_table
    mock_1.write(save_path+'mock_1f.hdf5', format='hdf5', path='data', overwrite=True)



    model_2a.populate_mock(halocat)
    mock_2 = model_2a.mock.galaxy_table
    mock_2.write(save_path+'mock_2a.hdf5', format='hdf5', path='data', overwrite=True)

    model_2b.populate_mock(halocat)
    mock_2 = model_2b.mock.galaxy_table
    mock_2.write(save_path+'mock_2b.hdf5', format='hdf5', path='data', overwrite=True)

    model_2c.populate_mock(halocat)
    mock_2 = model_2c.mock.galaxy_table
    mock_2.write(save_path+'mock_2c.hdf5', format='hdf5', path='data', overwrite=True)

    model_2d.populate_mock(halocat)
    mock_2 = model_2d.mock.galaxy_table
    mock_2.write(save_path+'mock_2d.hdf5', format='hdf5', path='data', overwrite=True)

    model_2e.populate_mock(halocat)
    mock_2 = model_2e.mock.galaxy_table
    mock_2.write(save_path+'mock_2e.hdf5', format='hdf5', path='data', overwrite=True)

    model_2f.populate_mock(halocat)
    mock_2 = model_2f.mock.galaxy_table
    mock_2.write(save_path+'mock_2f.hdf5', format='hdf5', path='data', overwrite=True)

if __name__ == '__main__':
    main()
