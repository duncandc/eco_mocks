"""
create and save all fiducial mocks
"""
import os
from fiducial_models import model_1, model_2
import numpy as np
from astropy.io import ascii

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
    model_1.populate_mock(halocat)
    mock_1 = model_1.mock.galaxy_table
    mock_1.write(save_path+'mock_1.hdf5', format='hdf5', path='data') 

    model_2.populate_mock(halocat)
    mock_2 = model_2.mock.galaxy_table
    mock_2.write(save_path+'mock_2.hdf5', format='hdf5', path='data') 

if __name__ == '__main__':
    main()