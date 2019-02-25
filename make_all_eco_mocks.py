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
    model_1a.param_dict['rho'] = 0.0
    model_1a.populate_mock(halocat)
    mock_1 = model_1a.mock.galaxy_table
    mock_1.write(save_path+'mock_1a_00.hdf5', format='hdf5', path='data', overwrite=True)

    model_1a.param_dict['rho'] = -0.1
    model_1a.populate_mock(halocat)
    mock_1 = model_1a.mock.galaxy_table
    mock_1.write(save_path+'mock_1a_-01.hdf5', format='hdf5', path='data', overwrite=True)

    model_1a.param_dict['rho'] = -0.2
    model_1a.populate_mock(halocat)
    mock_1 = model_1a.mock.galaxy_table
    mock_1.write(save_path+'mock_1a_-02.hdf5', format='hdf5', path='data', overwrite=True)

    model_1a.param_dict['rho'] = -0.3
    model_1a.populate_mock(halocat)
    mock_1 = model_1a.mock.galaxy_table
    mock_1.write(save_path+'mock_1a_-03.hdf5', format='hdf5', path='data', overwrite=True)

    model_1a.param_dict['rho'] = -0.4
    model_1a.populate_mock(halocat)
    mock_1 = model_1a.mock.galaxy_table
    mock_1.write(save_path+'mock_1a_-04.hdf5', format='hdf5', path='data', overwrite=True)

    model_1a.param_dict['rho'] = -0.5
    model_1a.populate_mock(halocat)
    mock_1 = model_1a.mock.galaxy_table
    mock_1.write(save_path+'mock_1a_-05.hdf5', format='hdf5', path='data', overwrite=True)

    model_1a.param_dict['rho'] = -0.6
    model_1a.populate_mock(halocat)
    mock_1 = model_1a.mock.galaxy_table
    mock_1.write(save_path+'mock_1a_-06.hdf5', format='hdf5', path='data', overwrite=True)

    model_1a.param_dict['rho'] = -0.7
    model_1a.populate_mock(halocat)
    mock_1 = model_1a.mock.galaxy_table
    mock_1.write(save_path+'mock_1a_-07.hdf5', format='hdf5', path='data', overwrite=True)

    model_1a.param_dict['rho'] = -0.8
    model_1a.populate_mock(halocat)
    mock_1 = model_1a.mock.galaxy_table
    mock_1.write(save_path+'mock_1a_-08.hdf5', format='hdf5', path='data', overwrite=True)

    model_1a.param_dict['rho'] = -0.9
    model_1a.populate_mock(halocat)
    mock_1 = model_1a.mock.galaxy_table
    mock_1.write(save_path+'mock_1a_-09.hdf5', format='hdf5', path='data', overwrite=True)

    model_1a.param_dict['rho'] = -1.0
    model_1a.populate_mock(halocat)
    mock_1 = model_1a.mock.galaxy_table
    mock_1.write(save_path+'mock_1a_-10.hdf5', format='hdf5', path='data', overwrite=True)



    model_1b.param_dict['rho'] = 0.0
    model_1b.populate_mock(halocat)
    mock_1 = model_1b.mock.galaxy_table
    mock_1.write(save_path+'mock_1b_00.hdf5', format='hdf5', path='data', overwrite=True)

    model_1b.param_dict['rho'] = 0.1
    model_1b.populate_mock(halocat)
    mock_1 = model_1b.mock.galaxy_table
    mock_1.write(save_path+'mock_1b_01.hdf5', format='hdf5', path='data', overwrite=True)

    model_1b.param_dict['rho'] = 0.2
    model_1b.populate_mock(halocat)
    mock_1 = model_1b.mock.galaxy_table
    mock_1.write(save_path+'mock_1b_02.hdf5', format='hdf5', path='data', overwrite=True)

    model_1b.param_dict['rho'] = 0.3
    model_1b.populate_mock(halocat)
    mock_1 = model_1b.mock.galaxy_table
    mock_1.write(save_path+'mock_1b_03.hdf5', format='hdf5', path='data', overwrite=True)

    model_1b.param_dict['rho'] = 0.4
    model_1b.populate_mock(halocat)
    mock_1 = model_1b.mock.galaxy_table
    mock_1.write(save_path+'mock_1b_04.hdf5', format='hdf5', path='data', overwrite=True)

    model_1b.param_dict['rho'] = 0.5
    model_1b.populate_mock(halocat)
    mock_1 = model_1b.mock.galaxy_table
    mock_1.write(save_path+'mock_1b_05.hdf5', format='hdf5', path='data', overwrite=True)

    model_1b.param_dict['rho'] = 0.6
    model_1b.populate_mock(halocat)
    mock_1 = model_1b.mock.galaxy_table
    mock_1.write(save_path+'mock_1b_06.hdf5', format='hdf5', path='data', overwrite=True)

    model_1b.param_dict['rho'] = 0.7
    model_1b.populate_mock(halocat)
    mock_1 = model_1b.mock.galaxy_table
    mock_1.write(save_path+'mock_1b_07.hdf5', format='hdf5', path='data', overwrite=True)

    model_1b.param_dict['rho'] = 0.8
    model_1b.populate_mock(halocat)
    mock_1 = model_1b.mock.galaxy_table
    mock_1.write(save_path+'mock_1b_08.hdf5', format='hdf5', path='data', overwrite=True)

    model_1b.param_dict['rho'] = 0.9
    model_1b.populate_mock(halocat)
    mock_1 = model_1b.mock.galaxy_table
    mock_1.write(save_path+'mock_1b_09.hdf5', format='hdf5', path='data', overwrite=True)

    model_1b.param_dict['rho'] = 1.0
    model_1b.populate_mock(halocat)
    mock_1 = model_1b.mock.galaxy_table
    mock_1.write(save_path+'mock_1b_10.hdf5', format='hdf5', path='data', overwrite=True)



    model_1c.param_dict['rho'] = 0.0
    model_1c.populate_mock(halocat)
    mock_1 = model_1c.mock.galaxy_table
    mock_1.write(save_path+'mock_1c_00.hdf5', format='hdf5', path='data', overwrite=True)

    model_1c.param_dict['rho'] = 0.1
    model_1c.populate_mock(halocat)
    mock_1 = model_1c.mock.galaxy_table
    mock_1.write(save_path+'mock_1c_01.hdf5', format='hdf5', path='data', overwrite=True)

    model_1c.param_dict['rho'] = 0.2
    model_1c.populate_mock(halocat)
    mock_1 = model_1c.mock.galaxy_table
    mock_1.write(save_path+'mock_1c_02.hdf5', format='hdf5', path='data', overwrite=True)

    model_1c.param_dict['rho'] = 0.3
    model_1c.populate_mock(halocat)
    mock_1 = model_1c.mock.galaxy_table
    mock_1.write(save_path+'mock_1c_03.hdf5', format='hdf5', path='data', overwrite=True)

    model_1c.param_dict['rho'] = 0.4
    model_1c.populate_mock(halocat)
    mock_1 = model_1c.mock.galaxy_table
    mock_1.write(save_path+'mock_1c_04.hdf5', format='hdf5', path='data', overwrite=True)

    model_1c.param_dict['rho'] = 0.5
    model_1c.populate_mock(halocat)
    mock_1 = model_1c.mock.galaxy_table
    mock_1.write(save_path+'mock_1c_05.hdf5', format='hdf5', path='data', overwrite=True)

    model_1c.param_dict['rho'] = 0.6
    model_1c.populate_mock(halocat)
    mock_1 = model_1c.mock.galaxy_table
    mock_1.write(save_path+'mock_1c_06.hdf5', format='hdf5', path='data', overwrite=True)

    model_1c.param_dict['rho'] = 0.7
    model_1c.populate_mock(halocat)
    mock_1 = model_1c.mock.galaxy_table
    mock_1.write(save_path+'mock_1c_07.hdf5', format='hdf5', path='data', overwrite=True)

    model_1c.param_dict['rho'] = 0.8
    model_1c.populate_mock(halocat)
    mock_1 = model_1c.mock.galaxy_table
    mock_1.write(save_path+'mock_1c_08.hdf5', format='hdf5', path='data', overwrite=True)

    model_1c.param_dict['rho'] = 0.9
    model_1c.populate_mock(halocat)
    mock_1 = model_1c.mock.galaxy_table
    mock_1.write(save_path+'mock_1c_09.hdf5', format='hdf5', path='data', overwrite=True)

    model_1c.param_dict['rho'] = 1.0
    model_1c.populate_mock(halocat)
    mock_1 = model_1c.mock.galaxy_table
    mock_1.write(save_path+'mock_1c_10.hdf5', format='hdf5', path='data', overwrite=True)



    model_1d.param_dict['rho'] = 0.0
    model_1d.populate_mock(halocat)
    mock_1 = model_1d.mock.galaxy_table
    mock_1.write(save_path+'mock_1d_00.hdf5', format='hdf5', path='data', overwrite=True)

    model_1d.param_dict['rho'] = -0.1
    model_1d.populate_mock(halocat)
    mock_1 = model_1d.mock.galaxy_table
    mock_1.write(save_path+'mock_1d_-01.hdf5', format='hdf5', path='data', overwrite=True)

    model_1d.param_dict['rho'] = -0.2
    model_1d.populate_mock(halocat)
    mock_1 = model_1d.mock.galaxy_table
    mock_1.write(save_path+'mock_1d_-02.hdf5', format='hdf5', path='data', overwrite=True)

    model_1d.param_dict['rho'] = -0.3
    model_1d.populate_mock(halocat)
    mock_1 = model_1d.mock.galaxy_table
    mock_1.write(save_path+'mock_1d_-03.hdf5', format='hdf5', path='data', overwrite=True)

    model_1d.param_dict['rho'] = -0.4
    model_1d.populate_mock(halocat)
    mock_1 = model_1d.mock.galaxy_table
    mock_1.write(save_path+'mock_1d_-04.hdf5', format='hdf5', path='data', overwrite=True)

    model_1d.param_dict['rho'] = -0.5
    model_1d.populate_mock(halocat)
    mock_1 = model_1d.mock.galaxy_table
    mock_1.write(save_path+'mock_1d_-05.hdf5', format='hdf5', path='data', overwrite=True)

    model_1d.param_dict['rho'] = -0.6
    model_1d.populate_mock(halocat)
    mock_1 = model_1d.mock.galaxy_table
    mock_1.write(save_path+'mock_1d_-06.hdf5', format='hdf5', path='data', overwrite=True)

    model_1d.param_dict['rho'] = -0.7
    model_1d.populate_mock(halocat)
    mock_1 = model_1d.mock.galaxy_table
    mock_1.write(save_path+'mock_1d_-07.hdf5', format='hdf5', path='data', overwrite=True)

    model_1d.param_dict['rho'] = -0.8
    model_1d.populate_mock(halocat)
    mock_1 = model_1d.mock.galaxy_table
    mock_1.write(save_path+'mock_1d_-08.hdf5', format='hdf5', path='data', overwrite=True)

    model_1d.param_dict['rho'] = -0.9
    model_1d.populate_mock(halocat)
    mock_1 = model_1d.mock.galaxy_table
    mock_1.write(save_path+'mock_1d_-09.hdf5', format='hdf5', path='data', overwrite=True)

    model_1d.param_dict['rho'] = -1.0
    model_1d.populate_mock(halocat)
    mock_1 = model_1d.mock.galaxy_table
    mock_1.write(save_path+'mock_1d_-10.hdf5', format='hdf5', path='data', overwrite=True)


    model_1e.param_dict['rho'] = 0.0
    model_1e.populate_mock(halocat)
    mock_1 = model_1e.mock.galaxy_table
    mock_1.write(save_path+'mock_1e_00.hdf5', format='hdf5', path='data', overwrite=True)

    model_1e.param_dict['rho'] = 0.1
    model_1e.populate_mock(halocat)
    mock_1 = model_1e.mock.galaxy_table
    mock_1.write(save_path+'mock_1e_01.hdf5', format='hdf5', path='data', overwrite=True)

    model_1e.param_dict['rho'] = 0.2
    model_1e.populate_mock(halocat)
    mock_1 = model_1e.mock.galaxy_table
    mock_1.write(save_path+'mock_1e_02.hdf5', format='hdf5', path='data', overwrite=True)

    model_1e.param_dict['rho'] = 0.3
    model_1e.populate_mock(halocat)
    mock_1 = model_1e.mock.galaxy_table
    mock_1.write(save_path+'mock_1e_03.hdf5', format='hdf5', path='data', overwrite=True)

    model_1e.param_dict['rho'] = 0.4
    model_1e.populate_mock(halocat)
    mock_1 = model_1e.mock.galaxy_table
    mock_1.write(save_path+'mock_1e_04.hdf5', format='hdf5', path='data', overwrite=True)

    model_1e.param_dict['rho'] = 0.5
    model_1e.populate_mock(halocat)
    mock_1 = model_1e.mock.galaxy_table
    mock_1.write(save_path+'mock_1e_05.hdf5', format='hdf5', path='data', overwrite=True)

    model_1e.param_dict['rho'] = 0.6
    model_1e.populate_mock(halocat)
    mock_1 = model_1e.mock.galaxy_table
    mock_1.write(save_path+'mock_1e_06.hdf5', format='hdf5', path='data', overwrite=True)

    model_1e.param_dict['rho'] = 0.7
    model_1e.populate_mock(halocat)
    mock_1 = model_1e.mock.galaxy_table
    mock_1.write(save_path+'mock_1e_07.hdf5', format='hdf5', path='data', overwrite=True)

    model_1e.param_dict['rho'] = 0.8
    model_1e.populate_mock(halocat)
    mock_1 = model_1e.mock.galaxy_table
    mock_1.write(save_path+'mock_1e_08.hdf5', format='hdf5', path='data', overwrite=True)

    model_1e.param_dict['rho'] = 0.9
    model_1e.populate_mock(halocat)
    mock_1 = model_1e.mock.galaxy_table
    mock_1.write(save_path+'mock_1e_09.hdf5', format='hdf5', path='data', overwrite=True)

    model_1e.param_dict['rho'] = 1.0
    model_1e.populate_mock(halocat)
    mock_1 = model_1e.mock.galaxy_table
    mock_1.write(save_path+'mock_1e_10.hdf5', format='hdf5', path='data', overwrite=True)


    model_1f.param_dict['rho'] = 0.0
    model_1f.populate_mock(halocat)
    mock_1 = model_1f.mock.galaxy_table
    mock_1.write(save_path+'mock_1f_00.hdf5', format='hdf5', path='data', overwrite=True)

    model_1f.param_dict['rho'] = 0.1
    model_1f.populate_mock(halocat)
    mock_1 = model_1f.mock.galaxy_table
    mock_1.write(save_path+'mock_1f_01.hdf5', format='hdf5', path='data', overwrite=True)

    model_1f.param_dict['rho'] = 0.2
    model_1f.populate_mock(halocat)
    mock_1 = model_1f.mock.galaxy_table
    mock_1.write(save_path+'mock_1f_02.hdf5', format='hdf5', path='data', overwrite=True)

    model_1f.param_dict['rho'] = 0.3
    model_1f.populate_mock(halocat)
    mock_1 = model_1f.mock.galaxy_table
    mock_1.write(save_path+'mock_1f_03.hdf5', format='hdf5', path='data', overwrite=True)

    model_1f.param_dict['rho'] = 0.4
    model_1f.populate_mock(halocat)
    mock_1 = model_1f.mock.galaxy_table
    mock_1.write(save_path+'mock_1f_04.hdf5', format='hdf5', path='data', overwrite=True)

    model_1f.param_dict['rho'] = 0.5
    model_1f.populate_mock(halocat)
    mock_1 = model_1f.mock.galaxy_table
    mock_1.write(save_path+'mock_1f_05.hdf5', format='hdf5', path='data', overwrite=True)

    model_1f.param_dict['rho'] = 0.6
    model_1f.populate_mock(halocat)
    mock_1 = model_1f.mock.galaxy_table
    mock_1.write(save_path+'mock_1f_06.hdf5', format='hdf5', path='data', overwrite=True)

    model_1f.param_dict['rho'] = 0.7
    model_1f.populate_mock(halocat)
    mock_1 = model_1f.mock.galaxy_table
    mock_1.write(save_path+'mock_1f_07.hdf5', format='hdf5', path='data', overwrite=True)

    model_1f.param_dict['rho'] = 0.8
    model_1f.populate_mock(halocat)
    mock_1 = model_1f.mock.galaxy_table
    mock_1.write(save_path+'mock_1f_08.hdf5', format='hdf5', path='data', overwrite=True)

    model_1f.param_dict['rho'] = 0.9
    model_1f.populate_mock(halocat)
    mock_1 = model_1f.mock.galaxy_table
    mock_1.write(save_path+'mock_1f_09.hdf5', format='hdf5', path='data', overwrite=True)

    model_1f.param_dict['rho'] = 1.0
    model_1f.populate_mock(halocat)
    mock_1 = model_1f.mock.galaxy_table
    mock_1.write(save_path+'mock_1f_10.hdf5', format='hdf5', path='data', overwrite=True)


    

    model_2a.param_dict['rho'] = 0.0
    model_2a.populate_mock(halocat)
    mock_2 = model_2a.mock.galaxy_table
    mock_2.write(save_path+'mock_2a_00.hdf5', format='hdf5', path='data', overwrite=True)

    model_2a.param_dict['rho'] = -0.1
    model_2a.populate_mock(halocat)
    mock_2 = model_2a.mock.galaxy_table
    mock_2.write(save_path+'mock_2a_-01.hdf5', format='hdf5', path='data', overwrite=True)

    model_2a.param_dict['rho'] = -0.2
    model_2a.populate_mock(halocat)
    mock_2 = model_2a.mock.galaxy_table
    mock_2.write(save_path+'mock_2a_-02.hdf5', format='hdf5', path='data', overwrite=True)

    model_2a.param_dict['rho'] = -0.3
    model_2a.populate_mock(halocat)
    mock_2 = model_2a.mock.galaxy_table
    mock_2.write(save_path+'mock_2a_-03.hdf5', format='hdf5', path='data', overwrite=True)

    model_2a.param_dict['rho'] = -0.4
    model_2a.populate_mock(halocat)
    mock_2 = model_2a.mock.galaxy_table
    mock_2.write(save_path+'mock_2a_-04.hdf5', format='hdf5', path='data', overwrite=True)

    model_2a.param_dict['rho'] = -0.5
    model_2a.populate_mock(halocat)
    mock_2 = model_2a.mock.galaxy_table
    mock_2.write(save_path+'mock_2a_-05.hdf5', format='hdf5', path='data', overwrite=True)

    model_2a.param_dict['rho'] = -0.6
    model_2a.populate_mock(halocat)
    mock_2 = model_2a.mock.galaxy_table
    mock_2.write(save_path+'mock_2a_-06.hdf5', format='hdf5', path='data', overwrite=True)

    model_2a.param_dict['rho'] = -0.7
    model_2a.populate_mock(halocat)
    mock_2 = model_2a.mock.galaxy_table
    mock_2.write(save_path+'mock_2a_-07.hdf5', format='hdf5', path='data', overwrite=True)

    model_2a.param_dict['rho'] = -0.8
    model_2a.populate_mock(halocat)
    mock_2 = model_2a.mock.galaxy_table
    mock_2.write(save_path+'mock_2a_-08.hdf5', format='hdf5', path='data', overwrite=True)

    model_2a.param_dict['rho'] = -0.9
    model_2a.populate_mock(halocat)
    mock_2 = model_2a.mock.galaxy_table
    mock_2.write(save_path+'mock_2a_-09.hdf5', format='hdf5', path='data', overwrite=True)

    model_2a.param_dict['rho'] = -1.0
    model_2a.populate_mock(halocat)
    mock_2 = model_2a.mock.galaxy_table
    mock_2.write(save_path+'mock_2a_-10.hdf5', format='hdf5', path='data', overwrite=True)



    model_2b.param_dict['rho'] = 0.0
    model_2b.populate_mock(halocat)
    mock_2 = model_2b.mock.galaxy_table
    mock_2.write(save_path+'mock_2b_00.hdf5', format='hdf5', path='data', overwrite=True)

    model_2b.param_dict['rho'] = 0.1
    model_2b.populate_mock(halocat)
    mock_2 = model_2b.mock.galaxy_table
    mock_2.write(save_path+'mock_2b_01.hdf5', format='hdf5', path='data', overwrite=True)

    model_2b.param_dict['rho'] = 0.2
    model_2b.populate_mock(halocat)
    mock_2 = model_2b.mock.galaxy_table
    mock_2.write(save_path+'mock_2b_02.hdf5', format='hdf5', path='data', overwrite=True)

    model_2b.param_dict['rho'] = 0.3
    model_2b.populate_mock(halocat)
    mock_2 = model_2b.mock.galaxy_table
    mock_2.write(save_path+'mock_2b_03.hdf5', format='hdf5', path='data', overwrite=True)

    model_2b.param_dict['rho'] = 0.4
    model_2b.populate_mock(halocat)
    mock_2 = model_2b.mock.galaxy_table
    mock_2.write(save_path+'mock_2b_04.hdf5', format='hdf5', path='data', overwrite=True)

    model_2b.param_dict['rho'] = 0.5
    model_2b.populate_mock(halocat)
    mock_2 = model_2b.mock.galaxy_table
    mock_2.write(save_path+'mock_2b_05.hdf5', format='hdf5', path='data', overwrite=True)

    model_2b.param_dict['rho'] = 0.6
    model_2b.populate_mock(halocat)
    mock_2 = model_2b.mock.galaxy_table
    mock_2.write(save_path+'mock_2b_06.hdf5', format='hdf5', path='data', overwrite=True)

    model_2b.param_dict['rho'] = 0.7
    model_2b.populate_mock(halocat)
    mock_2 = model_2b.mock.galaxy_table
    mock_2.write(save_path+'mock_2b_07.hdf5', format='hdf5', path='data', overwrite=True)

    model_2b.param_dict['rho'] = 0.8
    model_2b.populate_mock(halocat)
    mock_2 = model_2b.mock.galaxy_table
    mock_2.write(save_path+'mock_2b_08.hdf5', format='hdf5', path='data', overwrite=True)

    model_2b.param_dict['rho'] = 0.9
    model_2b.populate_mock(halocat)
    mock_2 = model_2b.mock.galaxy_table
    mock_2.write(save_path+'mock_2b_09.hdf5', format='hdf5', path='data', overwrite=True)

    model_2b.param_dict['rho'] = 1.0
    model_2b.populate_mock(halocat)
    mock_2 = model_2b.mock.galaxy_table
    mock_2.write(save_path+'mock_2b_10.hdf5', format='hdf5', path='data', overwrite=True)



    model_2c.param_dict['rho'] = 0.0
    model_2c.populate_mock(halocat)
    mock_2 = model_2c.mock.galaxy_table
    mock_2.write(save_path+'mock_2c_00.hdf5', format='hdf5', path='data', overwrite=True)

    model_2c.param_dict['rho'] = 0.1
    model_2c.populate_mock(halocat)
    mock_2 = model_2c.mock.galaxy_table
    mock_2.write(save_path+'mock_2c_01.hdf5', format='hdf5', path='data', overwrite=True)

    model_2c.param_dict['rho'] = 0.2
    model_2c.populate_mock(halocat)
    mock_2 = model_2c.mock.galaxy_table
    mock_2.write(save_path+'mock_2c_02.hdf5', format='hdf5', path='data', overwrite=True)

    model_2c.param_dict['rho'] = 0.3
    model_2c.populate_mock(halocat)
    mock_2 = model_2c.mock.galaxy_table
    mock_2.write(save_path+'mock_2c_03.hdf5', format='hdf5', path='data', overwrite=True)

    model_2c.param_dict['rho'] = 0.4
    model_2c.populate_mock(halocat)
    mock_2 = model_2c.mock.galaxy_table
    mock_2.write(save_path+'mock_2c_04.hdf5', format='hdf5', path='data', overwrite=True)

    model_2c.param_dict['rho'] = 0.5
    model_2c.populate_mock(halocat)
    mock_2 = model_2c.mock.galaxy_table
    mock_2.write(save_path+'mock_2c_05.hdf5', format='hdf5', path='data', overwrite=True)

    model_2c.param_dict['rho'] = 0.6
    model_2c.populate_mock(halocat)
    mock_2 = model_2c.mock.galaxy_table
    mock_2.write(save_path+'mock_2c_06.hdf5', format='hdf5', path='data', overwrite=True)

    model_2c.param_dict['rho'] = 0.7
    model_2c.populate_mock(halocat)
    mock_2 = model_2c.mock.galaxy_table
    mock_2.write(save_path+'mock_2c_07.hdf5', format='hdf5', path='data', overwrite=True)

    model_2c.param_dict['rho'] = 0.8
    model_2c.populate_mock(halocat)
    mock_2 = model_2c.mock.galaxy_table
    mock_2.write(save_path+'mock_2c_08.hdf5', format='hdf5', path='data', overwrite=True)

    model_2c.param_dict['rho'] = 0.9
    model_2c.populate_mock(halocat)
    mock_2 = model_2c.mock.galaxy_table
    mock_2.write(save_path+'mock_2c_09.hdf5', format='hdf5', path='data', overwrite=True)

    model_2c.param_dict['rho'] = 1.0
    model_2c.populate_mock(halocat)
    mock_2 = model_2c.mock.galaxy_table
    mock_2.write(save_path+'mock_2c_10.hdf5', format='hdf5', path='data', overwrite=True)



    model_2d.param_dict['rho'] = 0.0
    model_2d.populate_mock(halocat)
    mock_2 = model_2d.mock.galaxy_table
    mock_2.write(save_path+'mock_2d_00.hdf5', format='hdf5', path='data', overwrite=True)

    model_2d.param_dict['rho'] = -0.1
    model_2d.populate_mock(halocat)
    mock_2 = model_2d.mock.galaxy_table
    mock_2.write(save_path+'mock_2d_-01.hdf5', format='hdf5', path='data', overwrite=True)

    model_2d.param_dict['rho'] = -0.2
    model_2d.populate_mock(halocat)
    mock_2 = model_2d.mock.galaxy_table
    mock_2.write(save_path+'mock_2d_-02.hdf5', format='hdf5', path='data', overwrite=True)

    model_2d.param_dict['rho'] = -0.3
    model_2d.populate_mock(halocat)
    mock_2 = model_2d.mock.galaxy_table
    mock_2.write(save_path+'mock_2d_-03.hdf5', format='hdf5', path='data', overwrite=True)

    model_2d.param_dict['rho'] = -0.4
    model_2d.populate_mock(halocat)
    mock_2 = model_2d.mock.galaxy_table
    mock_2.write(save_path+'mock_2d_-04.hdf5', format='hdf5', path='data', overwrite=True)

    model_2d.param_dict['rho'] = -0.5
    model_2d.populate_mock(halocat)
    mock_2 = model_2d.mock.galaxy_table
    mock_2.write(save_path+'mock_2d_-05.hdf5', format='hdf5', path='data', overwrite=True)

    model_2d.param_dict['rho'] = -0.6
    model_2d.populate_mock(halocat)
    mock_2 = model_2d.mock.galaxy_table
    mock_2.write(save_path+'mock_2d_-06.hdf5', format='hdf5', path='data', overwrite=True)

    model_2d.param_dict['rho'] = -0.7
    model_2d.populate_mock(halocat)
    mock_2 = model_2d.mock.galaxy_table
    mock_2.write(save_path+'mock_2d_-07.hdf5', format='hdf5', path='data', overwrite=True)

    model_2d.param_dict['rho'] = -0.8
    model_2d.populate_mock(halocat)
    mock_2 = model_2d.mock.galaxy_table
    mock_2.write(save_path+'mock_2d_-08.hdf5', format='hdf5', path='data', overwrite=True)

    model_2d.param_dict['rho'] = -0.9
    model_2d.populate_mock(halocat)
    mock_2 = model_2d.mock.galaxy_table
    mock_2.write(save_path+'mock_2d_-09.hdf5', format='hdf5', path='data', overwrite=True)

    model_2d.param_dict['rho'] = -1.0
    model_2d.populate_mock(halocat)
    mock_2 = model_2d.mock.galaxy_table
    mock_2.write(save_path+'mock_2d_-10.hdf5', format='hdf5', path='data', overwrite=True)


    model_2e.param_dict['rho'] = 0.0
    model_2e.populate_mock(halocat)
    mock_2 = model_2e.mock.galaxy_table
    mock_2.write(save_path+'mock_2e_00.hdf5', format='hdf5', path='data', overwrite=True)

    model_2e.param_dict['rho'] = 0.1
    model_2e.populate_mock(halocat)
    mock_2 = model_2e.mock.galaxy_table
    mock_2.write(save_path+'mock_2e_01.hdf5', format='hdf5', path='data', overwrite=True)

    model_2e.param_dict['rho'] = 0.2
    model_2e.populate_mock(halocat)
    mock_2 = model_2e.mock.galaxy_table
    mock_2.write(save_path+'mock_2e_02.hdf5', format='hdf5', path='data', overwrite=True)

    model_2e.param_dict['rho'] = 0.3
    model_2e.populate_mock(halocat)
    mock_2 = model_2e.mock.galaxy_table
    mock_2.write(save_path+'mock_2e_03.hdf5', format='hdf5', path='data', overwrite=True)

    model_2e.param_dict['rho'] = 0.4
    model_2e.populate_mock(halocat)
    mock_2 = model_2e.mock.galaxy_table
    mock_2.write(save_path+'mock_2e_04.hdf5', format='hdf5', path='data', overwrite=True)

    model_2e.param_dict['rho'] = 0.5
    model_2e.populate_mock(halocat)
    mock_2 = model_2e.mock.galaxy_table
    mock_2.write(save_path+'mock_2e_05.hdf5', format='hdf5', path='data', overwrite=True)

    model_2e.param_dict['rho'] = 0.6
    model_2e.populate_mock(halocat)
    mock_2 = model_2e.mock.galaxy_table
    mock_2.write(save_path+'mock_2e_06.hdf5', format='hdf5', path='data', overwrite=True)

    model_2e.param_dict['rho'] = 0.7
    model_2e.populate_mock(halocat)
    mock_2 = model_2e.mock.galaxy_table
    mock_2.write(save_path+'mock_2e_07.hdf5', format='hdf5', path='data', overwrite=True)

    model_2e.param_dict['rho'] = 0.8
    model_2e.populate_mock(halocat)
    mock_2 = model_2e.mock.galaxy_table
    mock_2.write(save_path+'mock_2e_08.hdf5', format='hdf5', path='data', overwrite=True)

    model_2e.param_dict['rho'] = 0.9
    model_2e.populate_mock(halocat)
    mock_2 = model_2e.mock.galaxy_table
    mock_2.write(save_path+'mock_2e_09.hdf5', format='hdf5', path='data', overwrite=True)

    model_2e.param_dict['rho'] = 1.0
    model_2e.populate_mock(halocat)
    mock_2 = model_2e.mock.galaxy_table
    mock_2.write(save_path+'mock_2e_10.hdf5', format='hdf5', path='data', overwrite=True)


    model_2f.param_dict['rho'] = 0.0
    model_2f.populate_mock(halocat)
    mock_2 = model_2f.mock.galaxy_table
    mock_2.write(save_path+'mock_2f_00.hdf5', format='hdf5', path='data', overwrite=True)

    model_2f.param_dict['rho'] = 0.1
    model_2f.populate_mock(halocat)
    mock_2 = model_2f.mock.galaxy_table
    mock_2.write(save_path+'mock_2f_01.hdf5', format='hdf5', path='data', overwrite=True)

    model_2f.param_dict['rho'] = 0.2
    model_2f.populate_mock(halocat)
    mock_2 = model_2f.mock.galaxy_table
    mock_2.write(save_path+'mock_2f_02.hdf5', format='hdf5', path='data', overwrite=True)

    model_2f.param_dict['rho'] = 0.3
    model_2f.populate_mock(halocat)
    mock_2 = model_2f.mock.galaxy_table
    mock_2.write(save_path+'mock_2f_03.hdf5', format='hdf5', path='data', overwrite=True)

    model_2f.param_dict['rho'] = 0.4
    model_2f.populate_mock(halocat)
    mock_2 = model_2f.mock.galaxy_table
    mock_2.write(save_path+'mock_2f_04.hdf5', format='hdf5', path='data', overwrite=True)

    model_2f.param_dict['rho'] = 0.5
    model_2f.populate_mock(halocat)
    mock_2 = model_2f.mock.galaxy_table
    mock_2.write(save_path+'mock_2f_05.hdf5', format='hdf5', path='data', overwrite=True)

    model_2f.param_dict['rho'] = 0.6
    model_2f.populate_mock(halocat)
    mock_2 = model_2f.mock.galaxy_table
    mock_2.write(save_path+'mock_2f_06.hdf5', format='hdf5', path='data', overwrite=True)

    model_2f.param_dict['rho'] = 0.7
    model_2f.populate_mock(halocat)
    mock_2 = model_2f.mock.galaxy_table
    mock_2.write(save_path+'mock_2f_07.hdf5', format='hdf5', path='data', overwrite=True)

    model_2f.param_dict['rho'] = 0.8
    model_2f.populate_mock(halocat)
    mock_2 = model_2f.mock.galaxy_table
    mock_2.write(save_path+'mock_2f_08.hdf5', format='hdf5', path='data', overwrite=True)

    model_2f.param_dict['rho'] = 0.9
    model_2f.populate_mock(halocat)
    mock_2 = model_2f.mock.galaxy_table
    mock_2.write(save_path+'mock_2f_09.hdf5', format='hdf5', path='data', overwrite=True)

    model_2f.param_dict['rho'] = 1.0
    model_2f.populate_mock(halocat)
    mock_2 = model_2f.mock.galaxy_table
    mock_2.write(save_path+'mock_2f_10.hdf5', format='hdf5', path='data', overwrite=True)


if __name__ == '__main__':
    main()
