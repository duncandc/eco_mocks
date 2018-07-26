"""
process ROCKSTAR halo catalogues from Vishnu
"""

#load packages
from __future__ import print_function, division

import numpy as np
import re
import sys
import os
from halotools import sim_manager

# current directory
base_save_dir = os.path.dirname(os.path.realpath(__file__))

def main():
    
    #set the ascii file to process
    filepath = '../halo_catalogues/'
    if len(sys.argv)>1:
        filename = sys.argv[1]
    else:
        filename = 'Vishnu_z0p000.hlist.gz'
    
    
    if filename[-3:] == '.gz':
        was_zipped = True
        os.system("gunzip -f " + filepath + filename)
        filename = filename[:-3]
    else:
        was_zipped = False
    
    #set some properties
    simname = 'vishnu_130'
    version = 'custom'
    Lbox = 130.0
    particle_mass = 3.215*10**7.0
    halo_finder='Rockstar'
    
    #set the location and filename of the reduced catalogue
    savepath = base_save_dir+'Bolshoi/'
    savename = filename + '_' + version + '.hdf5'
    
    #extract the scale factor of the snapshot from the filename
    scale_factor = float(re.findall(r"[-+]?\d*\.\d+|\d+",filename)[0])
    redshift = 1.0/scale_factor -1.0
    
    columns_to_keep_dict = {'halo_id':              (1, 'i8'),
                            'halo_pid':             (5, 'i8'),
                            'halo_upid':            (6, 'i8'),
                            'halo_mvir':            (10, 'f4'),
                            'halo_rvir':            (11, 'f4'),
                            'halo_rs':              (12, 'f4'),
                            'halo_vmax':            (16, 'f4'),
                            'halo_x':               (17, 'f4'),
                            'halo_y':               (18, 'f4'),
                            'halo_z':               (19, 'f4'),
                            'halo_vx':              (20, 'f4'),
                            'halo_vy':              (21, 'f4'),
                            'halo_vz':              (22, 'f4'),
                            'halo_m200b':           (37, 'f4'),
                            'halo_m200c':           (38, 'f4'),
                            'halo_T/|U|':           (54, 'f4'),
                            'halo_macc':            (57, 'f4'),
                            'halo_mpeak':           (58, 'f4'),
                            'halo_vacc':            (59, 'f4'),
                            'halo_vpeak':           (60, 'f4'),
                            'halo_half_mass_scale': (61, 'f4'),
                            'halo_mpeak_scale':     (67, 'f4'),
                            'halo_acc_scale':       (68, 'f4'),
                            'halo_first_acc_scale': (69, 'f4'),
                            'halo_first_acc_mvir':  (70, 'f4'),
                            'halo_first_acc_vmax':  (71, 'f4'),
                            'halo_vmax_at_mpeak':   (72, 'f4'),
                            }
    
    columns_to_convert_from_kpc_to_mpc = ['halo_rvir','halo_rs']
    
    #apply cuts to catalogue
    row_cut_min_dict = {'halo_mpeak': particle_mass*50}
    processing_notes = ("all halos with halo_mpeak < 50 times the particle mass were \n"
                        "thrown out during the initial catalogue reduction.")
    
    #read in catalogue and save results
    reader = sim_manager.RockstarHlistReader(filepath+filename, columns_to_keep_dict,\
        savepath+savename, simname, halo_finder, redshift, version, Lbox, particle_mass,\
        row_cut_min_dict=row_cut_min_dict, processing_notes=processing_notes,\
        overwrite=True) 
    
    reader.read_halocat(columns_to_convert_from_kpc_to_mpc, write_to_disk = True, update_cache_log = True)
    
    if was_zipped:
        os.system("gzip " + filepath + filename)


if __name__ == '__main__':
    main()