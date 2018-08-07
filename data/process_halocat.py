"""
process ROCKSTAR halo catalogues from Vishnu
"""

# load packages
from __future__ import print_function, division
import sys
import os
from halotools import sim_manager

# current directory
base_save_dir = os.path.dirname(os.path.realpath(__file__)) + '/'

__author__ = ['Duncan Campbell']


def main():
    """
    Run this script to process ROCKSTAR halotools catalog.
    The script an be run by typing in the terminal:
        user$ python process_halocat.py
    """

    # set the ascii file to process
    filepath = './'
    if len(sys.argv) > 1:
        filename = sys.argv[1]
    else:
        filename = 'Vishnu_z0p000.hlist.gz'

    if filename[-3:] == '.gz':
        was_zipped = True
        os.system("gunzip -f " + filepath + filename)
        filename = filename[:-3]
    else:
        was_zipped = False

    # set some properties
    simname = 'vishnu_130'
    version = 'custom'
    Lbox = 130.0
    particle_mass = 3.215*10**7.0
    halo_finder = 'Rockstar'

    # set the location and filename of the reduced catalogue
    savepath = base_save_dir
    savename = filename + '_' + version + '.hdf5'

    # extract the scale factor of the snapshot from the filename
    # scale_factor = float(re.findall(r"[-+]?\d*\.\d+|\d+",filename)[0])
    scale_factor = 1.0
    redshift = 1.0/scale_factor - 1.0

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
                            'halo_m200b':           (39, 'f4'),
                            'halo_m200c':           (40, 'f4'),
                            'halo_T/|U|':           (56, 'f4'),
                            'halo_macc':            (60, 'f4'),
                            'halo_mpeak':           (61, 'f4'),
                            'halo_vacc':            (62, 'f4'),
                            'halo_vpeak':           (63, 'f4'),
                            'halo_half_mass_scale': (64, 'f4'),
                            'halo_acc_rate_inst':   (65, 'f4'),
                            'halo_acc_rate_100myr': (66, 'f4'),
                            'halo_acc_rate_1tdyn':  (67, 'f4'),
                            'halo_acc_rate_2tdyn':  (68, 'f4'),
                            'halo_acc_rate_mpeak':  (69, 'f4'),
                            'halo_mpeak_scale':     (70, 'f4'),
                            'halo_acc_scale':       (71, 'f4'),
                            'halo_first_acc_scale': (72, 'f4'),
                            'halo_first_acc_mvir':  (73, 'f4'),
                            'halo_first_acc_vmax':  (74, 'f4'),
                            'halo_vmax_at_mpeak':   (75, 'f4'),
                            }

    columns_to_convert_from_kpc_to_mpc = ['halo_rvir', 'halo_rs']

    # apply cuts to catalogue
    row_cut_min_dict = {'halo_mpeak': particle_mass*50}
    processing_notes = ("all halos with halo_mpeak < 50 times the particle \n"
                        "mass were discarded during the catalogue reduction.")

    # read in catalogue and save results
    reader = sim_manager.RockstarHlistReader(filepath+filename,
                                             columns_to_keep_dict,
                                             savepath+savename, simname,
                                             halo_finder, redshift, version,
                                             Lbox, particle_mass,
                                             row_cut_min_dict=row_cut_min_dict,
                                             processing_notes=processing_notes,
                                             overwrite=True)

    reader.read_halocat(columns_to_convert_from_kpc_to_mpc,
                        write_to_disk=True,
                        update_cache_log=True)

    # if the original hlist file was zipped, re-zip it.
    if was_zipped:
        os.system("gzip " + filepath + filename)


if __name__ == '__main__':
    main()
