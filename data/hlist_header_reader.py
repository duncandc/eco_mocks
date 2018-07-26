"""
read the header in the ROCKSTAR halo catalogue ascii files to make processing easier
"""

#load packages
from __future__ import print_function, division
import re
import sys

__all__ = ['read_header']

def main():
    """
    print out the column names and number
    """
    
    #set the ascii file to process
    if len(sys.argv)>1:
        fname = sys.argv[1]
    else:
        filepath = './'
        filename = 'Vishnu_z0p000.hlist'
        fname = filepath + filename
    
    col_names = read_header(fname)
    for i,name in enumerate(col_names):
        print(i, name)


def read_header(fname):
    """
    read the first line of ROCKSTAR hlist files and return a list of column names
    """
    
    with open(fname, 'r') as f:
        first_line = f.readline()
        cols_info = first_line.split(' ')
        col_names = []
        for col_info in cols_info:
            col_name = col_info.split('(')[0].strip('#')
            col_names.append(col_name)
    return col_names

if __name__ == '__main__':
    main()