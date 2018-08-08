# ECO/RESOLVE Mocks

This Repository contains tools to build mock galaxy catalogs for the RESOLVE and ECO [Surveys](https://resolve.astro.unc.edu).

## Dependancies
In order to run the code in this repository, you will need the following packages in your Python path:

* [Halotools](https://halotools.readthedocs.io/en/latest/)
* [AbundanceMatching](https://bitbucket.org/yymao/abundancematching)

## Making Mocks

First, follow the steps in the `README.md` located in the `/data` directory.  This will walk you through downloading and processing the halo catalog(s) and survey data needed to build mock catalogs. 

Next, you can create a set of fiducial mocks by running the `make_all_eco_mocks.py` Python script.  Mocks created using this script are stored in the `/mocks` directory as hdf5 files.

Example notebooks demonstrating the models are located in the `/notebooks` directory.  This directory also contains an example opening and plotting some data using the fiducial mocks.
 

Author: Duncan Campbell