# Data

* To download the halo catalog run the `download_data.sh` script.  This will save the unprocessed ROCKSTAR halo catalog for the z=0.0 snapshot of the Vishnu smulation to this directory (~500 Mb).

* To process the halo catalog into a Halotools formatted and cached catalog, run the `process_halocat.py` Python script.  This will save the  formatted (hdf5) catalog to this directoy (~500 Mb).  Halootls remembers the location of this file by storing it: `~/.astropy/cache/halotools/halo_table_cache_log.txt`.  To change the location where this file is stored, change the `base_save_dir` global variable in the 'process_halocat.py' script. 

* To download the ECO data, submit a SQL query [here](http://resolve.renci.org/eco/public/query.py/doit).  The script used to downlaod data for this project is located in `sql_query.txt`.  Note that a username and password is needed in order to access the query page.  The data should be saved as a .csv file called `eco_data.csv` and saved to this directory.