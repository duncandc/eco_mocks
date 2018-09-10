# ECO Mocks

This directory is used to store mocks.


## Available Fiducial Mocks

I have listed the mock filename names and the galaxy/halo properties used in the conditional abundance matching (CAM).
The format is:
(primary galaxy, primary halo, secondary galaxy, secondary halo, correlation stength)

- mock_1a (stellar mass, halo\_vpeak, u-r,    halo\_half\_mass\_scale, -1)
- mock_1b (stellar mass, halo\_vpeak, fsmgr,  halo\_half\_mass\_scale,  1)
- mock_1c (stellar mass, halo\_vpeak, fgas,   halo\_half\_mass\_scale,  1)
- mock_1d (stellar mass, halo\_vpeak, u-r,    halo\_acc\_rate\_100myr, -1)
- mock_1e (stellar mass, halo\_vpeak, fsmgr,  halo\_acc\_rate\_100myr,  1)
- mock_1f (stellar mass, halo\_vpeak, fgas,   halo\_acc\_rate\_100myr,  1)
- mock_2a (baryonic mass, halo\_vpeak, u-r,   halo\_half\_mass\_scale, -1)
- mock_2b (baryonic mass, halo\_vpeak, fsmgr, halo\_half\_mass\_scale,  1)
- mock_2c (baryonic mass, halo\_vpeak, fgas,  halo\_half\_mass\_scale,  1)
- mock_2d (baryonic mass, halo\_vpeak, u-r,   halo\_acc\_rate\_100myr, -1)
- mock_2e (baryonic mass, halo\_vpeak, fsmgr, halo\_acc\_rate\_100myr,  1)
- mock_2f (baryonic mass, halo\_vpeak, fgas,  halo\_acc\_rate\_100myr,  1)


## Columns

|column                       |description                  |units      |
|-----------------------------|-----------------------------|-----------|
|halo_id							|halo id                      |           |
|halo_upid						|halo uber id                 |           |
|halo_hostid						|halo id of host halo         |           |
|halo_x							|x-position of (sub-)halo     | h^-1 Mpc  |
|halo_y							|y-position of (sub-)halo     | h^-1 Mpc  |
|halo_z							|z-position of (sub-)halo     | h^-1 Mpc  |
|halo_vx							|x-velocity of (sub-)halo     | km/s      |
|halo_vy							|y-velocity of (sub-)halo     | km/s      |
|halo_vz							|z-velocity of (sub-)halo     | km/s      |
|halo_mvir						|instanateous (sub-)halo mass | h^-1 Msol |
|halo\_mvir\_host\_halo			|mass of host halo            | h^-1 Msol |
|halo_mpeak						|peak mass                    | h^-1 Msol |
|halo_vpeak						|peak maximum circular velcity| km/s      |
|halo_rvir						|virial radius                | h^-1 Mpc  |
|halo\_acc\_rate\_100my			|halo accretion rate          | ????      | 
|halo\_half\_mass\_scale		|formation time               |           |
|galid							   	|galaxy id                    |           |
|x									|x-position of galaxy         | h^-1 mpc  |
|y									|y-position of galaxy         | h^-1 mpc  |
|z									|z-position of galaxy         | h^-1 mpc  |
|vx									|x-velocity of galaxy         | km/s      |
|vy									|y-velocity of galaxy         | km/s      |
|vz									|z-velocity of galaxy         | km/s      |
|stellar_mass						|stellar mass                 | h^-2 Msol |
|baryonic_mass					|baryonic mass                | h^-2 Msol |
|gas_mass							|gas mass                     | h^-2 Msol |
|abs_umag							|absolute u-band magnitude    |           |
|abs_gmag							|absolute g-band magnitude    |           |
|abs_rmag							|absolute r-band magnitude    |           |
|g\_minus\_r						|g-r color                    |           |
|u\_minus\_r						|u-r color                    |           |
|fsmgr								|stellar mass growth rate     | ????      |
|b\_to\_a							|b/a axis ratio [0,1]         |           |
|bijected\_b\_to\_a				|remaped b/a axis ratio       |           |
|fgas								|gas fraction [0,1]           |           |
|bijected_fgas					|remapped gas fraction        |           |
|num\_morph\_type				|morphological type           |           |
|dgr								|g-r color gradient           |           |
|reference_idx					|index into eco catalog       |           |


## Notes

Tertiary galaxy properties (those other than the ones used in CAM) are assigned to mock galaxies by matching to the nearest neghbor based on the primary galaxy property into a reference ECO catalog, e.g. stellar mass.  This method preserves all the conditional distributions and correlations between galaxy properties as well as possible.  However, this does significanlty over-sample from the reference catalog.  

The gas fraction is used to assign baryonic mass and gas mass if stellar mass is assigined as the primary galaxy property.  Similarly, the gas fraction is used to assign stellar mass and gas mass if baryonic mass is assigined as the primary galaxy property.  this ensures that the baryonic mass is the sum of the stellar and gas masses for any individual mock galaxy.

For an example of how to load and access mocks, see the 'load_eco_mocks.ipynb' ipython notebook in the 'notebooks' directory.