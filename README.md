**dragn_hunter.py readme**

command line code for searching for pairs of radio components that satisfy specified criteria, and flagging likely Double Radio sources associated with Active Galactic Nuclei (DRAGNs). Additionally candidate hosts are identified by querying the CatWISE2020 catalog using the flux weighted centroid of the double as the reference position. Run as:

*python3 dragn_hunter.py radio_catalogue*

where *radio_catalogue* is the filename of the radio catalogue data file. By default this code assumes that a config file named *config.txt* exists but a different config file can be called by *--config=filename*.
The output consists of two files:

1) DRAGNs.fits is a table of unique radio pairs (candidate DRAGNs) and (if any) the closest CatWISE source to the flux weighted centroid of the pair.
2) host_candidates.fits is a table of all CatWISE sources within 10 arcsec of the flux weighted centroids of the radio pairs.

These two tables are joinable on the column name 'pair_name', and are detailed below. If importing the functionality from another code or using iPython in the command line, and outputting a file is not desired, then the option *--writepairs='False'* can be used to return an astropy table of the data instead.


**Config File Description**

parameter | description
----------|------------
name_col | name of the object name/id column in the radio catalogue file
ra_col | name of the RA column in the radio catalogue file (assumed to be in decimal degrees)
dec_col | name of the Dec column in the radio catalogue file (assumed to be in decimal degrees)
Speak_col | name of the peak flux column in the radio catalogue file
Stot_col | name of the total flux column in the radio catalogue file
maj_col | name of the major axis size column in the radio catalogue file
min_col | name of the minor axis size column in the radio catalogue file
pa_col | name of the position angle column in the radio catalogue file (assumed to be degrees)
min_brightness | minimum peak flux of radio sources to use in pair finding
min_size | minimum size of radio sources to use in pair finding
data_format | format of radio catalogue file


**Pair Data (DRAGNs.fits) Description**

the output data is a list of pairs where each **pair** appears only once - i.e. component_1 of a pair does not appear later on as component_2 of the pair. For each component of the pair, the columns specified in the config file will be returned in the pair data with _1 or _2 appended to the column name. In addition several derived columns are given and are described below.

parameter | description [units]
----------|------------
Sep_pair | the angular separation of the pair from the catalogue positions of the constituent components [arcsec]
PA_pair | the position angle East of North from component_1 to component_2 in the pair [deg]
dPA_n | difference in position angle of component_n to PA_pair [deg]
abs_dPA | magnitude of dPA_n
Tflux_ratio | ratio of total flux of component_1 to component_2
Pflux_ratio | ratio of peak flux of component_1 to component_2
medianRA | median (not flux weighted) RA of the pair [deg]
medianDEC | median (not flux weighted) DEC of the pair [deg]
cenRA | RA of the pair's flux weighted centroid [deg]
cenDEC | DEC of the pair's flux weighted centroid [deg]
d_2NN | distance to the second nearest neighbour of the pair [arcsec]
pair_name | identifier for the pair of the format Jhhmmss.ss±ddmmss.s
CWISE | name of the closest CatWISE2020 source to the flux weighted centroid
snrW1pm | signal/Noise in W1-band for CatWISE source
snrW2pm | signal/Noise in W2-band for CatWISE source
W1mproPM | magnitude (Vega) in W1 of CatWISE source [mag]
W2mproPM | magnitude (Vega) in W2 of CatWISE source [mag]
sep_from_radio_source | angular separation between the flux weighted centroid of the radio pair and the CatWISE source [arcsec]


**Host Data (host_candidates.fits) Description**

A list of **all** CatWISE2020 sources within 10 arcsec of the flux weighted centroid of radio pairs. CatWISE data are the default columns returned by querying https://cdsarc.unistra.fr/viz-bin/cat/II/365

parameter | description [units]
----------|------------
pair_name | identifier for the pair of the format Jhhmmss.ss±ddmmss.s
sep_from_radio_source | angular separation between the flux weighted centroid of the radio pair and the CatWISE source [arcsec]
RA_ICRS | Right Ascension of the CatWISE source [deg]
e_RA_ICRS | error in RA_ICRS [arcsec]
DE_ICRS | Declination of the CatWISE source [deg]
e_DE_ICRS | error in DEC_ICRS [arcsec]
CWISE | name of the closest CatWISE2020 source to the flux weighted centroid
nW1 | number of profile-fit flux measurements with S/N > 3 in W1-band
mW1 | number of profile-fit flux measurements in W1-band
nW2 | number of profile-fit flux measurements with S/N > 3 in W2-band
mW2 | number of profile-fit flux measurements in W2-band
RAPMdeg | Right Ascension at epoch=2015.4 [deg]
e_RAPMdeg | uncertainty in RAPMdeg [arcsec]
DEPMdeg | Declination at epoch=2015.4 [deg]
e_DEPMdeg | uncertainty in DEPMdeg [arcsec]
pmRA | proper motion in Right Ascension [arcsec/year]
e_pmRA | uncertainty in pmRA [arcsec/year]
pmDE | proper motion in Declination [arcsec/year]
e_pmDE | uncertainty in pmDE [arcsec/year]
snrW1pm | signal/Noise in W1-band for CatWISE source
snrW2pm | signal/Noise in W2-band for CatWISE source
FW1pm | raw flux in W1-band [count]
e_FW1pm | uncertainty in FW1pm [count]
FW2pm | raw flux in W2-band [count]
e_FW2pm | uncertainty in FW2pm [count]
W1mproPM | magnitude (Vega) in W1 of CatWISE source [mag]
W2mproPM | magnitude (Vega) in W2 of CatWISE source [mag]
pmQual | Quality of the proper motion solution
Dist | angular separation between apparitions
snrdELON | signal/noise for descending-ascending ecliptic longitude
snrdELAT | signal/noise for descending-ascending ecliptic latitude
chi2pmRA | chi-square for pmRA difference
chi2pmDE | chi-square for pmDE difference
ka | astrometry usage code
k1 | W1 photometry usage code
k2 | W2 photometry usage code
km | proper motion usage code
Sep | angular separation between CatWISE and AllWISE source if available [arcsec]
abf | AllWISE flag for W1 and W2
WISEA | AllWISE source within 3 arcsec


**Code Dependencies**
The following python packages are required to run this code (version used in development):
* numpy (1.17.4)
* pandas (1.0.3)
* argparse (1.1)
* astropy (3.2)
* astroquery (0.4.3) *note: will work with astroquery version 0.4.2 but not 0.4.1*
