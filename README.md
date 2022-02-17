**dragn_hunter.py readme**

command line code for searching for pairs of radio components that satisfy specified criteria, and flagging likely Double Radio sources associated with Active Galactic Nuclei (DRAGNs). Additionally candidate hosts are identified by querying the CatWISE2020 catalog using the flux weighted centroid of the double as the reference position. Run as:

    >python3 dragn_hunter.py radio_catalogue

where *radio_catalogue* is the filename of the radio catalogue data file. 

Additionally, there are a number of optional commands, these are [default values in square brackets]:\\
*--config* [config.txt], config file describing the input component table (see below);\\
*--findhosts* [True], tell the code whether to query the AllWISE source catalog for potential host candidates;\\
*--search_radius* [30arcsec], maximum radius to use when searching for radio cores (and AllWISE sources if *--findhosts=True*).


The output consists of two to three main files depending on whether *--findhosts=True*. In all cases the following two files are produced:
1) *DRAGNs.fits* is a table of unique radio pairs (candidate DRAGNs) and (if any) the closest CatWISE source to the flux weighted centroid of the pair
2) *candidate_cores.fits* is a table of all potential radio cores found for each pair
When *--findhosts==True* a third file is produced:
3) *host_candidates.fits* is a table of all the AllWISE sources within  the specified seach radius of the pair. Where a radio core has been found, the core position is used to search for AllWISE counterparts, otherwised the flux-weighted centroid of the pair is used.

These three tables are joinable on the column name *pair_name*, and are detailed below.


**Config File Description**

parameter | description
----------|------------
name_col | name of the object name/id column in the radio catalogue file
ra_col | name of the RA column in the radio catalogue file (assumed to be in decimal degrees)
dec_col | name of the Dec column in the radio catalogue file (assumed to be in decimal degrees)
Speak_col | name of the peak flux column in the radio catalogue file
Speak_err_col | name of the peak flux error column in the radio catalogue file
Stot_col | name of the total flux column in the radio catalogue file
Stot_err_col | name of the total flux error column in the radio catalogue file
maj_col | name of the major axis size column in the radio catalogue file
min_col | name of the minor axis size column in the radio catalogue file
pa_col | name of the position angle column in the radio catalogue file 
min_brightness | minimum peak flux of radio sources to use in pair finding (in same units as *Speak_col*)
min_size | minimum size of radio sources to use in pair finding (in same units as *maj_col*)
data_format | format of radio catalogue file


**Pair Data (DRAGNs.fits) Description**

the output data is a list of pairs where each **pair** appears only once - i.e. component_1 of a pair does not appear later on as component_2 of the pair. This data includes a number of derived parameters such as total source flux, size etc. Some of the precise column names are defined by the name of the column in the input radio catalogm and in such cases the column is listed below with the *config_file parameter* in the name. For instance if *Stot_col* is *Total_flux* then the column listed as *Stot_col_source* will actually be written to DRAGN.fits as *Total_flux_source*.

parameter | description [units]
----------|------------
pair_name | identifier for the pair of the format Jhhmmss.ss±ddmmss.s
RA_best | suggested central RA of the pair, if a core is found this is the core RA, otherwise this the flux-weighted RA of the pair [deg]
DE_best | suggested central Declination of the pair, if a core is found this is the core Dec, otherwise this the flux-weighted Dec of the pair [deg]
Sep_pair | the angular separation of the pair from the catalogue positions of the constituent components [arcsec]
Stot_col_source | Sum of all the component fluxes (both lobes and core), if no core is identified this is the same as *Stot_col_pair* [mJy]
E_Stot_col_source | quadrature sum of the uncertainties contributing to *Stot_col_source* [mJy]
CoreProm | ratio of core flux to *Stot_col_source* (if a core has been identified)
E_CoreProm | uncertainty in *CoreProm*
Tflux_ratio | ratio of total flux of name_col_1 to name_col_2 
E_Tflux_ratio | uncertainty in ratio of total flux of name_col_1 to name_col_2
Pflux_ratio | ratio of peak flux of name_col_1 to name_col_2
E_Pflux_ratio | uncertainty in ratio of peak flux of name_col_1 to name_col_2
Stot_col_pair | Sum of fluxes from the two lobes, if no core is identified this is the same as *Stot_col_source* [mJy]
E_Stot_col_pair | unceratainty in *Stot_col_pair* [mJy]
PA_pair | the position angle East of North from name_col_1 to name_col_2 in the pair [deg]
mean_misalignment | mean of the alignments of the two candidate lobes relative to the position angle of the pair [deg]
cenRA | RA of the flux weighted centroid of the pair  [deg]
cenDEC | DEC of the flux weighted centroid of the pair [deg]
medianRA | median (not flux weighted) RA of the pair [deg]
medianDEC | median (not flux weighted) DEC of the pair [deg]
ra_col_c | RA of core (if identified) [deg]
dec_col_c | | Dec of core (if identified) [deg]
sep_core_pcen | angular separation between the core position and the flux weighted centroid of the pair [arcsec]
d_2NN | distance to the second nearest neighbour of the pair [arcsec]
component_shares | product of the number of pairs identified that contain each component 
pref_pair | flag to identify the smallest separation pair if pair components are used in multiple pairs; 1=preferred, 2=not preferred 
candidate_core | the identifier of the candidate core component in the input catalog
name_col_1 | the identifier of the first lobe component from the pair in the input catalog
name_col_2 | the identifier of the second lobe component from the pair in the input catalog  
AllWISE | name of the closest AllWISE source to the position given by RA/DE_best if one found within search radius (only provided if *--findhosts=True*)
W1mag | magnitude (Vega) in W1 of AllWISE source (only provided if *--findhosts=True*) [mag] 
e_W1mag | uncertainty in *W1mag* of AllWISE source (only provided if *--findhosts=True*) [mag]
W2mag | magnitude (Vega) in W2 of AllWISE source (only provided if *--findhosts=True*) [mag]
e_W2mag | uncertainty in *W2mag* of AllWISE source (only provided if *--findhosts=True*) [mag]
W3mag | magnitude (Vega) in W3 of AllWISE source (only provided if *--findhosts=True*) [mag]
e_W3mag | uncertainty in *W3mag* of AllWISE source (only provided if *--findhosts=True*) [mag]
W3mag | magnitude (Vega) in W4 of AllWISE source (only provided if *--findhosts=True*) [mag]
e_W4mag | uncertainty in *W4mag* of AllWISE source (only provided if *--findhosts=True*) [mag]
sep_from_radio_source | angular separation between the flux weighted centroid of the radio pair and the AllWISE source (only provided if *--findhosts=True*) [arcsec]


**Core Data (candidate_cores.fits) Description**
A list of **all** the candidate cores found within the specified search radius or half the pair separation (whichever is the lesser) of the median position of the pair.

parameter | description [units]
----------|------------
pair_name | identifier for the pair of the format Jhhmmss.ss±ddmmss.s
candidate_core | the identifier of the candidate core component in the input catalog
sep_core_pmed | angular separation between the core position and the position of the pair [arcsec]
sep_core_pcen | angular separation between the core position and the flux weighted centroid of the pair [arcsec]


**Host Data (host_candidates.fits) Description**

A list of **all** AllWISE sources within the specified search radius of *RA/DE_best*. These data are only provided if *--findhosts=True* and are the default columns returned by querying https://cdsarc.cds.unistra.fr/viz-bin/cat/II/328

parameter | description [units]
----------|------------
pair_name | identifier for the pair of the format Jhhmmss.ss±ddmmss.s
sep_from_radio_source | angular separation between the flux weighted centroid of the radio pair and the AllWISE source [arcsec]
AllWISE | name of the closest AllWISE source to the position given by RA/DE_best if one found within search radius
RAJ2000 | RA of AllWISE source [deg]
DEJ2000 | Dec of AllWISE source [deg]
W1mag | magnitude (Vega) in W1 [mag] 
e_W1mag | uncertainty in *W1mag* [mag]
W2mag | magnitude (Vega) in W2 [mag]
e_W2mag | uncertainty in *W2mag* [mag]
W3mag | magnitude (Vega) in W3 [mag]
e_W3mag | uncertainty in *W3mag* [mag]
W3mag | magnitude (Vega) in W4 [mag]
e_W4mag | uncertainty in *W4mag* [mag]
Jmag | 2MASS J-band magnitude [mag]
e_Jmag | uncertainty in *Jmag* [mag]
Hmag | 2MASS H-band magnitude [mag]
e_Hmag | uncertainty in *Hmag* [mag] 
Kmag | 2MASS K-band magnitude [mag]
e_Kmag | uncertainty in *Kmag* [mag]
ccf | contamination and confucion flag
ex | extended source flag
var | variability flag
pmRA | proper motion in RA [mas/year]
e_pmRA | uncertainty in *pmRA* [mas/year]
pmDE | proper motion in Dec [mas/year]
e_pmDE | uncertainty in *pmDE* [mas/year]
qph | photometric quality flag
d2M | angular separation between AllWISE and 2MASS detections [arcsec]


**Code Dependencies**
The following python packages are required to run this code (version used in development):
* numpy (1.21.4)
* distutils (3.10.0)
* pandas (1.3.4)
* argparse (1.1)
* astropy (5.0)
* astroquery (0.4.5) *note: will work with astroquery version 0.4.2 but not 0.4.1*
