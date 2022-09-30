# dragn_hunter.py readme

command line code for searching for pairs of radio components that satisfy specified criteria, and flagging likely Double Radio sources associated with Active Galactic Nuclei (DRAGNs). Additionally candidate hosts are identified by querying the CatWISE2020 catalog using the flux weighted centroid of the double as the reference position. Run as:

    >python dragn_hunter.py radio_catalogue

where *radio_catalogue* is the filename of the radio catalogue data file. 

Additionally, there are a number of optional commands, these are [default values in square brackets]:\
*--config* [config.txt], config file describing the input component table (see below);\
*--find_hosts* [True], tell the code whether to query the AllWISE source catalog for potential host candidates;\
*--search_radius* [30arcsec], maximum radius to use when searching for radio cores (and AllWISE sources if *--find_hosts=True*);\
*--outdir* [output_files], directory in which to dump output files (will make directory if doesn't already exist).


## Output Data

The output of the code consists of a primary table (_dragns.fits_) and a number of suplementary tables described below. Note, the tables _sources.fits_ and _host_candidates.fits_ will only be produced if _--find_hosts_==True.

<br/>

### Table of DRAGNs (dragns.fits)
This is a list of DRAGNs selected from the parent sample of pairs of candidate lobes in _pairs.fits_. This is the primary output of _dragnhunter.py_.

parameter | description [units]
----------|------------
Name | Julian identifier (Jhhmmss.ss±ddmmss.s)
RA | suggested central R.A. of the DRAGN, if a core is found this is the core RA, otherwise this is the flux-weighted RA of the pair [deg]
DEC | suggested central Declination of the DRAGN, if a core is found this is the core Dec, otherwise this is the flux-weighted Dec of the pair [deg]
Flux | Total flux of the DRAGN, the sum of the fluxes of all associated components [mJy]
E_Flux | Uncertainty in _Flux_ [mJy]
Core_prom | Fraction of _Flux_ associated with the _Core_ (when identified)
E_Core_prom | Uncertainty in Core_prom
Lobe_flux_ratio | Ratio of the flux from _Lobe_1_ to the flux from _Lobe_2_
E_Lobe_flux_ratio | Uncertainty in _Lobe_flux_ratio_ 
LAS | Estimate of the Largest Angular Extent of the DRAGN [arcsec]
E_LAS | Uncertainty in LAS [arcsec]
Misalign_1 | Relative misalignment of _Lobe_1_ relative to the alignment of the DRAGN [deg]
E_Misalign_1 | Uncertainty in _Misalign_1_ [deg]
Misalign_2 | Relative misalignment of _Lobe_2_ relative to the alignment of the DRAGN [deg]
E_Misalign_2 | Uncertainty in _Misalign_2_ [deg]
Lobe_1 | Name of the component identified as the first lobe (lobe order is arbitrary)
Lobe_2 | Name of the component identified as the second lobe (lobe order is arbitrary)
Core | Name of the component identified as the radio core (if one was identified)
RA_core | R.A. of _Core_ (if identified) [deg]
DEC_core | Decl. of _Core_ (if identified) [deg]
RA_median | Median R.A. of two lobes [deg]
DEC_median | Median Decl. of two lobes [deg]
RA_fw | Flux weighted central R.A. of two lobes [deg]
DEC_fw | Flux weighted central Decl. of two lobes [deg]

<br/>

### Table of all identified pairs (pairs.fits)
This table describes the parent sample of 'pairs of candidate lobes' from which _dragns.fits_ is drawn. The table contains many of the columns from the input radio component catalog for each component identified as a canidate lobe or core (colname_1, colname_2, colname_c respectively). Additionaly some derived properties of the pairs are output with these column definitions below (some overlap with _dragns.fits_). These may be useful for refining _dragn_hunter.py_ for use on different data sets.

parameter | description [units]
----------|------------
pair_name | Julian identifier (Jhhmmss.ss±ddmmss.s)
Sep_pair | Angular separation between the pair of candidate lobes [arcsec]
E_Sep_pair | Uncertainty on _Sep_pair_ [arcsec]
PA_pair | Position angle (East of North) of the pair [deg]
E_PA_pair | Uncertainty on PA_pair [deg]
dPA_1 | difference between PA of candidate lobe 1 and _PA_pair_ [deg]
abs_dPA_1 | modulus of _dPA_1_ [deg]
E_dPA_1 | Uncertainty in _dPA_1_ [deg]
dPA_2 | difference between PA of candidate lobe 1 and _PA_pair_ [deg]
abs_dPA_2 | modulus of _dPA_1_ [deg]
E_dPA_2 | Uncertainty in _dPA_1_ [deg]
mean_misalignment | (_abs_dPA_1_ + _abs_dPA_2_)/ 2 [deg]
E_mean_misalignment | Uncertainty in _mean_misalignment_ [deg]
Tflux_ratio | flux_1(total)/flux_2(total)
E_Tflux_ratio | Uncertainty in _Tflux_ratio_
Pflux_ratio | flux_1(peak)/flux_2(peak)
E_Pflux_ratio | Uncertainty in _Pflux_ratio_
Flux_pair | flux_1(total)+flux_2(total); "_Flux_" will be replaced by the name of the total flux column in the input radio catalog [mJy]
E_Flux_pair | Uncertainty in _Flux_pair_; "_Flux_" will be replaced by the name of the total flux column in the input radio catalog [mJy]
medianRA | Median position of the two candidate lobes in R.A. [deg]
medianDEC | Median position of the two candidate lobes in Decl. [deg]
cenRA | Flux-weighted central position of the two candidate lobes in R.A. [deg]
cenDEC | Flux-weighted central position of the two candidate lobes in Decl. [deg]
d_2NN | Anglar distance to the second nearest candidate lobe [arcsec]
component_shares | Product of the number of pairs identified that contain each component
pref_pair | Flag to identify the smallest separation pair if pair components are used in multiple pairs; 1=preferred, 2=not preferred 
sep_core_pcen | Angular separation between the core position and the flux weighted centroid of the pair (if a core is identified) [arcsec]
Flux_source | Flux_c + _Flux_pair_; "_Flux_" will be replaced by the name of the total flux column in the input radio catalog [mJy]
E_Flux_source | Uncertainty in _Flux_source_; "_Flux_" will be replaced by the name of the total flux column in the input radio catalog [mJy]
CoreProm | Fraction of _Flux_source_ associated with the core (when identified)
E_CoreProm | Uncertainty in Core_prom
RA_best | 'Best' R.A. of pair. If a core is found this is the core RA, otherwise this is _cenRA_ [deg]
DE_best | 'Best' Decl. of pair. If a core is found this is the core Decl., otherwise this is _cenDEC_ [deg]

<br/>

### Core candidates (candidate_cores.fits)
A Table of all components identified as candidate core for the pairs listed in _pairs.fits_.

parameter | description [units]
----------|------------
pair_name | Name of the radio pair
core_name | Name of the candidate core
RA | R.A. of the source [deg]
sep_core_pmed | Angular separation between the candidate core and the median position of the pair [arcsec]
sep_core_pcen | Angular separation between the candidate core and the flux-weighted central position of the pair [arcsec]

<br/>

### Single component sources (single_comps.fits)
This is a table of all components (above the minimum flux limit set in the config file) not included as part of a DRAGN (lobe or core) in _dragns.fits_. These are treates as potential single-component sources for host finding purposes. Column names are the same as in the input radio component catalog file.

<br/>

### Source table (sources.fits)
A table of **all** sources (DRAGNs and single-component) for which host-finding is done. ONLY PRODUCED WHEN HOST_FINDING == True.

parameter | description [units]
----------|------------
Name | Julian identifier (Jhhmmss.ss±ddmmss.s)
RA | R.A. of the source [deg]
DEC | Decl. of the source [deg]
Flux | Total flux of the source [mJy]
E_Flux | Uncertainty in _Flux_ [mJy]
LAS | Estimate of the Largest Angular Extent of the source (deconvolved major axis for single-components) [arcsec]
E_LAS | Uncertainty in LAS [arcsec]
Type | Identifier on source type: 'S' = single-component, 'D' = DRAGN
AllWISE | Name of the closest AllWISE source
W1_mag | Vega magnitude of closest AllWISE source in the W1 band [mag]
e_W1_mag | Uncertainty in _W1_mag_ [mag]
W2_mag | Vega magnitude of closest AllWISE source in the W2 band [mag]
e_W2_mag | Uncertainty in _W2_mag_ [mag]
W3_mag | Vega magnitude of closest AllWISE source in the W3 band [mag]
e_W3_mag | Uncertainty in _W3_mag_ [mag]
W4_mag | Vega magnitude of closest AllWISE source in the W4 band [mag]
e_W4_mag | Uncertainty in _W4_mag_ [mag]
Sep_AllWISE | Angular separation between radio source and closest AllWISE source [arcsec]

<br/>

### All host candidates (host_candidates.fits)
A table of all AllWISE host candidates within the defined _search_radius_ of the radio source. ONLY PRODUCED WHEN HOST_FINDING == True. Columns are provided by the CDS VizieR version of the AllWISE catalog (https://cdsarc.cds.unistra.fr/viz-bin/cat/II/328) with the addition of:

parameter | description [units]
----------|------------
Name | Julian name of the radio source.
Sep_AllWISE | Angular separation between radio source and named AllWISE source [arcsec]

<br/>

## Configuration File
The config file should be a text file with two tab separated columns of data (an example, _config.txt_ is provided). The left-hand columns define the input data structure as listed below:

parameter | description
----------|------------
name_col | name of the object name/id column in the radio catalogue file
ra_col | name of the RA column in the radio catalogue file (assumed to be in decimal degrees)
dec_col | name of the Dec column in the radio catalogue file (assumed to be in decimal degrees)
ra_err_col | name of the RA uncertainty column in the radio catalogue file (assumed to be in decimal degrees)
dec_err_col | name of the Dec uncertainty column in the radio catalogue file (assumed to be in decimal degrees)
Speak_col | name of the peak flux column in the radio catalogue file
Speak_err_col | name of the peak flux error column in the radio catalogue file
Stot_col | name of the total flux column in the radio catalogue file
Stot_err_col | name of the total flux error column in the radio catalogue file
maj_col | name of the major axis size column in the radio catalogue file
min_col | name of the minor axis size column in the radio catalogue file
pa_col | name of the position angle column in the radio catalogue file
maj_err_col | name of the major axis size error column in the radio catalogue file 
min_err_col | name of the minor axis size error column in the radio catalogue file
pa_err_col | name of the position angle error column in the radio catalogue file
min_brightness | minimum peak flux of radio sources to use in pair finding (in same units as *Speak_col*)
min_size_lobe | minimum size of radio sources to use as candidate lobes in pair finding (in same units as *maj_col*)
data_format | format of radio catalogue file

<br/>

## Code Dependencies

The following python packages are required to run this code (version used in development):
* numpy (1.22.3)
* distutils (3.10.0)
* pandas (1.4.1)
* argparse (1.1)
* astropy (5.0.2)
* astroquery (0.4.6) *note: will work with astroquery version 0.4.2 but not 0.4.1*

