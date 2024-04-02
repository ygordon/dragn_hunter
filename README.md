# Code to find double sources in VLASS

This readme gives a brief, high level, overview of the CIRADA code to used to find DRAGNs and single-component radio sources in VLASS plus their AllWISE host IDs. If you find the code or its resultant data products helpful in your research we request that you cite the article _"A Quick Look at the 3GHz Radio Sky. II. Hunting for DRAGNs in the VLA Sky Survey"_ (Y. A. Gordon et al., 2023, submitted to ApJS, https://arxiv.org/abs/2303.12830), which describes the method and output data in detail. The code acts as a wrapper for two main algorithms:

1. Yjan Gordon's DRAGNhunter (https://github.com/ygordon/dragn_hunter), a code for finding likely double sources in component catalog data. This algorithm was designed using VLASS data and is also set up to find potential host candidates from AllWISE.
2. Leah Morabito's likelihood ratio host-identification code (https://github.com/lmorabit/likelihood_ratio).


In short this wrapper runs DRAGNhunter to identify likely doubles in VLASS and potential host candidates from AllWISE. The likelihood ratio code is then used to identify the most likely correct host for these doubles, as well as hosts for single-component sources. Where a double source has a radio core identified, in a few percent of cases the core and host identifications disagree. In such cases the core position is used to select the host over the lieklihood ratio identification. Run as 

    >python hunt_dragns_and_find_hosts.py radio_catalogue

where *radio_catalogue* is the filename of the radio catalogue data file.

By default this code will attempt to find redshifts for hosts it identifies. 
In some instances this step can fail, e.g. as a result of the CDS server bing down or a broken connection, etc.
In the event of this happening obtaining the redshifts can be reattempted without repeating the slow step of host finding by calling:

    >python fetch_z.py

There are a number of optional arguments to this _fetch_z.py_ which can be invoked (see help) but the default setup should be okay for most cases where this needs to be run. If still struggling with timeout issues, then the _--chunk_size_ argument should be used to query smaller data chunks at any one time.


A test input catalog has been provided in the directory _test_data_. This radio component data is a subset of the VLASS component catalog (https://ui.adsabs.harvard.edu/abs/2021ApJS..255...30G/abstract) covering the region $120^{\circ} < \alpha < 140^{\circ}$ and $-5^{\circ} < \delta +5^{\circ}$.



## Output Data

The output of the code consists of a table of double sources (_dragns.fits_) with properties and host IDs, and a table of all sources (single components and doubles) with host identifications and redshifts where available. These are output to the directory _output_files_. 
A number of suplementary provide more extensive meta data and are output to in the _output_files/supplementary_data_ folder.


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
Mean_misalign | Mean value of _Misalign_1_ and _Misalign_2_ [deg]
E_Mean_misalign | Uncertainty in _Mean_misalign_ [deg]
Lobe_1 | Name of the component identified as the first lobe (lobe order is arbitrary)
Lobe_2 | Name of the component identified as the second lobe (lobe order is arbitrary)
Core | Name of the component identified as the radio core (if one was identified)
RA_core | R.A. of _Core_ (if identified) [deg]
DEC_core | Decl. of _Core_ (if identified) [deg]
RA_median | Median R.A. of two lobes [deg]
DEC_median | Median Decl. of two lobes [deg]
RA_fw | Flux weighted central R.A. of two lobes [deg]
DEC_fw | Flux weighted central Decl. of two lobes [deg]
Source_flag | Source quality flag (>0 is suspect)
AllWISE | Name of the AllWISE host ID
Sep_AllWISE | Angular separation between radio source and AllWISE host ID [arcsec]
LR | Likelihood ratio of host ID
Rel | Probabilty that the host is correct from the likelihood ratio
Host_flag | Host ID flag (>0 is suspect)


<br/>

### Source table (sources.fits)
A table of **all** sources (DRAGNs and single-component) for which host-finding is done.

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
Source_flag | Source quality flag (>0 is suspect)
AllWISE | Name of the AllWISE host ID
RA_AllWISE | R.A. of the AllWISE host [deg]
DE_AllWISE | Decl. of the AllWISE host [deg]
Sep_AllWISE | Angular separation between radio source and AllWISE host ID [arcsec]
LR | Likelihood ratio of host ID
Rel | Probabilty that the host is correct from the likelihood ratio
Host_flag | Host ID flag (>0 is suspect) 
W1mag | Vega magnitude of AllWISE host in the W1 band [mag]
e_W1mag | Uncertainty in _W1_mag_ [mag]
W2mag | Vega magnitude of AllWISE host in the W2 band [mag]
e_W2mag | Uncertainty in _W2_mag_ [mag]
W3mag | Vega magnitude of AllWISE host in the W3 band [mag]
e_W3mag | Uncertainty in _W3_mag_ [mag]
W4mag | Vega magnitude of AllWISE host in the W4 band [mag]
e_W4mag | Uncertainty in _W4_mag_ [mag]
z | Host redshift
z_err | Uncertainty in _z_
z_type | Redshift type (spec/photo)
z_survey | Survey redshift was obtained from

<br/>

## Code Dependencies

This code was developed using Python 3.10 and the following packages are required to run (version used in development):
* applpy (2.1.0)
* argparse (1.1)
* astropy (5.0.2)
* astroquery (0.4.6)
* distutils (3.10.0) 
* matplotlib (3.5.1)
* numpy (1.22.3)
* pandas (1.4.1)
* scipy (1.8.0)

