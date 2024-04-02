import numpy as np
from astropy.table import Table, join, unique
from astropy import units as u
from collections import Counter


################################################################################
################################################################################
###parameters
###make command line args

las_scale = 0.3
point_sep = 1.8*u.arcsec

hostnamecol = 'AllWISE'
sourcenamecol = 'Name'
sepcol = 'Sep_AllWISE'
sizecol = 'LAS'
relcol = 'Rel'
flagcol = 'Host_flag'
ztype_col = 'z_type'

add_z = True

sfile = '../../../data/Nov2022_rerun/supplementary_data/sources_preclean.fits'
dfile = '../../../data/Nov2022_rerun/supplementary_data/dragns_preclean.fits'
zfile = '../../../data/Nov2022_rerun/supplementary_data/host_redshifts.fits'

dcols = ['Name', 'RA', 'DEC', 'Flux', 'E_Flux', 'Core_prom',
         'E_Core_prom', 'Lobe_flux_ratio', 'E_Lobe_flux_ratio',
         'LAS', 'E_LAS', 'Misalign_1', 'E_Misalign_1', 'Misalign_2',
         'E_Misalign_2', 'Mean_misalign', 'E_Mean_misalign',
         'Lobe_1', 'Lobe_2', 'Core', 'RA_core', 'DEC_core',
         'RA_median', 'DEC_median', 'RA_fw', 'DEC_fw']
         
hcols = ['AllWISE', 'Rel', 'Sep_AllWISE', 'Host_flag', 'W1mag', 'e_W1mag',
         'W2mag', 'e_W2mag', 'W3mag', 'e_W3mag', 'W4mag', 'e_W4mag']

################################################################################
################################################################################
###functions

def tidy_wise_masks(data, mask, maskcols, maskval=0):
    'tidy up masks for wise mags -- mask set to on joining'
    for col in maskcols:
        colmask = mask | (data[col]==maskval)
        data[col].mask = colmask
        
    return


def finalise_cats(sources, dragns, add_z=False, zdata=None,
                  hostnamecol='AllWISE', sourcenamecol='Name',
                  sepcol='Sep_AllWISE', sizecol='LAS', relcol='Rel',
                  flagcol='Host_flag', acol='RA', ztype_col='z_type',
                  z_col='z', zsurvey_col='z_survey', ze_col='z_err',
                  las_scale=0.3, point_sep=1.8*u.arcsec,
                  output_multi_source_host_info=False,
                  coremaskcol='Core_prom',
                  dcols = ['Name', 'RA', 'DEC', 'Flux', 'E_Flux',
                           'Core_prom','E_Core_prom', 'Lobe_flux_ratio',
                           'E_Lobe_flux_ratio', 'LAS', 'E_LAS',
                           'Misalign_1', 'E_Misalign_1',
                           'Misalign_2','E_Misalign_2',
                           'Mean_misalign', 'E_Mean_misalign',
                           'Lobe_1', 'Lobe_2', 'Core', 'RA_core',
                           'DEC_core','RA_median', 'DEC_median',
                           'RA_fw', 'DEC_fw'],
                  hcols = ['AllWISE', 'Rel', 'Sep_AllWISE',
                           'Host_flag', 'W1mag', 'e_W1mag',
                           'W2mag', 'e_W2mag', 'W3mag',
                           'e_W3mag', 'W4mag', 'e_W4mag'],
                  magerr_str='e_', magerr_maskval=0):
    'tidy up catalogs so that no hosts assigned to multiple hosts are included and LAS/psf limits obeyed in final catalog'
    
    ##1) filter out hosts outside sep limits: sep_allwise > MAX(las_scale*LAS, point_sep) ###act on sources file and join after
    hmask_bad = ((sources[sepcol]>point_sep)
                 & (sources[sepcol]>0.3*sources[sizecol]))
    hmask = sources[hostnamecol].mask | hmask_bad | sources[hostnamecol].mask
    for col in hcols:
        sources[col].mask = hmask
    
    ##2) if a host has been assigned to multiple sources don't trust
    hosts = sources[~sources['AllWISE'].mask] ##subset only sources with hosts
    hcount = Counter(hosts['AllWISE']) ##count number of times host is used
    hcount = Table({'AllWISE': list(hcount.keys()),
                    'n_sources': list(hcount.values())})
    hosts = join(hosts, hcount, keys='AllWISE', join_type='left') ###adds number of sources a host is assigned to to the host table
    
    if output_multi_source_host_info == True:
        multi = hosts[hosts['n_sources']>1]
    
    ##3) fix
    hosts = hosts[hosts['n_sources']==1]
    hosts.remove_column('n_sources')
    hosts = hosts[[sourcenamecol]+hcols]
    ##add in redshifts if provided
    if add_z==True:
        if zdata is not None:
            ###remove duplicates -- sorting puts masked first so reverse after sorting by z_type will prioritise spec_z
            zdata.sort(ztype_col)
            zdata.reverse()
            zdata = unique(zdata, hostnamecol)
            hosts = join(hosts, zdata, keys=hostnamecol, join_type='left')
        else:
            print('')
            print('Warning: No redshift data provided, cannot add redshift information')
    
    ###rejoin hosts with dragns and sources for tidy catalogs
    scols_only = [col for col in sources.colnames if col not in hcols]
    sources = sources[scols_only]
    dragns = join(dragns, hosts[[sourcenamecol, hostnamecol,
                                 relcol, sepcol, flagcol]],
                  keys=sourcenamecol, join_type='left')
    sources = join(sources, hosts, keys=sourcenamecol, join_type='left')
    
    ###sort by RA and clean up masks
    dragns.sort(acol)
    coremask = dragns[coremaskcol].mask
    corecols = [col for col in dragns.colnames if 'core' in col or 'Core' in col]
    for col in corecols:
        dragns[col].mask = coremask
    
    sources.sort(acol)
    tidy_wise_masks(data=sources, mask=sources[hostnamecol].mask,
                    maskcols=[col for col in hcols if magerr_str in col],
                    maskval=magerr_maskval)
    if z_col in sources.colnames:
        zmask = sources[z_col].mask
        for col in [z_col, ze_col, ztype_col, zsurvey_col]:
            sources[col].mask = zmask
    
    outdata = {'sources': sources, 'dragns': dragns}
    
    if output_multi_source_host_info==True:
        outdata['multi_source_hosts'] = multi
    
    return outdata


def make_flat_dragnsplus(dragns, sources, snamecol='Name',
                         hnamecol='AllWISE', sepcol='Sep_AllWISE',
                         relcol='Rel', flagcol='Host_flag'):
    'add all host info to dragns for a flat table (saves end-user joining with sources)'
    hosts = sources[~sources[hnamecol].mask]
    dragns.remove_columns(names=[hnamecol, sepcol, relcol, flagcol])
    dragns = join(dragns, hosts, keys=snamecol, join_type='left')
    
    return dragns

################################################################################
################################################################################
###main


sources = Table.read(sfile)

dragns = Table.read(dfile)
dragns = dragns[dcols]
zdata = Table.read(zfile)
    
outdata = finalise_cats(sources=sources, dragns=dragns,
                        add_z=True, zdata=zdata,
                        output_multi_source_host_info=True)

sources = outdata['sources']
dragns = outdata['dragns']

dplus = make_flat_dragnsplus(dragns=dragns.copy(), sources=sources)
    


###make command line callable? just import functions
###make function for flat_dragn+ cat (all host info not just key stuff)




####subset multiple source hosts for investigation
#multi = hosts[hosts['n_sources']>1]
####fail cases for multiple source hosts:
## i) bright sources with sidelobes
## ii) rejected doubles: too compact/misaligned, one or both are unresolved?
## iii) complex morphologies not described by double model
## very small fraction mask out as host not reliable and sources dubious.
