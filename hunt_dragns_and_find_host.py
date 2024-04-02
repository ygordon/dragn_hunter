###wrapper script for dragn_hunting and host_finding
###directories dragnhunter and LR_host_finding correspond to dragnhunter_v2 and LR_host_finding_v2 in main code -- oct 2022


import sys, argparse, os
sys.path.append('dragnhunter')
sys.path.append('LR_host_finding')
from dragn_hunter import *
from find_hosts_lr import *
from get_redshifts import *
from fetch_z import *
from clean_catalog import *


######################################################
######################################################
###parameters
min_p = 0.5
wise_psf = 1.8*u.arcsec
double_poscal = 6 ##calibration errors in arcsec
min_pairsep = 6 ##minimum pair sep for DRAGNs
sepmis_a = -96.01351351 ###eq1 in Gordon+2023 fit
sepmis_b = 225.32323928 ###eq1 in Gordon+2023 fit
host_las_lim = 0.3

#scols = ['Name', 'RA', 'DEC', 'Flux', 'E_Flux', 'LAS', 'E_LAS', 'Type']
#hcols = ['Name', 'AllWISE', 'W1_Rel', 'W1_separation']
#wcols = ['AllWISE', 'W1mag', 'e_W1mag', 'W2mag', 'e_W2mag',
#         'W3mag', 'e_W3mag', 'W4mag', 'e_W4mag']


hdescribe = {'AllWISE': 'AllWISE identifier of matched host',
             'Rel': 'Probability that match is the correct host',
             'LR': 'Likelihood ratio of match',
             'Sep_AllWISE': 'Angular separation between radio source and AllWISE host',
             'Host_flag': 'Host ID quality flag (>0 is suspect)',
             'Type': 'Is the source a DRAGN (D) or single-component (S) source?',
             'W1mag': 'Vega magnitude in WISE W1 band',
             'e_W1mag': 'Uncertainty in W1mag',
             'W2mag': 'Vega magnitude in WISE W2 band',
             'e_W2mag': 'Uncertainty in W2mag',
             'W3mag': 'Vega magnitude in WISE W3 band',
             'e_W3mag': 'Uncertainty in W3mag',
             'W4mag': 'Vega magnitude in WISE W4 band',
             'e_W4mag': 'Uncertainty in W4mag'}


sourcecols_out = ['Name', 'RA', 'DEC', 'Flux', 'E_Flux', 'LAS', 'E_LAS', 'Type',
                  'Source_flag', 'AllWISE', 'RA_AllWISE', 'DE_AllWISE',
                  'Sep_AllWISE', 'LR', 'Rel', 'Host_flag',
                  'W1mag', 'e_W1mag', 'W2mag', 'e_W2mag', 'W3mag', 'e_W3mag',
                  'W4mag', 'e_W4mag', 'z', 'z_err', 'z_type', 'z_survey']

dragncols_out = ['Name', 'RA', 'DEC', 'Flux', 'E_Flux', 'Core_prom',
                 'E_Core_prom', 'Lobe_flux_ratio', 'E_Lobe_flux_ratio',
                 'LAS', 'E_LAS', 'Misalign_1', 'E_Misalign_1', 'Misalign_2',
                 'E_Misalign_2', 'Mean_misalign', 'E_Mean_misalign',
                 'Lobe_1', 'Lobe_2', 'Core', 'RA_core', 'DEC_core',
                 'RA_median', 'DEC_median', 'RA_fw', 'DEC_fw', 'Source_flag',
                 'AllWISE', 'Sep_AllWISE', 'LR', 'Rel', 'Host_flag']

sourccols_desc = {'Name': 'Julian name',
                  'RA': 'R.A. (J2000) of the radio source',
                  'DEC': 'Decl. (J2000) of the radio source',
                  'Flux': 'Total flux density of the radio source',
                  'E_Flux': 'Uncertainty in Flux',
                  'LAS': 'Largest angular extent',
                  'E_LAS': 'Uncertainty in LAS',
                  'Type': 'Is the source a DRAGN (D) or single-component (S) source?',
                  'Source_flag': 'Source quality flag (>0 is suspect)',
                  'AllWISE': 'AllWISE Host ID',
                  'RA_AllWISE': 'Host R.A.',
                  'DE_AllWISE': 'Host Decl.',
                  'Sep_AllWISE': 'Angular separation between radio source and AllWISE host',
                  'LR': 'Likelihood ratio of the host ID match',
                  'Rel': 'Probability the host ID is correct',
                  'Host_flag': 'Host ID flag (>0 is suspect)',
                  'W1mag': 'Vega magnitude in WISE W1 band',
                  'e_W1mag': 'Uncertainty in W1mag',
                  'W2mag': 'Vega magnitude in WISE W2 band',
                  'e_W2mag': 'Uncertainty in W2mag',
                  'W3mag': 'Vega magnitude in WISE W3 band',
                  'e_W3mag': 'Uncertainty in W3mag',
                  'W4mag': 'Vega magnitude in WISE W2 band',
                  'e_W4mag': 'Uncertainty in W4mag',
                  'z': 'Host redshift',
                  'z_err': 'Uncertainty in z',
                  'z_type': 'Redshift type (spec/photo)',
                  'z_survey': 'Survey redshift was obtained from'}

dragncols_desc = {'Name': 'Julian name',
                  'RA': 'Best R.A. (J2000) of the DRAGN',
                  'DEC': 'Best Decl. (J2000) of the DRAGN',
                  'Flux': 'Total flux density of the radio source',
                  'E_Flux': 'Uncertainty in Flux',
                  'Core_prom': 'Fraction of Flux from Core',
                  'E_Core_prom': 'Uncertainty in Core_prom',
                  'Lobe_flux_ratio': 'Flux ratio of Lobe_1/Lobe_2',
                  'E_Lobe_flux_ratio':  'Uncertainty in Lobe_flux_ratio',
                  'LAS': 'Largest angular extent',
                  'E_LAS': 'Uncertainty in LAS',
                  'Misalign_1': 'Relative misalignment of Lobe_1',
                  'E_Misalign_1': 'Uncertainty in Misalign_1',
                  'Misalign_2': 'Relative misalignment of Lobe_2',
                  'E_Misalign_2': 'Uncertainty in Misalign_2',
                  'Mean_misalign': 'Mean misalignment of lobes',
                  'E_Mean_misalign': 'Uncertainty in Mean_misalign',
                  'Lobe_1': 'Component name of Lobe_1',
                  'Lobe_2': 'Component name of Lobe_2',
                  'Core': 'Component name of Core',
                  'RA_core': 'R.A. (J2000) of Core',
                  'DEC_core': 'Decl. (J2000) of Core',
                  'RA_median': 'Median R.A. (J2000) of lobes',
                  'DEC_median': 'Median Decl. (J2000) of lobes',
                  'RA_fw': 'Flux weighted central R.A. (J2000) of lobes',
                  'DEC_fw': 'Flux weighted central Decl. (J2000) of lobes',
                  'Source_flag': 'Source quality flag (>0 is suspect)',
                  'AllWISE': 'AllWISE Host ID',
                  'Sep_AllWISE': 'Angular separation between radio source and AllWISE host',
                  'LR': 'Likelihood ratio of the host ID match',
                  'Rel': 'Probability the host ID is correct',
                  'Host_flag': 'Host ID flag (>0 is suspect)'}

######################################################
######################################################
###functions

def prepocess_cirada_ql(filename, qflag='Quality_flag',
                        dflag='Duplicate_flag',
                        sflag='S_Code',
                        keepcols=['Component_name','RA','DEC',
                                  'E_RA','E_DEC','Total_flux',
                                  'E_Total_flux','Peak_flux',
                                  'E_Peak_flux','Isl_rms','DC_Maj',
                                  'E_DC_Maj','DC_Min','E_DC_Min',
                                  'DC_PA','E_DC_PA'],
                        outfile='components.fits'):
    'filter VLASS QL catalog to good quality data and required columns'
    data = Table.read(filename)
    
    print('')
    print('Selecting high quality components from CIRADA VLASSQL component catalog')
    
    ##filter good data
    good_filter = (data[sflag]!='E') & (data[dflag]<2) & ((data[qflag]==0) | (data[qflag]==4))
    good = data[good_filter]
    
    ##filter to required columns
    good = good[keepcols]
    
    ##write to file
    print(f'writing selection to temporary file {outfile}')
    good.write(outfile)
    
    return outfile


def write_dragn_data(data_dict, outdir='.',
                     inc_hosts=True,
                     dragnkey='dragns',
                     singlekey='single-comp',
                     sourcekey='sources',
                     hostkey='hosts'):
    'write output data from dragnhunter to file'
    data_dict[dragnkey].write('/'.join([outdir, 'dragns.fits']))
    data_dict[singlekey].write('/'.join([outdir, 'single_comps.fits']))
    if inc_hosts==True:
        data_dict[sourcekey].write('/'.join([outdir, 'sources.fits']))
        data_dict[hostkey].write('/'.join([outdir, 'hosts.fits']))
        
    return


def best_match(host_candidates, rname='Name', hname='AllWISE',
               pcol='Rel', lrcol='LR', sepcol='Sep_AllWISE',
               minp=0.5, maxp=1.0):
    'select the single host/source match with the highest probability and subset required columns'
    
    host_candidates.sort(pcol)
    host_candidates.reverse()
    hosts = unique(host_candidates, rname, keep='first')
    
    ##subset and rename columns
    hosts = hosts[((hosts[pcol]>=minp) & (hosts[pcol]<=maxp))]
    hosts = hosts[rname, hname, sepcol, lrcol, pcol]
                         
    return hosts
    
    
def join_hosts(sdata, hdata, joinkey='Name', racol='RA',
               sizecol='LAS', psf=wise_psf,
               sepcol='Sep_AllWISE', pcol='Rel', hname='AllWISE',
               lrcol='LR',
               magcols=['W1mag', 'e_W1mag', 'W2mag', 'e_W2mag',
                        'W3mag', 'e_W3mag', 'W4mag', 'e_W4mag']
               ):
    'join hosts and sources'
    ###retain unit info from separation col
    sepunit = hdata[sepcol].unit
    
    ##join tables and sort by ra
    matches = join(sdata, hdata, keys=joinkey, join_type='left')
    matches.sort(racol)
    
    ###set up mask (mask out anything not matched, or with sep>size/2 when size>psf
    mmask = matches[pcol].mask ##main mask
    
    ###size/sep/lr mask needs filling with na and converting to array first
    matches[sepcol].fill_value=np.nan
    matches[sepcol]=matches[sepcol].filled()
    
    matches[lrcol].mask = mmask
    
    sizes = np.array(matches[sizecol])
    seps = np.array(matches[sepcol])
    smask = ((seps>0.3*sizes) & (seps>psf.value))
    combmask = mmask | smask

    masknacols = [hname, pcol, sepcol, lrcol]+magcols
    for col in masknacols:
        mcol = MaskedColumn(data=np.array(matches[col]), mask=combmask)
        matches[col] = mcol
        
    ###add unit info back in
    matches[sepcol].unit = sepunit
    for col in magcols:
        matches[col].unit = u.mag
    
    return matches
    
    
def add_col_descriptions(data, info_dict):
    'add in descriptions to columns from dictionary'
    datacols = data.colnames
    for col in list(info_dict.keys()):
        if col in datacols:
            data[col].description = info_dict[col]
        
    return


def mask_bad_hosts(data, host_id='AllWISE',
                   sepcol='Sep_AllWISE', sizecol='LAS',
                   psf=wise_psf, las_scale=host_las_lim,
                   host_cols=['AllWISE', 'LR', 'Rel', 'Sep_AllWISE',
                              'Host_flag', 'W1mag', 'W2mag', 'W3mag',
                              'W4mag', 'e_W1mag', 'e_W2mag', 'e_W3mag',
                              'e_W4mag']):
    'mask out columns for bad hosts'
    hmask = data[host_id].mask
    bmask = (data[sepcol]>psf) & (data[sepcol]>las_scale*data[sizecol])
    use_mask = hmask | bmask
    for col in host_cols:
        if col in data.colnames:
            if type(data[col])!=MaskedColumn:
                data[col] = MaskedColumn(data[col])
            data[col].mask = use_mask
    
    return data

    
def finalise_cats(dragns, sources, wise, hosts, min_p=0.5, las_scale=host_las_lim,
                  sourcecols=['Name', 'RA', 'DEC', 'Flux', 'E_Flux',
                              'LAS', 'E_LAS', 'Type'],
                  wisecols=['AllWISE', 'RAJ2000', 'DEJ2000',
                            'W1mag', 'e_W1mag', 'W2mag',
                            'e_W2mag', 'W3mag', 'e_W3mag',
                            'W4mag', 'e_W4mag'],
                  radio_id='Name', wise_id='AllWISE', psf=wise_psf):
    'finalise the output catalogs for publication'
    ###subset to required columns -- dragnhunte finds nearest host by default, but we're replacing that info with LR matches
    sources = sources[sourcecols]
    if type(wise['AllWISE'])==MaskedColumn:
        wise = wise[~wise['AllWISE'].mask]
    wise = unique(wise[wisecols], wise_id)
    
    ###select best host for each source with p>min_p and join with WISE info
    hosts = best_match(hosts, minp=min_p)
    hosts = join(hosts, wise, keys=wise_id, join_type='left')
    
    ###join with sources and dragns
    sources = join_hosts(sdata=sources, hdata=hosts, psf=psf)
    dragns = join_hosts(sdata=dragns, hdata=hosts, psf=psf)
    
    ###ensure all hosts satisfy selection criteria
    sources = mask_bad_hosts(data=sources, host_id=wise_id,
                             sepcol='Sep_AllWISE', sizecol='LAS',
                             psf=psf, las_scale=las_scale)
    dragns = mask_bad_hosts(data=dragns, host_id=wise_id,
                            sepcol='Sep_AllWISE', sizecol='LAS',
                            psf=psf, las_scale=las_scale)
    
    ###add description to new columns in matches
    add_col_descriptions(data=sources, info_dict=hdescribe)
    add_col_descriptions(data=dragns, info_dict=hdescribe)
    
    ###clean output data
    outdata = {'dragns': dragns, 'host_ids': sources}
    for key in list(outdata.keys()):
        outdata[key].meta = {}
    
    return outdata


def finalise_output(dragns, sources, outdir='output_files',
                    supdir='supplementary_data',
                    lr_bin='lr_bin',
                    dcols=None, scols=None,
                    dcoldesc=None, scoldesc=None):
    ###mask out bad cores
    if type(dragns['Core'])!=MaskedColumn:
        dragns['Core'] = MaskedColumn(dragns['Core'])
    cmask = dragns['Core'].mask | (dragns['Core']=='0.0')
    dragns['Core'].mask = cmask
    
    ###set column order if specified
    if dcols is not None:
        dcols = [i for i in dcols if i in dragns.colnames]
        dragns = dragns[dcols]
        if dcoldesc is not None:
            add_col_descriptions(data=dragns, info_dict=dcoldesc)
    if scols is not None:
        scols = [i for i in scols if i in sources.colnames]
        sources = sources[scols]
        if scoldesc is not None:
            add_col_descriptions(data=sources, info_dict=scoldesc)
    
    ###write main data to file
    dragns.write('/'.join([outdir, 'dragns.fits']), format='fits', overwrite=True)
    sources.write('/'.join([outdir, 'sources.fits']), format='fits', overwrite=True)
    
    mainfiles = ['dragns.fits', 'sources.fits']
    
    ###create directory to move supplementary data to
    dlist = os.listdir(outdir)
    if supdir not in dlist:
        os.mkdir('/'.join([outdir, supdir]))
    for file in dlist:
        if file not in mainfiles:
            src = '/'.join([outdir, file])
            dst = '/'.join([outdir, supdir, file])
            shutil.move(src=src, dst=dst)
    
    ###move lr bin
    shutil.move(src=lr_bin, dst='/'.join([outdir, supdir, lr_bin]))
    
    return


def subset_dragns_lr(dragns, lrdata, radio_id='Name', wise_id='AllWISE'):
    'subset only the LR info for DRAGNs'
    lrdragns = join(dragns[[radio_id]], lrdata, keys=radio_id, join_type='left')
    lrdragns = lrdragns[~lrdragns[wise_id].mask]
    
    return lrdragns

    
def update_core_hosts(dragns, lrdata, wise, psf=wise_psf,
                      pcol='Rel', sepcol='Sep_AllWISE',
                      rid='Name', wid='AllWISE',
                      flagcol='Host_flag', sizecol='LAS',
                      p_min=0.5):
    'update the host IDs where needed in the presence of a core identification'
    ###accounts for low-prob counterparts being coincident with core
    ###subset require lr info (exclude where sep > LAS/2)
    lrinfo = join(lrdata, dragns[[rid, sizecol]], keys=rid, join_type='left')
    lrinfo = lrinfo[(lrinfo[sepcol]<lrinfo[sizecol]/2)]
    lrinfo.remove_column(sizecol)
    
    ##segregate dragns with and without core -- acount for 0.0 instead of mask
    if type(dragns['Core'])!=MaskedColumn:
        dragns['Core'] = MaskedColumn(dragns['Core'])
    
    dmask1 = dragns['Core'].mask
    dmask2 = dragns['Core'] == '0.0'
    dragns['Core'].mask = dmask1 | dmask2
    
    core = dragns[~dragns['Core'].mask]
    nocore = dragns[dragns['Core'].mask]
    
    ###add in number of WISE candidates per radio source in LR
    wisecount = Counter(lrinfo[rid])
    nwise = Table({rid: list(wisecount.keys()),
                   'n_wise': list(wisecount.values())})
    
    ###get lr data for all dragns with cores
    core_lr = join(core[[rid]], lrinfo, keys=rid, join_type='left')
    core_nohost = core_lr[core_lr[wid].mask] ##selects cores with NO candidates
    core_host = core_lr[~core_lr[wid].mask] ##removes those with no candidates
    core_host.sort(pcol)
    core_host.reverse()
    core_host = unique(core_host, rid, keep='first') ###keeps best match
    
    ###good IDs (core coincident with LR identified host), don't need changing
    ### includes missing (low-p so not selected)
    ##Flag = 0
    good_ids = core_host[(core_host[sepcol]<=psf)]
    good_ids[flagcol] = -2
    
    ###find incorrect IDs (core coincident with a lower Rel host ID)
    ##Flag = -1
    bad_ids = core_host[(core_host[sepcol]>psf)]
    bad_ids = join(bad_ids, nwise, keys=rid, join_type='left')
    
    mmatch = bad_ids[(bad_ids['n_wise']>1)]
    smatch = bad_ids[(bad_ids['n_wise']==1)]
    altmatch = join(mmatch[[rid]], lrinfo, keys=rid, join_type='left')
    
    altmatch = altmatch[(altmatch[sepcol]<psf)]
    ##keep best LR if multiple (probably shouldn't be)
    
    if len(altmatch)>0:
        altmatch.sort(pcol)
        altmatch.reverse()
        altmatch = unique(altmatch, rid, keep='first')
        altmatch[flagcol] = -1
    
    ###find uncertain IDs (core offset from host ID but not coincident with another candidate), bad_ids not in altmatch
    ##Flag = 1
    umlist = [id for id in list(bad_ids[rid]) if id not in list(altmatch[rid])]
    if len(umlist)>0:
        umdata = join(Table({rid: umlist}), lrinfo, keys=rid, join_type='left')
        umdata.sort(pcol)
        umdata.reverse()
        umdata = unique(umdata, rid, keep='first')
    
        umdata[flagcol] = 1 ###mark all these as untrustworthy
    
        ###vstack results
        newhosts = vstack([good_ids, altmatch, umdata, core_nohost])
    else:
        newhosts = vstack([good_ids, altmatch, core_nohost])
    ###need to mask AllWISE and image_id
    nhmask = newhosts[pcol].mask
    for col in newhosts.colnames:
        if type(newhosts[col])==MaskedColumn and col!=rid:
            newhosts[col].mask = nhmask
    
    ###fill masked flags with 0
    newhosts[flagcol].fill_value = 0
    newhosts[flagcol] = newhosts[flagcol].filled()
    
    ####need to rejoin with WISE and DRAGN info
    dcols = [rid] + [c for c in dragns.colnames if c not in lrinfo.colnames and c not in wise.colnames]
    lrcols = [c for c in lrinfo.colnames if c in dragns.colnames] + [flagcol]
    wcols = [c for c in wise.colnames if c in dragns.colnames and c not in [rid, sepcol]]
    
    newmatches = newhosts[~newhosts[wid].mask][lrcols]
    newmatches = join(newmatches, wise[wcols], keys=wid, join_type='left')
    ncdata = join(core[dcols], newmatches, keys=rid, join_type='left')
    
    ###stack with nocore data
    nhdata = vstack([ncdata, nocore])
    
    ###sort by name
    nhdata.sort(rid)
    
    ###fill masked flags with 0
    nhdata[flagcol].fill_value = 0
    nhdata[flagcol] = nhdata[flagcol].filled()
              
    ###remask AllWISE and Core (converts to unmasked when stacked!)
    hostmask = nhdata[pcol].mask
    nhdata[wid].mask = hostmask
    try:
        coremask = nhdata['Core_prom'].mask
        nhdata['Core'].mask = coremask
    except:
        print('No masking of core')
    
    return nhdata


def update_source_table(sources, new_dhosts, typecol='Type', rid='Name',
                        wid='AllWISE', pcol='Rel', hflagcol='Host_flag',
                        hcols=['Name', 'AllWISE', 'Rel', 'Sep_AllWISE',
                               'Host_flag', 'W1mag', 'e_W1mag', 'W2mag',
                               'e_W2mag', 'W3mag', 'e_W3mag', 'W4mag',
                               'e_W4mag']):
    'update source table after DRAGN hosts updated to account for cores'
    ####subset single/doubles sources. vstack after doubles updated
    d_sources = sources[(sources[typecol]=='D')]
    s_sources = sources[(sources[typecol]!='D')]
    
    ###set up columns from source table for joining
    scols = ['Name'] + [col for col in sources.colnames if col not in hcols]
    
    d_sources = join(d_sources[scols], new_dhosts[hcols],
                     keys=rid, join_type='left')
    
    new_sources = vstack([d_sources, s_sources])
    new_sources.sort('RA')
    
    ###fill masked host_flags
    new_sources[hflagcol].fill_value = 0
    new_sources[hflagcol] = new_sources[hflagcol].filled()
    
    ###need to remask AllWISE (bytes really doesnt like masks with vstack)
    new_sources[wid].mask = new_sources[pcol].mask
    
    return new_sources


def update_host_ids_with_core_info(dragns, lrdata, wise, sources,
                                   wise_id='AllWISE', las_scale=0.3,
                                   psf=wise_psf, sizecol='LAS',
                                   sepcol='Sep_AllWISE'):
    'use information that a DRAGN has a core to update its host ID from the LR if appropriate, and update dragns/sources to reflect this'
    new_dragns = update_core_hosts(dragns=dragns, lrdata=lrdata, wise=wise)
    new_sources = update_source_table(sources=sources, new_dhosts=new_dragns)
    
    ###mask bad hosts here
    new_dragns = mask_bad_hosts(data=new_dragns, host_id=wise_id, sepcol=sepcol,
                                 sizecol=sizecol, psf=psf, las_scale=las_scale)
    new_sources = mask_bad_hosts(data=new_sources, host_id=wise_id, sepcol=sepcol,
                                 sizecol=sizecol, psf=psf, las_scale=las_scale)
    
    updated_hosts = {'dragns': new_dragns, 'sources': new_sources}
    
    return updated_hosts


def z_targets(sources, hosts, namecol='AllWISE',
              acol='RAJ2000', dcol='DEJ2000', posunits=('deg', 'deg')):
    'code get AllWISE positions of sources to search for redshift data for'
    
    ###ensure units in positions
    if hosts[acol].unit is None:
        hosts[acol].unit = posunits[0]
    if hosts[dcol].unit is None:
        hosts[dcol].unit = posunits[1]
    
    ##ensure no duplicates in host candidates
    hosts = unique(hosts, namecol)
    
    ###select only sources with hosts
    sources = sources[~sources[namecol].mask]
    
    ###join tables
    targets = join(sources[[namecol]], hosts[[namecol, acol, dcol]], keys=namecol, join_type='inner')
    
    return targets


def parse_args():
    'parse command line arguments'
    parser = argparse.ArgumentParser(description="find radio pairs in catalogue data")
    parser.add_argument("radio_cat",
                        help="radio component catalogue to find pairs in")
    parser.add_argument("--config", action='store',
                        type=str, default='dragnhunter/config.txt',
                        help="config file")
    parser.add_argument("--search_radius", action='store',
                        type=str, default='30arcsec',
                        help="search radius to use for core and host candidates")
    parser.add_argument("--find_hosts", action='store', type=str, default='True',
                        help="set code to query AllWISE for potential host candidates")
    parser.add_argument("--outdir", action='store', type=str,
                        default='output_files',
                        help="directory to write files to")
    parser.add_argument("--image_meta", action='store',
                        type=str, default='LR_host_finding/AllWISE-image_atlas_metadata-lite.fits',
                        help="file containing image meta data for maskims")
    parser.add_argument("--get_redshifts", action='store', type=str, default='True',
                        help="set code to query for redshifts")
    parser.add_argument("--z_radius", action='store', type=str, default='1arcsec',
                        help="search radius for redshifts")
    parser.add_argument("--host_name", action='store', type=str, default='AllWISE',
                        help="name column for hosts")
    parser.add_argument("--host_ra", action='store', type=str, default='RAJ2000',
                        help="RA column for hosts")
    parser.add_argument("--host_dec", action='store', type=str, default='DEJ2000',
                        help="DEC column for hosts")
    parser.add_argument("--cirada_preprocess", action='store', type=str,
                        default='True',
                        help='select good quality data from full CIRADA QL catalog')
    args = parser.parse_args()
    
    ##make args.writepairs a bool and search_radius a quantity
    args.search_radius = u.Quantity(args.search_radius)
    args.z_radius = u.Quantity(args.z_radius)
    args.find_hosts = strtobool(args.find_hosts)
    args.get_redshifts = strtobool(args.get_redshifts)
    args.cirada_preprocess = strtobool(args.cirada_preprocess)
    
    return args



######################################################
######################################################
###main

if __name__ == '__main__':
    t0 = time.time()
    args = parse_args()
    if args.cirada_preprocess == True:
        temp_compfile = prepocess_cirada_ql(filename=args.radio_cat)
        args.radio_cat = temp_compfile
    
    cparams, cdict = config_parse(args.config)
    pairs = find_doubles(catfile=args.radio_cat, config_file=args.config,
                         search_rad=args.search_radius,
                         write_file=True, return_pairs=True,
                         outdir=args.outdir)
    data = select_sources_and_find_hosts(pairs=pairs,
                                         components=Table.read(args.radio_cat),
                                         minflux=cparams['flux_min'],
                                         cname=cdict['name'],
                                         cpeak=cdict['peak'],
                                         cflux=cdict['total'],
                                         cfluxerr=cdict['etot'],
                                         c_ra=cdict['ra'],
                                         c_dec=cdict['dec'],
                                         cpa=cdict['pa'],
                                         cmaj=cdict['maj'],
                                         emaj=cdict['emaj'],
                                         search_rad=args.search_radius,
                                         find_hosts=args.find_hosts,
                                         minsep=min_pairsep,
                                         sepmis_a=sepmis_a, sepmis_b=sepmis_b)
    ###account for missing hosts
    if args.find_hosts==True:
        if type(data['hosts']['AllWISE'])==MaskedColumn:
            data['hosts'] = data['hosts'][~data['hosts']['AllWISE'].mask]
    ##write data to file
    write_dragn_data(data_dict=data, outdir=args.outdir,
                     inc_hosts=args.find_hosts)
                     
    ###run LR host finding
    if args.find_hosts==True:
        image_info = Table.read(args.image_meta)
        host_matching(source_cat=data['sources'], host_cat=data['hosts'],
                      image_meta=image_info,
                      sra='RA', sdec='DEC', sid='Name', sflux='Flux',
                      sfluxerr='E_Flux', ssize='LAS',
                      iid='coadd_id', ira='ra', idec='dec',
                      hid='AllWISE', hra='RAJ2000', hdec='DEJ2000',
                      hmags=['W1mag'], hmagerrs=['e_W1mag'],
                      sepcol='Sep_AllWISE',
                      bin_dir='lr_bin', outdir=args.outdir,
                      assumed_psf=wise_psf, verbose=True,
                      radio_beam_size=3, cal_errors=double_poscal)
    
        ###need to add host info to sources for final table (and select core where available)
        lr_results=Table.read('/'.join([args.outdir, 'all_LR_matches.fits']))
        main_output = finalise_cats(dragns=data['dragns'], sources=data['sources'],
                                    wise=data['hosts'], hosts=lr_results, min_p=0.5,
                                    sourcecols=['Name', 'RA', 'DEC', 'Flux', 'E_Flux',
                                                'LAS', 'E_LAS', 'Type', 'Source_flag'],
                                    wisecols=['AllWISE', 'RAJ2000', 'DEJ2000',
                                              'W1mag', 'e_W1mag', 'W2mag',
                                              'e_W2mag', 'W3mag', 'e_W3mag', 'W4mag',
                                              'e_W4mag'],
                                    radio_id='Name', wise_id='AllWISE',
                                    psf=wise_psf)
    
        ###update hosts here
        lrdragns = subset_dragns_lr(dragns=main_output['dragns'], lrdata=lr_results)
        new_hosts = update_host_ids_with_core_info(dragns=main_output['dragns'],
                                                   lrdata=lrdragns,
                                                   wise=unique(data['hosts'], 'AllWISE'),
                                                   sources=main_output['host_ids'],
                                                   wise_id=args.host_name,
                                                   las_scale=host_las_lim,
                                                   psf=wise_psf, sizecol='LAS',
                                   sepcol='Sep_AllWISE')
        main_output['dragns'] = new_hosts['dragns']
        main_output['host_ids'] = new_hosts['sources']
        ##rename host pos cols and ensure masked if AllWISE masked
        print('')
        print('updating host position column names')
        main_output['host_ids'].rename_columns(names=['RAJ2000', 'DEJ2000'],
                                               new_names=['RA_AllWISE', 'DE_AllWISE'])
        if type(main_output['host_ids']['AllWISE'])==MaskedColumn:
            if type(main_output['host_ids']['RA_AllWISE'])!=MaskedColumn:
                main_output['host_ids']['RA_AllWISE'] = MaskedColumn(main_output['host_ids']['RA_AllWISE'])
            main_output['host_ids']['RA_AllWISE'].mask = main_output['host_ids']['AllWISE'].mask
            print('Host RA mask updated')
            if type(main_output['host_ids']['DE_AllWISE'])!=MaskedColumn:
                main_output['host_ids']['DE_AllWISE'] = MaskedColumn(main_output['host_ids']['DE_AllWISE'])
            main_output['host_ids']['DE_AllWISE'].mask = main_output['host_ids']['AllWISE'].mask
            print('Host Decl mask updated')
        else:
            print('mask not updated!')
        ###overwrite sources/dragns in output_files and move rest to new folder supplementary data
        finalise_output(dragns=main_output['dragns'],
                        sources=main_output['host_ids'],
                        outdir=args.outdir,
                        supdir='supplementary_data',
                        dcols=dragncols_out, scols=sourcecols_out,
                        dcoldesc=dragncols_desc, scoldesc=sourccols_desc)
    ###get redshifts?
    if args.get_redshifts==True and args.find_hosts==True:
        fetch_redshifts(sources=main_output['host_ids'], hosts=data['hosts'],
                        sourcecols_out=sourcecols_out,
                        sourccols_desc=sourccols_desc,
                        z_radius=args.z_radius,
                        host_name=args.host_name,
                        acol=args.host_ra,
                        dcol=args.host_dec,
                        outname='/'.join([args.outdir, 'sources.fits']))
    
    ###cleanup catalog
    cleanup(dragns=Table.read('/'.join([args.outdir, 'dragns.fits'])),
            sources=Table.read('/'.join([args.outdir, 'sources.fits'])),
            hosts=Table.read('/'.join([args.outdir, 'supplementary_data/hosts.fits'])),
            backup_data=False, hname='AllWISE', hmbackup='Sep_AllWISE',
            dfilename='/'.join([args.outdir, 'dragns.fits']),
            sfilename='/'.join([args.outdir, 'sources.fits']))
    
    ###print time taken -- for testing purposes
    if args.cirada_preprocess == True:
        if args.radio_cat == temp_compfile:
            print('')
            print(f'removing temporary file {args.radio_cat}')
            os.remove(args.radio_cat)
    t_elapsed = np.round(time.time()-t0, 2)
    print('')
    print('Time elapsed = ' + str(t_elapsed) + 's')
    

