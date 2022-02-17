###code for finding Double Radio AGN (DRAGNs) from radio catalogue data
###finds pairs of radio components based on flux and size criteria
###developed by Yjan Gordon (yjan.gordon@umanitba.ca) Jan 2021

import numpy as np, pandas as pd, argparse
from distutils.util import strtobool
from astropy.table import Table, hstack, vstack, unique, join, Column, MaskedColumn
from astropy.coordinates import SkyCoord, match_coordinates_sky
from astropy import units as u
from astroquery.vizier import Vizier
from collections import Counter

####use astroquery v 0.4.3 (0.4.1 breaks on multi object queries)

########################################
#time of code check - not important for code to run
import time
start_time = time.time()
########################################


##############################################################################
###functions

def source_name(ra, dec, aprec=2, dp=5, prefix=''):
    ###create source name to nearest arcsec based on ra/dec
    ###truncated to dp so compatible with IAU standards
    ra, dec = np.array(ra), np.array(dec)
    
    cat = SkyCoord(ra=ra, dec=dec, unit='deg')
    
    astring = cat.ra.to_string(sep='', unit='hour', precision=dp, pad=True)
    dstring = cat.dec.to_string(sep='', precision=dp, pad=True, alwayssign=True)
    
    ###truncation index
    tinda = aprec+7
    tindd = tinda
    if aprec==1:
        tindd = tindd-1
    
    if prefix == '':
        sname = ['J' + astring[i][:tinda] + dstring[i][:tindd] for i in range(len(astring))]
    else:
        sname = [prefix + ' J' + astring[i][:tinda] + dstring[i][:tindd] for i in
                 range(len(astring))]

    return(sname)


def round_cols(data, round_dict):
    'round columns in table to specified number of dp'
    ### data should be an astropy table
    ### round_dict should be of format {column_name: number of dp to round to}
    for col in data.colnames:
        if col in list(round_dict.keys()):
            colunit = data[col].unit
            data[col] = np.round(data[col], round_dict[col])
            data[col].unit = colunit
    
    return data
    
    
def fill_masks(data, coldict):
    'go through dict and fill mask values'
    for col in list(coldict.keys()):
        if col in data.colnames and type(data[col])==MaskedColumn:
            data[col].fill_value = coldict[col]
            data[col] = data[col].filled()
    
    return data
    

def identify_pair_shares(data, cid1='Component_name_1', cid2='Component_name_2',
                         pair_id='pair_name', sepcol='Sep_pair'):
    'identify cases where a component has been used in multiple pairs'
    ###MOVE THIS TO DRAGNhunter
    ##prioritise closest pair, but still flag
    
    ### 1) count usage of all component names
    compcount = Counter(np.array(list(data[cid1]) + list(data[cid2])))
    ccount = Table({'compid': list(compcount.keys()),
                    'n_pairs': list(compcount.values())})
    
    ### 2) qualify complexity -- setup another table to join with
    cq = data[[pair_id, cid1, cid2]]
    cq = join(cq, ccount, keys_left=cid1, keys_right='compid')
    cq = join(cq, ccount, keys_left=cid2, keys_right='compid')
    
    ##make column that is product of two n_pairs columns
    cq['component_shares'] = cq['n_pairs_1'] * cq['n_pairs_2']
    
    ###rejoin with data
    outdata = join(data, cq[[pair_id, 'component_shares']], keys=pair_id,
                   join_type='left')

    return outdata


def flag_closest_pair(data, cid1='Component_name_1', cid2='Component_name_2',
                      pair_id='pair_name', sort_on='Sep_pair', tokeep='first',
                      good_flag=1, bad_flag=2):
    'create flag for "closest" pair when component_shares>1'
    ###again move to dragnhunter.py
    c1data = data[[cid1, pair_id, sort_on]]
    c2data = data[[cid2, pair_id, sort_on]]
    
    ###rename columns for vstack
    c1data.rename_column(name=cid1, new_name='compid')
    c2data.rename_column(name=cid2, new_name='compid')
    
    comppair = vstack([c1data, c2data])
    ###sort and subset for rejoin, pairs to remove duplication then components
    comppair.sort(keys=sort_on)
    comppair = unique(comppair, keys=pair_id, keep=tokeep)
    
    comppair.sort(keys=sort_on)
    comppair = unique(comppair, keys='compid', keep=tokeep)
    
    ###add flag column
    comppair['pref_pair'] = good_flag
    
    ###rejoin with main data
    data = join(data, comppair[[pair_id, 'pref_pair']], keys=pair_id,
                join_type='left')
                
    ###fill masked values that are the not good_flags
    data['pref_pair'].fill_value = bad_flag
    data['pref_pair'] = data['pref_pair'].filled()
    
    return data
    

def find_cores(pairs, cand_cores, acol_p='cenRA', acol_med='medianRA',
               dcol_med='medianDEC', dcol_p='cenDEC', acol_c='RA',
               dcol_c='DEC', sizecol='Sep_pair', sepcolname='sep_from_cen',
               pair_id='pair_name', core_id='Component_name',
               corename_out='candidate_core', search_rad=30*u.arcsec):
    'find all candidate cores within search_radius and subset closest'
    
    ###setup poscats
    pcat = SkyCoord(ra=pairs[acol_med], dec=pairs[dcol_med])
    ccat = SkyCoord(ra=cand_cores[acol_c], dec=cand_cores[dcol_c])
    
    ###cross match
    xmatch = pcat.search_around_sky(ccat, seplimit=search_rad)
    
    ###create table of useful results
    xmtab = Table({pair_id: pairs[pair_id][xmatch[1]],
                   'pair_ra': pairs[acol_p][xmatch[1]],
                   'pair_dec': pairs[dcol_p][xmatch[1]],
                   core_id: cand_cores[core_id][xmatch[0]],
                   'core_ra': cand_cores[acol_c][xmatch[0]],
                   'core_dec': cand_cores[dcol_c][xmatch[0]],
                   'sep_core_pmed': xmatch[2].to(pairs[sizecol].unit)})
    
    ###add in sepcol to keep
    xmpos_pair = SkyCoord(ra=xmtab['pair_ra'], dec=xmtab['pair_dec'])
    xmpos_core = SkyCoord(ra=xmtab['core_ra'], dec=xmtab['core_dec'])
    xmtab[sepcolname] = xmpos_core.separation(xmpos_pair).to(search_rad.unit)
    
    ###remove any where sep_from_med > size/2
    xmtab = join(xmtab, pairs[[pair_id, sizecol]], keys=pair_id, join_type='left')
    xmtab = xmtab[(xmtab['sep_core_pmed']<xmtab[sizecol]/2)]
    
    ###remove extra columns from xmtab
    xmtab.remove_columns(names=['pair_ra', 'pair_dec', 'core_ra',
                                'core_dec', sizecol])
    
    ###subset best for joining with pairs
    xmtab.sort(sepcolname)
    bestmatch = unique(xmtab, keys=pair_id, keep='first')
    bestmatch.remove_column(name='sep_core_pmed')
    
    ###tidy up xmtab to save
    xmtab.sort(pair_id)
    xmtab.rename_column(name=core_id, new_name=corename_out)
    
    ###join best with pairs and cand_cores
    bestcores = join(bestmatch, cand_cores, keys=core_id, join_type='left')
    bestcores.rename_column(name=core_id, new_name=corename_out)
    pcdata = join(pairs, bestcores, keys=pair_id, join_type='left')
    
    return pcdata, xmtab


def find_allwise(data, acol='RA', dcol='DEC',
                 cwcols=['_q', 'AllWISE',
                         'W1mag', 'e_W1mag', 'W2mag', 'e_W2mag',
                         'W3mag', 'e_W3mag', 'W4mag', 'e_W4mag'],
                 sepcolname='sep_from_radio_source',
                 searchrad=30*u.arcsec):
    'query AllWISE via VizieR for host candidates using pairs data'
    
    ###creat skycoord list
    poscat = SkyCoord(ra=data[acol], dec=data[dcol])
    
    viz = Vizier(catalog='II/328', columns=['_q', '_r', '*']) ##setup vizier query and run
    viz.ROW_LIMIT = -1
    viz.TIMEOUT = 600
    awise = viz.query_region(poscat, radius=searchrad)[0] ##[0] returns astropy table
    
    ###join with data to create two outputs: all candidates and best match
    ###1) table of all returned candidates
    
    ##add index starting at 1 to data so as to join with _q from query
    data['_q'] = np.arange(len(data))+1
    
    ###create table of hosts
    hostcands = join(data[['pair_name', acol, dcol, '_q']], awise,
                     keys='_q', join_type='right')
    
    ###rename '_r' to sepcolname
    hostcands.rename_column(name='_r', new_name=sepcolname)
    
    ###2) data with added info about the closest CWISE host candidate
    ##subset useful columns
    bestmatch = hostcands[cwcols+[sepcolname]]
    
    ##sort by distance
    bestmatch.sort(keys=sepcolname)
    
    ##take only closest
    bestmatch = unique(bestmatch, keys='_q', keep='first')
    
    ###join with data
    dragns = join(data, bestmatch, keys='_q', join_type='left')

    
    ###remove superfluous columns from tables -- do this last
    hostcands.remove_columns(names=[acol, dcol, '_q'])
    dragns.remove_column(name='_q')
    data.remove_column(name='_q')
    
    ###clean metadata
    hostcands.meta = {}
    dragns.meta = {}
    ###alter description of sepcol
    hostcands[sepcolname].description = 'angular offset between radio and AllWISE positions'
    dragns[sepcolname].description = 'angular offset between radio and AllWISE positions'
    
    ###sort hosts by name and sep
    hostcands.sort(['pair_name', sepcolname])
                      
    return dragns, hostcands


def find_pairs(data, columndict = {'name': 'Component_name', 'ra': 'RA',
               'dec': 'DEC', 'peak': 'Peak_flux', 'total': 'Total_flux',
               'epeak': 'E_Peak_flux', 'etot': 'E_Total_flux', 'maj': 'DC_Maj',
               'min': 'DC_Min', 'pa': 'DC_PA'},
               only_unique_pairs=True):
    'find pair candidates, include basic neighbour info, sepatation, PA, and sep to 3rdNN'
    
    ###limit to only necessary columns - make generic
    #    usecols = ['Component_name', 'RA', 'DEC', 'Peak_flux', 'Total_flux',
    #               'DC_Maj', 'DC_Min', 'DC_PA']
    usecols = list(columndict.values())
    
    ###setup namecols
    name_1 = '_'.join([columndict['name'], '1'])
    name_2 = '_'.join([columndict['name'], '2'])
    
    ###setup paircols for hstack whilst at it
    p1cols = [i + '_1' for i in usecols]
    p2cols = [i + '_2' for i in usecols]
    
    pairdata = data[usecols]
    
    acol = columndict['ra']
    dcol = columndict['dec']
    spcol = columndict['peak']
    espcol = columndict['epeak']
    stcol = columndict['total']
    estcol = columndict['etot']
    majcol = columndict['maj']
    mincol = columndict['min']
    pacol = columndict['pa']
    
    ###add units to positions and pos angles - must be decimal degs (makes work if not fits/votable)
    pairdata[acol].unit = 'deg'
    pairdata[dcol].unit = 'deg'
    pairdata[pacol].unit = 'deg'
    
    ###sky position and find neighbours
    poscat = SkyCoord(ra=pairdata[acol], dec=pairdata[dcol])
    
    nn1 = match_coordinates_sky(poscat, poscat, nthneighbor=2)
    nn2 = match_coordinates_sky(poscat, poscat, nthneighbor=3)
    
    
    ###create table of pairs
    p1, p2 = pairdata.copy(), pairdata.copy()
    p1.rename_columns(names=usecols, new_names=p1cols)
    p2.rename_columns(names=usecols, new_names=p2cols)
    
    pairs = hstack([p1, p2[nn1[0]]])
    
    c1 = SkyCoord(ra=pairs[acol+'_1'], dec=pairs[dcol+'_1'])
    c2 = SkyCoord(ra=pairs[acol+'_2'], dec=pairs[dcol+'_2'])
    
    
    ###add in pair separation, position angle, and dist to 2nd NN, flux ratios, centroid
    ###calculate position angles
    pa_pair = c1.position_angle(c2).to(u.deg)
    sep_pair = c1.separation(c2)
    centroid = c1.directional_offset_by(pa_pair, sep_pair/2) ###need to add flux weighted
    flux_weighting = np.array(pairs[stcol+'_1']/(pairs[stcol+'_1'] + pairs[stcol+'_2']))
    ###use 1-flux weighting to bring FW centroid closer to brightest component
    fwcentroid = c1.directional_offset_by(pa_pair, sep_pair*(1-flux_weighting))
    
    
    dpa1 = np.array(pairs[pacol+'_1'] - pa_pair)
    dpa2 = np.array(pairs[pacol+'_2'] - pa_pair)
    dpa = [dpa1, dpa2]
    adpa = []
    for pa in dpa:
        pa[pa<-180] = pa[pa<-180]+360
        pa[pa<-90] = pa[pa<-90]+180
        pa[pa>90] = pa[pa>90]-180
        adpa.append(np.abs(pa))
    #        pa = pa*u.deg

    ###determine flux ratio errors
    tf_err = np.sqrt(pairs[estcol+'_1']**2 + pairs[estcol+'_2']**2)
    trat_rel_err = np.sqrt((pairs[estcol+'_1']/pairs[stcol+'_1'])**2
                           + (pairs[estcol+'_2']/pairs[stcol+'_2'])**2)
    prat_rel_err = np.sqrt((pairs[espcol+'_1']/pairs[spcol+'_1'])**2
                           + (pairs[espcol+'_2']/pairs[spcol+'_2'])**2)
    
    pairs['Sep_pair'] = nn1[1].to(u.arcsec)
    pairs['PA_pair'] = pa_pair
    pairs['dPA_1'] = dpa[0]*u.deg
    pairs['dPA_2'] = dpa[1]*u.deg
    pairs['abs_dPA_1'] = adpa[0]*u.deg
    pairs['abs_dPA_2'] = adpa[1]*u.deg
    pairs['mean_misalignment'] = (pairs['abs_dPA_1'] + pairs['abs_dPA_2'])/2
    pairs['Tflux_ratio'] = np.array(pairs[stcol+'_1']/pairs[stcol+'_2'])
    pairs['E_Tflux_ratio'] = pairs['Tflux_ratio']*trat_rel_err
    pairs['Pflux_ratio'] = np.array(pairs[spcol+'_1']/pairs[spcol+'_2'])
    pairs['E_Pflux_ratio'] = pairs['Pflux_ratio']*prat_rel_err
    pairs[stcol+'_pair'] = pairs[stcol+'_1'] + pairs[stcol+'_2']
    pairs[estcol+'_pair'] = tf_err
    pairs['medianRA'] = centroid.ra
    pairs['medianDEC'] = centroid.dec
    pairs['cenRA'] = fwcentroid.ra
    pairs['cenDEC'] = fwcentroid.dec
    pairs['d_2NN'] = nn2[1].to(u.arcsec)
    pairs['pair_name'] = source_name(ra=pairs['cenRA'], dec=pairs['cenDEC'])
    
    if only_unique_pairs==True:
        pairs = unique(pairs, 'pair_name')
        
    ###identify pairs where component is part of multiplt pairs, and flag smallest pair
    pairs = identify_pair_shares(data=pairs, cid1=name_1, cid2=name_2,
                                 pair_id='pair_name', sepcol='Sep_pair')
    pairs = flag_closest_pair(data=pairs, cid1=name_1, cid2=name_2, pair_id='pair_name',
                              sort_on='Sep_pair', tokeep='first')
    
    return(pairs)


def tidy_pairs(pairs, corefilcols, pair_outcols,
               pairflux_col, coreflux_col, sourceflux_col,
               pairflux_err, coreflux_err, sourceflux_err,
               acol_c, dcol_c, acol_p='cenRA', dcol_p='cenDEC',
               cpromcol='CoreProm', e_cpromcol='E_CoreProm',
               ra_out='RA_best', de_out='DE_best'):
    'tidy up data ready for output'
    ###add value best RA/DEC (core pos or flux weighted, source total flux+err, core prominence)
    ###fill masked float values
    filval=np.nan
    fill_dict = {}
    for col in corefilcols:
        fill_dict[col] = filval
    
    pairs = fill_masks(data=pairs, coldict=fill_dict)
    
    ###add value
    ##total flux of the source, create arrays then add column
    sunit = pairs[pairflux_col].unit
    raunit = pairs[acol_p].unit
    deunit = pairs[dcol_p].unit
    spair = np.array(pairs[pairflux_col])
    e_spair = np.array(pairs[pairflux_err])
    score = np.array(pairs[coreflux_col])
    e_score = np.array(pairs[coreflux_err])
    best_ra = np.array(pairs[acol_c])
    best_de = np.array(pairs[dcol_c])
    pair_ra = np.array(pairs[acol_p])
    pair_de = np.array(pairs[dcol_p])
    
    ##mask for no core
    nocore_mask = np.isnan(score)
    
    ###set no flux to zero for summing
    score[nocore_mask] = 0
    e_score[nocore_mask] = 0
    
    stot = (spair + score)
    e_stot = np.sqrt(e_spair**2 + e_score**2)
    coreprom = score/stot
    e_coreprom = coreprom*np.sqrt((e_stot/stot)**2 + (e_score/score)**2)
    coreprom[nocore_mask] = filval
    e_coreprom[nocore_mask] = filval
    
    ##best coordinates - set to pair coords if no core
    best_ra[nocore_mask] = pair_ra[nocore_mask]
    best_de[nocore_mask] = pair_de[nocore_mask]
    
    pairs[sourceflux_col] = stot*sunit
    pairs[sourceflux_err] = e_stot*sunit
    pairs[cpromcol] = coreprom
    pairs[e_cpromcol] = e_coreprom
    pairs[ra_out] = best_ra*raunit
    pairs[de_out] = best_de*deunit
    
    ###subset pair columns -- individual component info can be dropped, keep core position
    pairs = pairs[pair_outcols]
    
    return pairs


def preprocess_cat(catfile, flux_col, flux_min, size_col, size_min, fileformat='fits'):
    'preprocess radio catalogue to filter only components matching desired criteria'
    alldata = Table.read(catfile, format=fileformat)
    ldata = alldata[(alldata[flux_col]>=flux_min) & (alldata[size_col]>=size_min)]
    cdata = alldata[(alldata[flux_col]>=flux_min) & (alldata[size_col]<size_min)]
    
    return(ldata, cdata)


def hunt_dragns(catfile, config_file, write_file=True, find_hosts=False,
                host_chunk_size=5000, host_search_ra='RA_best', host_search_de='DE_best',
                search_rad=30*u.arcsec):
    'find pairs and flag candidate dragns'
    
    ##extract config params
    catparams, coldict = config_parse(config_file=config_file)
    
    ##load and preprocesses cat data
    catdata, coredata = preprocess_cat(catfile=catfile, flux_col=coldict['peak'],
                                       flux_min=catparams['flux_min'],
                                       size_col=coldict['maj'],
                                       size_min=catparams['lobesize_min'],
                                       fileformat=catparams['file_format'])
    
    ###remove nonessential columns from core data
    coredata = coredata[list(coldict.values())]
    
    ###find pairs and flag (need to add flagging)
    pairs = find_pairs(data=catdata, columndict=coldict)
    
    ###add in core data
    coresepcolname = 'sep_core_pcen'
    corename_col = 'candidate_core'
    pairs, candcores = find_cores(pairs=pairs, cand_cores=coredata, acol_p='cenRA',
                                  acol_med='medianRA', dcol_med='medianDEC',
                                  dcol_p='cenDEC', acol_c=coldict['ra'],
                                  dcol_c=coldict['dec'],
                                  sizecol='Sep_pair', sepcolname=coresepcolname,
                                  pair_id='pair_name', core_id=coldict['name'],
                                  corename_out=corename_col, search_rad=search_rad)
    
    ###add _c to core columns
    for col in pairs.colnames:
        if col in list(coldict.values()):
            pairs.rename_column(name=col, new_name='_'.join([col, 'c']))
    
    cfill_cols = [coresepcolname] + [i for i in pairs.colnames if i.rsplit('_')[-1]=='c']

    ###tidy up data
    best_acol = 'RA_best'
    best_dcol = 'DE_best'
    cpcol = 'CoreProm'
    e_cpcol = 'E_CoreProm'
    sflux_col = coldict['total']+'_source'
    sflux_err = coldict['etot']+'_source'
    
    poutcols = ['pair_name', best_acol, best_dcol,
                'Sep_pair', sflux_col, sflux_err,
                cpcol, e_cpcol, 'Tflux_ratio', 'E_Tflux_ratio',
                'Pflux_ratio', 'E_Pflux_ratio', coldict['total']+'_pair',
                coldict['etot']+'_pair', 'PA_pair', 'mean_misalignment',
                'cenRA', 'cenDEC', 'medianRA', 'medianDEC',
                coldict['ra']+'_c', coldict['dec']+'_c', coresepcolname,
                'd_2NN', 'component_shares', 'pref_pair',
                corename_col, coldict['name']+'_1', coldict['name']+'_2']
    
    
    pairs = tidy_pairs(pairs=pairs, corefilcols=cfill_cols, pair_outcols=poutcols,
                       pairflux_col=coldict['total']+'_pair',
                       coreflux_col=coldict['total']+'_c',
                       sourceflux_col=sflux_col,
                       pairflux_err=coldict['etot']+'_pair',
                       coreflux_err=coldict['etot']+'_c',
                       sourceflux_err=sflux_err, acol_c=coldict['ra']+'_c',
                       dcol_c=coldict['dec']+'_c',
                       acol_p='cenRA', dcol_p='cenDEC', cpromcol=cpcol,
                       e_cpromcol=e_cpcol, ra_out=best_acol, de_out=best_dcol)
    
    ##find host candidates - make this an argument
    ###need to break down large tables and do queries in sequence (e.g. if table > 10k rows... split into chunks then vstack results)
    if find_hosts == True:
        ###if len pairs > host_chunk_size split up table and run host finding in sections
        if len(pairs)>host_chunk_size:
            nrows = str(len(pairs)) + ' rows'
            chunksize = str(int(host_chunk_size)) + ' rows'
            print('')
            print('WARNING: pairs table is large (' +nrows+ '), host finding being performed in chunks of ' + chunksize)
            print('')
            dtabs, htabs = [], []
            steps = np.arange(0, len(pairs)+host_chunk_size, host_chunk_size)
            for i in range(len(steps)-1):
                subpairs = pairs[steps[i]: steps[i+1]]
                dtab, htab = find_allwise(data=subpairs, acol=host_search_ra,
                                          dcol=host_search_de, searchrad=search_rad)
                dtabs.append(dtab)
                htabs.append(htab)
            dragns = vstack(dtabs)
            hosts = vstack(htabs)
        else:
            dragns, hosts = find_allwise(data=pairs, acol=host_search_ra,
                                         dcol=host_search_de, searchrad=search_rad)
    else:
        print('PAIR INFO ONLY; HOST FINDING NOT PERFORMED')
        dragns = pairs

    ###round columns to sf
    cround_dict = {'sep_core_pmed': 2, coresepcolname: 2}

    dround_dict = {best_acol: 5, best_dcol: 5, 'Sep_pair': 2, sflux_col: 3,
                   sflux_err: 3, cpcol: 5, e_cpcol: 5, 'Tflux_ratio': 5,
                   'E_Tflux_ratio': 5, 'Pflux_ratio': 5, 'E_Pflux_ratio': 5,
                   coldict['total']+'_pair': 3, coldict['etot']+'_pair': 3,
                   'PA_pair': 2, 'mean_misalignment': 2, 'cenRA': 5, 'cenDEC': 5,
                   'medianRA': 5, 'medianDEC': 5, coldict['ra']+'_c': 5,
                   coldict['dec']+'_c': 5, coresepcolname: 2, 'd_2NN': 2,
                   'sep_from_radio_source': 2}

    dragns = round_cols(data=dragns, round_dict=dround_dict)
    candcores = round_cols(data=candcores, round_dict=cround_dict)
    
    candcores.write('candidate_cores.fits', format='fits')

    ###write to file or return data in code
    if write_file==True:
        ##create new file name based on old
        if find_hosts == True:
            dragns.write('DRAGNs.fits', format='fits')
            hosts.write('host_candidates.fits', format='fits')
        else:
            dragns.write('DRAGNs.fits', format='fits')
        return
    else:
        if find_hosts==True:
            return(dragns, hosts)
        else:
            return dragns


def config_parse(config_file):
    'extract required info from config file'
    cdata = pd.read_table(config_file, delim_whitespace=True).replace("'", "", regex=True)
    
    ##parameters for catalog preprocessing
    minflux = np.float(cdata['value'][np.where(cdata['parameter']=='min_brightness')[0][0]])
    minsizelobe = np.float(cdata['value'][np.where(cdata['parameter']=='min_size_lobe')[0][0]])
#    maxsizecore = np.float(cdata['value'][np.where(cdata['parameter']=='max_size_core')[0][0]])
    dformat = cdata['value'][np.where(cdata['parameter']=='data_format')[0][0]]
    
    catparams = {'flux_min': minflux, 'lobesize_min': minsizelobe,
                 'file_format':dformat}
    
    ##catalogue column names
    namecol = cdata['value'][np.where(cdata['parameter']=='name_col')[0][0]]
    racol = cdata['value'][np.where(cdata['parameter']=='ra_col')[0][0]]
    deccol = cdata['value'][np.where(cdata['parameter']=='dec_col')[0][0]]
    peakcol = cdata['value'][np.where(cdata['parameter']=='Speak_col')[0][0]]
    epeakcol = cdata['value'][np.where(cdata['parameter']=='Speak_err_col')[0][0]]
    totcol = cdata['value'][np.where(cdata['parameter']=='Stot_col')[0][0]]
    etotcol = cdata['value'][np.where(cdata['parameter']=='Stot_err_col')[0][0]]
    majcol = cdata['value'][np.where(cdata['parameter']=='maj_col')[0][0]]
    mincol = cdata['value'][np.where(cdata['parameter']=='min_col')[0][0]]
    pacol = cdata['value'][np.where(cdata['parameter']=='pa_col')[0][0]]
    
    coldict = {'name': namecol, 'ra': racol, 'dec': deccol, 'peak': peakcol,
               'epeak': epeakcol, 'total': totcol, 'etot': etotcol,'maj': majcol,
               'min': mincol, 'pa': pacol}
    
    return(catparams, coldict)


def parse_args():
    "parse input args, i.e. target and config file names"
    parser = argparse.ArgumentParser(description="find radio pairs in catalogue data")
    parser.add_argument("radio_cat", help="radio component catalogue to find pairs in")
    parser.add_argument("--config", action='store', type=str, default='config.txt',
                        help="config file")
    parser.add_argument("--findhosts", action='store', type=str, default='True',
                        help="if True look for WISE host candidates")
    parser.add_argument("--search_radius", action='store', type=str, default='30arcsec',
                        help="search radius to use for core and host candidates")
                                     
    args = parser.parse_args()
    
    ##make args.writepairs a bool and search_radius a quantity
    args.findhosts = strtobool(args.findhosts)
    args.search_radius = u.Quantity(args.search_radius)
    
    return args


##############################################################################
###main code

if __name__ == '__main__':
    args = parse_args()
    hunt_dragns(catfile=args.radio_cat, config_file=args.config, write_file=True,
                find_hosts=args.findhosts, search_rad=args.search_radius)

    #########################################################################
    print('END: ', np.round(time.time()-start_time, 2), 's')
