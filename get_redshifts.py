###code to get redshifts

import numpy as np, sys, warnings, argparse, os
from distutils.util import strtobool
from astropy import units as u
from astroquery.vizier import Vizier
from astroquery.utils.tap.core import TapPlus
from astroquery.xmatch import XMatch
from astropy.table import Table, join, unique, vstack, hstack, MaskedColumn
from astropy.coordinates import SkyCoord
from astropy.constants import c
from astropy.utils.exceptions import AstropyWarning
warnings.simplefilter('ignore', category=AstropyWarning) ###mutes mergeconflict warnings in printout


####sql queries to vizier have maximum upload size of 100k rows

########################################
#time of code check - not important for code to run
import time
start_time = time.time()
########################################

################################################################################
################################################################################
###parameters

sourcecols_out = ['Name', 'RA', 'DEC', 'Flux', 'E_Flux', 'LAS', 'E_LAS', 'Type',
                  'Source_flag', 'AllWISE', 'Sep_AllWISE', 'LR', 'Rel', 'Host_flag',
                  'W1mag', 'e_W1mag', 'W2mag', 'e_W2mag', 'W3mag', 'e_W3mag',
                  'W4mag', 'e_W4mag', 'z', 'z_err', 'z_type', 'z_survey']

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


################################################################################
################################################################################
###functions

def load_targets(source_file, host_file, namecol='AllWISE',
                 acol='RAJ2000', dcol='DEJ2000', posunits=('deg', 'deg')):
    'code get AllWISE positions of sources to search for redshift data for'
    ###load data
    sources = Table.read(source_file)
    hosts = Table.read(host_file)
    
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
    targets = join(sources[[namecol]], hosts[[namecol, acol, dcol]])
    
    return targets


def write_query(table_to_query='II/328/allwise',
                search_rad=30*u.arcsec, searchprec=0,
                upload_name='uploaded_data',
                upra='RA', updec='DEC',
                qra='RAJ2000', qdec='DEJ2000'):
    'construct SQL query string from input parameters'
    
    srad = str(np.round(search_rad.to('arcsec').value, searchprec))
    
    qstring = f"SELECT * \n FROM tap_upload.{upload_name} AS tup \n JOIN \"{table_to_query}\" AS db \n ON 1=CONTAINS(POINT('ICRS', tup.{upra}, tup.{updec}), CIRCLE('ICRS', db.{qra}, db.{qdec}, {srad}/3600.))"
    
    return qstring


def objectcols_to_str(data):
    'convert dtype==object cols to str to allow saving as fits files'
    ##alternative is to save as votable
    cols = data.colnames
    for col in cols:
        if data[col].dtype == object:
            data[col] = np.array(data[col]).astype(str)
    
    return
    

def tapq_vizier(data, acol='RAJ2000', dcol='DEJ2000',
                namecol='AllWISE', viztable='VII/292/north',
                search_radius=1*u.arcsec,
                vizra='RAJ2000', vizdec='DEJ2000',
                upload_name='uploaded_data',
                sep_col='ang_sep', verbose=True,
                qcols=None):
    'function to upload a table and perform positional cross match against a Vizier table'
    ##defaults to allwise if no vizier table input'
    
    ###mark time for testing
    t0 = time.time()
    
    ##construct query
    query = write_query(table_to_query=viztable,
                        search_rad=search_radius,
                        upload_name=upload_name,
                        upra=acol, updec=dcol,
                        qra=vizra, qdec=vizdec)
    
    ##setup tap and launch job
    tap = TapPlus(url='http://tapvizier.u-strasbg.fr/TAPVizieR/tap')
    job = tap.launch_job_async(query=query, upload_resource=data,
                               upload_table_name=upload_name,
                               verbose=verbose)
    
    ###get results
    outdata = job.get_data()
    
    ###rename viz pos cols if the same as input cols
    if vizra==acol:
        vizra = '_'.join([vizra, '2'])
    if vizdec==dcol:
        vizdec = '_'.join([vizdec, '2'])
    
    ###calculate angular separation
    posin = SkyCoord(ra=outdata[acol], dec=outdata[dcol])
    posviz = SkyCoord(ra=outdata[vizra], dec=outdata[vizdec])
    angsep = posin.separation(posviz).to(search_radius.unit)
    
    ###round to appropriate number of sig figs
    sf_dict = {u.arcsec: 2, u.arcmin: 6, u.deg:9}
    angsep = np.round(angsep, sf_dict[search_radius.unit])
    
    ##add to outdata
    outdata[sep_col] = angsep
    outdata[sep_col].unit = search_radius.unit
    
    ###subset columns to output
    if qcols is not None:
        outcols = [namecol] + qcols + [sep_col]
        outdata = outdata[outcols]
    
    ###mark time for testing
    if verbose==True:
        dt = time.time()-t0
        print('')
        print('time for host finding (N='+str(len(data))+') = ' + str(np.round(dt, 2)) +'s')
        print('')
    
    return outdata


def find_photozs(data, namecol='AllWISE', acol='RAJ2000', dcol='DEJ2000',
                 vizra='RAJ2000', vizdec='DEJ2000',
                 vizkey='VII/292/north',
                 zpcols=['zphot', 'e_zphot', 'fqual'],
                 sepcolname='ang_sep',
                 searchrad=1*u.arcsec,
                 chunk_size=40000):
    'query AllWISE using tap-vizier'
    ###query VizieR
    ###split into chunks if large input
    dlen = len(data)
    indata = data[[namecol, acol, dcol]]
    if dlen > chunk_size:
        job_results = []
        n_chunks = int(np.ceil(dlen/chunk_size))
        for i in range(n_chunks):
            print('')
            print('query chunk '+str(i+1)+'/'+str(n_chunks))
            print('')
            dchunk = indata[i*chunk_size: (i+1)*chunk_size]
            job = tapq_vizier(data=dchunk, acol=acol, dcol=dcol,
                              namecol=namecol, viztable=vizkey,
                              search_radius=searchrad,
                              qcols=zpcols,
                              sep_col=sepcolname, verbose=True)
            job_results.append(job)
        tap_results = vstack(job_results)
    else:
        tap_results = tapq_vizier(data=indata, acol=acol, dcol=dcol,
                                  namecol=namecol, viztable=vizkey,
                                  qcols=zpcols,
                                  search_radius=searchrad,
                                  sep_col=sepcolname, verbose=True)
    
    ###split into two tables (all host candidates, and source table with closest)
    ##remove acol, dcol from tap_results to prevent duplication
    tap_results.meta = {} ###removes meta for clean joins
    objectcols_to_str(data=tap_results) ###convert object cols to string for saving as fits
    tap_results.sort(sepcolname)
    bestmatch = unique(tap_results, namecol)
    
    ###convert bestmatch namecol dtype
    bestmatch[namecol] = np.array(bestmatch[namecol]).astype(str)
    bestmatch.sort(namecol)

    return bestmatch


def query_lsdr8_north_and_south(data, namecol='AllWISE', acol='RAJ2000',
                                dcol='DEJ2000', vizra='RAJ2000', vizdec='DEJ2000',
                                zpcols=['zphot', 'e_zphot', 'fqual'],
                                sepcolname='ang_sep', searchrad=1*u.arcsec,
                                chunk_size=20000, goodflag=True):
    'run find_photo_zs on both north and south tables and tidy results into single table'
    
    ###select subset in NGC for north query, everything else for south query
    ##determine b and use np arrays for filtering to avoid quantity issues
    ## NGC: b > 0 & d > 32.375deg
    skypos = SkyCoord(ra=data[acol], dec=data[dcol])
    galb = np.array(skypos.galactic.b.deg)
    decl = np.array(skypos.dec.deg)
    nfilt = (galb > 0) & (decl > 32.375)
    dnorth = data[nfilt]
    dsouth = data[~nfilt]
    
    ###account for possibility of data only in one or other fields and allow to still output
    if len(dnorth)>0:
        north = find_photozs(data=dnorth, namecol=namecol, acol=acol, dcol=dcol,
                             vizra=vizra, vizdec=vizdec, vizkey='VII/292/north',
                             zpcols=zpcols, sepcolname=sepcolname,
                             searchrad=searchrad, chunk_size=chunk_size)
    else:
        north = None
    
    if len(dsouth)>0:
        south = find_photozs(data=dsouth, namecol=namecol, acol=acol, dcol=dcol,
                             vizra=vizra, vizdec=vizdec, vizkey='VII/292/south',
                             zpcols=zpcols, sepcolname=sepcolname,
                             searchrad=searchrad, chunk_size=chunk_size)
    else:
        south = None
    
    table_list = [north, south]
    outtables = [t for t in table_list if t is not None]
                            
    output = vstack(outtables)
    
    ###keep only good flagged
    if goodflag==True:
        output = output[(output['fqual']==1)]
        output.remove_column('fqual') ##no need to output this column if using here
    
    ###remove duplicates, keep closest match
    output.sort(sepcolname)
    output = unique(output, namecol, keep='first')
    output.sort(namecol)
    
    return output


def cds_xmatch(data, racol='RAJ2000', decol='DEJ2000', maxsep=1*u.arcsec,
               catcols=['objID', 'zsp', 'e_zsp', 'Q', 'angDist'],
               namecol='AllWISE', colsuff=None,
               cat2='vizier:V/154/sdss16', zcol='zsp',
               timeout=1200):
    'use cds xmatch to query sdss'
    ###try to replace with async query (might not need to)
    xm = XMatch()
    xm.TIMEOUT = timeout
    
    xmatch = xm.query(cat1=data, cat2=cat2, max_distance=maxsep,
                      colRA1=racol, colDec1=decol)
    
    ###reduce to only spec data (and unique)
    outcols = [namecol] + catcols
    
    if type(xmatch[zcol]) == MaskedColumn:
        xmatch = xmatch[~xmatch[zcol].mask]
    xmatch.sort('angDist')
    
    if len(xmatch)>0:
        maxad = np.round(np.max(xmatch['angDist']), 3)
        print(f'max angDist = {maxad}')
        print(f'max search rad = {maxsep}')
        print('')
        
        xmatch = unique(xmatch, namecol) ###ensures no duplicates
    
    ###add in sep col units
    xmatch['angDist'] = np.round(xmatch['angDist'], 2)
    xmatch['angDist'].unit = maxsep.unit
    
    ###subset required cols and sort
    xmatch = xmatch[outcols]
    if colsuff is not None:
        newnames = ['_'.join([col, colsuff]) for col in catcols]
        xmatch.rename_columns(catcols, newnames)
    xmatch.sort(namecol)
    
    return xmatch


def get_spec_zs(targets, racol='RAJ2000', decol='DEJ2000', namecol='AllWISE',
                search_rad=1*u.arcsec, xmtimeout=1200,
                chunksize=50000):
    'obtain spectroscopic redshifts from: SDSS, GAMA, 2dF, 6dF, WiggleZ, 2MRS'
    
    ###subset to only name and position for querys to CDS xmatch service
    cat1 = targets[[namecol, racol, decol]]
    
    ######query each survey using cds xmatch -- can add queries to this later if new data becomes available
    print('')
    print('Querying spectroscopic catalogs in CDS XMatch servive')
    print('')
    if len(cat1)>chunksize:
        print('')
        print('Upload data is large, querying in chunks')
        
        n_chunks = int(np.ceil(len(cat1)/chunksize))
        
        ###SDSS
        print('SDSS DR16; Ahumada+ (2020, ApJS, 249, 3)')
        print('')
        sdss_results = []
        for i in range(n_chunks):
            print('')
            print('query chunk ' + str(i+1) + '/' + str(n_chunks))
            print('')
            data_chunk = cat1[i*chunksize: (i+1)*chunksize]
            job = cds_xmatch(data=data_chunk, racol=racol, decol=decol,
                             maxsep=search_rad,
                             catcols=['objID', 'zsp', 'e_zsp', 'Q', 'angDist'],
                             namecol=namecol, colsuff='SDSS-DR16',
                             cat2='vizier:V/154/sdss16', zcol='zsp',
                             timeout=xmtimeout)
            if len(job)>0:
                sdss_results.append(job)
        if len(sdss_results)>0:
            xm_sdss = vstack(sdss_results)
        else:
            xm_sdss = []
        
        ###GAMA
        print('GAMA DR3; Baldry+ (2018, MNRAS, 474, 3875)')
        print('')
        gama_results = []
        for i in range(n_chunks):
            print('')
            print('query chunk ' + str(i+1) + '/' + str(n_chunks))
            print('')
            data_chunk = cat1[i*chunksize: (i+1)*chunksize]
            job = cds_xmatch(data=data_chunk, racol=racol, decol=decol,
                             maxsep=search_rad,
                             catcols=['CataId', 'z', 'NQ', 'Survey', 'angDist'],
                             namecol=namecol, colsuff='GAMA',
                             cat2='vizier:J/MNRAS/474/3875/gamadr3', zcol='z',
                             timeout=xmtimeout)
            ##subset high q data and survey not in other survey searches
            job = job[((job['NQ_GAMA']>2) & (job['Survey_GAMA']!='SDSS')
                       & (job['Survey_GAMA']!='6dFGS')
                       & (job['Survey_GAMA']!='2dFGRS')
                       & (job['Survey_GAMA']!='WiggleZ'))]
            if len(job)>0:
                gama_results.append(job)
        if len(gama_results)>0:
            xm_gama = vstack(gama_results)
        else:
            xm_gama = []
    
        ###2dF
        print('2dFGRS; Colless+ (2001, MNRAS, 328, 1039)')
        print('')
        df2_results = []
        for i in range(n_chunks):
            print('')
            print('query chunk ' + str(i+1) + '/' + str(n_chunks))
            print('')
            data_chunk = cat1[i*chunksize: (i+1)*chunksize]
            job = cds_xmatch(data=data_chunk, racol=racol, decol=decol,
                             maxsep=search_rad,
                             catcols=['SeqNum', 'z', 'q_z', 'angDist'],
                             namecol=namecol, colsuff='2dFGRS',
                             cat2='vizier:VII/250/2dfgrs', zcol='z',
                             timeout=xmtimeout)
            ##subset hiqh q data
            job = job[(job['q_z_2dFGRS']>2)]
            if len(job)>0:
                df2_results.append(job)
        if len(df2_results)>0:
            xm_2df = vstack(df2_results)
        else:
            xm_2df = []
    
        ###6dF
        print('6dFGS; Jones+ (2009, MNRAS, 399, 683)')
        print('')
        df6_results = []
        for i in range(n_chunks):
            print('')
            print('query chunk ' + str(i+1) + '/' + str(n_chunks))
            print('')
            data_chunk = cat1[i*chunksize: (i+1)*chunksize]
            job = cds_xmatch(data=data_chunk, racol=racol, decol=decol,
                             maxsep=search_rad,
                             catcols=['6dFGS', 'cz', 'e_cz', 'q_cz', 'angDist'],
                             namecol=namecol, colsuff='6dFGS',
                             cat2='vizier:VII/259/6dfgs', zcol='cz',
                             timeout=xmtimeout)
            ##determine z from cz
            job['z_6dFGS'] = np.round(job['cz_6dFGS']/c.to('km/s').value, 5)
            job['e_z_6dFGS'] = np.round(job['e_cz_6dFGS']/c.to('km/s').value, 5)
            ##subset good quality
            job = job[(job['q_cz_6dFGS']>2)]
            if len(job)>0:
                df6_results.append(job)
        if len(df6_results)>0:
            xm_6df = vstack(df6_results)
        else:
            xm_6df = []
        
        ###WiggleZ
        print('WiggleZ; Drinkwater+ (2018, MNRAS, 474, 4151)')
        print('')
        wigglez_results = []
        for i in range(n_chunks):
            print('')
            print('query chunk ' + str(i+1) + '/' + str(n_chunks))
            print('')
            data_chunk = cat1[i*chunksize: (i+1)*chunksize]
            job = cds_xmatch(data=data_chunk, racol=racol, decol=decol,
                             maxsep=search_rad,
                             catcols=['Name', 'z', 'q_z', 'angDist'],
                             namecol=namecol, colsuff='WiggleZ',
                             cat2='vizier:J/MNRAS/474/4151/wigglez', zcol='z',
                             timeout=xmtimeout)
            ##subset high q data
            job = job[(job['q_z_WiggleZ']>2)]
            if len(job)>0:
                wigglez_results.append(job)
        if len(wigglez_results)>0:
            xm_wigglez = vstack(wigglez_results)
        else:
            xm_wigglez = []
        
        ###2MRS
        print('2MRS; Huchra+ (2012, ApJS, 199, 26)')
        print('')
        mrs2_results = []
        for i in range(n_chunks):
            print('')
            print('query chunk ' + str(i+1) + '/' + str(n_chunks))
            print('')
            data_chunk = cat1[i*chunksize: (i+1)*chunksize]
            job = cds_xmatch(data=data_chunk, racol=racol, decol=decol,
                             maxsep=search_rad,
                             catcols=['ID', 'cz', 'CAT', 'angDist'],
                             namecol=namecol, colsuff='2MRS',
                             cat2='vizier:J/ApJS/199/26/table3', zcol='cz',
                             timeout=xmtimeout)
            ##subset high q data
            job['z_2MRS'] = np.round(job['cz_2MRS']/c.to('km/s').value, 5)
            if len(job)>0:
                mrs2_results.append(job)
        if len(mrs2_results)>0:
            xm_2mr = vstack(mrs2_results)
        else:
            xm_2mr = []

    else:
        ###SDSS
        print('SDSS DR16; Ahumada+ (2020, ApJS, 249, 3)')
        print('')
        xm_sdss = cds_xmatch(data=cat1, racol=racol, decol=decol, maxsep=search_rad,
                             catcols=['objID', 'zsp', 'e_zsp', 'Q', 'angDist'],
                             namecol=namecol, colsuff='SDSS-DR16',
                             cat2='vizier:V/154/sdss16', zcol='zsp',
                             timeout=xmtimeout)
    
        ###GAMA
        print('GAMA DR3; Baldry+ (2018, MNRAS, 474, 3875)')
        print('')
        xm_gama = cds_xmatch(data=cat1, racol=racol, decol=decol, maxsep=search_rad,
                             catcols=['CataId', 'z', 'NQ', 'Survey', 'angDist'],
                             namecol=namecol, colsuff='GAMA',
                             cat2='vizier:J/MNRAS/474/3875/gamadr3', zcol='z',
                             timeout=xmtimeout)
        ##subset high q data and survey not in other survey searches
        xm_gama = xm_gama[((xm_gama['NQ_GAMA']>2) & (xm_gama['Survey_GAMA']!='SDSS')
                           & (xm_gama['Survey_GAMA']!='6dFGS')
                           & (xm_gama['Survey_GAMA']!='2dFGRS')
                           & (xm_gama['Survey_GAMA']!='WiggleZ'))]

        ###2dF
        print('2dFGRS; Colless+ (2001, MNRAS, 328, 1039)')
        print('')
        xm_2df = cds_xmatch(data=cat1, racol=racol, decol=decol, maxsep=search_rad,
                            catcols=['SeqNum', 'z', 'q_z', 'angDist'],
                            namecol=namecol, colsuff='2dFGRS',
                            cat2='vizier:VII/250/2dfgrs', zcol='z',
                            timeout=xmtimeout)
        ##subset hiqh q data
        xm_2df = xm_2df[(xm_2df['q_z_2dFGRS']>2)]
    
        ###6dF
        print('6dFGS; Jones+ (2009, MNRAS, 399, 683)')
        print('')
        xm_6df = cds_xmatch(data=cat1, racol=racol, decol=decol, maxsep=search_rad,
                            catcols=['6dFGS', 'cz', 'e_cz', 'q_cz', 'angDist'],
                            namecol=namecol, colsuff='6dFGS',
                            cat2='vizier:VII/259/6dfgs', zcol='cz',
                            timeout=xmtimeout)
        ##determine z from cz
        xm_6df['z_6dFGS'] = np.round(xm_6df['cz_6dFGS']/c.to('km/s').value, 5)
        xm_6df['e_z_6dFGS'] = np.round(xm_6df['e_cz_6dFGS']/c.to('km/s').value, 5)
        ##subset good quality
        xm_6df = xm_6df[(xm_6df['q_cz_6dFGS']>2)]
    
        ###WiggleZ
        print('WiggleZ; Drinkwater+ (2018, MNRAS, 474, 4151)')
        print('')
        xm_wigglez = cds_xmatch(data=cat1, racol=racol, decol=decol, maxsep=search_rad,
                                catcols=['Name', 'z', 'q_z', 'angDist'],
                                namecol=namecol, colsuff='WiggleZ',
                                cat2='vizier:J/MNRAS/474/4151/wigglez', zcol='z',
                                timeout=xmtimeout)
        ##subset high q data
        xm_wigglez = xm_wigglez[(xm_wigglez['q_z_WiggleZ']>2)]
    
        ###2MRS
        print('2MRS; Huchra+ (2012, ApJS, 199, 26)')
        print('')
        xm_2mr = cds_xmatch(data=cat1, racol=racol, decol=decol, maxsep=search_rad,
                            catcols=['ID', 'cz', 'CAT', 'angDist'],
                            namecol=namecol, colsuff='2MRS',
                            cat2='vizier:J/ApJS/199/26/table3', zcol='cz',
                            timeout=xmtimeout)
        ##subset high q data
        xm_2mr['z_2MRS'] = np.round(xm_2mr['cz_2MRS']/c.to('km/s').value, 5)

    ####join all tables -- allow for no matches
    sz = cat1.copy()
    if len(xm_sdss)>0:
        sz = join(sz, xm_sdss, keys=namecol, join_type='left')
    if len(xm_gama)>0:
        sz = join(sz, xm_gama, keys=namecol, join_type='left')
    if len(xm_2df)>0:
        sz = join(sz, xm_2df, keys=namecol, join_type='left')
    if len(xm_wigglez)>0:
        sz = join(sz, xm_wigglez, keys=namecol, join_type='left')
    if len(xm_6df)>0:
        sz = join(sz, xm_6df, keys=namecol, join_type='left')
    if len(xm_2mr)>0:
        sz = join(sz, xm_2mr, keys=namecol, join_type='left')
    
    ###break here if no spec-zs found (use collist len in sz as proxy)
    if len(sz.colnames) <= len(cat1.colnames):
        print('no spec-zs found')
        return
    
    else:
        ####select best z data for targets
        ##dictionary of required survey info (in order of priority) -- will need to add to this if adding additional survey catalogs
        spec_surveys = {'SDSS-DR16': {'z': 'zsp_SDSS-DR16', 'ez': 'e_zsp_SDSS-DR16'},
                        'GAMA': {'z': 'z_GAMA', 'ez': None},
                        '2dFGRS': {'z': 'z_2dFGRS', 'ez': None},
                        '6dFGS': {'z': 'z_6dFGS', 'ez': 'e_z_6dFGS'},
                        'WiggleZ': {'z': 'z_WiggleZ', 'ez': None},
                        '2MRS': {'z': 'z_2MRS', 'ez': None}}
    
        ###set up empty arrays for z, e_z and survey name
        svy = np.empty(len(sz), dtype='<U'+str(int(np.max([len(i) for i in list(spec_surveys.keys())]))))
        z, ez = np.nan*(np.zeros(len(sz))), np.nan*(np.zeros(len(sz)))
    
        ###iterate through surveys and add data to z, e_z and survey name
        for survey in list(spec_surveys.keys()):
            ##apply only if in sz
            if spec_surveys[survey]['z'] in sz.colnames:
                ##determine masks
                if type(sz[spec_surveys[survey]['z']]) != MaskedColumn: ##ensure mask available
                    sz[spec_surveys[survey]['z']] = MaskedColumn(sz[spec_surveys[survey]['z']])
                zmask = sz[spec_surveys[survey]['z']].mask ##no need to assign to dict here, just use immediately!
                ##if z available replace nan in z array (and error if possible) + survey
                zupdate = ~zmask & np.isnan(z)
                z[zupdate] = sz[spec_surveys[survey]['z']][zupdate]
                svy[zupdate] = survey
                if spec_surveys[survey]['ez'] is not None:
                    ez[zupdate] = sz[spec_surveys[survey]['ez']][zupdate]
    
        ###create final table
        zdat = Table({'z': z, 'z_err': ez, 'Survey': svy})
        zdat = hstack([cat1[namecol], zdat])
    
        return zdat


def combine_spec_photo_zs(speczs, photozs, namecol='AllWISE',
                          zscol='z', zpcol='zphot', ezscol='z_err',
                          ezpcol='e_zphot', zcol_out='z',
                          ezcol_out='z_err', ztout='z_type',
                          zs_surveycol='Survey', surveycol='z_survey',
                          photo_survey='LS-DR8', zdp=5):
    'combine spec and photo zs'
    
    spec_comb = join(speczs, photozs, keys=namecol, join_type='left')
    
    ###mask fill z data if needed
    for col in [zscol, ezscol, zpcol, ezpcol]:
        if type(spec_comb[col])==MaskedColumn:
            spec_comb[col].fill_value = np.nan
            spec_comb[col] = spec_comb[col].filled()
    
    ###set up arrays (zspec as default data)
    z = np.array(spec_comb[zscol])
    ez = np.array(spec_comb[ezscol])
    zp = np.array(spec_comb[zpcol])
    ezp = np.array(spec_comb[ezpcol])
    survey = np.array(spec_comb[zs_surveycol])
    
    z_type = np.empty(len(z), dtype='<U5')
    z_type[~np.isnan(z)] = 'spec'
    
    ##filter for no z-spec
    zspec_filt = np.isnan(z) & ~np.isnan(zp)
    
    ##replace with photo if needed
    z[zspec_filt] = zp[zspec_filt]
    ez[zspec_filt] = ezp[zspec_filt]
    z_type[zspec_filt] = 'photo'
    survey[zspec_filt] = photo_survey
    
    zdata = Table({zcol_out: np.round(z, zdp), ezcol_out: np.round(ez, zdp),
                   ztout: z_type, surveycol: survey})
    zdata = hstack([spec_comb[[namecol]], zdata])
    
    return zdata


def add_col_descriptions(data, info_dict):
    'add in descriptions to columns from dictionary'
    datacols = data.colnames
    for col in list(info_dict.keys()):
        if col in datacols:
            data[col].description = info_dict[col]
        
    return


def parse_args():
    "parse input args, i.e. target and config file names"
    parser = argparse.ArgumentParser(description="find redshifts for input list of sources")
    parser.add_argument("sources",
                        help="table of radio sources with host ids")
    parser.add_argument("hosts",
                        help="table of host id positions")
    parser.add_argument("--namecol", action='store', type=str, default='AllWISE',
                        help="name of id col in hosts")
    parser.add_argument("--racol", action='store', type=str, default='RAJ2000',
                        help="name of RA col in hosts")
    parser.add_argument("--decol", action='store', type=str, default='DEJ2000',
                        help="name of Decl. col in hosts")
    parser.add_argument("--outdir", action='store', type=str, default='.',
                        help="directory to write files to")
    parser.add_argument("--search_radius", action='store', type=str,
                        default='1arcsec',
                        help="search radius to use for finding redshifts")
    parser.add_argument("--update_sources", action='store', type=str,
                        default='True',
                        help="directory to write files to")
                        
    args = parser.parse_args()
    
    ##make args.search_radius a quantity
    args.search_radius = u.Quantity(args.search_radius)
    args.update_sources = strtobool(args.update_sources)
    
    return args


################################################################################
################################################################################
###main -- use fetch_z.py instead

#if __name__ == '__main__':
#    t0 = time.time()
#    args = parse_args()
#    targets = load_targets(source_file=args.sources, host_file=args.hosts,
#                           namecol=args.namecol, acol=args.racol, dcol=args.decol)
#    zspec = get_spec_zs(targets, racol=args.racol, decol=args.decol,
#                        namecol=args.namecol, search_rad=args.search_radius)
#    zphot = query_lsdr8_north_and_south(targets, namecol=args.namecol,
#                                        acol=args.racol, dcol=args.decol,
#                                        searchrad=args.search_radius)
#    zdata = combine_spec_photo_zs(speczs=zspec, photozs=zphot,
#                                  namecol=args.namecol)
#
#    ###write to file
#    dlist = os.listdir(args.outdir)
#    if 'supplementary_data' in dlist:
#        filename = '/'.join([args.outdir, 'supplementary_data',
#                             'host_redshifts.fits'])
#    else:
#        filename = '/'.join([args.outdir, 'host_redshifts.fits'])
#    zdata.write(filename, format='fits')
#
#    ###join with source data
#    if args.update_sources==True:
#        print('')
#        print('Updating source table to include redshifts')
#        sources = Table.read(args.sources)
#        skeys = sources[~sources[args.namecol].mask][['Name', args.namecol]]
#        szdat = join(skeys, zdata, keys=args.namecol, join_type='inner')
#        szdat = unique(szdat, 'Name')
#        szdat.remove_column(args.namecol)
#        szdat = join(sources, szdat, keys='Name', join_type='left')
#        zmask = np.isnan(szdat['z'])
#        if type(szdat['z'])==MaskedColumn:
#            zmask = zmask | szdat['z'].mask
#        for zcol in ['z', 'z_err', 'z_type', 'z_survey']:
#            if type(szdat[zcol])==MaskedColumn:
#                szdat[zcol].mask = zmask
#
#        ###write to file
#        zfilename = '/'.join([args.outdir, 'sources.fits'])
#        szdat.sort('RA')
#        szcols = [i for i in sourcecols_out if i in szdat.colnames]
#        szdat = szdat[szcols]
#        add_col_descriptions(data=szdat, info_dict=sourccols_desc)
#        szdat.write(zfilename, format='fits', overwrite=True)
#
#    ###print time taken -- for testing purposes
#    t_elapsed = np.round(time.time()-t0, 2)
#    print('')
#    print('Time elapsed = ' + str(t_elapsed) + 's')
    
    
