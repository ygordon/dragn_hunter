###query vizier using TAP
import numpy as np
from astroquery.utils.tap.core import TapPlus
from astropy.table import Table
from astropy.io.votable.tree import VOTableFile
from astropy import units as u
from astropy.coordinates import SkyCoord
from io import BytesIO
import time


#################################################################
#################################################################
###parameters

test_dat = Table({'Name': ['VLASS1QLCIR J100926.89-015746.9',
                           'VLASS1QLCIR J071435.04-013628.0',
                           'VLASS1QLCIR J090114.17-040540.4',
                           'VLASS1QLCIR J071650.91-002057.8',
                           'VLASS1QLCIR J082418.14-003309.6'],
                  'RA': np.array([152.36206125, 108.64602113,
                                  135.30904595,109.21213064,
                                  126.07559758])*u.deg,
                  'DEC': np.array([-1.96305167, -1.60779641,
                                   -4.09457912, -0.34939501,
                                   -0.55268124])*u.deg})

#################################################################
#################################################################
###functions

def write_query(table_to_query='II/328/allwise',
                search_rad=30*u.arcsec, searchprec=0,
                upload_name='uploaded_data',
                upra='RA', updec='DEC',
                qra='RAJ2000', qdec='DEJ2000'):
    'construct SQL query string from input parameters'
    
    srad = str(np.round(search_rad.to('arcsec').value, searchprec))
    
    qstring = f"SELECT * \n FROM tap_upload.{upload_name} AS tup \n JOIN \"{table_to_query}\" AS db \n ON 1=CONTAINS(POINT('ICRS', tup.{upra}, tup.{updec}), CIRCLE('ICRS', db.{qra}, db.{qdec}, {srad}/3600.))"
    
    return qstring


def tapq_vizier(data, acol='RA', dcol='DEC',
                viztable='II/328/allwise',
                search_radius=30*u.arcsec,
                vizra='RAJ2000', vizdec='DEJ2000',
                upload_name='uploaded_data',
                sep_col='ang_sep', verbose=True):
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
    
    ###mark time for testing
    dt = time.time()-t0
    print('')
    print('time for host finding (N='+str(len(data))+') = ' + str(np.round(dt, 2)) +'s')
    print('')
    
    return outdata



    
