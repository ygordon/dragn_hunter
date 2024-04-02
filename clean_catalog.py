###code to fix host positions in sources.fits
##to be run after the fact by pointing to WISE data; long term solution requires fix to dragnhunter
##also finds any cases where host has been assigned to multiple radio sources (edge of AWISE images; LR matching done image by image then stacked)
###any lingering shared component doubles also cleaned up here
##add units to host position

import numpy as np, argparse
from datetime import datetime
from distutils.util import strtobool
from collections import Counter
from astropy.table import Table, join, unique, MaskedColumn, setdiff
from astropy.coordinates import SkyCoord

###run as
## >python clean_catalog.py output_files/sources.fits output_files/dragns.fits output_files/supplementary_data/hosts.fits

##replace with correct path to sources, dragns and hosts data files

################################################################################
################################################################################
###parameters -- args

################################################################################
################################################################################
###functions

def lighten_hosts(data, namecol='AllWISE',
                  keepcols=['AllWISE', 'RAJ2000',
                            'DEJ2000', 'W1mag',
                            'e_W1mag', 'W2mag',
                            'e_W2mag', 'W3mag',
                            'e_W3mag']):
    'select lite version of data -- no duplicates'
    udata = unique(data, namecol)
    if len(udata)<len(data):
        print('reducing host data')
        udata = udata[keepcols]
        udata.sort(namecol)
        return udata
    else:
        return data


def find_shared(ddata, namecol='Name', sortcol='LAS',
                l1col='Lobe_1', l2col='Lobe_2'):
    'find dragns that share components with another dragn and return sources to remove (keep smallest)'
    ##count components
    c1 = np.array(ddata[l1col])
    c2 = np.array(ddata[l2col])
    c12 = np.append(c1, c2)
    ccount = Counter(c12)
    ccount = Table({'Component_name': list(ccount.keys()),
                    'n_dragns': list(ccount.values())})
    ###identify dragns that share components
    sharedcomps = ccount[ccount['n_dragns']>1]
    mdragns, comps = [], []
    for i in range(len(sharedcomps)):
        cname = sharedcomps['Component_name'][i]
        dset = ddata[(ddata['Lobe_1']==cname) | (ddata['Lobe_2']==cname)]
        for j in range(len(dset)):
            comps.append(cname)
            mdragns.append(dset[namecol][j])
    mdragns = Table({'Component_name': comps, namecol: mdragns})
    
    ###identify dragns to remove--keep smallest
    if len(mdragns)>0:
        mdragns = join(mdragns, ddata, keys=namecol, join_type='left')
        mdragns.sort(sortcol)
        bestdragns = unique(mdragns, 'Component_name', keep='first')
    else:
        for col in ddata.colnames:
            if col not in mdragns.colnames:
                mdragns.add_column([], name=col)
        bestdragns = mdragns
    
    ###need to list those sources to be removed
    removed = [name for name in mdragns[namecol] if name not in bestdragns[namecol]]
    removed = Table({namecol: removed})
    
    print('')
    print(f'{len(removed)} sources identified for removal after checking component shares')
    print('')
    
    return removed


def remove_bad(data, bad, namecol='Name', dataname='data'):
    'remove sources identified as bad from table'
    datadiff = setdiff(data[[namecol]], bad)
    newdata = join(data, datadiff, keys=namecol, join_type='inner')
    
    n_removed = len(data) - len(newdata)
    print('')
    print(f'{n_removed} rows removed from {dataname}')
    print(f'old len({dataname}) = {len(data)}')
    print(f'new len({dataname}) = {len(newdata)}')
    print('')
    
    return newdata


def test_positions(sdata, hdata, joincol='AllWISE',
                   scols=['Name', 'AllWISE', 'RA_AllWISE',
                          'DE_AllWISE'],
                   hcols=['AllWISE', 'RAJ2000', 'DEJ2000'],
                   pa1='RA_AllWISE', pd1='DE_AllWISE',
                   pa2='RAJ2000', pd2='DEJ2000',
                   pu1='deg', pu2='deg',
                   sepunit='arcsec', dp=2, tolerance=0,
                   return_n=True):
    'test if host positions are the same in sources and wise data'
    shdat = sdata[~sdata[joincol].mask]
    test = join(shdat[scols], hdata[hcols], keys=joincol, join_type='left')
    p1 = SkyCoord(ra=test[pa1], dec=test[pd1], unit=pu1)
    p2 = SkyCoord(ra=test[pa2], dec=test[pd2], unit=pu2)
    seps = p1.separation(p2)
    test['sep'] = seps.to(sepunit)
    maxsep = np.round(np.max(test['sep']), dp)
    n_wrong = len(test[test['sep']>tolerance])
    
    print('')
    print(f'{n_wrong} hosts have reported positions > {tolerance} {sepunit} from real position')
    print(f'Maximum offset between real and reported host positions: {maxsep} {sepunit}')
    print('')
    
    if return_n == True:
        return n_wrong
    else:
        return


def fix_host_positions(sdata, hdata,
                       hostcol='AllWISE',
                       namecol='Name',
                       acol_s = 'RA_AllWISE',
                       dcol_s = 'DE_AllWISE',
                       acol_h = 'RAJ2000',
                       dcol_h = 'DEJ2000',
                       aunit='deg', dunit='deg',
                       adesc='Host R.A.', ddesc='Host Decl.'):
    'bring in correct host positions and update'
    ###set up necessary info for joins, updates
    hcols = [hostcol, acol_h, dcol_h]
    shjoincols = [namecol, acol_h, dcol_h]
    ra_desc = sdata[acol_s].description
    ra_unit = sdata[acol_s].unit
    de_desc = sdata[dcol_s].description
    de_unit = sdata[dcol_s].unit
    if ra_desc is None:
        ra_desc = adesc
    if de_desc is None:
        de_desc = ddesc
    if ra_unit is None:
        ra_unit = aunit
    if de_unit is None:
        de_unit = dunit
    
    ###join sources with hosts with host data, then join result with sdata
    shosts = sdata[~sdata[hostcol].mask]
    snew = join(shosts, hdata[hcols], keys=hostcol, join_type='left')
    snew = join(sdata, snew[shjoincols], keys=namecol, join_type='left')
    
    ###replace cold poscols
    snew[acol_s] = snew[acol_h]
    snew[dcol_s] = snew[dcol_h]
    snew.remove_columns(names=[acol_h, dcol_h])
    
    ###update meta
    snew.meta = {}
    snew[acol_s].unit = ra_unit
    snew[acol_s].description = ra_desc
    snew[dcol_s].unit = de_unit
    snew[dcol_s].description = de_desc
    
    return snew


def count_multiradio_hosts(sdata, hostcol='AllWISE'):
    'count hosts assigned to multiple radio sources'
    shosts = sdata[~sdata[hostcol].mask]
    hcount = Counter(shosts[hostcol])
    hcount = Table({hostcol: list(hcount.keys()),
                    'n_radio': list(hcount.values())})
                    
    ###number of multi-source hosts
    n_multi = len(hcount[hcount['n_radio']>1])
    
    print('')
    print(f'Number of hosts assigned to multiple radio sources: {n_multi}')
    print('')
    
    return n_multi


def fix_multisource_hosts(sdata, namecol='Name',
                          hostcol='AllWISE',
                          lrcol='LR', racol='RA',
                          zcol='z', ztcol='z_type',
                          zscol='z_survey',
                          hcols=['AllWISE', 'RA_AllWISE',
                                 'DE_AllWISE', 'Sep_AllWISE',
                                 'LR', 'Rel', 'Host_flag',
                                 'W1mag', 'e_W1mag', 'W2mag',
                                 'e_W2mag', 'W3mag', 'e_W3mag',
                                 'W4mag', 'e_W4mag', 'z',
                                 'z_err', 'z_type', 'z_survey']):
    'select best host/source match for hosts assigned to multiple sources'
    
    shosts = sdata[~sdata[hostcol].mask]
    shosts.sort(lrcol)
    shosts.reverse()
    bhosts = unique(shosts, hostcol, keep='first')
    bhosts.sort(racol)
    bhosts = bhosts[[namecol]+hcols]
    ###update z masks so consistent
    if zcol in hcols:
        zmask = bhosts[zcol].mask
        if ztcol in hcols:
            bhosts[ztcol].mask = zmask
        if zscol in hcols:
            bhosts[zscol].mask = zmask
    
    sdata.remove_columns(hcols) ##remove columns that will be readded through join
    sdata = join(sdata, bhosts, keys=namecol, join_type='left')
    
    return sdata


def update_dragns(ddata, sdata, namecol='Name',
                  hcols=['AllWISE', 'Sep_AllWISE',
                         'LR', 'Rel', 'Host_flag']):
    'update host info in dragns table from sources table'
    ddata.remove_columns(hcols)
    ddata = join(ddata, sdata[[namecol]+hcols], keys=namecol,
                 join_type='left')
    
    return ddata


def date_filename(fname, suffix='backup'):
    'make a backup filename including date'
    
    fnsplit = fname.rsplit('.', 1)
    
    now = datetime.now()
    datestr = f'{now.day}-{now.month}-{now.year}'
    timestr = f'{now.hour}h{now.minute}m{now.second}s'
    
    newname = '_'.join([fnsplit[0], suffix, datestr, timestr])
    newname = '.'.join([newname, fnsplit[-1]])
    
    return newname


def cleanup(dragns, sources, hosts,
            backup_data=True, hname='AllWISE', hmbackup='Sep_AllWISE',
            dfilename='output_files/dragns.fits',
            sfilename='output_files/sources.fits'):
    'function to do all the cleaning'
    
    ###make sure host name is a masked column in sources/dragns
    if type(sources[hname])!=MaskedColumn:
        print(f"sources['{hname}'] is not a Masked Column, taking mask from sources['{hmbackup}']")
        sources[hname] = MaskedColumn(sources[hname])
        sources[hname].mask = sources[hmbackup].mask
    if type(dragns[hname])!=MaskedColumn:
        print(f"dragns['{hname}'] is not a Masked Column, taking mask from dragns['{hmbackup}']")
        dragns[hname] = MaskedColumn(dragns[hname])
        dragns[hname].mask = dragns[hmbackup].mask
    
    ###lighten hosts
    hosts = lighten_hosts(data=hosts)
    
    ###check if data needs fixing
    ##find mult-dragn comps
    to_remove = find_shared(ddata=dragns)
    if len(to_remove)>0:
        sources = remove_bad(data=sources, bad=to_remove,
                             dataname='sources')
        dragns = remove_bad(data=dragns, bad=to_remove,
                             dataname='dragns')
    ##find wrong host pos
    n_wrongpos = test_positions(sdata=sources, hdata=hosts)
    ##find hosts assigned to multiple radio sources
    n_multi = count_multiradio_hosts(sdata=sources)
    if n_wrongpos>0 or n_multi>0 or len(to_remove)>0:
        ###backup
        if backup_data==True:
            print('writing back up files')
            print('')
            dragns.write(date_filename(fname=dfilename, suffix='backup'))
            sources.write(date_filename(fname=sfilename, suffix='backup'))
        ###fix data if needed
        if n_wrongpos>0:
            print('updating host coordinates')
            need_to_update_dragns = True
            sources = fix_host_positions(sdata=sources, hdata=hosts)
            n_wrong2 = test_positions(sdata=sources, hdata=hosts)
        if n_multi>0:
            print('removing duplicate hosts')
            need_to_update_dragns = True
            sources = fix_multisource_hosts(sdata=sources)
            nm2 = count_multiradio_hosts(sdata=sources)
        ##update dragns
        dragns = update_dragns(ddata=dragns, sdata=sources)
        ###sort by RA
        sources.sort('RA')
        dragns.sort('RA')
        ###make sure host name is a masked column in sources/dragns
        if len(sources[~sources[hname].mask])!=len(sources[~sources[hmbackup].mask]):
            print(f"applying mask from sources['{hmbackup}'] to sources['{hname}']")
            sources[hname].mask = sources[hmbackup].mask
        if len(dragns[~dragns[hname].mask])!=len(dragns[~dragns[hmbackup].mask]):
            print(f"applying mask from dragns['{hmbackup}'] to dragns['{hname}']")
            dragns[hname].mask = dragns[hmbackup].mask
        print('writing fixed data to file')
        sources.write(sfilename, overwrite=True)
        dragns.write(dfilename, overwrite=True)
    
    return


def parse_args():
    'parse command line arguments'
    parser = argparse.ArgumentParser(description="update hosts in data")
    parser.add_argument("sources",
                        help="source file")
    parser.add_argument("dragns",
                        help="dragn file")
    parser.add_argument("hosts",
                        help="source file")
    parser.add_argument("--backup", action='store', type=str,
                        default='True',
                        help="create backups of data before fixing")
    
    args = parser.parse_args()
    
    ###make strings quantities where needed
    args.backup = strtobool(args.backup)

    return args

################################################################################
################################################################################
###main ###make a function to import to hunt_dragns

if __name__ == '__main__':
    args = parse_args()
    ##load data
    dragns = Table.read(args.dragns)
    sources = Table.read(args.sources)
    hosts = Table.read(args.hosts)
    
    cleanup(dragns=dragns, sources=sources, hosts=hosts,
            backup_data=True, hname='AllWISE', hmbackup='Sep_AllWISE',
            dfilename='output_files_bl_clean/dragns.fits',
            sfilename='output_files_bl_clean/sources.fits')
            
