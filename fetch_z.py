###script to run in case hunt_dragns timesout at get_redshifts stage

from hunt_dragns_and_find_host import *
import time

######################################################
######################################################
###parameters


######################################################
######################################################
###functions


def parse_args():
    'parse command line arguments'
    parser = argparse.ArgumentParser(description="get redshifts for sources found by hunt_dragns")
    parser.add_argument("--sources", action='store',
                        type=str, default='sources.fits',
                        help="source file")
    parser.add_argument("--hosts", action='store',
                        type=str, default='supplementary_data/hosts.fits',
                        help="hosts file")
    parser.add_argument("--host_name", action='store', type=str, default='AllWISE',
                        help="name column for hosts")
    parser.add_argument("--host_ra", action='store', type=str, default='RAJ2000',
                        help="RA column for hosts")
    parser.add_argument("--host_dec", action='store', type=str, default='DEJ2000',
                        help="DEC column for hosts")
    parser.add_argument("--z_radius", action='store', type=str, default='1arcsec',
                        help="search radius for redshifts")
    parser.add_argument("--outdir", action='store', type=str,
                        default='output_files',
                        help="directory to write files to")
    parser.add_argument("--chunksize", action='store', type=int,
                        default=10000,
                        help="query chunk size")
    
    args = parser.parse_args()
    
    ###make strings quantities where needed
    args.z_radius = u.Quantity(args.z_radius)

    return args
    
    
def fetch_redshifts(sources, hosts,
                    sourcecols_out=None,
                    sourccols_desc=None,
                    host_name='AllWISE',
                    acol='RAJ2000', dcol='DEJ2000',
                    z_radius=1*u.arcmin,
                    outname='sources.fits',
                    chunk_size=10000,
                    sortcol='RA'):
    'get redshifts for sources'
    print('')
    print('Searching for available redshift measurements')
    
    targets = z_targets(sources=sources, hosts=hosts,
                        namecol=host_name, acol=acol,
                        dcol=dcol)
    zspec = get_spec_zs(targets, racol=acol, decol=dcol,
                        namecol=host_name, search_rad=z_radius,
                        chunksize=chunk_size)
    zphot = query_lsdr8_north_and_south(targets, namecol=host_name,
                                        acol=acol, dcol=dcol,
                                        searchrad=z_radius,
                                        chunk_size=chunk_size)
    zdata = combine_spec_photo_zs(speczs=zspec, photozs=zphot,
                                      namecol=host_name)

    ###join with source data
#    sources = main_output['host_ids']
    skeys = sources[~sources[host_name].mask][['Name', host_name]]
    szdat = join(skeys, zdata, keys=host_name, join_type='inner')
    szdat = unique(szdat, 'Name')
    szdat.remove_column(host_name)
    szdat = join(sources, szdat, keys='Name', join_type='left')
    zmask = np.isnan(szdat['z'])
    if type(szdat['z'])==MaskedColumn:
        zmask = zmask | szdat['z'].mask
    for zcol in ['z', 'z_err', 'z_type', 'z_survey']:
        if type(szdat[zcol])==MaskedColumn:
            szdat[zcol].mask = zmask

    ###write to file
    if sortcol not in szdat.colnames:
        sortcols = acol
    szdat.sort(sortcol)
    
    if sourcecols_out is not None:
        szcols = [i for i in sourcecols_out if i in szdat.colnames]
        szdat = szdat[szcols]
    else:
        print('WARNING: no list of output columns provided')
        
    if sourccols_desc is not None:
        add_col_descriptions(data=szdat, info_dict=sourccols_desc)
    else:
        print('WARNING: no column descriptions provided')
        
    szdat.write(outname, format='fits', overwrite=True)

    return

######################################################
######################################################
###main


if __name__ == '__main__':
    t0 = time.time()
    args = parse_args()
    source_file = '/'.join([args.outdir, args.sources])
    host_file = '/'.join([args.outdir, args.hosts])
    sources = Table.read(source_file)
    hosts = Table.read(host_file)
    fetch_redshifts(sources=sources, hosts=hosts,
                    sourcecols_out=sourcecols_out,
                    sourccols_desc=sourccols_desc,
                    host_name=args.host_name,
                    acol=args.host_ra,
                    dcol=args.host_dec,
                    outname=source_file)
    time_elapsed = np.round(time.time()-t0, 2)
    print('Time elapsed: ' + str(time_elapsed) + 's')
