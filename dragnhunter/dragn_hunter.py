###wrapper script for pair_finding, select_dragns and lr_host_matching
from pair_finding import *
from select_dragns import *


######################################################
######################################################
###parameters (make args where appropriate)

######################################################
######################################################
###functions

def parse_args():
    'parse command line arguments'
    parser = argparse.ArgumentParser(description="find radio pairs in catalogue data")
    parser.add_argument("radio_cat",
                        help="radio component catalogue to find pairs in")
    parser.add_argument("--config", action='store',
                        type=str, default='config.txt',
                        help="config file")
    parser.add_argument("--search_radius", action='store',
                        type=str, default='30arcsec',
                        help="search radius to use for core and host candidates")
    parser.add_argument("--find_hosts", action='store', type=str, default='True',
                        help="set code to query AllWISE for potential host candidates")
    parser.add_argument("--outdir", action='store', type=str,
                        default='output_files',
                        help="directory to write files to")
    args = parser.parse_args()
    
    ##make args.writepairs a bool and search_radius a quantity
    args.search_radius = u.Quantity(args.search_radius)
    args.find_hosts = strtobool(args.find_hosts)
    
    return args

######################################################
######################################################
###main

if __name__ == '__main__':
    args = parse_args()
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
                                         find_hosts=args.find_hosts)
    ##write data to file
    data['dragns'].write('/'.join([args.outdir, 'dragns.fits']))
    data['single-comp'].write('/'.join([args.outdir, 'single_comps.fits']))
    if args.find_hosts==True:
        data['sources'].write('/'.join([args.outdir, 'sources.fits']))
        data['hosts'].write('/'.join([args.outdir, 'host_candidates.fits']))

