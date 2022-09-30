###select dragns from pairs (using misalign/sep; output plot?) and tidy (rounded, output columns defined etc, add in LAS)
###input pair list and component list
###output two tables: DRAGNs, and single component sources (all components not used in dragns)
###host candidate finding (as an option): i.e. all AWISE sources within, e.g. 30''; may want to make up to 1' to account for largest (~2') doubles (i.e. largest sep/2).

import numpy as np, os, argparse
from astropy.table import Table, join, unique, vstack
from astropy import units as u
from astropy.coordinates import SkyCoord
from astroquery.vizier import Vizier
from astropy.io.votable import parse
from io import BytesIO
from distutils.util import strtobool

#plt.interactive(True)

###############################################################
###############################################################
###parameters

##dictionary of colnames to output in dragns table {'input_name': output_name'}
##currently hardcoded in for VLASS, make configurable in future
colname_dict = {'pair_name': {'name': 'Name',
                              'description': 'Julian name'},
                'RA_best': {'name': 'RA',
                            'description': 'Best R.A. (J2000) of DRAGN'},
                'DE_best': {'name': 'DEC',
                            'description': 'Best Decl. (J2000) of DRAGN'},
                'Total_flux_source': {'name': 'Flux',
                                      'description': 'Sum of the total flux of all associated components'},
                'E_Total_flux_source': {'name': 'E_Flux',
                                        'description': 'Uncertainty in Flux'},
                'CoreProm': {'name': 'Core_prom',
                             'description': 'Fraction of Flux from Core'},
                'E_CoreProm': {'name': 'E_Core_prom',
                               'description': 'Uncertainty in Core_prom'},
                'Tflux_ratio': {'name': 'Lobe_flux_ratio',
                                'description': 'Ratio of flux from Lobe_1 to flux from Lobe_2'},
                'E_Tflux_ratio': {'name': 'E_Lobe_flux_ratio',
                                  'description': 'Uncertainty in Lobe_flux_ratio'},
                'LAS': {'name': 'LAS',
                        'description': 'Largest Angular Extent'},
                'E_LAS': {'name': 'E_LAS',
                          'description': 'Uncertainty in LAS'},
                'abs_dPA_1': {'name': 'Misalign_1',
                              'description': 'Relative misalignment of Lobe_1'},
                'E_dPA_1': {'name': 'E_Misalign_1',
                            'description': 'Uncertainty in Misalign_1'},
                'abs_dPA_2': {'name': 'Misalign_2',
                              'description': 'Relative misalignment of Lobe_2'},
                'E_dPA_2': {'name': 'E_Misalign_2',
                            'description': 'Uncertainty in Misalign_2'},
                'mean_misalignment': {'name': 'Mean_misalign',
                                      'description': 'Mean misalignment of lobes'},
                'E_mean_misalignment': {'name': 'E_Mean_misalign',
                                        'description': 'Uncertainty in Mean_misalign'},
                'Component_name_1': {'name': 'Lobe_1',
                                     'description': 'Component name of Lobe_1'},
                'Component_name_2': {'name': 'Lobe_2',
                                     'description': 'Component_name of Lobe_2'},
                'candidate_core': {'name': 'Core',
                                   'description': 'Component_name of Core'},
                'RA_c': {'name': 'RA_core',
                         'description': 'R.A. (J2000) of Core'},
                'DEC_c': {'name': 'DEC_core',
                          'description': 'Decl. (J2000) of Core'},
                'medianRA': {'name': 'RA_median',
                             'description': 'Median R.A,. (J2000) of lobes'},
                'medianDEC': {'name': 'DEC_median',
                              'description': 'Median Decl. (J2000) of lobes'},
                'cenRA': {'name': 'RA_fw',
                          'description': 'Flux weighted central R.A. (J2000) of lobes'},
                'cenDEC': {'name': 'DEC_fw',
                           'description': 'Flux weighted central Decl. (J2000) of Lobes'}}
                



###############################################################
###############################################################
###functions

def select_dragns(data,
                  sepcol='Sep_pair',
                  misacol='mean_misalignment',
                  prefpaircol='pref_pair',
                  minsep=10,
                  sepmis_a=-96.01351351,
                  sepmis_b=225.32323928):
    'select DRAGNs from pairs'
    ###set up in first instance to use relation determined for paper, upgrade later to determine this automatically
    
    ###subset preferred pairs
    dragns = data[data[prefpaircol]==1]
    
    ##subset minumum size
    dragns = data[data[sepcol]>minsep]
    
    ##subset on sep v misalign -- v1 coefficients are input, but later versions may upgrade to determine this automatically
    dragns = dragns[(dragns[misacol]<(sepmis_a*np.log10(dragns[sepcol])+sepmis_b))]
    
    return dragns


def proxy_LAS(ddat, racol='RA', decol='DEC', pname='pair_name',
              pacol='DC_PA', majcol='DC_Maj',
              esep_col='E_Sep_pair', emajcol='E_DC_Maj',
              lasunit=u.arcsec, lascolname='LAS'):
    'join doubles and components to estimate LAS using component geometry'
    
    ###alter position and size cols to '_1/2'
    a1 = '_'.join([racol, '1'])
    d1 = '_'.join([decol, '1'])
    p1 = '_'.join([pacol, '1'])
    m1 = '_'.join([majcol, '1'])
    a2 = '_'.join([racol, '2'])
    d2 = '_'.join([decol, '2'])
    p2 = '_'.join([pacol, '2'])
    m2 = '_'.join([majcol, '2'])
    em1 = '_'.join([emajcol, '1'])
    em2 = '_'.join([emajcol, '2'])
    
    ###determine offset positions - do pos and neg offset for each component
    ###need to make angle obj not columns for sep and pa
    pa1 = np.array(ddat[p1])*ddat[p1].unit
    smaj1 = np.array(ddat[m1]/2)*ddat[m1].unit
    pa2 = np.array(ddat[p2])*ddat[p2].unit
    smaj2 = np.array(ddat[m2]/2)*ddat[m2].unit
    
    pos1 = SkyCoord(ra=ddat[a1], dec=ddat[d1])
    pos2 = SkyCoord(ra=ddat[a2], dec=ddat[d2])
    
    off1_pve = pos1.directional_offset_by(position_angle=pa1,
                                          separation=smaj1)
    off1_nve = pos1.directional_offset_by(position_angle=pa1,
                                          separation=-smaj1)
    off2_pve = pos2.directional_offset_by(position_angle=pa2,
                                          separation=smaj2)
    off2_nve = pos2.directional_offset_by(position_angle=pa2,
                                          separation=-smaj2)
                                          
    sp1_p2 = off1_pve.separation(off2_pve)
    sp1_n2 = off1_pve.separation(off2_nve)
    sn1_p2 = off1_nve.separation(off2_pve)
    sn1_n2 = off1_nve.separation(off2_nve)
    
    ###proxy LAS is max of (sp1_p2, sp1_n2, sn1_p1, sn1_n2) for each pair
    cplas_a = np.array([sp1_p2, sp1_n2, sn1_p2, sn1_n2])
    las_proxy = np.round(np.max(cplas_a*sp1_p2.unit, axis=0).to(lasunit), 2)
    
    eunit = ddat[esep_col].unit
    elas = np.round(np.sqrt(ddat[esep_col]**2 + ddat[em1]**2 + ddat[em2]**2), 2)

    ddat[lascolname] = las_proxy
    ddat['E_'+lascolname] = elas
    ddat['E_'+lascolname].unit = eunit
    
    return ddat


def tidy_dragns(data, col_info=colname_dict, namekey='name',
                desckey='description'):
    'reduce output of dragns table to useful columns'
    
    if col_info is not None:
        oldnames = list(col_info.keys())
        newnames = [col_info[i][namekey] for i in oldnames]
        descriptions = [col_info[i][desckey] for i in oldnames]
        data = data[oldnames]
        for i in range(len(oldnames)):
            oname = oldnames[i]
            nname = newnames[i]
            descript = descriptions[i]
            data.rename_column(name=oname, new_name=nname)
            data[nname].info.description = descript
    
    data.meta = {}
    
    return data


def find_unused(dragns, components,
                l1col='Lobe_1', l2col='Lobe_2',
                corecol='Core', compcol='Component_name',
                flux_col='Peak_flux', minflux=3):
    'find and return a list of components not used in the DRAGNs selection (not Lobe_1, Lobe_2, or Core)'
    
    ###subset only components brighter than minflux
    cdata = components[components[flux_col]>minflux]
    
    ###list all components in dragns
    comps_used = list(dragns[l1col]) + list(dragns[l2col]) + list(dragns[corecol][~dragns[corecol].mask])
    comps_used.sort()
    comps_used = Table({compcol: comps_used})
    comps_used['used'] = 1 ##dummy column to be masked when joining with parent sample
    
    ###join with component data
    unused = join(cdata, comps_used, keys=compcol, join_type='left')
    unused = unused[unused['used'].mask]
    unused.remove_column('used') ##output is now same format as component table listing those components above minflux not used in dragns selection
    
    return unused
    
    
def make_source_list(dragns, single_comps,
                     dcols=['Name', 'RA', 'DEC', 'Flux',
                            'E_Flux', 'LAS', 'E_LAS'],
                     scols=['Component_name', 'RA', 'DEC',
                            'Total_flux', 'E_Total_flux',
                            'DC_Maj', 'E_DC_Maj'],
                     typecol='Type', sortcol='RA'):
    'make a single list of sources (dragns and unused single_components) for host finding'
    
    ##subset to required columns, rename single_comp to be same as dragns
    ddat = dragns[dcols]
    sdat = single_comps[scols]
    sdat.rename_columns(names=scols, new_names=dcols)
    
    ###add in type (S=single-component, D=DRAGN)
    ddat[typecol] = 'D'
    sdat[typecol] = 'S'
    
    ##stack and sort
    sources = vstack([ddat, sdat])
    sources.sort(keys=sortcol)
    
    return sources


def vobytes_to_table(data):
    'reads in votable in bytes (i.e. direct from query not xml file) and outputs astropy table'
    bytes = BytesIO(data) ##make a BytesIO object of the data (input as bytes)
    parsed_bytes = parse(bytes) ##parse into astropy VOTableFile object
    outdata = parsed_bytes.get_first_table().to_table() ##returns first table in VOTableFile; wont work for multiple
    
    return outdata


def find_allwise(data, acol='RA', dcol='DEC', namecol='Name',
                 cwcols=['_q', 'AllWISE',
                         'W1mag', 'e_W1mag', 'W2mag', 'e_W2mag',
                         'W3mag', 'e_W3mag', 'W4mag', 'e_W4mag'],
                 sepcolname='Sep_AllWISE',
                 searchrad=30*u.arcsec,
                 q_async=True,
                 maxtimeout=600):
    'query AllWISE via VizieR for host candidates using pairs data'
    ###async query waits for job to finish before sending response, no benefit here -- maybe better with SQL query to IRSA?
    ###creat skycoord list
    poscat = SkyCoord(ra=data[acol], dec=data[dcol])

    viz = Vizier(catalog='II/328', columns=['_q', '_r', '*']) ##setup vizier query and run
    viz.ROW_LIMIT = -1
    viz.TIMEOUT = maxtimeout
    
    if q_async==True:
        qres = viz.query_region_async(poscat, radius=searchrad)
        awise = vobytes_to_table(qres.content)
    else:
        qres = viz.query_region(poscat, radius=searchrad) ##returns table list
        awise = qres[0] ##returns astropy table

    ###join with data to create two outputs: all candidates and best match
    ###1) table of all returned candidates

    ##add index starting at 1 to data so as to join with _q from query
    data['_q'] = np.arange(len(data))+1

    ###create table of hosts
    hostcands = join(data[[namecol, acol, dcol, '_q']], awise,
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
    hostcands.sort([namecol, sepcolname])

    return dragns, hostcands


def select_sources_and_find_hosts(pairs, components,
                                  minflux=3, pair_sep_col='Sep_pair',
                                  misalignment_col='mean_misalignment',
                                  prefpaircol='pref_pair', minsep=10,
                                  sepmis_a=-96.01351351, sepmis_b=225.32323928,
                                  p_ra='RA', p_dec='DEC', p_name='pair_name',
                                  cname='Component_name', cpeak='Peak_flux',
                                  cflux='Total_flux', cfluxerr='E_Total_flux',
                                  c_ra='RA', c_dec='DEC',
                                  cpa='DC_PA', cmaj='DC_Maj', esep='E_Sep_pair',
                                  emaj='E_DC_Maj', sizeunit=u.arcsec,
                                  find_hosts=True, search_rad=30*u.arcsec):
    'select DRAGNs and single component sources, find hosts if requested'
#    ###load data -- moved outside of function
#    pairs = Table.read(pairs_file)
#    components = Table.read(component_file)
    
    ###select dragns, add LAS and tidy
    dragns = select_dragns(data=pairs, sepcol=pair_sep_col,
                           misacol=misalignment_col,
                           prefpaircol=prefpaircol, minsep=minsep,
                           sepmis_a=sepmis_a,sepmis_b=sepmis_b)
    dragns = proxy_LAS(ddat=dragns, racol=p_ra, decol=p_dec, pname=p_name,
                       pacol=cpa, majcol=cmaj, esep_col=esep,
                       emajcol=emaj, lasunit=sizeunit)
    dragns = tidy_dragns(dragns)
    
    ###single component sources
    single_comps = find_unused(dragns=dragns, components=components,
                               compcol=cname, flux_col=cpeak,
                               minflux=minflux)
    
    ###make source list
    dragncols = ['Name', p_ra, p_dec, 'Flux', 'E_Flux', 'LAS', 'E_LAS']
    sccols = [cname, c_ra, c_dec, cflux, cfluxerr, cmaj, emaj]
    source_list = make_source_list(dragns=dragns, single_comps=single_comps,
                                   dcols=['Name', p_ra, p_dec, 'Flux',
                                          'E_Flux', 'LAS', 'E_LAS'],
                                   scols=[cname, c_ra, c_dec, cflux,
                                          cfluxerr, cmaj, emaj],
                                   sortcol=p_ra)
    
    outdata = {'dragns': dragns, 'single-comp': single_comps}
    
    ###add in host finding
    if find_hosts==True:
        host_data = find_allwise(source_list, acol=p_ra, dcol=p_dec,
                                 searchrad=search_rad)
        outdata['sources'] = host_data[0]
        outdata['hosts'] = host_data[1]
    
    return outdata



def parse_args():
    "parse input args, i.e. target and config file names"
    parser = argparse.ArgumentParser(description="select sources and find potential hosts")
    parser.add_argument("double_table",
                        help="catalog of radio pairs for selecting DRAGNs from")
    parser.add_argument("component_table",
                        help="radio component catalogue for identifyinf single component sources")
    parser.add_argument("--outdir", action='store', type=str, default='.',
                        help="directory to write files to")
    parser.add_argument("--min_flux", action='store', type=int, default=3,
                        help="minimum peak_flux to use in selecting single component sources")
    parser.add_argument("--search_radius", action='store', type=str, default='30arcsec',
                        help="search radius to use for core and host candidates")
    parser.add_argument("--find_hosts", action='store', type=str, default='True',
                        help="set code to query AllWISE for potential host candidates")
                                     
    args = parser.parse_args()
    
    ##make args.writepairs a bool and search_radius a quantity
    args.search_radius = u.Quantity(args.search_radius)
    args.find_hosts = strtobool(args.find_hosts)
    
    return args

###############################################################
###############################################################
###main

#if __name__ == '__main__':
#    args = parse_args()
#    pairs = Table.read(pairs_file)
#    components = Table.read(component_file)
#    data = select_sources_and_find_hosts(pairs_data=pairs,
#                                         component_data=components,
#                                         minflux=args.min_flux)
#    data['dragns'].write('/'.join([args.outdir, 'dragns.fits']))
#    data['single-comp'].write('/'.join([args.outdir, 'single_comps.fits']))
#    if args.find_hosts==True:
#        data['sources'].write('/'.join([args.outdir, 'sources.fits']))
#        data['hosts'].write('/'.join([args.outdir, 'host_candidates.fits']))
