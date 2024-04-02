###script to run likelihood ratio on doubles with host candidates and no core
###is 30'' radius enough to determine background? may need to download more sources out to e.g. 1', 5' etc radius


import numpy as np, sys, os, shutil, argparse
#from scipy.stats import chisquare
#sys.path.append('/Users/yjangordon/Documents/science/science_code/astro_tools/')
#from astro_tools import *
from astropy.table import Table, join, hstack, vstack, unique
#from matplotlib.colors import LogNorm
#from matplotlib.ticker import MultipleLocator, FormatStrFormatter,AutoMinorLocator
from astropy import cosmology as cos, units as u
from astropy.coordinates import SkyCoord, match_coordinates_sky
from collections import Counter
from astropy.wcs import WCS
#from astropy.nddata import Cutout2D
from likelihood_ratio_matching import *

########################################
#time of code check - not important for code to run
import time
start_time = time.time()
########################################

################################################################################
################################################################################
### parameters
sizelim = 0.3


################################################################################
################################################################################
### functions

def quick_cleanup(keep_files=['likelihood_ratio_matching.py',
                              'matches_bin', '__pycache__',
                              'find_hosts.py']):
    dlist = os.listdir()
    for file in dlist:
        if file not in keep_files:
            os.remove(file)
    return


def trim_slash(string):
    'remove forward slash from end of string, e.g. allows functions here to work if given directory name of format dir/ rather than dir'
    if string[-1] == '/':
        string = string[:-1]
    
    return string
    
    
def subset_required_wise_images(catalog, wise_meta, acol='RA', dcol='DEC',
                                wra='ra', wdec='dec', wid='coadd_id'):
    'function to subset required list of WISE images for host finding'
    catpos = SkyCoord(ra=catalog[acol], dec=catalog[dcol])
    wpos = SkyCoord(ra=wise_meta[wra], dec=wise_meta[wdec])
    
    ###cross match positions
    xmatch = catpos.match_to_catalog_sky(wpos)
    
    ###unique list of indeces of nearest image to each catalog position
    need_ims_idx = np.unique(xmatch[0])
    
    ###Table of metadata for only required images
    need_ims = wise_meta[need_ims_idx]
    need_ims.sort(wra)
    
    return need_ims


def headinfo_from_meta(metadata, image_id, idcol='coadd_id', naxis='naxis',
                       naxis1='naxis1', naxis2='naxis2', ctype1='ctype1',
                       ctype2='ctype2', crval1='ra', crval2='dec',
                       crpix1='crpix1', crpix2='crpix2', cdelt1='cdelt1',
                       cdelt2='cdelt2', cunit1='deg', cunit2='deg'):
    'obtain fits header info from wise meta'
    
    image_meta = metadata[(metadata[idcol]==image_id)][0]
    
    wcs_dict = {'NAXIS': image_meta[naxis],
                'NAXIS1': image_meta[naxis1],
                'NAXIS2': image_meta[naxis2],
                'CTYPE1': image_meta[ctype1],
                'CTYPE2': image_meta[ctype2],
                'CRVAL1': image_meta[crval1],
                'CRVAL2': image_meta[crval2],
                'CRPIX1': image_meta[crpix1],
                'CRPIX2': image_meta[crpix2],
                'CDELT1': image_meta[cdelt1],
                'CDELT2': image_meta[cdelt2],
                'CUNIT1': cunit1,
                'CUNIT2': cunit2}
    
    return wcs_dict


def fakeim_from_header_info(header_dict):
    'make a fake (2D) image (all zero values) with WCS info from a dictionary with header info in'
    
    ###build primary hdu
    hdu = fits.PrimaryHDU(np.zeros(shape=(header_dict['NAXIS1'], header_dict['NAXIS2'])))
    hdu.header.update(header_dict)
    hdulist = fits.HDUList(hdu)

    return hdulist


def minmax_coords_of_image(header):
    'calculate the min/max ra/dec covered by an image'
    wcs = WCS(header)
    
    coords = wcs.calc_footprint().transpose()
    
    extreme_coords = {'ra_min': np.min(coords[0]),
                      'ra_max': np.max(coords[0]),
                      'de_min': np.min(coords[1]),
                      'de_max': np.max(coords[1])}
    
    return extreme_coords
    
    
def extract_matches_and_tidy(prefix, target_dir='lr_output', outdir='matches',
                             matchfile='W1_LR_matches.dat'):
    'get match info and move to clean folder then remove target folder'
    target_dir = trim_slash(target_dir)
    outdir = trim_slash(outdir)
    
    ###setup filename
    matchfile_in = '/'.join([target_dir, matchfile])
    matchfile_out = '/'.join([outdir, '_'.join([prefix, matchfile])])
    
    dlist = os.listdir()
    if outdir not in dlist:
        os.mkdir(outdir)
        
    ###move file
    if target_dir in dlist:
        t_list = os.listdir(target_dir)
        if matchfile in t_list:
            shutil.move(src=matchfile_in, dst=matchfile_out)
        ###remove superfluous LR output
        shutil.rmtree(target_dir)
    
    return
    

def concat_results(target_directory, extension='dat', format='ascii',
                   oldcolnames=['radio_ID', 'W1_ID', 'W1_LR',
                                'W1_Rel', 'W1_n_cont', 'W1_separation'],
                   newcolnames=['Name', 'AllWISE', 'LR', 'Rel', 'N_cont',
                                'Sep_AllWISE'], round=True,
                   outname='all_lr_matches.fits', outformat='fits',
                   file_identifier='_W1_LR_matches.dat',
                   lrcol='W1_LR', ridcol='radio_ID', oidcol='W1_ID',
                   ncontcol='W1_n_cont', relcol='W1_Rel', sepcol='W1_separation',
                   add_original_file_identifier=True,
                   sepunit=u.arcsec):
    'stack all tables in target directory'
    target_directory = trim_slash(target_directory)
    
    ###list available files
    file_list = os.listdir(target_directory)
    ###limit to those with correct exension
    file_list = [file for file in file_list if extension in file]
    
    ###open data and append tables to list
    table_list = []
    for file in file_list:
        fname = '/'.join([target_directory, file])
        data = Table.read(fname, format=format)
        if len(data)>0:
            if add_original_file_identifier==True:
                file_id = file.replace(file_identifier, '')
                data['image_id'] = file_id
            table_list.append(data)
        
    ####need to remove empty tables before stacking
    table_list = [t for t in table_list if len(t)>0]
        
    ###stack tables
    matchdata = vstack(table_list)
    
    ###remove cases where the same match is duplicated due to image edge effects (keep highest n_cont/lowest LR)
    matchdata.sort(lrcol)
    matchname = [matchdata[ridcol][i]+'---'+matchdata[oidcol][i] for i in range(len(matchdata))]
    matchdata['match_id'] = matchname
    matchdata = unique(matchdata, keys='match_id', keep='first')
    matchdata.remove_column('match_id')
    matchdata.sort(ridcol)
    
    ###round data
    if round==True:
        rounddict={lrcol: 6, relcol: 6, ncontcol: 6, sepcol: 3}
        for key in list(rounddict.keys()):
            matchdata[key] = np.round(matchdata[key], rounddict[key])
    
    ###add units to sep
    matchdata[sepcol].unit = sepunit
    
    ###rename columns
    matchdata.rename_columns(names=oldcolnames, new_names=newcolnames)
    
    ###write to file
    matchdata.write(outname, format=outformat)
    
    return


def prep_cats(radio, iropt, radname='Name', rada='RA', radd='DEC',
              radflux='Flux', radfluxerr='E_Flux', radsize='LAS',
              irname='AllWISE', ira='RAJ2000', ird='DEJ2000',
              irmags=['W1mag'], irmagerrs=['e_W1mag'], irsep='Sep_AllWISE',
              minsize=3*u.arcsec, fillval=np.nan):
    'prepare catalogs for likelihood ratio matching'
    ###radio needs name, ra, dec, flux, flux_err
    ###ir/opt needs name, ra, dec, mag(s), mag_err(s), stargal, mask
    
    ###iropt need size info for masking
    iroptdata = join(iropt, radio[radname, radsize], keys=radname, join_type='left')
    
    ####need to remove those without s/n > 5 in AT LEAST 1 band, fill errors with nans
#    ((1.08574/hosts['e_W1mag'])>5)
    band_sn = []
    for err in irmagerrs:
        iroptdata[err].fill_value = fillval
        iroptdata[err] = iroptdata[err].filled()
        sn = np.array(1.08574/iroptdata[err])
        snflag = np.zeros(len(iroptdata)).astype(int)
        snflag[sn>5] = 1
        band_sn.append(snflag)
    band_sn = np.array(band_sn)
    sn_filter = np.max(band_sn.transpose(), axis=1).astype(bool)
    iroptdata = iroptdata[sn_filter]
    
    ###add in mask and stargal cols
    maskcol = np.ones(len(iroptdata)).astype(int)
    maskfilt = ((iroptdata[irsep]<iroptdata[radsize]*sizelim) | (iroptdata[irsep]<minsize))
    maskcol[maskfilt] = 0
    
    iroptdata['Mask'] = maskcol
    iroptdata['StarGal'] = 1
    
    ircols = [irname, ira, ird]
    for i in range(len(irmags)):
        ircols.append(irmags[i])
        ircols.append(irmagerrs[i])
    ircols = ircols + ['Mask', 'StarGal']
    
    icatout = iroptdata[ircols]
    ###if same candidate is within 30'' of two sources it may be duplicated in the iropt catalog, remove duplicates as this screws with LR code.
    icatout = unique(icatout, irname)

    ###radio just needs subset of columns -- rename radname to 'Source_id'
    rcatout = radio[[radname, rada, radd, radflux, radfluxerr]]
    if radname!= 'Source_id':
        rcatout.rename_column(name=radname, new_name='Source_id')
    
    lrcats = {'radio': rcatout, 'mw': icatout}
    
    return lrcats


def setup_config(beam_size=3, outdir='.', bands='W1',
                 id_col='AllWISE', ra_col='RAJ2000', dec_col='DEJ2000',
                 mask_col='Mask', sg_col='StarGal', mag_col='Xmag',
                 mag_err_col='e_Xmag', flux_col='Flux',
                 flux_err_col='E_Flux', cal_errors=7.5,
                 write=False, outname='lr_config.txt'):
    'setup config file for likelihood ratio'
    
    ###dictionary of config params'
    config_dict = {'parameter': 'value',
                   'outdir': outdir,
                   'bands': bands,
                   'id_col': id_col,
                   'ra_col': ra_col,
                   'dec_col': dec_col,
                   'mask_col': mask_col,
                   'sg_col': sg_col,
                   'mag_col': mag_col,
                   'mag_err_col': mag_err_col,
                   'flux_col': flux_col,
                   'flux_err_col': flux_err_col,
                   'beam_size': str(beam_size),
                   'cal_errors': str(cal_errors)}
    
    if write==True:
        line_list = []
        for k, v in config_dict.items():
            line_list.append('\t'.join([k, v]))
        outstring = '\n'.join(line_list)
        outfile = open(outname, 'w')
        outfile.write(outstring)
        outfile.close()
        return
    else:
        return config_dict


def run_lr(rfile, mwfile, imfile, region_name, config='lr_config.txt',
           outdir='lr_bin', matchfile='W1_LR_matches.dat'):
    'create inputs, run LR host finding, and tidy up for single image'
    
    ###list directory contents -- need for cleaning up later
    dlist = os.listdir()
    
    ###setup dictionary for LR arguments
    arg_dict = {'multiwave_cat': mwfile, 'radio_cat': rfile,
                'mask_image': imfile,'config_file': config,
                'overwrite': False, 'snr_cut': 5.0,
                'LR_thresh': 0.8}
    
    ###run LR host_finding
    lr_main(multiwave_cat=arg_dict['multiwave_cat'], radio_cat=arg_dict['radio_cat'],
            mask_image=arg_dict['mask_image'], config_file=arg_dict['config_file'],
            overwrite=arg_dict['overwrite'], snr_cut=arg_dict['snr_cut'],
            LR_threshold=arg_dict['LR_thresh'])
    
    ###tidy up output files
    cleandir(dlist)
    extract_matches_and_tidy(prefix=region_name, outdir=outdir,
                             matchfile=matchfile) ###these two functions have some redundancy; cleandir is generic from LR code, extract tidies further
    
    ###remove created input files
    os.remove(mwfile)
    os.remove(rfile)
    os.remove(imfile)
    os.remove(config)
    os.remove('/'.join([outdir, 'radio_masked.fits']))
    
    return


def does_file_exist(filename, directory):
    'check if file exists'
    ##ensure filename doesn't contain path
    filename = filename.split('/')[-1]
    
    ###list directory contents
    dlist = os.listdir(directory)
    
    ##is file in there?
    in_dir = filename in dlist
    
    return in_dir


def image_has_lr_results(image_id, directory,
                         resfile_contains='LR_matches'):
    'check if an image has LR results in a directory'
    ###list directory contents
    dlist = os.listdir(directory)
    
    ###subset with image_id and resfile_contains in filename
    image_results = [f for f in dlist if str(image_id) in f and resfile_contains in f]
    
    image_has_results = len(image_results)>0
    
    return image_has_results
    


###under construction -- needs to ID required images then loop through subsetting required data for host finding for each image, and concat results
def host_matching(source_cat, host_cat, image_meta,
                  sra='RA', sdec='DEC', sid='Name', sflux='Flux',
                  sfluxerr='E_Flux', ssize='LAS',
                  iid='coadd_id', ira='ra', idec='dec',
                  hid='AllWISE', hra='RAJ2000', hdec='DEJ2000',
                  hmags=['W1mag'], hmagerrs=['e_W1mag'],
                  sepcol='Sep_AllWISE',
                  bin_dir='lr_bin', outdir='output_files',
                  assumed_psf=6.1*u.arcsec, verbose=True,
                  radio_beam_size=3, cal_errors=7.5, maxraspread=10):
    'iterate through images performing LR host finding'
    
    t0 = time.time() ##function start time
    ##0)remove meta for clean joins, ensure bindir and ourdir exist
    source_cat.meta={}
    host_cat.meta={}
    dlist = os.listdir()
    if bin_dir != '.' and bin_dir not in dlist:
        os.mkdir(bin_dir)
    if outdir != '.' and outdir not in dlist:
        os.mkdir(outdir)
    
    ##1) list required images
    use_images = subset_required_wise_images(catalog=source_cat, wise_meta=image_meta,
                                             acol=sra, dcol=sdec, wra=ira, wdec=idec,
                                             wid=iid)

    ##2) loop though
#    use_images = use_images[use_images[iid]=='1364m379_ac51'] #test bad image
    N_iters = len(use_images)
#    N_iters = 3 ##test small subset
    for i in range(N_iters):
        ##define image specific values to use
        im_id = use_images[iid][i]
        
        if verbose == True:
            im_n = i+1
            print('')
            print('Likelihood ratio host finding for ' + str(im_id) + ', image ' + str(im_n) + '/' + str(N_iters))
            print('')
        
        ###check if image results exist already
        already_done = image_has_lr_results(image_id=im_id, directory=bin_dir)
        
        if already_done == True:
            print('   LR_matches already exist, moving to next image')
            print('')
        else:
            ##fake image and subset catalogs appropriately
            header_info = headinfo_from_meta(metadata=use_images, image_id=im_id)
            fake_im = fakeim_from_header_info(header_info)
            coordlims = minmax_coords_of_image(header_info)
            
            ###filter input catalogs for position
            sfilter = ((source_cat[sra]>coordlims['ra_min'])
                       & (source_cat[sra]<coordlims['ra_max'])
                       & (source_cat[sdec]>coordlims['de_min'])
                       & (source_cat[sdec]<coordlims['de_max']))
            
            hfilter = ((host_cat[hra]>coordlims['ra_min'])
                       & (host_cat[hra]<coordlims['ra_max'])
                       & (host_cat[hdec]>coordlims['de_min'])
                       & (host_cat[hdec]<coordlims['de_max']))
            ###account for ra wrapping -- ra_min could be > ra_max at RA~0
            ra_spread = np.abs(coordlims['ra_max']-coordlims['ra_min'])
            
            if ra_spread > maxraspread:
                if coordlims['ra_min'] + 360 - coordlims['ra_max'] < maxraspread:
                    sfilter = (((source_cat[sra]>coordlims['ra_max'])
                                | (source_cat[sra]<coordlims['ra_min']))
                               & (source_cat[sdec]>coordlims['de_min'])
                               & (source_cat[sdec]<coordlims['de_max']))
            
                    hfilter = (((host_cat[hra]>coordlims['ra_max'])
                                | (host_cat[hra]<coordlims['ra_min']))
                               & (host_cat[hdec]>coordlims['de_min'])
                               & (host_cat[hdec]<coordlims['de_max']))
            
            s_data = source_cat[sfilter]
            h_data = host_cat[hfilter]
        
            ###need to prep catalog data for likelihood_ratio_matching.py
            if len(s_data)>0 and len(h_data)>0: ##ensures data in both sub catalogs
                catdata = prep_cats(radio=s_data, iropt=h_data, radname=sid,
                                    rada=sra, radd=sdec, radflux=sflux,
                                    radfluxerr=sfluxerr, radsize=ssize,
                                    irname=hid, ira=hra, ird=hdec,
                                    irmags=hmags, irmagerrs=hmagerrs,
                                    irsep=sepcol, minsize=assumed_psf)
                ###add in if statement here to abandon if all wise sources masked
                if len(catdata['mw'][(catdata['mw']['Mask']==0)])>0 and len(catdata['radio'])>0:
                
                    ###write to bin
                    rfile = '/'.join([bin_dir, 'radio.fits'])
                    hfile = '/'.join([bin_dir, 'wise.fits'])
                    ifile = '/'.join([bin_dir, 'maskim.fits'])
        
                    catdata['radio'].write(rfile)
                    catdata['mw'].write(hfile)
                    fake_im.writeto(ifile)
        
                    ###set up config file
                    band_list = [i.replace('mag', '') for i in hmags]
                    band_list = ', '.join(band_list)
                    setup_config(beam_size=radio_beam_size, outdir='.', bands='W1',
                                 id_col=hid, ra_col=hra, dec_col=hdec,
                                 mask_col='Mask', sg_col='StarGal', mag_col='Xmag',
                                 mag_err_col='e_Xmag', flux_col=sflux,
                                 flux_err_col=sfluxerr, cal_errors=cal_errors,
                                 write=True, outname='lr_config.txt')
        
                    ####run LR
                    run_lr(rfile=rfile, mwfile=hfile, imfile=ifile,
                           region_name=im_id, config='lr_config.txt',
                           outdir='lr_bin', matchfile='W1_LR_matches.dat')
                
                else:
                    print('No matches in image!')
        
            else:
                print('   No sources in image!')
                print('')
    
    ###stack results
    stacked_file = '/'.join([outdir, 'all_LR_matches.fits'])
    sf_exists = does_file_exist(filename=stacked_file, directory=outdir)
    vnumber = 1
    while sf_exists==True:
        stacked_file = 'all_LR_matches_v' + str(vnumber) + '.fits'
        stacked_file = '/'.join([outdir, stacked_file])
        vnumber = vnumber + 1
        sf_exists = does_file_exist(filename=stacked_file, directory=outdir)
    
    concat_results(target_directory=bin_dir, extension='dat', format='ascii',
                   oldcolnames=['radio_ID', 'W1_ID', 'W1_LR',
                                'W1_Rel', 'W1_n_cont', 'W1_separation'],
                   newcolnames=[sid, hid, 'LR', 'Rel', 'N_cont', 'Sep_AllWISE'],
                   round=True, outname=stacked_file,
                   outformat='fits')
    print('')
    print('Concatenated results saved to ' + stacked_file)
    
    if verbose == True:
        telapsed = np.round(time.time()-t0, 2)
        print('')
        print('Time to run = ' + str(telapsed))
        print('')
    
    return


def parse_args():
    'parse command line arguments'
    parser = argparse.ArgumentParser(description="find radio pairs in catalogue data")
    parser.add_argument("source_cat",
                        help="catalog of radio sources")
    parser.add_argument("host_cat",
                        help="catalog of candidate hosts")
    parser.add_argument("--image_meta", action='store',
                        type=str, default='AllWISE-image_atlas_metadata-lite.fits',
                        help="file containing image meta data for maskims")
    
    args = parser.parse_args()
    
    return args


################################################################################
################################################################################
### main

if __name__ == '__main__':
    args = parse_args()
    sources = Table.read(args.source_cat)
    hosts = Table.read(args.host_cat)
    image_info = Table.read(args.image_meta)
    
    host_matching(source_cat=sources, host_cat=hosts, image_meta=image_info,
                  sra='RA', sdec='DEC', sid='Name', sflux='Flux',
                  sfluxerr='E_Flux', ssize='LAS',
                  iid='coadd_id', ira='ra', idec='dec',
                  hid='AllWISE', hra='RAJ2000', hdec='DEJ2000',
                  hmags=['W1mag'], hmagerrs=['e_W1mag'],
                  sepcol='Sep_AllWISE',
                  bin_dir='lr_bin', outdir='output_files',
                  assumed_psf=6.1*u.arcsec, verbose=True,
                  radio_beam_size=3)
    ###make a config file for this column definitions and additional arguments in
    

###need to post-process stacked matches
### -- in some cases same match found twice and Rel is halved as a result (i.e. 0.49 instead of 0.98) -- duplicate allwise sources in input handled as separate by LR
### other cases same match found multiple times because of image overlaps (Rel is correct in these cases, just pick one with highest n_cont)


################################################################################
###testing

####refactor such that code finds out which WISE images are required then iterates through, doing host finding on an image by image basis
####subset rwquired catalog info for each image
###mask out host candidates that are more than LAS > 2 away from position, where LAS/2 > 3''(accounts for WISE and VLASS PSF).


###PSFs: 3'' VLASS, 6.1'' AWISE (W1), use 6.1'' for masking purposes
