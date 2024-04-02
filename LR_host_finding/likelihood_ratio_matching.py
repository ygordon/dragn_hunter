import argparse
import numpy as np
import pandas
import os
from astropy.io import fits
from astropy.wcs import WCS
from astropy.table import Table, Column
from scipy.optimize import curve_fit
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt
import aplpy
import shutil

########################################
#time of code check - not important for code to run
import time
start_time = time.time()
########################################


def cenang( a1, d1, a2, d2 ):

        # convert to radians
        a1 = a1 * np.pi / 180.
        d1 = d1 * np.pi / 180.
        a2 = a2 * np.pi / 180.
        d2 = d2 * np.pi / 180.

        ## angular distance on a great circle
        # vincenty method
        numerator = np.sqrt( np.power( np.cos(d2)*np.sin(a1-a2), 2. ) + np.power( (np.cos(d1)*np.sin(d2)) - (np.sin(d1)*np.cos(d2)*np.cos(a1-a2)) , 2. )  )
        denominator = ( np.sin(d1)*np.sin(d2) ) + ( np.cos(d1)*np.cos(d2)*np.cos(a1-a2) )
        theta = np.arctan2(numerator,denominator)

        # convert back to deg
        theta = theta * 180. / np.pi
        return theta

def apply_mask( catalogue, mask_image, ra_col='RA', dec_col='DEC', overwrite=True ):

    ## open the catalogue
    mycat = fits.open( catalogue )
    data = mycat[1].data
    ## generate pairs of RA, Dec
    ra_dec_vals = list( zip( data[ra_col], data[dec_col] ) )
    
    ## open the mask WCS (this should exist in the header
    mymask_WCS = WCS( mask_image )
    
    ## get the pixel coordinates
    pix_coords = mymask_WCS.wcs_world2pix( ra_dec_vals, 0 )  ## the second argument is the origin (1 for FITS)

    ## get the pixel values
    mymask = fits.open( mask_image )
    mask_vals = mymask[0].data
    naxis1 = mymask[0].header['NAXIS1']
    naxis2 = mymask[0].header['NAXIS2']

    ## find the flags -- this can take a while if the catalogue is large
    flags = []
    for x, y in pix_coords:
        ## check that it is inside the area
        x_check = x >= 0 and x < naxis1
        y_check = y >= 0 and y < naxis2
        if x_check and y_check:
            flags.append( mask_vals[int(y),int(x)] )
        else:
            flags.append(1.0)

    ## make new data where unmasked
    unmasked_idx = np.where( np.array(flags) == 0.0 )
    new_data = data[unmasked_idx]

    ## write a new file
    new_header = mycat[1].header
    new_header['NAXIS2'] = len(unmasked_idx[0])
    new_table = fits.BinTableHDU( data=new_data, header=new_header )
    new_name = catalogue.replace('.fits', '_masked.fits' )
    new_table.writeto( new_name, overwrite=overwrite )
    
    return( new_name )

def get_unmasked_area( mask_image, overwrite=True ):

    if os.path.isfile( 'unmasked_area.dat' ) and not overwrite:
        ## get the area from the file
        with open( 'unmasked_area.dat', 'r' ) as f:
            lines = f.readlines()
        f.close()
        area_asec = np.float( lines[0] )
    else:
        ## open the mask image
        mymask = fits.open( mask_image )
        ## get the pixel scale
        xpix_asec = mymask[0].header['CDELT1'] * 60. * 60.
        ypix_asec = mymask[0].header['CDELT2'] * 60. * 60.
        pix_area_asec = np.abs(xpix_asec) * np.abs(ypix_asec) 
        ## and the number of unmasked pixels
        npix_unmasked = np.count_nonzero(mymask[0].data==0.0)
        ## find the unmasked area
        area_asec = pix_area_asec * npix_unmasked
        mymask.close()
        with open( 'unmasked_area.dat', 'w' ) as f:
            f.write(str(area_asec))
        f.close()
	    
    print( 'Area is %s (in arcsec^2)'%str(area_asec) )
    return( area_asec )

def make_master_cat( multiwave_cat, overwrite=False, outdir='.', my_bands='J,H,Ks', id_col='UID', ra_col='RA', dec_col='DEC', mask_col='Mask', sg_col='SGsep', mag_col='mag_X', mag_err_col='mag_err_X' ):

    ## check if output file exists
    outtmp = multiwave_cat.split('/')[-1]
    outcat = os.path.join( outdir, outtmp )
    outcat = outcat.replace('.fits', '_master.fits' )

    ## read in catalogue
    mw_hdul = fits.open( multiwave_cat )
    mw_dat = mw_hdul[1].data
    mw_hdr = mw_hdul[1].header
    ## convert to table
    mw_tab = Table( mw_dat )
    mw_hdul.close()

    if os.path.isfile( outcat ) and not overwrite:
        print( 'File already exists and overwrite is not set to true, exiting.' )
    else:
        # columns to keep
        basic_cols = [id_col, ra_col, dec_col, mask_col, sg_col]

        # magnitude columns
        keep_bands = [ mag_col.replace('X', my_band) for my_band in my_bands ]
        keep_bands_err = [ mag_err_col.replace('X', my_band) for my_band in my_bands ]
       
        ## make a smaller table with only the columns we want
        keep_cols = basic_cols + keep_bands + keep_bands_err
        new_tab = mw_tab[keep_cols]

        ## select non-halo objects
        print( ' ... removing objects where '+mask_col+' is 1' )
        non_halo_idx = np.where(new_tab[mask_col] == 0)[0]
        new_tab = new_tab[non_halo_idx]

        ## select only galaxies
        print( ' ... selecting galaxies ('+sg_col+' is 1)' )
        galaxy_idx = np.where(new_tab[sg_col] == 1)[0]
        new_tab = new_tab[galaxy_idx]
   
        ## add signal to noise columns
        snr_cols = [ keep_band + '_SNR' for keep_band in keep_bands ]
        for keep_band_err, snr_col in zip(keep_bands_err, snr_cols):
            new_tab[snr_col] = 1.08574 / new_tab[keep_band_err]

        ## write to file
        new_fits = fits.BinTableHDU( data=new_tab )
        new_fits.writeto( outcat, overwrite=overwrite )

        print( 'Master catalogue %s created.'%outcat )

    return( outcat )

def calculate_positional_errors( data, beam_size=5., cal_errors=0.1, flux_col='Total_flux', flux_err_col='E_Total_flux' ):

    ## positional errors using Ivison's method
    pos_noise = 0.6 * beam_size / ( data[flux_col] / data[flux_err_col] )
    ## add in calibration errors
    sigma_pos = np.sqrt( np.power( cal_errors, 2. ) + np.power( pos_noise, 2. ) )
    r_max = 5. * np.max( np.sqrt( np.power( sigma_pos, 2. ) ) )
    print( 'r_max is %f arcsec'%r_max )
    
    return( r_max, sigma_pos )

def make_mag_bins( magnitudes ):

    mag_bins = np.arange( np.floor(np.min(magnitudes)), np.ceil(np.max(magnitudes)), 0.4 )
    ## check that it spans the entire range
    
    if len(mag_bins)!=0: ###added to account for lack of bin range in poor data sampling, return statement added in main to account for len=0
        if np.ceil(np.max(magnitudes)) > np.max( mag_bins ):
            mag_bins = np.append( mag_bins, np.max(mag_bins)+0.4 )
        if np.floor(np.min(magnitudes)) < np.min( mag_bins ):
            mag_bins = np.insert( mag_bins, 0, np.min( mag_bins )-0.4 )
    return( mag_bins )

def make_match_magnitudes( band, band_data, radio_data, outfile='', ra_col='RA', dec_col='DEC', mag='', r_max=0.0, overwrite=True ):
    
    ## get number of radio sources
    n_radio_sources = radio_data.shape[0]

    ## convert r_max to deg
    r_max_deg = r_max / 60. / 60.

    ## only run if the outfile does not exist and overwrite is False
    if os.path.isfile( outfile ) and not overwrite:
        print( 'File exits, reading from disk.' )
        m_mags = Table.read( outfile, format='ascii' )
    else:
        print( 'Finding magnitudes of all sources within r_max.' )
        match_mags = []
        radio_ids = []
        for x in np.arange( n_radio_sources ):
            distances = cenang( radio_data['RA'][x], radio_data['DEC'][x], band_data[ra_col], band_data[dec_col] )
            candidate_idx = np.where( distances <= r_max_deg )[0]
            match_mags = match_mags + band_data[mag][candidate_idx].tolist()
            radio_ids = radio_ids + np.repeat( radio_data['Source_id'][x], len(candidate_idx) ).tolist()
        ## make a table and write the file
        m_mags = Table()
        m_mags['radio_id'] = Column( radio_ids )
        m_mags['matched_mag'] = Column( match_mags )
        m_mags.write( outfile, format='ascii', overwrite=True )

    return( m_mags )

def find_number_no_counterparts( radio_dat, band_dat, radii, ra_col='RA', dec_col='DEC' ):

    ## what is the number of sources with no possible counterparts
    n_counterparts = np.zeros( (radio_dat.shape[0], len(radii)) )
    for xx in np.arange( radio_dat.shape[0] ):
        ## calculate the distance from the source to all other sources, convert to asec
        distances = cenang( radio_dat['RA'][xx], radio_dat['DEC'][xx], band_dat[ra_col], band_dat[dec_col] ) * 60. * 60. 
        ## loop through radii to find number of counterparts
        for yy in np.arange( len(radii) ):
            n_counterparts[xx,yy] = len( np.where(distances <= radii[yy])[0] )

    ## find number of sources with no counterpart as function of radius
    n_blanks = np.zeros( len(radii) )
    for xx in np.arange( len(radii) ):
        n_blanks[xx] = np.count_nonzero( n_counterparts[:,xx] == 0 )

    return( n_blanks )

def q0_function( x, q0, sig ):
    return 1. - q0 * ( 1. - np.exp( -np.power( x, 2. ) / ( 2. * np.power( sig, 2. ) ) ) )

def random_points_on_a_sphere( npoints, ra_limits, dec_limits ):

    ## convert to radians first
    ra_limits = ra_limits * np.pi / 180.
    dec_limits = dec_limits * np.pi / 180.

    ## RA is uniform in azimuthal angle
    np.random.seed(30)
    random_ra = np.random.uniform( ra_limits[0], ra_limits[1], npoints )
    
    ### - YG, account for RA wrapping here! - need to account in RA spread too
    if ra_limits[1] - ra_limits[0] > np.pi:
        amin_new = ra_limits[1] - 2*np.pi
        random_ra = np.random.uniform(amin_new, ra_limits[0], npoints)
        random_ra[(random_ra<0)] = random_ra[(random_ra<0)] + 2*np.pi

    ## Dec is not
#    np.random.seed(20)
#    random_angles = np.random.uniform( 0, 1, npoints ) * ( dec_limits[1] - dec_limits[0] ) + dec_limits[0]
#    random_dec = np.arcsin( random_angles )####why does this line exist (polar clustering), breaks code?

    np.random.seed(20)
    random_angles = np.random.uniform( 0, 1, npoints )
    random_angles = np.arcsin(random_angles)
    random_dec = random_angles * ( dec_limits[1] - dec_limits[0] ) + dec_limits[0]



    ## convert back to degrees
    random_ra = random_ra * 180. / np.pi
    random_dec = random_dec * 180. / np.pi


    return random_ra, random_dec


def find_Q0_fleuren( band, radio_dat, band_dat, radii, mask_image, ra_col='RA', dec_col='DEC', overwrite=True ):


    ## first check if the number of no counterparts has already been found
    f_file = band + '_Fleuren_no_counterparts.dat'

    if os.path.isfile( f_file ) and not overwrite:
        print( 'Blanks already found for this band, reading file.' )
        t1 = Table.read( f_file, format='ascii' )
    else:
        ## first find the number of radio sources with no counterparts
        ##  (as a function of radius)
        REAL_no_counterparts = find_number_no_counterparts( radio_dat, band_dat, radii, ra_col=ra_col, dec_col=dec_col )

        ## check if a random radio catalogue already exists
        n_srcs = radio_dat.shape[0]
        f_rand_file = band + '_Fleuren_random_catalogue.fits'
        if os.path.isfile( f_rand_file ) and not overwrite:
            print( 'Random radio catalogue already exists, will use.' )
        else:
            print( "Random radio catalogue does not exist or overwrite=True." )


            ###YG edits - account for RA wrap at 360
            ## boundaries for random catalogue -- doesn't account for RA wrapping
            ## expand the area and apply the mask afterwards
            min_RA = np.min( radio_dat['RA'] )
            max_RA = np.max( radio_dat['RA'] )
            min_DEC = np.min( radio_dat['DEC'] )
            max_DEC = np.max( radio_dat['DEC'] )
            RA_spread = max_RA - min_RA
            
            ###YG - alter RA_spread for RA wrapping instances - move above padding
            if RA_spread > 180:
                a_array = np.array(radio_dat['RA'])
                new_max = max(a_array[(a_array<180)])
                new_min = min(a_array[(a_array>180)])
                RA_spread = new_max - new_min + 360

                max_RA = new_min
                min_RA = new_max
            
                min_RA = min_RA + RA_spread * 0.05
                max_RA = max_RA - RA_spread * 0.05
            else:
                min_RA = min_RA - RA_spread * 0.05
                max_RA = max_RA + RA_spread * 0.05
            
            DEC_spread = max_DEC - min_DEC
            min_DEC = min_DEC - DEC_spread * 0.05
            max_DEC = max_DEC + DEC_spread * 0.05


            ## randomly generate RA/DEC pairs - npoints to n_src*10
            rand_RA, rand_DEC = random_points_on_a_sphere( n_srcs*20, np.array([min_RA,max_RA]), np.array([min_DEC,max_DEC]) )
            ## create a table and write a fits catalogue (for apply_mask)
            t = Table()
            t['RA'] = rand_RA
            t['DEC'] = rand_DEC
            rand_table = fits.BinTableHDU( data=t )
            rand_table.writeto( f_rand_file, overwrite=overwrite )

        masked_rand_cat = apply_mask( f_rand_file, mask_image, overwrite=True )
        ## read in the masked radio catalogue 
        masked_hdul = fits.open( masked_rand_cat )
        masked_dat = masked_hdul[1].data
        if masked_dat.shape[0] >= n_srcs:
            ## trim it down to the right size
            masked_dat = masked_dat[np.arange(n_srcs)]
        else:
            print( 'Not enough data points in the catalogue, try increasing the initial number of random sources.' )
            return(1)

        ## find the number of random radio sources with no counterparts
        ##  (as a function of radius)
        RANDOM_no_counterparts = find_number_no_counterparts( masked_dat, band_dat, radii, ra_col=ra_col, dec_col=dec_col )

        ## take the ratio of real to random no counterparts
        no_counterpart_ratio = REAL_no_counterparts / RANDOM_no_counterparts

        ## YG: replace nan with 0
        no_counterpart_ratio[np.isnan(no_counterpart_ratio)] = 0

        ## make a table to write
        t1 = Table()
        t1['Radius'] = radii
        t1['Real'] = REAL_no_counterparts
        t1['Random'] = RANDOM_no_counterparts
        t1['Ratio'] = no_counterpart_ratio

        t1.write( f_file, format='ascii', overwrite=overwrite )


    ## fit the model to the data
    init_coeffs = np.array([1.3, 1.0])
    try: ##allows code to proceed if can't fit with default number of attempts (maxfev=600)
        coeff, covariance = curve_fit( q0_function, t1['Radius'], t1['Ratio'], init_coeffs )
    except:
        coeff, covariance = curve_fit( q0_function, t1['Radius'], t1['Ratio'], init_coeffs , maxfev=1200)
    coeff_err = np.sqrt(np.diag(covariance))
    print( '... q_0 = %s +/- %s'%(str(coeff[0]),str(coeff_err[0])) )

    ## output the result
    t2 = Table()
    t2['parameter'] = np.array(['Q0','sig'])
    t2['value'] = coeff
    t2['1sig'] = coeff_err
    t2.write( band+'_Q0_estimates.dat', format='ascii', overwrite=True )

    ## and plot
    n_srcs = radio_dat.shape[0]
    plt.plot( t1['Radius'], t1['Random']/n_srcs, linewidth=2, color='0.5', label='Random' )
    plt.plot( t1['Radius'], t1['Real']/n_srcs, linewidth=2, color='black', label='Real' )
    yvals = q0_function( t1['Radius'], coeff[0], coeff[1] )
    plt.plot( t1['Radius'], yvals, linewidth=2, color='red', label='Fit' )
    plt.plot( t1['Radius'], t1['Ratio'], marker='o', color='blue', label='Ratio' )
    plt.legend()
    plt.savefig( band + '_q0_estimate.png' )
    plt.close()

    return( coeff[0], coeff_err[0] )

def LR_and_reliability( band, band_dat, radio_dat, qm_nm, sigma_pos, mag_bins, r_max, q0, LR_threshold=0.8, ra_col='RA', dec_col='DEC', mag_col='', id_col='' ):

    ## housekeeping
    band_col = mag_col.replace('X', band)
    sig_sq = np.power( sigma_pos, 2. )

    ## initialize empty lists
    radio_id = []
    band_id = []
    lr_value = []
    lr_rel = []
    n_cont = []
    separation = []
    count = 0

    for xx in np.arange(radio_dat.shape[0]):

        ## calculate f(r)
        distances = cenang( radio_dat['RA'][xx], radio_dat['DEC'][xx], band_dat[ra_col], band_dat[dec_col] ) * 60. * 60.
        candidate_idx = np.where( distances <= r_max )[0]
        n_cand = len( candidate_idx )
        
        if n_cand > 0:
            ## select the data
            tmp_dat = band_dat[candidate_idx]
            ## get the magnitudes
            band_mags = tmp_dat[band_col]
            ## calculate the radial probability distribution for the candidates
            f_r = 1. / ( 2. * np.pi * sig_sq[xx] ) * np.exp( - np.power( distances[candidate_idx], 2. ) / ( 2. * sig_sq[xx] ) ) 
            ## loop through the candidates
            LR = []
            for yy in np.arange( len(candidate_idx) ):
                try:
                    qm_nm_tmp = qm_nm[ np.max( np.where( mag_bins <= band_mags[yy] )[0] ) ]
                except: ###accounts for very very rare case when only match has magnitude exactly at the faintest bin edge
                    m_nm_tmp = qm_nm[ np.max( np.where( mag_bins < band_mags[yy] )[0] ) ]
                tmpval = qm_nm_tmp*f_r[yy]
                LR.append(tmpval)
            LR = np.array(LR)

            ## calculate the reliability
            LR_reliability = LR / ( np.sum( LR ) + 1. - q0 )

            ## save information to lists            
            radio_id = radio_id + np.repeat( radio_dat['Source_id'][xx], n_cand ).tolist()
            band_id = band_id + tmp_dat[id_col].tolist()
            lr_value = lr_value + LR.tolist()
            lr_rel = lr_rel + LR_reliability.tolist()
            # contaminants ?
            lr_cont = np.sum( 1. - LR_reliability[ np.where( LR_reliability > LR_threshold )[0] ] ).tolist()
            n_cont = n_cont + np.repeat( lr_cont, n_cand ).tolist()
            separation = separation + distances[candidate_idx].tolist()
            count = count + n_cand

    ## write results for the band
    t = Table()
    t['radio_ID'] = radio_id
    band_name_id = band + '_ID'
    t[band_name_id] = band_id
    band_name_LR = band + '_LR'
    t[band_name_LR] = lr_value
    band_name_Rel = band + '_Rel'
    t[band_name_Rel] = lr_rel
    band_name_n_cont = band + '_n_cont'
    t[band_name_n_cont] = n_cont
    band_name_sep = band + '_separation'
    t[band_name_sep] = separation

    print( 'Found a total of %s candidates'%str(n_cand) )
    outfile = band + '_LR_matches.dat'
    print( 'Saving matches to %s'%outfile )
    t.write( outfile, format='ascii' )
    
    return( outfile )


def lr_main( multiwave_cat, radio_cat, mask_image, config_file='lr_config.txt', overwrite=False, snr_cut=5.0, LR_threshold=0.8 ):

    ## read the configuration file and parse
    print( 'Reading in configuration file.' )
    config_params = pandas.read_table( config_file, delim_whitespace=True ).replace("'","",regex=True)
    outdir = config_params['value'][np.where( config_params['parameter'] == 'outdir' )[0][0]]
    bands = config_params['value'][np.where( config_params['parameter'] == 'bands' )[0][0]]
    my_bands = bands.split(',')
    id_col = config_params['value'][np.where( config_params['parameter'] == 'id_col' )[0][0]]
    ra_col = config_params['value'][np.where( config_params['parameter'] == 'ra_col' )[0][0]]
    dec_col = config_params['value'][np.where( config_params['parameter'] == 'dec_col' )[0][0]]
    mask_col = config_params['value'][np.where( config_params['parameter'] == 'mask_col' )[0][0]]
    sg_col = config_params['value'][np.where( config_params['parameter'] == 'sg_col' )[0][0]]
    mag_col = config_params['value'][np.where( config_params['parameter'] == 'mag_col' )[0][0]]
    mag_err_col = config_params['value'][np.where( config_params['parameter'] == 'mag_err_col' )[0][0]]
    flux_col = config_params['value'][np.where( config_params['parameter'] == 'flux_col' )[0][0]]
    flux_err_col = config_params['value'][np.where( config_params['parameter'] == 'flux_err_col' )[0][0]]
    beam_size = np.float(config_params['value'][np.where( config_params['parameter'] == 'beam_size' )[0][0]])
    cal_errors = np.float(config_params['value'][np.where( config_params['parameter'] == 'cal_errors' )[0][0]])


    ## From the multiwavelength catalogue, make a master catalogue with some cuts and keeping only some of the catalogues
    print( 'Making a master catalogue.'  )
    master_cat = make_master_cat( multiwave_cat, overwrite=overwrite, outdir=outdir, my_bands=my_bands, id_col=id_col, ra_col=ra_col, 
            dec_col=dec_col, mask_col=mask_col, sg_col=sg_col, mag_col=mag_col, mag_err_col=mag_err_col )

    ## mask the master cat just to be sure it's right
    masked_master = apply_mask( master_cat, mask_image, ra_col=ra_col, dec_col=dec_col, overwrite=overwrite )
    master_cat = masked_master

    ## read in the master catalogue
    master_hdul = fits.open( master_cat )
    master_hdr = master_hdul[1].header
    master_dat = master_hdul[1].data
    
    ## get the un-masked area
    area_asec = get_unmasked_area( mask_image, overwrite=overwrite )

    ## mask the radio data
    masked_radio_cat = apply_mask( radio_cat, mask_image, ra_col='RA', dec_col='DEC', overwrite=overwrite )

    ## read in the masked radio data
    radio_hdul = fits.open( masked_radio_cat )
    radio_dat = radio_hdul[1].data
    ## count the radio sources
    n_radio_sources = radio_dat.shape[0]
    
    ####put in break (return) if either of the masked catalogs are empty? --YG: allows function to fail without crashing parent code calling it
    if len(master_dat)<1 or n_radio_sources<1:
        return
    
    ## find the positional error
    r_max, sigma_pos = calculate_positional_errors( radio_dat, beam_size=beam_size, cal_errors=cal_errors, flux_col=flux_col, flux_err_col=flux_err_col )
    ## where sigma_pos is less than 0.5 arcsec, make it 0.5 arcsec
    sigma_idx = np.where( sigma_pos < 0.5 )[0]
    sigma_pos[sigma_idx] = 0.5

    ## make a plot of the sky coverage
    ## first get a colormap
    cmap = matplotlib.cm.get_cmap('magma')
    im = aplpy.FITSFigure( mask_image )
    im.show_markers( master_dat[ra_col], master_dat[dec_col], marker=',', facecolor='0.4', edgecolor='0.4' )
    im.show_markers( radio_dat['RA'], radio_dat['DEC'], marker='*', facecolor=cmap(0.75), edgecolor=cmap(0.15), linewidth=0.3, s=70  )
    im.save('Sky_coverage.png')
    im.close()
    

    for my_band in my_bands:

        ## column names
        band_col = mag_col.replace( 'X', my_band )
        band_col_err = mag_err_col.replace( 'X', my_band )
        snr_col = band_col + '_SNR'

        ## select based on SNR
        snr_idx = np.where( master_dat[snr_col] >= snr_cut )
        band_dat = master_dat[snr_idx]

        ## get rid of negative things
        pos_idx = np.where( band_dat[band_col] > 0 )
        band_dat = band_dat[pos_idx]

        ## make magnitude bins
        mag_bins = make_mag_bins( band_dat[band_col] )
        if len(mag_bins)==0:
            print('')
            print('WARNING: magnitude bins could not be determined, check multiwavelength sample size')
            print('') ###accounts for code aborting if mag bins can't be estimated --YG
            return

        ## calculate the density of background sources: n(m)
        nmhist = np.histogram( band_dat[band_col], bins=mag_bins )
        nm_counts = nmhist[0]
        nm = nm_counts / area_asec

        ## find the matched magnitudes -- this takes a long time so it first checks if the file already exists
        mag_file = my_band + '_matched_mags_r' + str( round( r_max, ndigits=2 ) ) + '_SNR' + str( snr_cut ) + '.dat'
        m_mags = make_match_magnitudes( my_band, band_dat, radio_dat, outfile=mag_file, ra_col=ra_col, dec_col=dec_col, mag=band_col, r_max=r_max, overwrite=overwrite )
        match_magnitudes = np.array(m_mags['matched_mag'].tolist())
        match_radio_sources = np.array(m_mags['radio_id'].tolist())

        ## find q0 (Fleuren+ 2012 method)
        radii = np.arange( 1., beam_size, 0.2 )
        Q0, Q0_err = find_Q0_fleuren( my_band, radio_dat, band_dat, radii, mask_image, ra_col=ra_col, dec_col=dec_col, overwrite=overwrite )      
        
        ## find the expected distribution of true counterparts -- q(m)
        tmhist = np.histogram( match_magnitudes, bins=mag_bins )
        total_m = np.cumsum( tmhist[0] )
        background = nm * n_radio_sources * np.pi * np.power( r_max, 2. )
        real_m = total_m - background
        qm = real_m / np.sum( real_m ) * Q0

        ## find ratio of q(m)/n(m) -- the expected divided by the background
	## if there are values of zero then there will be issues
        qm_nm = [ qq/nn if not  qq == 0 and not nn == 0 else 0. for qq, nn in zip( qm, nm ) ]
        qm_nm = np.array( qm_nm )

        log_total = [ np.log10( tt ) if not tt==0 else 0. for tt in total_m  ]
        log_total = np.array( log_total )

        log_real = [ np.log10( rr ) if not rr<=0 else 0. for rr in real_m ]
        log_real = np.array( log_real )

        log_bkg = [ np.log10( bb ) if not bb <= 0 else 0. for bb in background ]
        log_bkg = np.array( log_bkg )

        log_qm_nm = [ np.log10( qn ) if not qn <= 0 else 0. for qn in qm_nm ]
        log_qm_nm = np.array( log_qm_nm )

        log_qm = [ np.log10( qq ) if not qq <= 0 else 0. for qq in qm ]
        log_qm = np.array( log_qm )
	

        ## make some plots
        mag_bin_mids = mag_bins[1:] + 0.2
        fig, axs = plt.subplots( 3, sharex=True, gridspec_kw={'hspace': 0})
	## plot 1
        axs[0].step( mag_bin_mids, log_total, where='mid', label='Total', color='0.5', linewidth=2 )
        axs[0].step( mag_bin_mids, log_real, where='mid', label='Real', linewidth=2, linestyle='dashed', color='black' )
        axs[0].step( mag_bin_mids, log_bkg, where='mid', label='Background', linewidth=2, linestyle='dotted', color='0.4' )
	## plot 2
        axs[1].step( mag_bin_mids, log_qm_nm, where='mid' )
	## plot 3
        axs[2].step( mag_bin_mids, log_qm, where='mid' )
	## axis labels
        axs[0].set( ylabel='log(N(cparts))' )
        axs[1].set( ylabel='log(P(m))' )
        axs[2].set( xlabel=my_band+' magnitudes', ylabel='log(q(m))' )
	## add legend
        axs[0].legend()
        ## adjust ylimits
        fig.savefig( my_band + '_magnitude_distributions.png' )
        fig.clear()
        
        ## Now calculate the LR and reliability
        final_file = LR_and_reliability( my_band, band_dat, radio_dat, qm_nm, sigma_pos, mag_bins, r_max, Q0, LR_threshold=LR_threshold, ra_col=ra_col, dec_col=dec_col, mag_col=mag_col, id_col=id_col )

        ## make another plot -- LR vs. separation
        final_matches = Table.read( final_file, format='ascii' )
        ## plot the LR values where the reliability is above the threshold
        rel_col = my_band + '_Rel'
        lr_thresh_idx = np.where( final_matches[rel_col] >= LR_threshold )[0]
        ## for ease of plotting
        sep_col = my_band + '_separation'
        xvals = final_matches[sep_col]
        fig, axs = plt.subplots( 2, sharex=True, gridspec_kw={'hspace':0})
        ## reliability
        axs[0].scatter( xvals, final_matches[rel_col], marker='.', color='0.5' )
        axs[0].scatter( xvals[lr_thresh_idx], final_matches[rel_col][lr_thresh_idx], marker='.', color='blue' )
        
        ###YG this breaks if no matches found
#        xr = np.arange( np.max( xvals ) )
        lr_col = my_band + '_LR'
        if len(xvals) == 0:
            xr = np.arange(0)
            ax1ylim = 1
        else:
            xr = np.arange(np.max( xvals ))
            ax1ylim = 3.*np.median(final_matches[lr_col])
        axs[0].plot( xr, np.repeat( LR_threshold, len(xr) ), color='black' )
        axs[1].scatter( xvals, final_matches[lr_col], marker='.', color='0.5' )
        axs[1].scatter( xvals[lr_thresh_idx], final_matches[lr_col][lr_thresh_idx], marker='.', color='blue' )
        ## axis limits
        axs[0].set( ylim=(0,1.1))
        axs[1].set( ylim=(0,ax1ylim) )
        ## axis labels
        axs[0].set( ylabel='Reliability' )
        axs[1].set( xlabel='separation [arcsec]', ylabel='LR value' )
        fig.savefig( my_band + '_LR_values.png' )
        fig.clear()


####directory clean up functions added by Yjan
def cleandir(original_files, outdir='lr_output'):
    newdlist = os.listdir()
    os.mkdir(outdir)
    for file in newdlist:
        if file not in original_files:
            dst = '/'.join([outdir, file])
            os.rename(file, dst)
    return


if __name__ == "__main__":

    dlist = os.listdir() ###added by YG for cleanup
    parser = argparse.ArgumentParser()
    parser.add_argument('multiwave_cat', type=str, help='input catalogue from which to generate a master catalogue' )
    parser.add_argument('radio_cat', type=str, help='input radio catalogue' )
    parser.add_argument('mask_image', type=str, help='Mask image where 1 is masked, 0 is unmasked' )
    parser.add_argument('--config_file', type=str, help='configuration file (defaukt lf_config.txt)', default='lr_config.txt' )
    parser.add_argument('--overwrite', action='store_true', help='set if you want to overwite a catalogue with the same name' )
    parser.add_argument('--snr_cut', type=float, default=5.0, help='SNR cut (applied individually to each band, default 5)' )
    parser.add_argument('--LR_thresh', type=float, default=0.8, help='threshold for LR (default 0.8)' )
    args = parser.parse_args()

    lr_main( args.multiwave_cat, args.radio_cat, args.mask_image, config_file=args.config_file, overwrite=args.overwrite, snr_cut=args.snr_cut, LR_threshold=args.LR_thresh )

    cleandir(dlist)

    #########################################
    #time taken for code to run
    print('END: ', ("--- %s seconds ---" % (time.time() - start_time)))
