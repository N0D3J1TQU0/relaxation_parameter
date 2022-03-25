# relaxation_parameter
Relaxation parameter (Wen &amp; Han et al. 2013) for SPT-CLJ2344-4224
This script generates the smoothed optical map, calculate the Asymmetry factor, the Ridge flatness, the Normalized deviation, and finally, the relaxation parameter from Wen &amp; Han et al. 2013. The example cluster SPT-CLJ2344-4224 has public data in https://deslabs.ncsa.illinois.edu/index.html.
INPUTS:
cluster redshift = 0.282384
SZ center (deg) = (356.1481, -42.41)
Cosmology from Bocquet et al. 2015
cluster m* i-band = 18.5887
the following lines shows the columns on each table data used in the code
SPT-CLJ2344-4224_redsequence.des
['COADD_OBJECT_ID',
 'RA',
 'DEC',
 'MAG_AUTO_G',
 'MAGERR_AUTO_G',
 'MAG_AUTO_R',
 'MAGERR_AUTO_R',
 'MAG_AUTO_I',
 'MAGERR_AUTO_I',
 'MAG_AUTO_Z',
 'MAGERR_AUTO_Z',
 'MOF_PSF_MAG_G',
 'MOF_PSF_MAG_ERR_G',
 'MOF_PSF_MAG_R',
 'MOF_PSF_MAG_ERR_R',
 'MOF_PSF_MAG_I',
 'MOF_PSF_MAG_ERR_I',
 'MOF_PSF_MAG_Z',
 'MOF_PSF_MAG_ERR_Z',
 'MOF_BDF_MAG_G_CORRECTED',
 'MOF_BDF_MAG_R_CORRECTED',
 'MOF_BDF_MAG_I_CORRECTED',
 'MOF_BDF_MAG_Z_CORRECTED',
 'FLAGS_FOOTPRINT',
 'FLAGS_BADREGIONS',
 'FLAGS_GOLD',
 'FLAGS_FOREGROUND',
 'EXT_MASH',
 'A_IMAGE',
 'B_IMAGE',
 'THETA_J2000',
 'DNF_ZMC_SOF',
 'DNF_ZMEAN_SOF',
 'DNF_ZSIGMA_SOF']
OUTPUT:
reltab.csv:
['SPT-CL', 'alpha', 'beta', 'delta', 'gamma']

