from marvin.tools import Maps
from marvin.tools import Cube
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
import numpy as np
import numpy.ma as ma
import streamlit as st

from marvin import config
config.switchSasUrl(sasmode='mirror')
n = '7961-12704'
n = st.text_input("Type in a Marvin Galaxy ID: ", '7961-12704')
maps = Maps(n) #Pull Marvin map for a galaxy
cube = Cube(n)

snrdict = {}
masks, fig, axes = maps.get_bpt(snr_min=snrdict) #get bpt for selected galaxy MAP

ha = maps['emline_gflux_ha_6564'] #pulls ha flux
hb = maps['emline_gflux_hb_4862'] #pulls hb flux

fig, ax = ha.plot(cblabel='H-alpha')
st.pyplot(fig)

masks = maps.get_bpt(show_plot=False, return_figure=False) #show that BPT, baby

#make a mask of spaxels designated "unusable for science" and "non-starforming"
mask_non_sf = ~masks['sf']['global'] * ha.pixmask.labels_to_value('DONOTUSE')
# Do a bitwise OR between DAP mask and non-star-forming mask.
hamasked = (ha.mask | mask_non_sf)



bmask_non_sf = ~masks['sf']['global'] * ha.pixmask.labels_to_value('DONOTUSE') #mask out all non-starforming spaxels

# Do a bitwise OR between DAP mask and non-star-forming mask.
hbmasked = (hb.mask | bmask_non_sf)

quot = ha/hb

quotmask_non_sf = ~masks['sf']['global'] * ha.pixmask.labels_to_value('DONOTUSE') #mask out all non-starforming spaxels

# Do a bitwise OR between DAP mask and non-star-forming mask.
quotmasked = (quot.mask | quotmask_non_sf)

#dust_ha = hamasked * 10**(2.468*(0.934*np.log((hamasked) / 2.86))) #Dust corrected halpha flux
dust_ha = ha * 10**(2.468*(0.934*(np.log(quot.value) / 2.86))) #Dust corrected halpha flux


from astropy.io import fits

drpall = fits.open('drpall-v3_1_1.fits')  #assumes you are in the same directory as the DRPall file
tbdata = drpall[1].data

print('redshift =', tbdata['nsa_z'][0])   #printing the redshift of this index corresponding to our galaxy

c = 299792  # speed of light [km/s]
H0 = 70  # [km s^-1 Mpc^-1]
D = c * tbdata['nsa_z'][0] / H0  # approx. distance to galaxy [Mpc]
dist = 174000 #kpc from online cosmology calc. find a better way

lum_ha = dust_ha * 4 * np.pi * ((D*10**24)**2) #Dust-corrected luminosity of Halpha

import sys
np.set_printoptions(threshold=sys.maxsize)
sfr = 5.5*(10**(-42))*lum_ha * 10**-17 #sfr from Kennicutt's, adjusted for 10^-17

erad = maps['spx_ellcoo_elliptical_radius'] #elliptical radius
r = maps.header['reff'] #effective radius

roverr = erad/float(r) #elliptical radius / effective radius

spaxel_size = 0.5  # [arcsec]
redshift = tbdata['nsa_z'][0]

scale = 1 / 206265 * D * 1e6  # 1 radian = 206265 arcsec [pc / arcsec]
spaxel_area = (scale * spaxel_size)**2  # [pc^2]

sfrdlog = np.log10(10**sfr.value / spaxel_area)  # [sfr / pc^2]
sfrd = sfr/spaxel_area

#plt.imshow(sfr)
#plt.colorbar()
#plt.show()

fig, ax = sfrd.plot(mask=quotmasked, cblabel='Star Formation Rate')
st.pyplot(fig)