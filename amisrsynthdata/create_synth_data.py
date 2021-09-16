# create_synth_data.py
# Create an AMISR data file with synthetic data

from .Ionosphere import Ionosphere
from .Radar import Radar
import numpy as np
import datetime as dt
import pymap3d as pm
import h5py
try:
    import ConfigParser as configparser
except ImportError:
    import configparser
import matplotlib.pyplot as plt

from argparse import ArgumentParser, RawDescriptionHelpFormatter

config_file_help = 'Some help string'

# Build the argument parser tree
parser = ArgumentParser(description=config_file_help,
                        formatter_class=RawDescriptionHelpFormatter)
arg = parser.add_argument('synth_config_file',help='Configuration file for synthetic data set.')
args = vars(parser.parse_args())

config = configparser.ConfigParser()
config.read(args['synth_config_file'])
starttime = dt.datetime.fromisoformat(config['GENERAL']['STARTTIME'])
endtime = dt.datetime.fromisoformat(config['GENERAL']['ENDTIME'])
output_filename = config['GENERAL']['OUTPUT_FILENAME']

# generate ionosphere object
iono = Ionosphere(args['synth_config_file'])

# generate radar object
radar = Radar(args['synth_config_file'])

# caculate scalar ionosphere parameters at each fitted radar bin
ne = iono.density(radar.lat, radar.lon, radar.alt)
te = iono.etemp(radar.lat, radar.lon, radar.alt)
ti = iono.itemp(radar.lat, radar.lon, radar.alt)

# calculate LoS velocity for each bin by taking the dot product of the radar kvector and the velocity field
kvec = radar.kvec_all_gates()
Vvec = iono.velocity(radar.lat, radar.lon, radar.alt)
Vlos = np.einsum('...i,...i->...',kvec, Vvec)

# calculate density in ACF bins
ne_nb = iono.density(radar.lat_nb, radar.lon_nb, radar.alt_nb)


# create unix time array
ust = (starttime-dt.datetime.utcfromtimestamp(0)).total_seconds()
num_tstep = int((endtime-starttime).total_seconds()/radar.integration_period)
utime = np.array([[ust+t*radar.integration_period, ust+(t+1)*radar.integration_period] for t in range(num_tstep)])


# create fit and error arrays that match the shape of whats in the processed fitted files
# Fit Array: Nrecords x Nbeams x Nranges x Nions+1 x 4 (fraction, temperature, coll. freq., LoS speed)
# assume only O+, but include fields for other parameters so array is a general shape
s = (utime.shape[0],)+radar.fit_slant_range.shape
fit_array = np.full(s+(len(iono.ion_mass)+1,4), np.nan)
fit_array[:,:,:,0,1] = ti
fit_array[:,:,:,-1,1] = te
fit_array[:,:,:,0,3] = Vlos
fit_array[:,:,:,-1,3] = Vlos
fit_array[:,:,:,:,0] = np.zeros(s+(len(iono.ion_mass)+1,))
fit_array[:,:,:,0,0] = np.ones(s)
fit_array[:,:,:,-1,0] = np.ones(s)

ne_array = np.full(s, np.nan)
ne_array[:,:,:] = ne

# empty error arrays
err_array = np.full(fit_array.shape,10.)
dne = np.full(ne.shape, 1.e10)
noise = np.full(s+(3,), np.nan)

# fit parameters
chi2 = np.full(s, 1.0)
dof = np.full(s, 26)
fitcode = np.full(s, 1)
nfev = np.full(s, 0)


snr = np.full(radar.slant_range.shape, np.nan)
dnefrac = np.full(ne_nb.shape, np.nan)


# create output hdf5 file
with h5py.File(output_filename, mode='w') as h5:

    h5.create_dataset('BeamCodes', data=radar.beam_codes)

    h5.create_group('/FittedParams')
    h5.create_dataset('/FittedParams/Altitude', data=radar.alt)
    h5.create_dataset('/FittedParams/Errors', data=err_array)
    h5.create_dataset('/FittedParams/Fits', data=fit_array)
    h5.create_dataset('/FittedParams/IonMass', data=iono.ion_mass)
    h5.create_dataset('/FittedParams/Ne', data=ne_array)
    h5.create_dataset('/FittedParams/Noise', data=noise)
    h5.create_dataset('/FittedParams/Range', data=radar.fit_slant_range)
    h5.create_dataset('/FittedParams/dNe', data=dne)

    h5.create_group('/FittedParams/FitInfo')
    h5.create_dataset('/FittedParams/FitInfo/chi2', data=chi2)
    h5.create_dataset('/FittedParams/FitInfo/dof', data=dof)
    h5.create_dataset('/FittedParams/FitInfo/fitcode', data=fitcode)
    h5.create_dataset('/FittedParams/FitInfo/nfev', data=nfev)

    h5.create_group('/Calibration')

    # need to call IGRF
    h5.create_group('/Geomag')
    h5.create_dataset('/Geomag/Latitude', data=radar.lat)
    h5.create_dataset('/Geomag/Longitude', data=radar.lon)
    h5.create_dataset('/Geomag/Altitude', data=radar.alt)
    h5.create_dataset('/Geomag/ke', data=kvec[:,:,0])
    h5.create_dataset('/Geomag/kn', data=kvec[:,:,1])
    h5.create_dataset('/Geomag/kz', data=kvec[:,:,2])

    # need to call MSIS-E
    h5.create_group('/MSIS')

    h5.create_group('/NeFromPower')
    h5.create_dataset('/NeFromPower/Altitude', data=radar.alt_nb)
    h5.create_dataset('/NeFromPower/Ne_Mod', data=ne_nb)
    h5.create_dataset('/NeFromPower/Ne_NoTr', data=ne_nb)
    h5.create_dataset('/NeFromPower/Range', data=radar.slant_range)
    h5.create_dataset('/NeFromPower/SNR', data=snr)
    h5.create_dataset('/NeFromPower/dNeFrac', data=dnefrac)

    h5.create_group('/ProcessingParams')

    h5.create_group('Site')
    h5.create_dataset('/Site/Latitude', data=radar.site_lat)
    h5.create_dataset('/Site/Longitude', data=radar.site_lon)
    h5.create_dataset('/Site/Altitude', data=radar.site_alt)
    h5.create_dataset('/Site/Code', data=0)
    h5.create_dataset('/Site/Name', data=radar.radar_name, dtype=h5py.string_dtype(encoding='utf-8',length=len(radar.radar_name.encode('utf-8'))))
    h5.create_dataset('/Site/MagneticLatitude', data=radar.site_lat)
    h5.create_dataset('/Site/MagneticLongitude', data=radar.site_lon)
    h5.create_dataset('/Site/MagneticLocalTimeMidnight', data=radar.site_alt)

    # all these fields can be calclulated from time array
    h5.create_group('Time')
    h5.create_dataset('/Time/Day', data=utime)
    h5.create_dataset('/Time/MagneticLocalTimeSite', data=utime)
    h5.create_dataset('/Time/Month', data=utime)
    h5.create_dataset('/Time/UnixTime', data=utime)
    h5.create_dataset('/Time/Year', data=utime)
    h5.create_dataset('/Time/doy', data=utime)
    h5.create_dataset('/Time/dtime', data=utime)



# summary plot of output
glat, glon, galt = np.meshgrid(np.arange(radar.site_lat,radar.site_lat+10., 1.), np.arange(radar.site_lon-20.,radar.site_lon+20, 2.), np.arange(200000., 500000., 200000.))

ne0 = iono.density(glat, glon, galt)
te0 = iono.etemp(glat, glon, galt)
ti0 = iono.itemp(glat, glon, galt)
ve = iono.velocity(glat, glon, galt)

# scaling/rotation of vector to plot in cartopy
# https://github.com/SciTools/cartopy/issues/1179
es = ve[:,:,:,0]/np.cos(glat*np.pi/180.)
ns = ve[:,:,:,1]
sf = np.sqrt(ve[:,:,:,0]**2+ve[:,:,:,1]**2)/np.sqrt(es**2+ns**2)
e = es*sf
n = ns*sf


import cartopy.crs as ccrs
import matplotlib.gridspec as gridspec

proj = ccrs.AzimuthalEquidistant(central_latitude=radar.site_lat, central_longitude=radar.site_lon)

fig = plt.figure(figsize=(12,6))
gs = gridspec.GridSpec(3,4)

for j in range(2):
    ax = fig.add_subplot(gs[j,0], projection=proj)
    ax.coastlines()
    ax.scatter(glon[:,:,j], glat[:,:,j], c=ne0[:,:,j], vmin=0., vmax=4e11, transform=ccrs.Geodetic())

    ax = fig.add_subplot(gs[j,1], projection=proj)
    ax.coastlines()
    ax.quiver(glon[:,:,j], glat[:,:,j], e[:,:,j], n[:,:,j], color='blue', transform=ccrs.PlateCarree())

    ax = fig.add_subplot(gs[j,2], projection=proj)
    ax.coastlines()
    ax.scatter(glon[:,:,j], glat[:,:,j], c=te0[:,:,j], vmin=0., vmax=5e3, transform=ccrs.Geodetic())

    ax = fig.add_subplot(gs[j,3], projection=proj)
    ax.coastlines()
    ax.scatter(glon[:,:,j], glat[:,:,j], c=ti0[:,:,j], vmin=0., vmax=3e3, transform=ccrs.Geodetic())

x, y, z = pm.geodetic2ecef(radar.lat, radar.lon, radar.alt)
fp = np.isfinite(radar.alt)

ax = fig.add_subplot(gs[-1,0], projection='3d')
c = ax.scatter(x[fp], y[fp], z[fp], c=ne[fp], vmin=0., vmax=4e11)
ax.xaxis.set_ticklabels([])
ax.yaxis.set_ticklabels([])
ax.zaxis.set_ticklabels([])
ax.set_box_aspect((1,1,1))
fig.colorbar(c)

ax = fig.add_subplot(gs[-1,1], projection='3d')
c = ax.scatter(x[fp], y[fp], z[fp], c=Vlos[fp], vmin=-500., vmax=500.)
ax.xaxis.set_ticklabels([])
ax.yaxis.set_ticklabels([])
ax.zaxis.set_ticklabels([])
ax.set_box_aspect((1,1,1))
fig.colorbar(c)

ax = fig.add_subplot(gs[-1,2], projection='3d')
c = ax.scatter(x[fp], y[fp], z[fp], c=te[fp], vmin=0., vmax=5e3)
ax.xaxis.set_ticklabels([])
ax.yaxis.set_ticklabels([])
ax.zaxis.set_ticklabels([])
ax.set_box_aspect((1,1,1))
fig.colorbar(c)

ax = fig.add_subplot(gs[-1,3], projection='3d')
c = ax.scatter(x[fp], y[fp], z[fp], c=ti[fp], vmin=0., vmax=3e3)
ax.xaxis.set_ticklabels([])
ax.yaxis.set_ticklabels([])
ax.zaxis.set_ticklabels([])
ax.set_box_aspect((1,1,1))
fig.colorbar(c)

plt.show()
