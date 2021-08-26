# create_synth_data.py
# Create an AMISR data file with synthetic data

from Field import Field
from Radar import Radar
import numpy as np
import datetime as dt
import h5py
import matplotlib.pyplot as plt

from argparse import ArgumentParser, RawDescriptionHelpFormatter

config_file_help = 'Some help string'

# Build the argument parser tree
parser = ArgumentParser(description=config_file_help,
                        formatter_class=RawDescriptionHelpFormatter)
arg = parser.add_argument('synth_config_file',help='Configuration file for synthetic data set.')
# arg = parser.add_argument('vvels_config_file',help='Vvels config file.')
args = vars(parser.parse_args())


# generate field object
field = Field(args['synth_config_file'])
# field.plot_ionosphere()

# generate radar object
radar = Radar(args['synth_config_file'])

print(radar.X.shape)

ne = field.density(radar.X, radar.Y, radar.Z)
print(ne.shape)

# fig = plt.figure(figsize=(10,10))
# ax = fig.add_subplot(111,projection='3d')
# c = ax.scatter(radar.X, radar.Y, radar.Z, c=ne)
# fig.colorbar(c)
# plt.show()


# # interpolate the field to the radar bin locations
# Vx = field.interpVx(np.array([radar.X, radar.Y, radar.Z]).T)
# Vy = field.interpVy(np.array([radar.X, radar.Y, radar.Z]).T)
# Vz = field.interpVz(np.array([radar.X, radar.Y, radar.Z]).T)
# Vvec = np.array([Vx, Vy, Vz]).T
# print(Vvec.shape)

Vvec = field.velocity(radar.X, radar.Y, radar.Z)
print(Vvec.shape)
#
#
# calculate LoS velocity for each bin by taking the dot product of the radar kvector and the interpolated field
# Vlos = np.tile(np.einsum('...i,...i->...',radar.kvec, Vvec), (len(self.times),1,1))
Vlos = np.einsum('...i,...i->...',radar.kvec, Vvec)
# assume constant error
dVlos = np.full(Vlos.shape, radar.vel_error)

print(Vlos.shape)

# fig = plt.figure(figsize=(10,10))
# ax = fig.add_subplot(111,projection='3d')
# c = ax.scatter(radar.X, radar.Y, radar.Z, c=Vlos)
# fig.colorbar(c)
# plt.show()

# create unix time array
# Time array is 10 time steps (of integration period defined by the radar portion of the config file) after midnight
#   on January 1 of the year defined in the apex_year portion of the config file.  Because the field is defined manually
#   and not based on some empirical model, the time really doesn't matter and is mostly included to be consistent
#   with the read data file format.
time0 = (dt.datetime(field.apex_year,1,1)-dt.datetime.utcfromtimestamp(0)).total_seconds()
utime = np.array([[time0+t*radar.integration_period, time0+(t+1)*radar.integration_period] for t in range(10)])


# create fit and error arrays that match the shape of whats in the processed fitted files
s = (utime.shape[0],)+radar.X.shape
fit_array = np.full(s+(6,4), np.nan)
fit_array[:,:,:,0,3] = Vlos
err_array = np.full(s+(6,4),np.nan)
err_array[:,:,:,0,3] = dVlos

chi2 = np.full(s, 1.0)
fitcode = np.full(s, 1)

print(fit_array.shape, chi2.shape, fitcode.shape)

# # generate dummy density array
# self.ne = np.full(s, 1e11)

# create output hdf5 file
with h5py.File(radar.output_filename, mode='w') as h5:
    h5.create_group('FittedParams')
    h5.create_group('Geomag')
    h5.create_group('Time')
    h5.create_group('Site')

    h5.create_group('/FittedParams/FitInfo')

    h5.create_dataset('BeamCodes', data=radar.beam_codes)

    h5.create_dataset('/FittedParams/Fits', data=fit_array)
    h5.create_dataset('/FittedParams/Errors', data=err_array)
    h5.create_dataset('/FittedParams/Ne', data=ne)

    h5.create_dataset('/FittedParams/FitInfo/chi2', data=chi2)
    h5.create_dataset('/FittedParams/FitInfo/fitcode', data=fitcode)

    h5.create_dataset('/Geomag/Latitude', data=radar.lat)
    h5.create_dataset('/Geomag/Longitude', data=radar.lon)
    h5.create_dataset('/Geomag/Altitude', data=radar.alt)

    h5.create_dataset('/Geomag/ke', data=radar.ke)
    h5.create_dataset('/Geomag/kn', data=radar.kn)
    h5.create_dataset('/Geomag/kz', data=radar.ku)

    h5.create_dataset('/Time/UnixTime', data=utime)

    h5.create_dataset('/Site/Latitude', data=radar.site_coords[0])
    h5.create_dataset('/Site/Longitude', data=radar.site_coords[1])
    h5.create_dataset('/Site/Altitude', data=radar.site_coords[2])
