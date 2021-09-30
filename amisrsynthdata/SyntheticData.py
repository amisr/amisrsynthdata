# SyntheticData.py

from .Ionosphere import Ionosphere
from .Radar import Radar

import numpy as np
import datetime as dt
import configparser
import h5py
import pymap3d as pm

import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import matplotlib.gridspec as gridspec

class SyntheticData(object):

    def __init__(self, configfile):

        # parse items from config file
        config = configparser.ConfigParser()
        config.read(configfile)
        starttime = dt.datetime.fromisoformat(config['GENERAL']['STARTTIME'])
        endtime = dt.datetime.fromisoformat(config['GENERAL']['ENDTIME'])
        output_filename = config['GENERAL']['OUTPUT_FILENAME']

        # generate ionosphere object
        self.iono = Ionosphere(configfile)

        # generate radar object
        self.radar = Radar(configfile)

        self.generate_time_array(starttime, endtime)
        self.calc_radar_measurements()
        self.save_hdf5_output(output_filename)

    def generate_time_array(self, starttime, endtime):
        # create unix time array
        ust = (starttime-dt.datetime.utcfromtimestamp(0)).total_seconds()
        num_tstep = int((endtime-starttime).total_seconds()/self.radar.integration_period)
        self.utime = np.array([[ust+t*self.radar.integration_period, ust+(t+1)*self.radar.integration_period] for t in range(num_tstep)])


    def calc_radar_measurements(self):
        # caculate scalar ionosphere parameters at each fitted radar bin
        self.ne = self.iono.density(self.radar.lat, self.radar.lon, self.radar.alt)
        self.te = self.iono.etemp(self.radar.lat, self.radar.lon, self.radar.alt)
        self.ti = self.iono.itemp(self.radar.lat, self.radar.lon, self.radar.alt)

        # calculate LoS velocity for each bin by taking the dot product of the radar kvector and the velocity field
        self.kvec = self.radar.kvec_all_gates()
        Vvec = self.iono.velocity(self.radar.lat, self.radar.lon, self.radar.alt)
        self.Vlos = np.einsum('...i,...i->...',self.kvec, Vvec)

        # calculate density in ACF bins
        ne_nb_ot = self.iono.density(self.radar.lat_nb, self.radar.lon_nb, self.radar.alt_nb)
        self.ne_nb = np.full((self.utime.shape[0],)+ne_nb_ot.shape, np.nan)
        self.ne_nb[:,:,:] = ne_nb_ot


        # create fit and error arrays that match the shape of whats in the processed fitted files
        # Fit Array: Nrecords x Nbeams x Nranges x Nions+1 x 4 (fraction, temperature, coll. freq., LoS speed)
        # assume only O+, but include fields for other parameters so array is a general shape
        s = (self.utime.shape[0],)+self.radar.fit_slant_range.shape
        self.fit_array = np.full(s+(len(self.iono.ion_mass)+1,4), np.nan)
        self.fit_array[:,:,:,0,1] = self.ti
        self.fit_array[:,:,:,-1,1] = self.te
        self.fit_array[:,:,:,0,3] = self.Vlos
        self.fit_array[:,:,:,-1,3] = self.Vlos
        self.fit_array[:,:,:,:,0] = np.zeros(s+(len(self.iono.ion_mass)+1,))
        self.fit_array[:,:,:,0,0] = np.ones(s)
        self.fit_array[:,:,:,-1,0] = np.ones(s)

        self.ne_array = np.full(s, np.nan)
        self.ne_array[:,:,:] = self.ne

        # empty error arrays
        self.err_array = np.full(self.fit_array.shape,10.)
        self.dne = np.full(self.ne.shape, 1.e10)
        self.noise = np.full(s+(3,), np.nan)

        # fit parameters
        self.chi2 = np.full(s, 1.0)
        self.dof = np.full(s, 26)
        self.fitcode = np.full(s, 1)
        self.nfev = np.full(s, 0)


        self.snr = np.full(self.radar.slant_range.shape, np.nan)
        self.dnefrac = np.full(self.ne_nb.shape, np.nan)


    def save_hdf5_output(self, outfilename):
        # create output hdf5 file
        with h5py.File(outfilename, mode='w') as h5:

            h5.create_dataset('BeamCodes', data=self.radar.beam_codes)

            h5.create_group('/FittedParams')
            h5.create_dataset('/FittedParams/Altitude', data=self.radar.alt)
            h5.create_dataset('/FittedParams/Errors', data=self.err_array)
            h5.create_dataset('/FittedParams/Fits', data=self.fit_array)
            h5.create_dataset('/FittedParams/IonMass', data=self.iono.ion_mass)
            h5.create_dataset('/FittedParams/Ne', data=self.ne_array)
            h5.create_dataset('/FittedParams/Noise', data=self.noise)
            h5.create_dataset('/FittedParams/Range', data=self.radar.fit_slant_range)
            h5.create_dataset('/FittedParams/dNe', data=self.dne)

            h5.create_group('/FittedParams/FitInfo')
            h5.create_dataset('/FittedParams/FitInfo/chi2', data=self.chi2)
            h5.create_dataset('/FittedParams/FitInfo/dof', data=self.dof)
            h5.create_dataset('/FittedParams/FitInfo/fitcode', data=self.fitcode)
            h5.create_dataset('/FittedParams/FitInfo/nfev', data=self.nfev)

            h5.create_group('/Calibration')

            # need to call IGRF
            h5.create_group('/Geomag')
            h5.create_dataset('/Geomag/Latitude', data=self.radar.lat)
            h5.create_dataset('/Geomag/Longitude', data=self.radar.lon)
            h5.create_dataset('/Geomag/Altitude', data=self.radar.alt)
            h5.create_dataset('/Geomag/ke', data=self.kvec[:,:,0])
            h5.create_dataset('/Geomag/kn', data=self.kvec[:,:,1])
            h5.create_dataset('/Geomag/kz', data=self.kvec[:,:,2])

            # need to call MSIS-E
            h5.create_group('/MSIS')

            h5.create_group('/NeFromPower')
            h5.create_dataset('/NeFromPower/Altitude', data=self.radar.alt_nb)
            h5.create_dataset('/NeFromPower/Ne_Mod', data=self.ne_nb)
            h5.create_dataset('/NeFromPower/Ne_NoTr', data=self.ne_nb)
            h5.create_dataset('/NeFromPower/Range', data=self.radar.slant_range)
            h5.create_dataset('/NeFromPower/SNR', data=self.snr)
            h5.create_dataset('/NeFromPower/dNeFrac', data=self.dnefrac)

            h5.create_group('/ProcessingParams')

            h5.create_group('Site')
            h5.create_dataset('/Site/Latitude', data=self.radar.site_lat)
            h5.create_dataset('/Site/Longitude', data=self.radar.site_lon)
            h5.create_dataset('/Site/Altitude', data=self.radar.site_alt)
            h5.create_dataset('/Site/Code', data=0)
            h5.create_dataset('/Site/Name', data=self.radar.radar_name, dtype=h5py.string_dtype(encoding='utf-8',length=len(self.radar.radar_name.encode('utf-8'))))
            h5.create_dataset('/Site/MagneticLatitude', data=self.radar.site_lat)
            h5.create_dataset('/Site/MagneticLongitude', data=self.radar.site_lon)
            h5.create_dataset('/Site/MagneticLocalTimeMidnight', data=self.radar.site_alt)

            # all these fields can be calclulated from time array
            h5.create_group('Time')
            h5.create_dataset('/Time/Day', data=self.utime)
            h5.create_dataset('/Time/MagneticLocalTimeSite', data=self.utime)
            h5.create_dataset('/Time/Month', data=self.utime)
            h5.create_dataset('/Time/UnixTime', data=self.utime)
            h5.create_dataset('/Time/Year', data=self.utime)
            h5.create_dataset('/Time/doy', data=self.utime)
            h5.create_dataset('/Time/dtime', data=self.utime)

    def summary_plot(self):

        # summary plot of output
        alt_layers = np.arange(100.,500.,100.)
        e, n, u = np.meshgrid(np.arange(-500.,500.,50.)*1000., np.arange(0.,500.,50.)*1000., alt_layers*1000.)
        glat, glon, galt = pm.enu2geodetic(e, n, u, self.radar.site_lat, self.radar.site_lon, 0.)

        # glat, glon, galt = np.meshgrid(np.arange(radar.site_lat,radar.site_lat+10., 1.), np.arange(radar.site_lon-20.,radar.site_lon+20, 2.), np.arange(200000., 500000., 100000.))

        ne0 = self.iono.density(glat, glon, galt)
        te0 = self.iono.etemp(glat, glon, galt)
        ti0 = self.iono.itemp(glat, glon, galt)
        ve = self.iono.velocity(glat, glon, galt)

        # scaling/rotation of vector to plot in cartopy
        # https://github.com/SciTools/cartopy/issues/1179
        es = ve[:,:,:,0]/np.cos(glat*np.pi/180.)
        ns = ve[:,:,:,1]
        sf = np.sqrt(ve[:,:,:,0]**2+ve[:,:,:,1]**2)/np.sqrt(es**2+ns**2)
        e = es*sf
        n = ns*sf



        proj = ccrs.AzimuthalEquidistant(central_latitude=self.radar.site_lat, central_longitude=self.radar.site_lon)

        fig = plt.figure(figsize=(12,12))
        gs = gridspec.GridSpec(len(alt_layers)+1,4)

        for j in range(len(alt_layers)):
            ax = fig.add_subplot(gs[j,0], projection=proj)
            ax.coastlines()
            # ax.scatter(glon[:,:,j], glat[:,:,j], c=ne0[:,:,j], vmin=0., vmax=4e11, transform=ccrs.Geodetic())
            ax.contourf(glon[:,:,j], glat[:,:,j], ne0[:,:,j], vmin=0., vmax=4e11, cmap='viridis', transform=ccrs.PlateCarree())
            ax.set_title('{} km'.format(alt_layers[j]))

            ax = fig.add_subplot(gs[j,1], projection=proj)
            ax.coastlines()
            ax.quiver(glon[:,:,j], glat[:,:,j], e[:,:,j], n[:,:,j], color='blue', transform=ccrs.PlateCarree())
            ax.set_title('{} km'.format(alt_layers[j]))

            ax = fig.add_subplot(gs[j,2], projection=proj)
            ax.coastlines()
            # ax.scatter(glon[:,:,j], glat[:,:,j], c=te0[:,:,j], vmin=0., vmax=5e3, transform=ccrs.Geodetic())
            ax.contourf(glon[:,:,j], glat[:,:,j], te0[:,:,j], vmin=0., vmax=5e3, cmap='inferno', transform=ccrs.PlateCarree())
            ax.set_title('{} km'.format(alt_layers[j]))

            ax = fig.add_subplot(gs[j,3], projection=proj)
            ax.coastlines()
            # ax.scatter(glon[:,:,j], glat[:,:,j], c=ti0[:,:,j], vmin=0., vmax=3e3, transform=ccrs.Geodetic())
            ax.contourf(glon[:,:,j], glat[:,:,j], ti0[:,:,j], vmin=0., vmax=3e3, cmap='magma', transform=ccrs.PlateCarree())
            ax.set_title('{} km'.format(alt_layers[j]))

        x, y, z = pm.geodetic2ecef(self.radar.lat, self.radar.lon, self.radar.alt)
        fp = np.isfinite(self.radar.alt)

        ax = fig.add_subplot(gs[-1,0], projection='3d')
        c = ax.scatter(x[fp], y[fp], z[fp], c=self.ne[fp], vmin=0., vmax=4e11)
        ax.xaxis.set_ticklabels([])
        ax.yaxis.set_ticklabels([])
        ax.zaxis.set_ticklabels([])
        ax.set_box_aspect((1,1,1))
        fig.colorbar(c)

        ax = fig.add_subplot(gs[-1,1], projection='3d')
        c = ax.scatter(x[fp], y[fp], z[fp], c=self.Vlos[fp], vmin=-500., vmax=500.)
        ax.xaxis.set_ticklabels([])
        ax.yaxis.set_ticklabels([])
        ax.zaxis.set_ticklabels([])
        ax.set_box_aspect((1,1,1))
        fig.colorbar(c)

        ax = fig.add_subplot(gs[-1,2], projection='3d')
        c = ax.scatter(x[fp], y[fp], z[fp], c=self.te[fp], vmin=0., vmax=5e3)
        ax.xaxis.set_ticklabels([])
        ax.yaxis.set_ticklabels([])
        ax.zaxis.set_ticklabels([])
        ax.set_box_aspect((1,1,1))
        fig.colorbar(c)

        ax = fig.add_subplot(gs[-1,3], projection='3d')
        c = ax.scatter(x[fp], y[fp], z[fp], c=self.ti[fp], vmin=0., vmax=3e3)
        ax.xaxis.set_ticklabels([])
        ax.yaxis.set_ticklabels([])
        ax.zaxis.set_ticklabels([])
        ax.set_box_aspect((1,1,1))
        fig.colorbar(c)

        plt.show()
