# SyntheticData.py

from .ionosphere import Ionosphere
from .radar import Radar

import numpy as np
import datetime as dt
import h5py
from argparse import ArgumentParser, RawDescriptionHelpFormatter
import yaml

class SyntheticData(object):

    def __init__(self, config):

        starttime = config['GENERAL']['starttime']
        endtime = config['GENERAL']['endtime']
        output_filename = config['GENERAL']['output_filename']
        err_coef = config['GENERAL']['err_coef']
        add_noise = config['GENERAL']['noise']

        # generate ionosphere object
        self.iono = Ionosphere(config)

        # generate radar object
        self.radar = Radar(config)

        self.generate_time_array(starttime, endtime)
        self.generate_geomag()
        self.generate_site(starttime)
        self.calc_radar_measurements()
        self.calc_errors(err_coef)
        if add_noise:
            self.add_measurment_noise()
        self.save_hdf5_output(output_filename)
        if 'SUMMARY_PLOT' in config:
            self.summary_plot(config['SUMMARY_PLOT'])

    def generate_time_array(self, starttime, endtime):
        # create time arrays
        ust = (starttime-dt.datetime.utcfromtimestamp(0)).total_seconds()
        num_tstep = int((endtime-starttime).total_seconds()/self.radar.integration_period)
        self.utime = np.array([ust+np.arange(0,num_tstep)*self.radar.integration_period, ust+np.arange(1,num_tstep+1)*self.radar.integration_period]).T
        time = np.array([[dt.datetime.utcfromtimestamp(t) for t in ut] for ut in self.utime])

        self.Time = {'UnixTime':self.utime}

        self.Time['Day'] = np.array([[t.day for t in ip] for ip in time])
        self.Time['Month'] = np.array([[t.month for t in ip] for ip in time])
        self.Time['Year'] = np.array([[t.year for t in ip] for ip in time])
        self.Time['doy'] = np.array([[t.timetuple().tm_yday for t in ip] for ip in time])
        self.Time['dtime'] = np.array([[(t-dt.datetime(t.year,t.month,t.day,0,0,0)).total_seconds()/3600. for t in ip] for ip in time])

        _, mlt0 = self.iono.apex.convert(self.radar.site_lat, self.radar.site_lon, 'geo', 'mlt', height=self.radar.site_alt/1000., datetime=time[:,0])
        _, mlt1 = self.iono.apex.convert(self.radar.site_lat, self.radar.site_lon, 'geo', 'mlt', height=self.radar.site_alt/1000., datetime=time[:,1])
        self.Time['MagneticLocalTimeSite'] = np.array([mlt0,mlt1]).T

    def generate_geomag(self):
        # generate Geomag array
        # Reqires running IGRF to get all fields
        self.Geomag = {'Latitude':self.radar.lat, 'Longitude':self.radar.lon, 'Altitude':self.radar.alt}


    def generate_site(self, st):
        self.Site = {'Latitude':self.radar.site_lat, 'Longitude':self.radar.site_lon, 'Altitude':self.radar.site_alt, 'Code':0}
        mlat, mlon = self.iono.apex.geo2apex(self.radar.site_lat, self.radar.site_lon, height=self.radar.site_alt/1000.)
        _, mlt = self.iono.apex.convert(self.radar.site_lat, self.radar.site_lon, 'geo', 'mlt', height=self.radar.site_alt/1000., datetime=dt.datetime(st.year,st.month,st.day,0,0,0))
        self.Site.update(MagneticLatitude=mlat, MagneticLongitude=mlon, MagneticLocalTimeMidnight=mlt)
        # self.Site.update(MagneticLatitude=0., MagneticLongitude=0., MagneticLocalTimeMidnight=0.)


    def calc_radar_measurements(self):
        # caculate scalar ionosphere parameters at each fitted radar bin
        self.ne = self.iono.density(self.utime[:,0], self.radar.lat, self.radar.lon, self.radar.alt)
        self.te = self.iono.etemp(self.utime[:,0], self.radar.lat, self.radar.lon, self.radar.alt)
        self.ti = self.iono.itemp(self.utime[:,0], self.radar.lat, self.radar.lon, self.radar.alt)

        # calculate LoS velocity for each bin by taking the dot product of the radar kvector and the velocity field
        kvec = self.radar.kvec_all_gates()
        Vvec = self.iono.velocity(self.utime[:,0], self.radar.lat, self.radar.lon, self.radar.alt)
        self.Vlos = np.einsum('...i,k...i->k...', kvec, Vvec)

        self.Geomag.update(ke=kvec[:,:,0], kn=kvec[:,:,1], kz=kvec[:,:,2])


        self.FittedParams = {'Altitude':self.radar.alt, 'IonMass':self.iono.ion_mass, 'Range':self.radar.slant_range}


        # create fit and error arrays that match the shape of whats in the processed fitted files
        # Fit Array: Nrecords x Nbeams x Nranges x Nions+1 x 4 (fraction, temperature, coll. freq., LoS speed)
        # assume only O+, but include fields for other parameters so array is a general shape
        s = (self.utime.shape[0],)+self.radar.slant_range.shape
        self.FittedParams['Fits'] = np.full(s+(len(self.iono.ion_mass)+1,4), np.nan)
        self.FittedParams['Fits'][:,:,:,0,1] = self.ti
        self.FittedParams['Fits'][:,:,:,-1,1] = self.te
        self.FittedParams['Fits'][:,:,:,0,3] = self.Vlos
        self.FittedParams['Fits'][:,:,:,-1,3] = self.Vlos
        self.FittedParams['Fits'][:,:,:,:,0] = np.zeros(s+(len(self.iono.ion_mass)+1,))
        self.FittedParams['Fits'][:,:,:,0,0] = np.ones(s)
        self.FittedParams['Fits'][:,:,:,-1,0] = np.ones(s)

        # do this kind of stuff with explicit broadcasting
        # self.FittedParams['Ne'] = np.full(s, np.nan)
        # self.FittedParams['Ne'] = np.broadcast_to(self.ne, s)
        self.FittedParams['Ne'] = self.ne

        self.FittedParams['Noise'] = np.full(s+(3,), np.nan)

        # fit info
        self.FitInfo = {'chi2':np.full(s, 1.0), 'dof':np.full(s, 26), 'fitcode':np.full(s, 1), 'nfev':np.full(s, 0)}


        # calculate density in ACF bins
        self.NeFromPower = {'Altitude':self.radar.alt_p, 'Range':self.radar.slant_range_p}

        ne_p = self.iono.density(self.utime[:,0], self.radar.lat_p, self.radar.lon_p, self.radar.alt_p)
        # self.NeFromPower['Ne_Mod'] = np.broadcast_to(ne_p, (self.utime.shape[0],)+ne_p.shape)
        # self.NeFromPower['Ne_NoTr'] = np.broadcast_to(ne_p, (self.utime.shape[0],)+ne_p.shape)
        self.NeFromPower['Ne_Mod'] = ne_p
        self.NeFromPower['Ne_NoTr'] = ne_p
        self.NeFromPower['SNR'] = np.full(self.radar.slant_range_p.shape, np.nan)
        self.NeFromPower['dNeFrac'] = np.full(self.NeFromPower['Ne_NoTr'].shape, np.nan)

    def calc_errors(self, err_coef):

        ne_err = err_coef[0] * self.radar.slant_range**2
        ve_err = err_coef[1] * self.radar.slant_range**2
        te_err = err_coef[2] * self.radar.slant_range**2
        ti_err = err_coef[3] * self.radar.slant_range**2

        # fill out hdf5 arrays
        self.FittedParams['dNe'] = np.broadcast_to(ne_err, self.FittedParams['Ne'].shape)
        self.FittedParams['Errors'] = np.full(self.FittedParams['Fits'].shape, np.nan)
        self.FittedParams['Errors'][:,:,:,0,1] = ti_err
        self.FittedParams['Errors'][:,:,:,-1,1] = te_err
        self.FittedParams['Errors'][:,:,:,0,3] = ve_err
        self.FittedParams['Errors'][:,:,:,-1,3] = ve_err


    def add_measurment_noise(self):

        # fitted parameters
        self.FittedParams['Ne'] = np.random.normal(loc=self.FittedParams['Ne'], scale=self.FittedParams['dNe'])

        # ACF parameters
        self.NeFromPower['Ne_Mod'] = np.random.normal(loc=self.NeFromPower['Ne_Mod'], scale=self.NeFromPower['dNeFrac'])
        self.NeFromPower['Ne_NoTr'] = np.random.normal(loc=self.NeFromPower['Ne_NoTr'], scale=self.NeFromPower['dNeFrac'])

    def save_hdf5_output(self, outfilename):
        # create output hdf5 file
        with h5py.File(outfilename, mode='w') as h5:

            h5.create_dataset('BeamCodes', data=self.radar.beam_codes)

            h5.create_group('/FittedParams')
            for k, v in self.FittedParams.items():
                h5.create_dataset('/FittedParams/{}'.format(k), data=v)

            h5.create_group('/FittedParams/FitInfo')
            for k, v in self.FitInfo.items():
                h5.create_dataset('/FittedParams/FitInfo/{}'.format(k), data=v)

            h5.create_group('/Calibration')

            h5.create_group('/Geomag')
            for k, v in self.Geomag.items():
                h5.create_dataset('/Geomag/{}'.format(k), data=v)

            # need to call MSIS-E
            h5.create_group('/MSIS')

            h5.create_group('/NeFromPower')
            for k, v in self.NeFromPower.items():
                h5.create_dataset('/NeFromPower/{}'.format(k), data=v)

            h5.create_group('/ProcessingParams')

            h5.create_group('Site')
            for k, v in self.Site.items():
                h5.create_dataset('/Site/{}'.format(k), data=v)
            h5.create_dataset('/Site/Name', data=self.radar.radar_name, dtype=h5py.string_dtype(encoding='utf-8',length=len(self.radar.radar_name.encode('utf-8'))))

            h5.create_group('Time')
            for k, v in self.Time.items():
                h5.create_dataset('/Time/{}'.format(k), data=v)


    def summary_plot(self, config):

        # optional imports used ONLY for creating summary plots
        # matplotlib and cartopy are not listed in the package requirments
        try:
            import pymap3d as pm
            import matplotlib.pyplot as plt
            import matplotlib.gridspec as gridspec
            import cartopy.crs as ccrs
        except ModuleNotFoundError as e:
            raise ModuleNotFoundError('In order to create summary plots, the optional modules matplotlib (https://matplotlib.org/) '
                                      'and cartopy (https://scitools.org.uk/cartopy/docs/latest/) must be installed.') from e


        # form grid of coordinates for plotting
        # Make customizable in config file
        alt_layers = np.array(config['alt_slices'])
        e, n = np.meshgrid(np.arange(*config['slice_lon_rng']), np.arange(*config['slice_lat_rng']))
        glat, glon, galt = pm.enu2geodetic(e, n, 0., self.radar.site_lat, self.radar.site_lon, 0.)

        glat = np.tile(glat, alt_layers.shape).reshape(glat.shape+alt_layers.shape, order='F')
        glon = np.repeat(glon, alt_layers.shape).reshape(glon.shape+alt_layers.shape, order='C')
        galt = np.broadcast_to(alt_layers, galt.shape+alt_layers.shape)

        tidx = np.argmin(np.abs((config['time']-dt.datetime.utcfromtimestamp(0)).total_seconds()-self.Time['UnixTime'][:,0]))

        ne0 = np.squeeze(self.iono.density(self.Time['UnixTime'][tidx,0], glat, glon, galt))
        te0 = np.squeeze(self.iono.etemp(self.Time['UnixTime'][tidx,0], glat, glon, galt))
        ti0 = np.squeeze(self.iono.itemp(self.Time['UnixTime'][tidx,0], glat, glon, galt))
        ve = np.squeeze(self.iono.velocity(self.Time['UnixTime'][tidx, 0], glat, glon, galt))


        # scaling/rotation of vector to plot in cartopy
        # https://github.com/SciTools/cartopy/issues/1179
        es = ve[:,:,:,0]/np.cos(glat*np.pi/180.)
        ns = ve[:,:,:,1]
        sf = np.sqrt(ve[:,:,:,0]**2+ve[:,:,:,1]**2)/np.sqrt(es**2+ns**2)
        e = es*sf
        n = ns*sf


        plotting_params = [dict(synthdata=self.ne, param=ne0, cparam=config['dens_colors'], label=r'Ne (m$^{-3}$)', title='Electron Density', output=config['output_prefix']+'ne.png'),
                           dict(synthdata=self.te, param=te0, cparam=config['etemp_colors'], label=r'Te (K)', title='Electron Temperature', output=config['output_prefix']+'te.png'),
                           dict(synthdata=self.ti, param=ti0, cparam=config['itemp_colors'], label=r'Ti (K)', title='Ion Temperature', output=config['output_prefix']+'ti.png'),
                           dict(synthdata=self.Vlos, param=[e,n], cparam=config['vlos_colors'], label=r'Vlos (m/s)', title='Plasma Velocity', output=config['output_prefix']+'vi.png')]

        proj = ccrs.AzimuthalEquidistant(central_latitude=self.radar.site_lat, central_longitude=self.radar.site_lon)
        Nalt = len(alt_layers)
        wr = [1]*Nalt
        wr.append(Nalt/2.)
        gs = gridspec.GridSpec(2, Nalt+1, width_ratios=wr, left=0.1, right=0.95, bottom=0.1, top=0.95, wspace=0.1, hspace=0.1)
        bidx = np.where(self.radar.beam_codes[:,0]==config['beam'])[0][0]

        for p in plotting_params:

            fig = plt.figure(figsize=(15,7))
            fig.suptitle(p['title'])


            # Create slice plots
            for j in range(len(alt_layers)):
    
                ax = fig.add_subplot(gs[0,j], projection=proj)
                ax.coastlines()
                ax.set_title('{} km'.format(alt_layers[j]/1000.))

                if p['title'] == 'Plasma Velocity':
                    ax.quiver(glon[:,:,j], glat[:,:,j], p['param'][0][:,:,j], p['param'][1][:,:,j], color='blue', transform=ccrs.PlateCarree())

                else:
                    cs = ax.contourf(glon[:,:,j], glat[:,:,j], p['param'][:,:,j], **p['cparam'], transform=ccrs.PlateCarree())
                    cs.cmap.set_over('white')
                    cs.cmap.set_under('grey')
                    cs.changed()

                # Add beam positions
                aidx = np.nanargmin(np.abs(self.radar.alt-alt_layers[j]),axis=1)
                aidx0 = np.array([aidx]).T
                slice_lat = np.take_along_axis(self.radar.lat, aidx0, axis=1)
                slice_lon = np.take_along_axis(self.radar.lon, aidx0, axis=1)
                ax.scatter(slice_lon, slice_lat, facecolors='none', edgecolors='k', transform=ccrs.Geodetic())
                ax.scatter(slice_lon[bidx], slice_lat[bidx], facecolors='none', edgecolors='magenta', transform=ccrs.Geodetic())



            # Create RTI
            ax = fig.add_subplot(gs[1,:-1])
            time = self.utime[:,0].astype('datetime64[s]')
            alt = self.radar.alt[bidx,:]
            c = ax.pcolormesh(time, alt[np.isfinite(alt)], p['synthdata'][:,bidx,np.isfinite(alt)].T, **p['cparam'])
            ax.axvline(x=config['time'], color='magenta')
            ax.set_xlabel('Universal Time')
            ax.set_ylabel('Altitude (m)')
            ax.set_title('Beam Number: {:.0f} (az:{:.1f}, el:{:.1f})'.format(self.radar.beam_codes[bidx,0], self.radar.beam_codes[bidx,1], self.radar.beam_codes[bidx,2]))


            # Create 3D FoV plot
            x, y, z = pm.geodetic2enu(self.radar.lat, self.radar.lon, self.radar.alt, self.radar.site_lat, self.radar.site_lon, self.radar.site_alt)
            fp = np.isfinite(self.radar.alt)
    
            ax = fig.add_subplot(gs[:,-1], projection='3d')
            c = ax.scatter(x[fp], y[fp], z[fp], c=p['synthdata'][tidx,fp], **p['cparam'])
            #ax.xaxis.set_ticklabels([])
            #ax.yaxis.set_ticklabels([])
            #ax.zaxis.set_ticklabels([])
            ax.set_box_aspect((1,1,1))
            ax.set_box_aspect([ub - lb for lb, ub in (getattr(ax, f'get_{a}lim')() for a in 'xyz')])
            fig.colorbar(c, label=p['label'])
    
            fig.savefig(p['output'])



def main():
    help_string = 'Program to generate synthetic data for multibeam, AMISR like incoherent scatter radars.'

    # Build the argument parser tree
    parser = ArgumentParser(description=help_string,
                            formatter_class=RawDescriptionHelpFormatter)
    arg = parser.add_argument('synth_config_file', help='Configuration file for synthetic data set.')
    args = vars(parser.parse_args())


    with open(args['synth_config_file'], 'r') as cf:
        config = yaml.safe_load(cf)

    sd = SyntheticData(config)
