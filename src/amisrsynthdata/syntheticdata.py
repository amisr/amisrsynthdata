# syntheticdata.py

from .ionosphere import Ionosphere
from .radar import Radar

import numpy as np
import datetime as dt
import h5py
from argparse import ArgumentParser, RawDescriptionHelpFormatter
import yaml


class SyntheticData(object):

    def __init__(self, config):

        starttime = config['GENERAL']['starttime']  # init input
        endtime = config['GENERAL']['endtime']      # init input
        err_coef = config['GENERAL']['err_coef']

        # generate ionosphere object
        self.iono = Ionosphere(config)

        # generate radar object
        self.radar = Radar(config)

        self.utime, self.time = self.generate_time_array(starttime, endtime)

        (
            self.ne,
            self.ti,
            self.te,
            self.vlos,
            self.ne_notr
        ) = self.generate_radar_measurements()

        (
            self.ne_err,
            self.ti_err,
            self.te_err,
            self.vlos_err,
            self.ne_notr_err
        ) = self.generate_errors(err_coef)

    def generate_time_array(self, starttime, endtime):
        # create time arrays
        ust = (starttime - dt.datetime.utcfromtimestamp(0)).total_seconds()
        num_tstep = int(
            (endtime - starttime).total_seconds() /
            self.radar.integration_period)
        utime = np.array([
            ust + np.arange(0, num_tstep) * self.radar.integration_period,
            ust + np.arange(1, num_tstep + 1) * self.radar.integration_period
            ]).T
        time = np.array([
            [dt.datetime.utcfromtimestamp(t) for t in ut] for ut in utime
            ])
        return utime, time

    def generate_radar_measurements(self):
        # caculate scalar ionosphere parameters at each fitted radar bin
        ne = self.iono.density(self.utime[:, 0], self.radar.lat,
                               self.radar.lon, self.radar.alt)
        ti = self.iono.itemp(self.utime[:, 0], self.radar.lat,
                             self.radar.lon, self.radar.alt)
        te = self.iono.etemp(self.utime[:, 0], self.radar.lat,
                             self.radar.lon, self.radar.alt)

        # calculate LoS velocity for each bin by taking the dot product of the
        # radar kvector and the velocity field
        kvec = self.radar.kvec_all_gates()
        Vvec = self.iono.velocity(
            self.utime[:, 0], self.radar.lat, self.radar.lon, self.radar.alt)
        vlos = np.einsum('...i,k...i->k...', kvec, Vvec)

        # calculate ACF density
        ne_notr = self.iono.density(self.utime[:, 0], self.radar.acf_lat,
                                    self.radar.acf_lon, self.radar.acf_alt)
        return ne, ti, te, vlos, ne_notr

    def generate_errors(self, err_coef):

        # Need to make this more rigerous
        ne_err = err_coef[0] * self.radar.slant_range**2
        ti_err = err_coef[1] * self.radar.slant_range**2
        te_err = err_coef[2] * self.radar.slant_range**2
        vlos_err = err_coef[3] * self.radar.slant_range**2
        ne_notr_err = err_coef[0] * self.radar.acf_slant_range**2

        return ne_err, ti_err, te_err, vlos_err, ne_notr_err

    def noisy_measurements(self):
        ne = np.random.normal(loc=self.ne, scale=self.ne_err)
        ti = np.random.normal(loc=self.ti, scale=self.ti_err)
        te = np.random.normal(loc=self.te, scale=self.te_err)
        vlos = np.random.normal(loc=self.vlos, scale=self.vlos_err)
        ne_notr = np.random.normal(loc=self.ne_notr, scale=self.ne_notr_err)

        return ne, ti, te, vlos, ne_notr

    def generate_beamcodes(self):

        beamcodes = np.array([self.radar.beam_codes,
                              self.radar.beam_azimuth,
                              self.radar.beam_elevation,
                              self.radar.beam_ksys]).T
        return beamcodes

    def generate_fitted_params(self):

        FittedParams = {
            'Altitude': self.radar.alt,
            'IonMass': self.iono.ion_mass,
            'Range': self.radar.slant_range}

        # create fit and error arrays that match the shape of whats in the
        # processed fitted files
        # Fit Array: Nrecords x Nbeams x Nranges x Nions+1 x 4
        #   (fraction, temperature, coll. freq., LoS speed)
        # assume only O+, but include fields for other parameters so array is a
        # general shape
        s = (self.utime.shape[0],) + self.radar.slant_range.shape
        FittedParams['Fits'] = np.full(
            s + (len(self.iono.ion_mass) + 1, 4), np.nan)
        FittedParams['Fits'][:, :, :, 0, 1] = self.ti
        FittedParams['Fits'][:, :, :, -1, 1] = self.te
        FittedParams['Fits'][:, :, :, 0, 3] = self.vlos
        FittedParams['Fits'][:, :, :, -1, 3] = self.vlos
        FittedParams['Fits'][:, :, :, :, 0] = np.zeros(
            s + (len(self.iono.ion_mass) + 1,))
        FittedParams['Fits'][:, :, :, 0, 0] = np.ones(s)
        FittedParams['Fits'][:, :, :, -1, 0] = np.ones(s)

        FittedParams['Ne'] = self.ne

        FittedParams['Noise'] = np.full(s + (3,), np.nan)

        FittedParams['dNe'] = np.broadcast_to(self.ne_err, s)
        FittedParams['Errors'] = np.full(
            s + (len(self.iono.ion_mass) + 1, 4), np.nan)
        FittedParams['Errors'][:, :, :, 0, 1] = self.ti_err
        FittedParams['Errors'][:, :, :, -1, 1] = self.te_err
        FittedParams['Errors'][:, :, :, 0, 3] = self.vlos_err
        FittedParams['Errors'][:, :, :, -1, 3] = self.vlos_err

        # fit info
        FitInfo = {
            'chi2': np.full(s, 1.0),
            'dof': np.full(s, 26),
            'fitcode': np.full(s, 1),
            'nfev': np.full(s, 0)
            }

        # calculate density in ACF bins
        NeFromPower = {
            'Altitude': self.radar.acf_alt,
            'Range': self.radar.acf_slant_range
            }

        NeFromPower['Ne_Mod'] = self.ne_notr
        NeFromPower['Ne_NoTr'] = self.ne_notr
        NeFromPower['SNR'] = np.full(self.radar.acf_slant_range.shape, np.nan)
        NeFromPower['dNeFrac'] = self.ne_notr_err / self.ne_notr

        return FittedParams, FitInfo, NeFromPower

    def generate_time(self):
        # Everything below this can be a seperate "create time arrays" function
        # - only useful for output

        Time = {'UnixTime': self.utime}

        Time['Day'] = np.array([[t.day for t in ip] for ip in self.time])
        Time['Month'] = np.array([[t.month for t in ip] for ip in self.time])
        Time['Year'] = np.array([[t.year for t in ip] for ip in self.time])
        Time['doy'] = np.array(
            [[t.timetuple().tm_yday for t in ip] for ip in self.time])
        Time['dtime'] = np.array([
            [(t - dt.datetime(t.year, t.month, t.day, 0, 0, 0)
              ).total_seconds() / 3600. for t in ip]
            for ip in self.time])

        _, mlt0 = self.iono.apex.convert(
                self.radar.site_lat, self.radar.site_lon, 'geo', 'mlt',
                height=self.radar.site_alt / 1000., datetime=self.time[:, 0]
                )
        _, mlt1 = self.iono.apex.convert(
                self.radar.site_lat, self.radar.site_lon, 'geo', 'mlt',
                height=self.radar.site_alt / 1000., datetime=self.time[:, 1]
                )
        Time['MagneticLocalTimeSite'] = np.array([mlt0, mlt1]).T

        return Time

    def generate_geomag(self):
        # generate Geomag array
        # Reqires running IGRF to get all fields
        Geomag = {
            'Latitude': self.radar.lat,
            'Longitude': self.radar.lon,
            'Altitude': self.radar.alt}
        kvec = self.radar.kvec_all_gates()
        Geomag.update(ke=kvec[:, :, 0], kn=kvec[:, :, 1], kz=kvec[:, :, 2])

        return Geomag

    def generate_site(self):
        Site = {
            'Latitude': self.radar.site_lat,
            'Longitude': self.radar.site_lon,
            'Altitude': self.radar.site_alt,
            'Code': 0}
        mlat, mlon = self.iono.apex.geo2apex(
            self.radar.site_lat, self.radar.site_lon,
            height=self.radar.site_alt / 1000.)
        _, mlt = self.iono.apex.convert(
            self.radar.site_lat, self.radar.site_lon, 'geo', 'mlt',
            height=self.radar.site_alt / 1000., datetime=self.time[0, 0])
        Site.update(
            MagneticLatitude=mlat,
            MagneticLongitude=mlon,
            MagneticLocalTimeMidnight=mlt)

        return Site

    def create_hdf5_output(self, outfilename):

        beamcodes = self.generate_beamcodes()

        # create output hdf5 file
        with h5py.File(outfilename, mode='w') as h5:

            h5.create_dataset('BeamCodes', data=beamcodes)

            FittedParams, FitInfo, NeFromPower = self.generate_fitted_params()
            h5.create_group('/FittedParams')
            for k, v in FittedParams.items():
                h5.create_dataset('/FittedParams/{}'.format(k), data=v)

            h5.create_group('/FittedParams/FitInfo')
            for k, v in FitInfo.items():
                h5.create_dataset('/FittedParams/FitInfo/{}'.format(k), data=v)

            h5.create_group('/NeFromPower')
            for k, v in NeFromPower.items():
                h5.create_dataset('/NeFromPower/{}'.format(k), data=v)

            h5.create_group('/Calibration')

            Geomag = self.generate_geomag()
            h5.create_group('/Geomag')
            for k, v in Geomag.items():
                h5.create_dataset('/Geomag/{}'.format(k), data=v)

            # need to call MSIS-E
            h5.create_group('/MSIS')

            h5.create_group('/ProcessingParams')

            Site = self.generate_site()
            h5.create_group('Site')
            for k, v in Site.items():
                h5.create_dataset('/Site/{}'.format(k), data=v)
            h5.create_dataset(
                '/Site/Name',
                data=self.radar.radar_name,
                dtype=h5py.string_dtype(
                    encoding='utf-8',
                    length=len(
                        self.radar.radar_name.encode('utf-8'))))

            Time = self.generate_time()
            h5.create_group('Time')
            for k, v in Time.items():
                h5.create_dataset('/Time/{}'.format(k), data=v)

    def create_summary_plots(self,
                             plot_time=None,
                             plot_beam=None,
                             output_prefix='synthdata_summary_',
                             alt_slices=[100000.,
                                         200000.,
                                         300000.,
                                         400000.,
                                         500000.],
                             slice_xrng=[-500000.,
                                         500000.,
                                         10000.],
                             slice_yrng=[-200000.,
                                         800000.,
                                         10000.],
                             dens_colors={'vmin': 0.,
                                          'vmax': 4.e11,
                                          'cmap': 'viridis'},
                             itemp_colors={'vmin': 0.,
                                           'vmax': 3000.,
                                           'cmap': 'magma'},
                             etemp_colors={'vmin': 0.,
                                           'vmax': 5000.,
                                           'cmap': 'inferno'},
                             vlos_colors={'vmin': -500.,
                                          'vmax': 500.,
                                          'cmap': 'coolwarm'}):

        # optional imports used ONLY for creating summary plots
        # matplotlib and cartopy are not listed in the package requirments
        try:
            import pymap3d as pm
            import matplotlib.pyplot as plt
            import matplotlib.gridspec as gridspec
            import cartopy.crs as ccrs
        except ModuleNotFoundError as e:
            raise ModuleNotFoundError(
                'In order to create summary plots, the optional modules '
                'matplotlib (https://matplotlib.org/) and cartopy '
                '(https://scitools.org.uk/cartopy/docs/latest/) must be '
                'installed.') from e

        # form grid of coordinates for plotting
        # Make customizable in config file
        alt_layers = np.array(alt_slices)
        e, n = np.meshgrid(np.arange(*slice_xrng), np.arange(*slice_yrng))
        glat, glon, galt = pm.enu2geodetic(
            e, n, 0., self.radar.site_lat, self.radar.site_lon, 0.)

        glat = np.tile(
            glat,
            alt_layers.shape).reshape(
            glat.shape +
            alt_layers.shape,
            order='F')
        glon = np.repeat(
            glon,
            alt_layers.shape).reshape(
            glon.shape +
            alt_layers.shape,
            order='C')
        galt = np.broadcast_to(alt_layers, galt.shape + alt_layers.shape)

        tidx = np.argmin(np.abs(
            (plot_time - dt.datetime.utcfromtimestamp(0)).total_seconds() -
            self.utime[:, 0]))

        ne0 = np.squeeze(self.iono.density(
            self.utime[tidx, 0], glat, glon, galt))
        te0 = np.squeeze(self.iono.etemp(
            self.utime[tidx, 0], glat, glon, galt))
        ti0 = np.squeeze(self.iono.itemp(
            self.utime[tidx, 0], glat, glon, galt))
        ve = np.squeeze(self.iono.velocity(
            self.utime[tidx, 0], glat, glon, galt))

        # scaling/rotation of vector to plot in cartopy
        # https://github.com/SciTools/cartopy/issues/1179
        es = ve[:, :, :, 0] / np.cos(glat * np.pi / 180.)
        ns = ve[:, :, :, 1]
        sf = np.sqrt(ve[:, :, :, 0]**2 + ve[:, :, :, 1]**2) / \
            np.sqrt(es**2 + ns**2)
        e = es * sf
        n = ns * sf

        plotting_params = [
            dict(
                synthdata=self.ne,
                param=ne0,
                cparam=dens_colors,
                label=r'Ne (m$^{-3}$)',
                title='Electron Density',
                output=output_prefix + 'ne.png'),
            dict(
                synthdata=self.ti,
                param=ti0,
                cparam=itemp_colors,
                label=r'Ti (K)',
                title='Ion Temperature',
                output=output_prefix + 'ti.png'),
            dict(
                synthdata=self.te,
                param=te0,
                cparam=etemp_colors,
                label=r'Te (K)',
                title='Electron Temperature',
                output=output_prefix + 'te.png'),
            dict(
                synthdata=self.vlos,
                param=[
                    e,
                    n],
                cparam=vlos_colors,
                label=r'Vlos (m/s)',
                title='Plasma Velocity',
                output=output_prefix + 'vlos.png')]

        proj = ccrs.AzimuthalEquidistant(
            central_latitude=self.radar.site_lat,
            central_longitude=self.radar.site_lon)
        Nalt = len(alt_layers)
        wr = [1] * Nalt
        wr.append(Nalt / 2.)
        gs = gridspec.GridSpec(
            2,
            Nalt + 1,
            width_ratios=wr,
            left=0.1,
            right=0.95,
            bottom=0.1,
            top=0.95,
            wspace=0.1,
            hspace=0.1)
        bidx = np.where(self.radar.beam_codes == plot_beam)[0][0]

        for p in plotting_params:

            fig = plt.figure(figsize=(15, 7))
            fig.suptitle(p['title'])

            # Create slice plots
            for j in range(len(alt_layers)):

                ax = fig.add_subplot(gs[0, j], projection=proj)
                ax.coastlines()
                ax.set_title('{} km'.format(alt_layers[j] / 1000.))

                if p['title'] == 'Plasma Velocity':
                    ax.quiver(glon[:, :, j], glat[:, :, j],
                              p['param'][0][:, :, j], p['param'][1][:, :, j],
                              color='blue', transform=ccrs.PlateCarree())

                else:
                    cs = ax.contourf(glon[:, :, j], glat[:, :, j],
                                     p['param'][:, :, j], **p['cparam'],
                                     transform=ccrs.PlateCarree())
                    cs.cmap.set_over('white')
                    cs.cmap.set_under('grey')
                    cs.changed()

                # Add beam positions
                aidx = np.nanargmin(
                    np.abs(self.radar.alt - alt_layers[j]),
                    axis=1)
                aidx0 = np.array([aidx]).T
                slice_lat = np.take_along_axis(self.radar.lat, aidx0, axis=1)
                slice_lon = np.take_along_axis(self.radar.lon, aidx0, axis=1)
                ax.scatter(
                    slice_lon,
                    slice_lat,
                    facecolors='none',
                    edgecolors='k',
                    transform=ccrs.Geodetic())
                ax.scatter(
                    slice_lon[bidx],
                    slice_lat[bidx],
                    facecolors='none',
                    edgecolors='magenta',
                    transform=ccrs.Geodetic())

            # Create RTI
            ax = fig.add_subplot(gs[1, :-1])
            time = self.utime[:, 0].astype('datetime64[s]')
            alt = self.radar.alt[bidx, :]
            c = ax.pcolormesh(time,
                              alt[np.isfinite(alt)],
                              p['synthdata'][:,
                                             bidx,
                                             np.isfinite(alt)].T,
                              **p['cparam'])
            ax.axvline(x=plot_time, color='magenta')
            ax.set_xlabel('Universal Time')
            ax.set_ylabel('Altitude (m)')
            ax.set_title(
                'Beam Number: {:.0f} (az:{:.1f}, el:{:.1f})'.format(
                    self.radar.beam_codes[bidx],
                    self.radar.beam_azimuth[bidx],
                    self.radar.beam_elevation[bidx]))

            # Create 3D FoV plot
            x, y, z = pm.geodetic2enu(
                self.radar.lat, self.radar.lon, self.radar.alt,
                self.radar.site_lat, self.radar.site_lon, self.radar.site_alt)
            fp = np.isfinite(self.radar.alt)

            ax = fig.add_subplot(gs[:, -1], projection='3d')
            c = ax.scatter(x[fp], y[fp], z[fp], c=p['synthdata']
                           [tidx, fp], **p['cparam'])
            # ax.xaxis.set_ticklabels([])
            # ax.yaxis.set_ticklabels([])
            # ax.zaxis.set_ticklabels([])
            ax.set_box_aspect((1, 1, 1))
            ax.set_box_aspect(
                [ub - lb for lb, ub in (getattr(ax, f'get_{a}lim')()
                 for a in 'xyz')])
            fig.colorbar(c, label=p['label'])

            fig.savefig(p['output'])


def main():
    help_string = 'Program to generate synthetic data for multibeam, AMISR-'
    'like incoherent scatter radars.'

    # Build the argument parser tree
    parser = ArgumentParser(description=help_string,
                            formatter_class=RawDescriptionHelpFormatter)
    parser.add_argument(
        'synth_config_file',
        help='Configuration file for synthetic data set.')
    args = vars(parser.parse_args())

    with open(args['synth_config_file'], 'r') as cf:
        config = yaml.safe_load(cf)

    sd = SyntheticData(config)
    sd.create_hdf5_output(config['GENERAL']['output_filename'])
    sd.create_summary_plots(**config['SUMMARY_PLOT'])
