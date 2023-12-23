# syntheticdata.py

from .ionosphere import Ionosphere
from .radar import Radar

import numpy as np
import datetime as dt
import h5py
from argparse import ArgumentParser, RawDescriptionHelpFormatter
import yaml


class SyntheticData(object):
    """
    Class for creating synthetic amisr data based on config options.

    Parameters
    ----------
    config : :obj:`yaml`
        Configuration parameters
    """

    def __init__(self, config):

        starttime = config['GENERAL']['starttime']  # init input
        endtime = config['GENERAL']['endtime']      # init input
        rel_err = config['GENERAL']['rel_err']
        err_ref_rng = config['GENERAL']['err_ref_rng']
        self.include_noise = config['GENERAL']['noise']

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
        ) = self.generate_errors(rel_err, err_ref_rng)

    def generate_time_array(self, starttime, endtime):
        """
        Generate time arrays.

        Parameters
        ----------
        starttime : :obj:`datetime.datetime`
            Start time for synthetic data
        endtime : :obj:`datetime.datetime`
            End time for synthetic data

        Returns
        -------
        utime : np.ndarray of floats
            Array of unix timestamps
        time : np.ndarray of :obj:`datetime.datetime`
            Array of datetime.datetime objects
        """
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
        """
        Generate the basic radar measurmeents at each point in the
        field-of-view.  These include:
        * Elecron density (Ne)
        * Ion temperature (Ti)
        * Electron temperature (Te)
        * Line-of-Site plasma velocity (Vlos)

        Returns
        -------
        ne : np.ndarray
            Electron density at the fitted range gates
        ti : np.ndarray
            Ion temperature at the fitted range gates
        te : np.ndarray
            Electron temperature at the fitted range gates
        vlos : np.ndarray
            Line-of-Site velocity at the fitted range gates
        ne_notr : np.ndarray
            Electron density at the ACF range gates
        """

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

    def generate_errors(self, rel_err, err_ref_rng):
        """
        Generate error values associated with each radar measurment. This is
        a simple approximation where errors are proportional to the distance
        to the radar squared.

            err = C*r^2

        Here, C is `rel_err` and r is the radar slant range divided by
        `err_ref_rng`, such that the ouput will have the specified relative
        error at the specified reference range.

        Parameters
        ----------
        rel_err : float
            Relative error
        err_ref_rng: float
            Error reference range, or the range at which the output
            error equals the specified `rel_err`

        Returns
        -------
        ne_err : np.ndarray
            Electron density error at the fitted range gates
        ti_err : np.ndarray
            Ion temperatrue error at the fitted range gates
        te_err : np.ndarray
            Electron temperatrue error at the fitted range gates
        vlos_err : np.ndarray
            LoS velocity error at the fitted range gates
        ne_notr_err : np.ndarray
            Electron density error at the ACF range gates
        """

        # Need to make this more rigerous
        r = self.radar.slant_range / err_ref_rng
        ne_err = rel_err * r**2 * self.ne
        ti_err = rel_err * r**2 * self.ti
        te_err = rel_err * r**2 * self.te
        vlos_err = rel_err * r**2 * np.abs(self.vlos)
        r_sr = self.radar.acf_slant_range / err_ref_rng
        ne_notr_err = rel_err * r_sr**2 * self.ne_notr

        return ne_err, ti_err, te_err, vlos_err, ne_notr_err

    def noisy_measurements(self):
        """
        Return the radar measurements with gaussian random noise proportional
        to their error measurements.

        Returns
        -------
        ne : np.ndarray
            Noisy electron density at the fitted range gates
        ti : np.ndarray
            Noisy ion temperature at the fitted range gates
        te : np.ndarray
            Noisy electron temperature at the fitted range gates
        vlos : np.ndarray
            Noisy line-of-Site velocity at the fitted range gates
        ne_notr : np.ndarray
            Noisy electron density at the ACF range gates
        """

        ne = np.random.normal(loc=self.ne, scale=self.ne_err)
        ti = np.random.normal(loc=self.ti, scale=self.ti_err)
        te = np.random.normal(loc=self.te, scale=self.te_err)
        vlos = np.random.normal(loc=self.vlos, scale=self.vlos_err)
        ne_notr = np.random.normal(loc=self.ne_notr, scale=self.ne_notr_err)

        return ne, ti, te, vlos, ne_notr

    def generate_beamcodes(self):
        """
        Generate beamcode array.  This is a Nbeam x 4 array where the four
        columns are beamcode, azimuth, elevation, and ksys.

        Returns
        -------
        beamcodes : np.ndarray
            Radar beamcode array
        """

        beamcodes = np.array([self.radar.beam_codes,
                              self.radar.beam_azimuth,
                              self.radar.beam_elevation,
                              self.radar.beam_ksys]).T
        return beamcodes

    def generate_fitted_params(self):
        """
        Generate FittedParams, FitInfo, and NeFromPower output dictionaries.
        These are the main outputs that contain actual measurements and
        information about the measurements.

        Returns
        -------
        FittedParams : dict
            Main fitted meausrement output
        FitInfo : dict
            Information about the quality of fit
        NeFromPower : dict
            Measruements from the ACF range gates (electron density only)
        """

        FittedParams = {
            'Altitude': self.radar.alt,
            'IonMass': self.iono.ion_mass,
            'Range': self.radar.slant_range}

        # Add noise to measured parameters
        if self.include_noise:
            ne, ti, te, vlos, ne_notr = self.noisy_measurements()
        else:
            ne, ti, te, vlos, ne_notr = (
                    self.ne, self.ti, self.te, self.vlos, self.ne_notr)

        # create fit and error arrays that match the shape of whats in the
        # processed fitted files
        # Fit Array: Nrecords x Nbeams x Nranges x Nions+1 x 4
        #   (fraction, temperature, coll. freq., LoS speed)
        # assume only O+, but include fields for other parameters so array is a
        # general shape
        s = (self.utime.shape[0],) + self.radar.slant_range.shape
        FittedParams['Fits'] = np.full(
            s + (len(self.iono.ion_mass) + 1, 4), np.nan)
        FittedParams['Fits'][:, :, :, 0, 1] = ti
        FittedParams['Fits'][:, :, :, -1, 1] = te
        FittedParams['Fits'][:, :, :, 0, 3] = vlos
        FittedParams['Fits'][:, :, :, -1, 3] = vlos
        FittedParams['Fits'][:, :, :, :, 0] = np.zeros(
            s + (len(self.iono.ion_mass) + 1,))
        FittedParams['Fits'][:, :, :, 0, 0] = np.ones(s)
        FittedParams['Fits'][:, :, :, -1, 0] = np.ones(s)

        FittedParams['Ne'] = ne

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

        NeFromPower['Ne_Mod'] = ne_notr
        NeFromPower['Ne_NoTr'] = ne_notr
        NeFromPower['SNR'] = np.full(self.radar.acf_slant_range.shape, np.nan)
        NeFromPower['dNeFrac'] = self.ne_notr_err / ne_notr

        return FittedParams, FitInfo, NeFromPower

    def generate_time(self):
        """
        Generate dictionary of the various time formats included in the
        standard output.

        Returns
        -------
        Time : dict
           Various formats for output timestamps
        """

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
        """
        Generate dictionary with fitted range gate coordinates and information.

        Returns
        -------
        Geomag : dict
            Coordinates of fitted range gates
        """

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
        """
        Generate Site dictionary containing information about the site
        coordinates.

        Returns
        -------
        Site : dict
            Coordinates for the site in various systems
        """

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
        """
        Generate output parameter arrays and save them to an hdf5 file.

        Parameters
        ----------
        outfilename : str
            Name of output file
        """

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
        """
        Create basic summary plots of the synthetic dataset.  These are
        thumbnail plots intended to confirm you set up your dataset correctly
        rather than for any kind of analysis.  They include altitude slices
        of the ionosphere at different altitudes, a RTI plot of a particular
        beam, and a 3D plot of the full radar FoV.  A seperate png will be
        created for each of the four ISR parameters.

        Parameters
        ----------
        plot_time : :obj:`datetime.datetime`
            Time to plot in the aliude slice and 3D FoV plots
        plot_beam : float
            Beamcode of beam to plot for the RTI
        output_prefix : str
            Start of the name of the saved output plots, including path
        alt_slices : list of floats
            Altitudes at which to create the slice plots
        slice_xrng : list of floats
            Definition of the horizontal grid in the x-drection
            (start, stop, step)
        slice_yrng : list of floats
            Definition of the horizontal grid in the y-drection
            (start, stop, step)
        dens_colors : dict
            Dictionary listing plotting parameters for electron density,
            specifically ``vmin``, ``vmax``, and ``cmap``
        itemp_colors : dict
            Dictionary listing plotting parameters for ion temperature,
            specifically ``vmin``, ``vmax``, and ``cmap``
        etemp_colors : dict
            Dictionary listing plotting parameters for electron temperature,
            specifically ``vmin``, ``vmax``, and ``cmap``
        vlos_colors : dict
            Dictionary listing plotting parameters for line-of-sight velocity,
            specifically ``vmin``, ``vmax``, and ``cmap``
        """

        # optional imports used ONLY for creating summary plots
        # matplotlib and cartopy are not listed in the package requirments
        try:
            import pymap3d as pm
            import matplotlib as mpl
            import matplotlib.pyplot as plt
            import matplotlib.gridspec as gridspec
            import matplotlib.dates as mdates
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

        if self.include_noise:
            ne, ti, te, vlos, _ = self.noisy_measurements()
        else:
            ne, ti, te, vlos = (self.ne, self.ti, self.te, self.vlos)

        plotting_params = [
            dict(
                synthdata=ne,
                param=ne0,
                cparam=dens_colors,
                label=r'Ne (m$^{-3}$)',
                title='Electron Density',
                output=output_prefix + 'ne.png'),
            dict(
                synthdata=ti,
                param=ti0,
                cparam=itemp_colors,
                label=r'Ti (K)',
                title='Ion Temperature',
                output=output_prefix + 'ti.png'),
            dict(
                synthdata=te,
                param=te0,
                cparam=etemp_colors,
                label=r'Te (K)',
                title='Electron Temperature',
                output=output_prefix + 'te.png'),
            dict(
                synthdata=vlos,
                param=[e, n],
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
            fig.suptitle(p['title'], fontsize=20, fontweight=3)

            cmap = mpl.cm.get_cmap(p['cparam']['cmap']).copy()
            cmap.set_over('white')
            cmap.set_under('grey')

            norm = mpl.colors.Normalize(vmin=p['cparam']['vmin'],
                                        vmax=p['cparam']['vmax'])

            # Create slice plots
            for j in range(len(alt_layers)):

                ax = fig.add_subplot(gs[0, j], projection=proj)
                ax.coastlines()
                ax.set_title('{} km'.format(alt_layers[j] / 1000.))

                if p['title'] == 'Plasma Velocity':
                    # Only plot a subset of the vector grid to keep the plot
                    #  readable
                    s = [int(N/10)+1 for N in glon[:, :, j].shape]
                    q = ax.quiver(glon[::s[0], ::s[1], j],
                                  glat[::s[0], ::s[1], j],
                                  p['param'][0][::s[0], ::s[1], j],
                                  p['param'][1][::s[0], ::s[1], j],
                                  color='blue', transform=ccrs.PlateCarree())
                    if gs[0, j].is_first_col():
                        u = p['cparam']['vmax']
                        ax.quiverkey(q, 0.1, -0.1, u, f'{u} m/s', labelpos='E')

                else:
                    ax.contourf(glon[:, :, j], glat[:, :, j],
                                p['param'][:, :, j], cmap=cmap, norm=norm,
                                transform=ccrs.PlateCarree())

                # Add site location
                ax.scatter(self.radar.site_lon, self.radar.site_lat,
                           marker='^', color='k',
                           transform=ccrs.PlateCarree())

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
            alt = self.radar.alt[bidx, :] / 1000.
            c = ax.pcolormesh(time,
                              alt[np.isfinite(alt)],
                              p['synthdata'][:, bidx,
                                             np.isfinite(alt)].T,
                              cmap=cmap, norm=norm)
            ax.axvline(x=plot_time, color='magenta')
            ax.text(0.0, -0.15, np.datetime_as_string(time[0], unit='D'),
                    transform=ax.transAxes)
            ax.xaxis.set_major_formatter(mdates.DateFormatter('%H:%M'))
            ax.set_xlabel('Universal Time')
            ax.set_ylabel('Altitude (km)')
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
            c = ax.scatter(x[fp] / 1000., y[fp] / 1000., z[fp] / 1000.,
                           c=p['synthdata'][tidx, fp], cmap=cmap, norm=norm)
            ax.scatter(0., 0., 0., marker='^', color='k')
            ax.set_xlabel('East (km)')
            ax.set_ylabel('North (km)')
            ax.set_zlabel('Altitude (km)')
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
