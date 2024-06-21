import numpy as np
import pymap3d as pm
import datetime as dt
from .utils import output_shape, gemini_helper


class Density(object):
    def __init__(self, type, params, utime0):

        # set density function
        self.Ne_function = getattr(self, type)
        self.utime0 = utime0

        # set parameters as class attributes
        for param, value in params.items():
            setattr(self, param, value)

    def __call__(self, utime, glat, glon, galt):
        return self.Ne_function(utime, glat, glon, galt)

    def uniform(self, utime, glat, glon, galt):
        """
        Uniform electron density at all points.

        Parameters
        ----------
        value: float
            The value to assign at all points (m-3).
        """

        s = output_shape(utime, galt)
        if not s:
            if np.isnan(galt):
                Ne0 = np.nan
            else:
                Ne0 = self.value
        else:
            Ne0 = np.full(s, self.value)
            Ne0[..., np.isnan(galt)] = np.nan

        return Ne0

    def chapman(self, utime, glat, glon, galt):
        """
        Standard Chapman profile.

        z' = (z-z0)/H

        Ne = N0 * exp(0.5(1-z'-exp(-z'))/cos(sza))

        From Schunk and Nagy, 2009; Equation 11.57

        Parameters
        ----------
        N0: float
            Peak electron density (m-3)
        z0: float
            Peak altitude (m)
        H: float
            Scale height (m)
        sza: float
            Solar zenith angle (degrees)
        """

        # From Schunk and Nagy, 2009; eqn 11.57
        zp = (galt - self.z0) / self.H
        Ne = self.N0 * \
            np.exp(0.5 * (1 - zp - np.exp(-zp) / np.cos(
                self.sza * np.pi / 180.)))

        s = output_shape(utime, galt)
        if not s:
            Ne0 = Ne
        else:
            Ne0 = np.broadcast_to(Ne, s)

        return Ne0

    def epstein(self, utime, glat, glon, galt):
        """
        Standard Epstein profile.

        z' = (z-z0)/H

        Ne = N0 / cosh^2(z')

        From Themens et al., 2019; Equation 7

        Parameters
        ----------
        N0: float
            Peak electron density (m-3)
        z0: float
            Peak altitude (m)
        H: float
            Scale height (m)
        """

        # From Themens et al., 2019; eqn 7
        zp = (galt - self.z0) / self.H
        Ne = self.N0 / np.cosh(zp)**2

        s = output_shape(utime, galt)
        if not s:
            Ne0 = Ne
        else:
            Ne0 = np.broadcast_to(Ne, s)

        return Ne0

    def gradient(self, utime, glat, glon, galt):
        """
        Single horizontal gradient crated with a hyperbolic tangent function.

        Parameters
        ----------
        N0: float
            Electron density at high side of gradient (m-3)
        L: float
            Gradient scale length (m)
        cent_lat: float
            Geodetic latitude of point to center the gradient around (deg)
        cent_lon: float
            Geodetic longitude of point to center the gradient around (deg)
        az: float
            Azimuth of gradient normal direction from geodetic north (deg)
        """

        # ECEF vector to the center point
        center_vec = np.array(
            pm.geodetic2ecef(
                self.cent_lat,
                self.cent_lon,
                0.))

        # define norm vector and array of point vectors in ECEF
        norm_vec = np.array(
            pm.aer2ecef(
                self.az,
                0.,
                1.,
                self.cent_lat,
                self.cent_lon,
                0.)) - center_vec
        point_vec = np.moveaxis(
            np.array(
                pm.geodetic2ecef(
                    glat, glon, galt)), 0, -1) - center_vec

        # calculate distance between each point and the plane
        r = np.einsum('...i,i->...', point_vec, norm_vec)

        # apply hyperbolic tangent function to create gradient
        Ne = self.N0 * (np.tanh(r / self.L) + 1)

        # s = (utime.shape[0],)+galt.shape
        # Ne0 = np.broadcast_to(Ne, s)
        #
        s = output_shape(utime, galt)
        if not s:
            Ne0 = Ne
        else:
            Ne0 = np.broadcast_to(Ne, s)

        return Ne0

    def tubular_patch(self, utime, glat, glon, galt):
        """
        "Infintite tubular patch" of the variety that has traditionally been
        used to model polar cap patchs in numerical models.  Formed from two
        hyperbolic tangients in the horizontal direction and a Gaussian in the
        vertical direction.  A non-zero velocity will make the patch move in
        time.

        Parameters
        ----------
        N0: float
            Peak electron density within the patch (m-3)
        L: float
            Gradient scale length (m)
        width: float
            Patch width (m)
        height: float
            Patch height (m)
        cent_lat: float
            Geodetic latitude of point to center the gradient around (deg)
        cent_lon: float
            Geodetic longitude of point to center the gradient around (deg)
        cent_alt: float
            Geodetic altitude of point to center the gradient around (m)
        az: float
            Azimuth of gradient normal direction from geodetic north (deg)
        velocity: list
            Patch velocity [E, N, U] (m/s)
        """

        w = self.width / 2.
        h = self.height / 2.

        s = output_shape(utime, galt)
        if utime.shape:
            # Some fancy array axis manipulation to get t in the proper shape
            # for brodcasting
            t = np.moveaxis(np.broadcast_to(
                utime - self.utime0, s[1:] + (s[0],)), -1, 0)
        else:
            t = utime - self.utime0

        # define array of point vectors in ENU
        e, n, u = pm.geodetic2enu(
            glat, glon, galt, self.cent_lat, self.cent_lon, self.cent_alt)
        point_vec = np.moveaxis(np.array([e, n, u]), 0, -1)

        # define norm vector in ENU
        norm_vec = np.array(pm.aer2enu(self.az, 0., 1.))

        # create wave with sine function
        Vn = np.einsum('i,...i->...', np.array(self.velocity), norm_vec)
        pn = np.einsum('...i,i->...', point_vec, norm_vec)

        Ne0 = self.N0 / 2. * (np.tanh((pn - Vn * t + w) / self.L) - np.tanh(
            (pn - Vn * t - w) / self.L)) * np.exp(
                    -0.5 * (galt - self.cent_alt)**2 / h**2)

        return Ne0

    def circle_patch(self, utime, glat, glon, galt):
        """
        Polar cap patch formed from a 3D Gaussian.

        Parameters
        ----------
        N0: float
            Electron density at patch peak (m-3)
        width: float
            Patch width - Horizontal Gaussian FWHM (m)
        height: float
            Patch height - Horizontal Gaussian FWHM (m)
        cent_lat: float
            Geodetic latitude of center of patch (deg)
        cent_lon: float
            Geodetic longitude of center of patch (deg)
        cent_alt: float
            Geodetic altitude of center of patch (m)
        """

        r = self.width / 2.
        h = self.height / 2.

        e, n, u = pm.geodetic2enu(
            glat, glon, galt, self.cent_lat, self.cent_lon, self.cent_alt)
        Ne = self.N0 * np.exp(-0.5 * (e**2 / r**2 + n**2 / r**2 + u**2 / h**2))

        # s = (utime.shape[0],)+galt.shape
        # Ne0 = np.broadcast_to(Ne, s)
        #
        s = output_shape(utime, galt)
        if not s:
            Ne0 = Ne
        else:
            Ne0 = np.broadcast_to(Ne, s)

        return Ne0

    # add wave fluctuation functions - L. Goodwin code
    def wave(self, utime, glat, glon, galt):
        """
        A wave-like structure, such as traveling ionospheric distrubances
        (TID) or gravity waves (GW). The wave has the functional form:

        Ne = A0 * sin( 2pi * (k*r + t/P) ) * exp( (z-wave_alt)^2 / H^2 )

        This function produces wave purturbations around zero with amplitude
        decaying as H^-2 moving away from the specified peak wave altitude.

        Parameters
        ----------
        A0: float
            Wave amplitude (m^-3)
        k: list
            Wavevector [E, N, U] (m^-1)
        P: float
            Wave period (s)
        orig_lat: float
            Geodetic latitude wavevector is defined from (deg)
        orig_lon: float
            Geodetic longitude wavevector is defined from (deg)
        orig_alt: float
            Geodetic altiude wavevector is defined from (m)
        H: float
            Scale height for wave amplitude decay (m)
        wave_alt: float
            Peak altitude of wave that it should decay away from (m)

        Notes
        -----
        If this function is not used in conjunction with a suitable background,
        it WILL generate negative values for electron density.  If a uniform
        background density is selected, make sure A0 is less than the uniform
        background value.  If a Chapman layer background is desired, the
        easiest way to ensure positive definite density is to follow these
        guidlines:

        N0_chap > A0

        z0_chap = wave_alt

        H_chap > H

        """

        s = output_shape(utime, galt)
        if utime.shape:
            # Some fancy array axis manipulation to get t in the proper shape
            # for brodcasting
            t = np.moveaxis(np.broadcast_to(
                utime - self.utime0, s[1:] + (s[0],)), -1, 0)
        else:
            t = utime - self.utime0

        # define array of point vectors in ECEF
        e, n, u = pm.geodetic2enu(
            glat, glon, galt, self.orig_lat, self.orig_lon, self.orig_alt)
        point_vec = np.moveaxis(np.array([e, n, u]), 0, -1)

        # create wave with sine function
        kr = np.einsum('i,...i->...', np.array(self.k), point_vec)
        Ne0 = self.A0 * np.sin(2 * np.pi * (kr + t / self.P)) * \
            np.exp(-0.5 * (galt - self.wave_alt)**2 / self.H**2)

        return Ne0

    def gemini(self, utime, glat, glon, galt):
        """
        Electron density output from GEMINI model.

        Parameters
        ----------
        gemini_output_dir: string
            Path to directory of GEMINI output files
        """

        gh = gemini_helper(self.gemini_output_dir, glat, glon, galt)

        if not utime.shape:
            Ne0 = gh.query_model(
                dt.datetime.utcfromtimestamp(utime),
                'ne')

        else:
            s = output_shape(utime, galt)
            Ne0 = np.empty(s)
            for i, ut in enumerate(utime):
                Ne0[i] = gh.query_model(
                    dt.datetime.utcfromtimestamp(ut), 'ne')

        return Ne0
