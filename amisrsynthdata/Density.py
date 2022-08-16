# Density.py
import numpy as np
import pymap3d as pm

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

        s = (utime.shape[0],)+galt.shape
        Ne0 = np.full(s, self.value)

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
        zp = (galt-self.z0)/self.H
        Ne = self.N0*np.exp(0.5*(1-zp-np.exp(-zp)/np.cos(self.sza*np.pi/180.)))

        s = (utime.shape[0],)+galt.shape
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
        center_vec = np.array(pm.geodetic2ecef(self.cent_lat, self.cent_lon, 0.))

        # define norm vector and array of point vectors in ECEF
        norm_vec = np.array(pm.aer2ecef(self.az, 0., 1., self.cent_lat, self.cent_lon, 0.))-center_vec
        point_vec = np.moveaxis(np.array(pm.geodetic2ecef(glat, glon, galt)), 0, -1)-center_vec

        # calculate distance between each point and the plane
        r = np.einsum('...i,i->...', point_vec, norm_vec)

        # apply hyperbolic tangent function to create gradient
        Ne = self.N0*(np.tanh(r/self.L)+1)

        s = (utime.shape[0],)+galt.shape
        Ne0 = np.broadcast_to(Ne, s)

        return Ne0

    def tubular_patch(self, utime, glat, glon, galt):
        """
        "Infintite tubular patch" of the variety that has traditionally been used to model polar cap patchs
        in numerical models.  Formed from two hyperbolic tangients in the horizontal direction and a
        Gaussian in the vertical direction.  A non-zero velocity will make the patch move in time.

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

        w = self.width/2.
        h = self.height/2.

        s = (utime.shape[0],)+galt.shape
        Ne0 = np.empty(s)

        for i in range(len(utime)):

            # Progress center point to new location
            t = utime[i,0]-self.utime0
            cent_lat, cent_lon, cent_alt = pm.enu2geodetic(self.velocity[0]*t, self.velocity[1]*t, self.velocity[2]*t, self.orig_lat, self.orig_lon, self.orig_alt)

            # ECEF vector to the center point
            center_vec = np.array(pm.geodetic2ecef(cent_lat, cent_lon, cent_alt))

            # define norm vector and array of point vectors in ECEF
            norm_vec = np.array(pm.aer2ecef(self.az, 0., 1., cent_lat, cent_lon, cent_alt))-center_vec

            # print(norm_vec.shape)
            point_vec = np.moveaxis(np.array(pm.geodetic2ecef(glat, glon, galt)), 0, -1)-center_vec

            # calculate distance between each point and the plane
            r = np.einsum('...i,i->...', point_vec, norm_vec)

            # apply hyperbolic tangent function to create gradient
            Ne = self.N0/2.*(np.tanh((r+w)/self.L)-np.tanh((r-w)/self.L))*np.exp(-0.5*(galt-cent_alt)**2/h**2)
            # Ne = N0/2.*(np.tanh((r+w)/L)-np.tanh((r-w)/L))

            Ne0[i] = Ne

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

        r = self.width/2.
        h = self.height/2.

        e, n, u  = pm.geodetic2enu(glat, glon, galt, self.cent_lat, self.cent_lon, self.cent_alt)
        Ne = self.N0*np.exp(-0.5*(e**2/r**2 + n**2/r**2 + u**2/h**2))

        s = (utime.shape[0],)+galt.shape
        Ne0 = np.broadcast_to(Ne, s)

        return Ne0

    # add wave fluctuation functions - L. Goodwin code
