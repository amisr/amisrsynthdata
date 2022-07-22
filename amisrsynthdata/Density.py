# Density.py
import numpy as np
import pymap3d as pm

class Density(object):
    def __init__(self, utime0, config_params):
        # set density function
        self.Ne_function = getattr(self, config_params['type'])
        # set starttime
        self.utime0 = utime0

        # assign remaining config options to parameters to be handled by each function
        config_params.pop('type')
        self.params = config_params


    def __call__(self, utime, glat, glon, galt):
        return self.Ne_function(utime, glat, glon, galt)


    def uniform_density(self, utime, glat, glon, galt):
        Ne = float(self.params['value'])

        s = (utime.shape[0],)+galt.shape
        Ne0 = np.full(s, Ne)

        return Ne0


    def chapman(self, utime, glat, glon, galt):
        N0 = float(self.params['N0'])
        H = float(self.params['H'])
        z0 = float(self.params['z0'])
        sza = float(self.params['sza'])*np.pi/180.

        # From Schunk and Nagy, 2009; eqn 11.57
        zp = (galt-z0)/H
        Ne = N0*np.exp(0.5*(1-zp-np.exp(-zp)/np.cos(sza)))

        s = (utime.shape[0],)+galt.shape
        Ne0 = np.broadcast_to(Ne, s)

        return Ne0

    def gradient(self, utime, glat, glon, galt):

        cent_lat = float(self.params['cent_lat'])
        cent_lon = float(self.params['cent_lon'])
        N0 = float(self.params['N0'])
        L = float(self.params['L'])
        az = float(self.params['az'])

        # ECEF vector to the center point
        center_vec = np.array(pm.geodetic2ecef(cent_lat, cent_lon, 0.))

        # define norm vector and array of point vectors in ECEF
        norm_vec = np.array(pm.aer2ecef(az, 0., 1., cent_lat, cent_lon, 0.))-center_vec
        point_vec = np.moveaxis(np.array(pm.geodetic2ecef(glat, glon, galt)), 0, -1)-center_vec

        # calculate distance between each point and the plane
        r = np.einsum('...i,i->...', point_vec, norm_vec)

        # apply hyperbolic tangent function to create gradient
        Ne = N0*(np.tanh(r/L)+1)

        s = (utime.shape[0],)+galt.shape
        Ne0 = np.broadcast_to(Ne, s)

        return Ne0

    def tubular_patch(self, utime, glat, glon, galt):

        lat0 = float(self.params['cent_lat'])
        lon0 = float(self.params['cent_lon'])
        alt0 = float(self.params['cent_alt'])
        N0 = float(self.params['N0'])
        L = float(self.params['L'])
        w = float(self.params['width'])/2.
        az = float(self.params['az'])
        h = float(self.params['height'])/2.
        V = np.array(self.params['velocity'])

        s = (utime.shape[0],)+galt.shape
        Ne0 = np.empty(s)

        for i in range(len(utime)):

            # Progress center point to new location
            t = utime[i,0]-self.utime0
            cent_lat, cent_lon, cent_alt = pm.enu2geodetic(V[0]*t, V[1]*t, V[2]*t, lat0, lon0, alt0)

            # ECEF vector to the center point
            center_vec = np.array(pm.geodetic2ecef(cent_lat, cent_lon, cent_alt))

            # define norm vector and array of point vectors in ECEF
            norm_vec = np.array(pm.aer2ecef(az, 0., 1., cent_lat, cent_lon, cent_alt))-center_vec

            # print(norm_vec.shape)
            point_vec = np.moveaxis(np.array(pm.geodetic2ecef(glat, glon, galt)), 0, -1)-center_vec

            # calculate distance between each point and the plane
            r = np.einsum('...i,i->...', point_vec, norm_vec)

            # apply hyperbolic tangent function to create gradient
            Ne = N0/2.*(np.tanh((r+w)/L)-np.tanh((r-w)/L))*np.exp(-0.5*(galt-cent_alt)**2/h**2)
            # Ne = N0/2.*(np.tanh((r+w)/L)-np.tanh((r-w)/L))

            Ne0[i] = Ne

        return Ne0


    def circle_patch(self, utime, glat, glon, galt):
        # circular gaussian polar cap patch
        cent_lat = float(self.params['cent_lat'])
        cent_lon = float(self.params['cent_lon'])
        cent_alt = float(self.params['cent_alt'])
        N0 = float(self.params['n0'])
        r = float(self.params['width'])/2.
        h = float(self.params['height'])/2.

        e, n, u  = pm.geodetic2enu(glat, glon, galt, cent_lat, cent_lon, cent_alt)
        Ne = N0*np.exp(-0.5*(e**2/r**2 + n**2/r**2 + u**2/h**2))

        s = (utime.shape[0],)+galt.shape
        Ne0 = np.broadcast_to(Ne, s)

        return Ne0

    # add wave fluctuation functions - L. Goodwin code
