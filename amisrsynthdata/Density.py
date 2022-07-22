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

        s = (utime.shape[0],)+galt.shape
        Ne0 = np.full(s, self.value)

        return Ne0


    def chapman(self, utime, glat, glon, galt):

        # From Schunk and Nagy, 2009; eqn 11.57
        zp = (galt-self.z0)/self.H
        Ne = self.N0*np.exp(0.5*(1-zp-np.exp(-zp)/np.cos(self.sza*np.pi/180.)))

        s = (utime.shape[0],)+galt.shape
        Ne0 = np.broadcast_to(Ne, s)

        return Ne0

    def gradient(self, utime, glat, glon, galt):

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

        r = self.width/2.
        h = self.height/2.

        e, n, u  = pm.geodetic2enu(glat, glon, galt, self.cent_lat, self.cent_lon, self.cent_alt)
        Ne = self.N0*np.exp(-0.5*(e**2/r**2 + n**2/r**2 + u**2/h**2))

        s = (utime.shape[0],)+galt.shape
        Ne0 = np.broadcast_to(Ne, s)

        return Ne0

    # add wave fluctuation functions - L. Goodwin code
