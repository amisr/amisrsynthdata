# Velocity.py
import numpy as np
import pymap3d as pm
from apexpy import Apex

class Velocity(object):
    def __init__(self, utime0, config_params, apex=None):
        # set density function
        self.Vi_function = getattr(self, config_params['TYPE'])
        # set starttime
        self.utime0 = utime0

        # initialize Apex object for all mapping
        self.apex = apex

        # assign remaining config options to parameters to be handled by each function
        config_params.pop('type')
        self.params = config_params


    def __call__(self, utime, glat, glon, galt):
        return self.Vi_function(utime, glat, glon, galt)


    def uniform(self, utime, glat, glon, galt):

        alat, alon = self.apex.geo2apex(glat.ravel(), glon.ravel(), galt.ravel()/1000.)
        map_glat, map_glon, _ = self.apex.apex2geo(alat, alon, 300.)

        V = np.array([float(i) for i in self.params['value'].split(',')])

        # Find ECEF velocity components for given geodetic velocity at center of points
        u, v, w = pm.enu2uvw(V[0], V[1], V[2], np.nanmean(map_glat), np.nanmean(map_glon))
        # Find ENU components for same velosity translated to all mapped locations
        e, n, u = pm.uvw2enu(u, v, w, map_glat, map_glon)
        u = np.zeros(u.shape)   # set up component to zero
        V_map = np.array([e,n,u])
        # rescale velocities so that they all have the same magnitude as the original vector
        V_scale = V_map*np.linalg.norm(V)/np.linalg.norm(V_map, axis=0)
        # map velocities along field lines back to the original heights
        V0 = self.apex.map_V_to_height(alat, alon, 300., galt.ravel()/1000., V_scale)
        # reform original array shape
        VE = V0.T.reshape(galt.shape+(3,))

        s = (utime.shape[0],)+galt.shape+(3,)
        VE0 = np.broadcast_to(VE, s)

        return VE0

    def uniform_mlat_aligned(self, utime, glat, glon, galt):
        Ve1, Ve2, Ve3 = [float(i) for i in self.params['value'].split(',')]

        # Find base vector at each location
        _, _, _, _, _, _, _, _, _, e1, e2, e3 = self.apex.basevectors_apex(glat.ravel(), glon.ravel(), galt.ravel()/1000.)
        # reshape basevector arrays to match the original input
        e1 = e1.T.reshape(glat.shape+(3,))
        e2 = e2.T.reshape(glat.shape+(3,))
        e3 = e3.T.reshape(glat.shape+(3,))

        # calculate V in geodetic coordinates
        VE = Ve1*e1 + Ve2*e2 + Ve3*e3

        s = (utime.shape[0],)+galt.shape+(3,)
        VE0 = np.broadcast_to(VE, s)

        return VE0


    def uniform_glat_aligned(self, utime, glat, glon, galt):
        VE = [float(i) for i in self.params['value'].split(',')]

        s = (utime.shape[0],)+galt.shape+(3,)
        VE0 = np.broadcast_to(VE, s)

        return VE0
