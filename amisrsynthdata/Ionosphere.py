# Ionosphere.py
import numpy as np
import datetime as dt
from apexpy import Apex
import pymap3d as pm


try:
    import ConfigParser as configparser
except ImportError:
    import configparser

# Figure out how to specify multiple species sensibly
# this can probably wait as a feature - unlikely to come up often
# may be nice to be able to have some kind of placeholder just to make the arrays the full shape though

class Ionosphere(object):

    def __init__(self, *args, **kwargs):

        if len(args) == 1:
            self.read_config(args[0])
        else:
            self.apex_year = kwargs['apex_year']
            # self.field_coords = np.array(kwargs['field_coords'])
            # self.field_values = np.array(kwargs['field_values'])

        # initialize Apex object
        self.apex = Apex(date=self.apex_year)


    def read_config(self, config_file):

        config = configparser.ConfigParser()
        config.read(config_file)

        self.density = getattr(self, config['DENSITY']['TYPE'])
        self.density_params = dict(config['DENSITY'])
        self.density_params.pop('type')

        self.velocity = getattr(self, config['VELOCITY']['TYPE'])
        self.velocity_params = dict(config['VELOCITY'])
        self.velocity_params.pop('type')

        self.etemp = getattr(self, config['ETEMP']['TYPE'])
        self.etemp_params = dict(config['ETEMP'])
        self.etemp_params.pop('type')

        self.itemp = getattr(self, config['ITEMP']['TYPE'])
        self.itemp_params = dict(config['ITEMP'])
        self.itemp_params.pop('type')

        self.apex_year = dt.datetime.fromisoformat(config['GENERAL']['STARTTIME'])
        self.ion_mass = np.array([float(im) for im in config['GENERAL']['ION_MASS'].split(',')])
        # self.field_coords = np.array(eval(config.get('FIELD', 'field_coords')))
        # self.field_values = np.array(eval(config.get('FIELD', 'field_values')))



    # add different kinds of velocity fields

    def uniform_velocity(self, glat, glon, galt):

        alat, alon = self.apex.geo2apex(glat.ravel(), glon.ravel(), galt.ravel()/1000.)
        map_glat, map_glon, _ = self.apex.apex2geo(alat, alon, 300.)

        V = np.array([float(i) for i in self.velocity_params['value'].split(',')])

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
        V0 = V0.T.reshape(galt.shape+(3,))

        return V0

    def uniform_mlat_aligned(self, glat, glon, galt):
        Ve1, Ve2, Ve3 = [float(i) for i in self.velocity_params['value'].split(',')]

        # Find base vector at each location
        _, _, _, _, _, _, _, _, _, e1, e2, e3 = self.apex.basevectors_apex(glat.ravel(), glon.ravel(), galt.ravel()/1000.)
        # reshape basevector arrays to match the original input
        e1 = e1.T.reshape(glat.shape+(3,))
        e2 = e2.T.reshape(glat.shape+(3,))
        e3 = e3.T.reshape(glat.shape+(3,))

        # calculate V in geodetic coordinates
        V0 = Ve1*e1 + Ve2*e2 + Ve3*e3

        return V0

    def chapman(self, glat, glon, galt):
        N0 = float(self.density_params['n0'])
        H = float(self.density_params['h'])
        z0 = float(self.density_params['z0'])
        sza = float(self.density_params['sza'])*np.pi/180.

        # From Schunk and Nagy, 2009; eqn 11.57
        zp = (galt-z0)/H
        Ne = N0*np.exp(0.5*(1-zp-np.exp(-zp)/np.cos(sza)))

        return Ne


    def uniform_density(self, glat, glon, galt):
        Ne = float(self.density_params['value'])
        return np.full(galt.shape, Ne)


    def circle_patch(self, glat, glon, galt):
        # circular gaussian polar cap patch
        cent_lat = float(self.density_params['cent_lat'])
        cent_lon = float(self.density_params['cent_lon'])
        cent_alt = float(self.density_params['cent_alt'])
        N0 = float(self.density_params['n0'])
        r = float(self.density_params['width'])/2.
        h = float(self.density_params['height'])/2.

        e, n, u  = pm.geodetic2enu(glat, glon, galt, cent_lat, cent_lon, cent_alt)
        Ne = N0*np.exp(-0.5*(e**2/r**2 + n**2/r**2 + u**2/h**2))

        return Ne

    # add wave fluctuation functions - L. Goodwin code


    def uniform_Te(self, glat, glon, galt):
        Te = float(self.etemp_params['value'])
        return np.full(galt.shape, Te)

    def uniform_Ti(self, glat, glon, galt):
        Ti = float(self.itemp_params['value'])
        return np.full(galt.shape, Ti)
