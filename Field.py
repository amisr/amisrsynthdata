# Field.py
import numpy as np
from apexpy import Apex
from scipy import interpolate
import pymap3d as pm
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

try:
    import ConfigParser as configparser
except ImportError:
    import configparser

# Figure out how to specify multiple species sensibly
# this can probably wait as a feature - unlikely to come up often
# may be nice to be able to have some kind of placeholder just to make the arrays the full shape though

class Field(object):

    def __init__(self, *args, **kwargs):

        if len(args) == 1:
            self.read_config(args[0])
        else:
            self.apex_year = kwargs['apex_year']
            self.field_coords = np.array(kwargs['field_coords'])
            self.field_values = np.array(kwargs['field_values'])

        # initialize Apex object
        self.apex = Apex(date=self.apex_year)

        self.map_velocity_field(self.field_coords, self.field_values)
        self.convert_to_ECEF()
        self.create_interpolators()

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

        self.apex_year = config.getint('FIELD', 'apex_year')
        self.field_coords = np.array(eval(config.get('FIELD', 'field_coords')))
        self.field_values = np.array(eval(config.get('FIELD', 'field_values')))


    def map_velocity_field(self, coords, field):
        # coords - array (N,3) of geodetic lat, lon, alt
        # field - array (N,3) of geodetic E, N, U components of the velocity field at each position

        # define output altitudes
        altitude = np.arange(50., 1000., 50.)
        # create array in proper shape to be applied to every input coordinate
        self.altitude = np.repeat(altitude,coords.shape[-1])

        # map to diffent altitudes manually - the current expected input/output arrays of apexpy.map_to_height makes this function difficult to use for this purpose
        alat, alon = self.apex.geo2apex(coords[0], coords[1], coords[2])
        # find positions at each altitude
        self.latitude, self.longitude, __ = self.apex.apex2geo(np.tile(alat,len(altitude)), np.tile(alon,len(altitude)), self.altitude)

        # map field to each altitude
        f = np.array([self.apex.map_V_to_height(alat, alon, coords[2], a, field.T).T for a in altitude])
        self.field = f.reshape(-1,f.shape[-1])


    def convert_to_ECEF(self):

        self.X, self.Y, self.Z = pm.geodetic2ecef(self.latitude, self.longitude, self.altitude*1000.)
        self.Vx, self.Vy, self.Vz = pm.enu2uvw(self.field[:,0], self.field[:,1], self.field[:,2], self.latitude, self.longitude)

    def create_interpolators(self):

        self.interpVx = interpolate.LinearNDInterpolator(np.array([self.X, self.Y, self.Z]).T, self.Vx)
        self.interpVy = interpolate.LinearNDInterpolator(np.array([self.X, self.Y, self.Z]).T, self.Vy)
        self.interpVz = interpolate.LinearNDInterpolator(np.array([self.X, self.Y, self.Z]).T, self.Vz)

    # add different kinds of velocity fields

    def uniform_velocity(self, x, y, z):
        glat, glon, galt = pm.ecef2geodetic(x.flatten(), y.flatten(), z.flatten())
        alat, alon = self.apex.geo2apex(glat, glon, galt/1000.)
        map_glat, map_glon, _ = self.apex.apex2geo(alat, alon, 300.)

        V = np.array([float(i) for i in self.velocity_params['value'].split(',')])

        # V = [500., 700., 0.]
        # Find ECEF velocity components for given geodetic velocity at center of points
        u, v, w = pm.enu2uvw(V[0], V[1], V[2], np.nanmean(map_glat), np.nanmean(map_glon))
        # Find ENU components for same velosity translated to all mapped locations
        e, n, u = pm.uvw2enu(u, v, w, map_glat, map_glon)
        u = np.zeros(u.shape)   # set up component to zero
        V_map = np.array([e,n,u])
        # rescale velocities so that they all have the same magnitude as the original vector
        V_scale = V_map*np.linalg.norm(V)/np.linalg.norm(V_map, axis=0)
        # map velocities along field lines back to the original heights
        V0 = self.apex.map_V_to_height(alat, alon, 300., galt/1000., V_scale)
        # convert to ECEF components
        u, v, w = pm.enu2uvw(V0[0], V0[1], V0[2], glat, glon)
        # reform original array shape
        V0 = np.array([u, v, w]).T.reshape(x.shape+(3,))
        # return np.tile([500,0,0], x.shape+(1,))
        return V0

    def chapman(self, x, y, z):
        N0 = float(self.density_params['n0'])
        H = float(self.density_params['h'])
        z0 = float(self.density_params['z0'])
        glat, glon, galt = pm.ecef2geodetic(x, y, z)
        return N0*np.exp(1-(galt/1000.-z0)/H-np.exp(-(galt/1000.-z0)/H))

        # L. Goodwin's Chapman Layer function - ask about this
          # zprime = ((altdata[i,j]-zm0)/scaleheight)
          # y0 =nem0*math.exp(0.5*(1-zprime-abs(1/math.cos(SolarZenData[t,i,j]*math.pi/180))*math.exp(-1*zprime)))


    def uniform_density(self, x, y, z):
        Ne = float(self.density_params['value'])
        return np.full(x.shape, Ne)

    # add polar cap patch function

    # add wave fluctuation functions - L. Goodwin code


    def uniform_Te(self, x, y, z):
        Te = float(self.etemp_params['value'])
        return np.full(x.shape, Te)

    def uniform_Ti(self, x, y, z):
        Ti = float(self.itemp_params['value'])
        return np.full(x.shape, Ti)

    # def plot_ionosphere(self):
    #
    #     fig = plt.figure(figsize=(10,10))
    #     ax = fig.add_subplot(111,projection='3d')
    #
    #     for x, y, z, vx, vy, vz in zip(self.X, self.Y, self.Z, self.Vx, self.Vy, self.Vz):
    #         ax.quiver(x, y, z, vx, vy, vz, length=0.4*np.sqrt(vx**2+vy**2+vz**2), color='green')
    #
    #     plt.show()
