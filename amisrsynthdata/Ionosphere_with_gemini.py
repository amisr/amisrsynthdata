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

        self.utime0 = (dt.datetime.fromisoformat(config['GENERAL']['STARTTIME'])-dt.datetime.utcfromtimestamp(0)).total_seconds()


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



    def time_interp():
        # this morphs between two states
        N1 = self.density1(glat, glon, galt)
        N2 = self.density2(glat, glon, galt)

        # interpolate in time between N1 and N2

    # Interpolating between two states is NOT the same a structre propigating
    # interpolating will just result in the patch slowly fading in one location
    # while reappearing in another
    # easiest/most flexible way to do this is to just build it into the existing state functions
    # all state functions take time and return full array shape
    # if state is static, just repeteat the array as is done in Synthetic Data

    def uniform_velocity(self, utime, glat, glon, galt):

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
        VE = V0.T.reshape(galt.shape+(3,))

        s = (utime.shape[0],)+galt.shape+(3,)
        VE0 = np.broadcast_to(VE, s)

        return VE0

    def uniform_mlat_aligned(self, utime, glat, glon, galt):
        Ve1, Ve2, Ve3 = [float(i) for i in self.velocity_params['value'].split(',')]

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
        VE = [float(i) for i in self.velocity_params['value'].split(',')]

        s = (utime.shape[0],)+galt.shape+(3,)
        VE0 = np.broadcast_to(VE, s)

        return VE0

    def gemini_Vi(self, utime, glat, glon, galt):
        gemini_output_dir = self.velocity_params['gemini_output_dir']

        from gemini3d.grid.gridmodeldata import model2pointsgeogcoords
        import gemini3d.read as read

        cfg = read.config(gemini_output_dir)
        xg = read.grid(gemini_output_dir)

        s = (utime.shape[0],)+galt.shape+(3,)
        Vi0 = np.empty(s)

        for i in range(len(utime)):

            dat = read.frame(gemini_output_dir, dt.datetime.utcfromtimestamp(utime[i,0]), var='v1')
            V1 = model2pointsgeogcoords(xg, dat['v1'], galt, glon, glat)
            dat = read.frame(gemini_output_dir, dt.datetime.utcfromtimestamp(utime[i,0]), var='v2')
            V2 = model2pointsgeogcoords(xg, dat['v2'], galt, glon, glat)
            dat = read.frame(gemini_output_dir, dt.datetime.utcfromtimestamp(utime[i,0]), var='v3')
            V3 = model2pointsgeogcoords(xg, dat['v3'], galt, glon, glat)

# [egalt,eglon,eglat]=unitvecs_geographic(xg)
# #^ returns a set of geographic unit vectors on xg; these are in ECEF geomag comps
# #    like all other unit vectors in xg
#
# # each of the components in models basis projected onto geographic unit vectors
# vgalt=( np.sum(xg["e1"]*egalt,3)*dat["v1"] + np.sum(xg["e2"]*egalt,3)*dat["v2"] +
#     np.sum(xg["e3"]*egalt,3)*dat["v3"] )
# vglat=( np.sum(xg["e1"]*eglat,3)*dat["v1"] + np.sum(xg["e2"]*eglat,3)*dat["v2"] +
#     np.sum(xg["e3"]*eglat,3)*dat["v3"] )
# vglon=( np.sum(xg["e1"]*eglon,3)*dat["v1"] + np.sum(xg["e2"]*eglon,3)*dat["v2"] +
#     np.sum(xg["e3"]*eglon,3)*dat["v3"] )



            V = np.array([V1, V2, V3]).T
            print(galt.shape, V1.shape, Vi0.shape, V.shape)

            Vi0[i] = V.reshape(galt.shape+(3,))
            # Vi0[i,:,1] = V2.reshape(galt.shape)
            # Vi0[i,:,2] = V3.reshape(galt.shape)

        return Vi0



    def chapman(self, utime, glat, glon, galt):
        N0 = float(self.density_params['n0'])
        H = float(self.density_params['h'])
        z0 = float(self.density_params['z0'])
        sza = float(self.density_params['sza'])*np.pi/180.

        # From Schunk and Nagy, 2009; eqn 11.57
        zp = (galt-z0)/H
        Ne = N0*np.exp(0.5*(1-zp-np.exp(-zp)/np.cos(sza)))

        s = (utime.shape[0],)+galt.shape
        Ne0 = np.broadcast_to(Ne, s)

        return Ne0

    def gradient(self, utime, glat, glon, galt):

        cent_lat = float(self.density_params['cent_lat'])
        cent_lon = float(self.density_params['cent_lon'])
        N0 = float(self.density_params['n0'])
        L = float(self.density_params['l'])
        az = float(self.density_params['az'])

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

        lat0 = float(self.density_params['cent_lat'])
        lon0 = float(self.density_params['cent_lon'])
        alt0 = float(self.density_params['cent_alt'])
        N0 = float(self.density_params['n0'])
        L = float(self.density_params['l'])
        w = float(self.density_params['width'])/2.
        az = float(self.density_params['az'])
        h = float(self.density_params['height'])/2.
        V = np.array([float(i) for i in self.density_params['velocity'].split(',')])

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


    def uniform_density(self, utime, glat, glon, galt):
        Ne = float(self.density_params['value'])

        s = (utime.shape[0],)+galt.shape
        Ne0 = np.full(s, Ne)

        return Ne0


    def circle_patch(self, utime, glat, glon, galt):
        # circular gaussian polar cap patch
        cent_lat = float(self.density_params['cent_lat'])
        cent_lon = float(self.density_params['cent_lon'])
        cent_alt = float(self.density_params['cent_alt'])
        N0 = float(self.density_params['n0'])
        r = float(self.density_params['width'])/2.
        h = float(self.density_params['height'])/2.

        e, n, u  = pm.geodetic2enu(glat, glon, galt, cent_lat, cent_lon, cent_alt)
        Ne = N0*np.exp(-0.5*(e**2/r**2 + n**2/r**2 + u**2/h**2))

        s = (utime.shape[0],)+galt.shape
        Ne0 = np.broadcast_to(Ne, s)

        return Ne0

    def gemini(self, utime, glat, glon, galt):

        gemini_output_dir = self.density_params['gemini_output_dir']

        from gemini3d.grid.gridmodeldata import model2pointsgeogcoords
        import gemini3d.read as read

        cfg = read.config(gemini_output_dir)
        xg = read.grid(gemini_output_dir)

        s = (utime.shape[0],)+galt.shape
        Ne0 = np.empty(s)

        for i in range(len(utime)):

            dat = read.frame(gemini_output_dir, dt.datetime.utcfromtimestamp(utime[i,0]), var='ne')

            # Call pygemini to queary gemini results?
            Ne = model2pointsgeogcoords(xg, dat['ne'], galt, glon, glat)

            Ne0[i] = Ne.reshape(galt.shape)

        return Ne0



    # add wave fluctuation functions - L. Goodwin code


    def uniform_Te(self, utime, glat, glon, galt):
        Te = float(self.etemp_params['value'])

        s = (utime.shape[0],)+galt.shape
        Te0 = np.full(s, Te)

        return Te0

    def hypertan_Te(self, utime, glat, glon, galt):
        maxTe = float(self.etemp_params['maxte'])
        scale_height = float(self.etemp_params['scale_height'])

        Te = maxTe*np.tanh(galt/scale_height)

        s = (utime.shape[0],)+galt.shape
        Te0 = np.full(s, Te)

        return Te0

    def gemini_Te(self, utime, glat, glon, galt):

        gemini_output_dir = self.density_params['gemini_output_dir']

        from gemini3d.grid.gridmodeldata import model2pointsgeogcoords
        import gemini3d.read as read

        cfg = read.config(gemini_output_dir)
        xg = read.grid(gemini_output_dir)

        s = (utime.shape[0],)+galt.shape
        Te0 = np.empty(s)

        for i in range(len(utime)):

            dat = read.frame(gemini_output_dir, dt.datetime.utcfromtimestamp(utime[i,0]), var='Te')

            # Call pygemini to queary gemini results?
            Te = model2pointsgeogcoords(xg, dat['Te'], galt, glon, glat)

            Te0[i] = Te.reshape(galt.shape)

        return Te0


    def uniform_Ti(self, utime, glat, glon, galt):
        Ti = float(self.itemp_params['value'])

        s = (utime.shape[0],)+galt.shape
        Ti0 = np.full(s, Ti)

        return Ti0

    def hypertan_Ti(self, utime, glat, glon, galt):
        maxTi = float(self.itemp_params['maxti'])
        scale_height = float(self.itemp_params['scale_height'])

        Ti = maxTi*np.tanh(galt/scale_height)

        s = (utime.shape[0],)+galt.shape
        Ti0 = np.full(s, Ti)

        return Ti0


    def gemini_Ti(self, utime, glat, glon, galt):

        gemini_output_dir = self.density_params['gemini_output_dir']

        from gemini3d.grid.gridmodeldata import model2pointsgeogcoords
        import gemini3d.read as read

        cfg = read.config(gemini_output_dir)
        xg = read.grid(gemini_output_dir)

        s = (utime.shape[0],)+galt.shape
        Ti0 = np.empty(s)

        for i in range(len(utime)):

            dat = read.frame(gemini_output_dir, dt.datetime.utcfromtimestamp(utime[i,0]), var='Ti')

            # Call pygemini to queary gemini results?
            Ti = model2pointsgeogcoords(xg, dat['Ti'], galt, glon, glat)

            Ti0[i] = Ti.reshape(galt.shape)

        return Ti0