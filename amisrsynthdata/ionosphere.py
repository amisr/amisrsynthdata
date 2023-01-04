# Ionosphere.py
import numpy as np
import datetime as dt
import pymap3d as pm
from apexpy import Apex
try:
    from gemini3d.grid.gridmodeldata import model2pointsgeogcoords
    import gemini3d.read as read
except ImportError:
    print('WARNING: pygemini is not installed.  GEMINI functionality will not be available.')





# Figure out how to specify multiple species sensibly
# this can probably wait as a feature - unlikely to come up often
# may be nice to be able to have some kind of placeholder just to make the arrays the full shape though
def output_shape(ut, x):
    # determine the appropriate output shape for the given time and position inputs
    if (np.isscalar(ut) and np.isscalar(x)):
        s = None
    elif np.isscalar(ut):
        s = x.shape
    elif np.isscalar(x):
        s = ut.shape
    else:
        s = ut.shape+x.shape
    return s



class Ionosphere(object):

    def __init__(self, config):

        starttime = config['GENERAL']['starttime']
        start_utime = (starttime-dt.datetime.utcfromtimestamp(0)).total_seconds()

        # initialize Apex object
        self.apex = Apex(date=starttime)

        self.ion_mass = config['GENERAL']['ion_mass']

        # create lists of all the functions that will be used for each ionospheric state
        self.density_functions = [Density(type, params, start_utime) for type, params in config['DENSITY'].items()]
        self.velocity_functions = [Velocity(type, params, start_utime, apex=self.apex) for type, params in config['VELOCITY'].items()]
        self.etemp_functions = [Temperature(type, params, start_utime) for type, params in config['ETEMP'].items()]
        self.itemp_functions = [Temperature(type, params, start_utime) for type, params in config['ITEMP'].items()]

    def zero_array(self, ut, x, vec=False):
        s = output_shape(ut, x)
        if vec:
            if not s:
                s = (3,)
            else:
                s = s + (3,)
            return np.zeros(s)
        if not s:
            return 0.0
        else:
            return np.zeros(s)


    def density(self, utime, glat, glon, galt):
        # call each density instance and sum the results
        dens = self.zero_array(utime, galt)
        # dens = np.zeros((utime.shape[0],)+galt.shape)
        for fun in self.density_functions:
            dens = dens + fun(utime, glat, glon, galt)
        return dens

    def velocity(self, utime, glat, glon, galt):
        # vel = np.zeros((utime.shape[0],)+galt.shape+(3,))
        vel = self.zero_array(utime, galt, vec=True)
        for fun in self.velocity_functions:
            vel = vel + fun(utime, glat, glon, galt)
        return vel

    def etemp(self, utime, glat, glon, galt):
        # call each electron temperature instance and sum the results
        # etemp = np.zeros((utime.shape[0],)+galt.shape)
        etemp = self.zero_array(utime, galt)
        for fun in self.etemp_functions:
            etemp = etemp + fun(utime, glat, glon, galt)
        return etemp

    def itemp(self, utime, glat, glon, galt):
        # call each ion temperature instance and sum the results
        # itemp = np.zeros((utime.shape[0],)+galt.shape)
        itemp = self.zero_array(utime, galt)
        for fun in self.itemp_functions:
            itemp = itemp + fun(utime, glat, glon, galt)
        return itemp



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
            Ne0 = self.value
        else:
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
        center_vec = np.array(pm.geodetic2ecef(self.cent_lat, self.cent_lon, 0.))

        # define norm vector and array of point vectors in ECEF
        norm_vec = np.array(pm.aer2ecef(self.az, 0., 1., self.cent_lat, self.cent_lon, 0.))-center_vec
        point_vec = np.moveaxis(np.array(pm.geodetic2ecef(glat, glon, galt)), 0, -1)-center_vec

        # calculate distance between each point and the plane
        r = np.einsum('...i,i->...', point_vec, norm_vec)

        # apply hyperbolic tangent function to create gradient
        Ne = self.N0*(np.tanh(r/self.L)+1)

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

    def gemini(self, utime, glat, glon, galt):
        """
        Density output from GEMINI model.

        Parameters
        ----------
        gemini_output_dir: string
            Path to directory of GEMINI output files
        """

        xg = read.grid(self.gemini_output_dir)

        if not utime.shape:
            dat = read.frame(self.gemini_output_dir, dt.datetime.utcfromtimestamp(utime), var='ne')
            Ne = model2pointsgeogcoords(xg, dat['ne'], galt, glon, glat)
            Ne0 = Ne.reshape(galt.shape)

        else:
            s = output_shape(utime, galt)
            Ne0 = np.empty(s)
            for i, ut in enumerate(utime):
                dat = read.frame(self.gemini_output_dir, dt.datetime.utcfromtimestamp(ut), var='ne')
                Ne = model2pointsgeogcoords(xg, dat['ne'], galt, glon, glat)
                Ne0[i] = Ne.reshape(galt.shape)

        print(Ne0.shape)
        return Ne0



# Because ion and electron temperature functions tend to be similar, they are both caputred
#   within this function.  Ions and Electrons DO NOT have to use the same background function
#   or parameters, these are specified seperately in the config file.
class Temperature(object):
    def __init__(self, type, params, utime0):

        # set density function
        self.Ts_function = getattr(self, type)
        # set starttime
        self.utime0 = utime0

        # set parameters as class attributes
        for param, value in params.items():
            setattr(self, param, value)


    def __call__(self, utime, glat, glon, galt):
        return self.Ts_function(utime, glat, glon, galt)


    def uniform(self, utime, glat, glon, galt):
        """
        Uniform temperature at all points.

        Parameters
        ----------
        value: float
            The value to assign at all points (K)
        """

        # s = (utime.shape[0],)+galt.shape
        # Ts0 = np.full(s, self.value)

        s = output_shape(utime, galt)
        if not s:
            Ts0 = self.value
        else:
            Ts0 = np.full(s, self.value)

        return Ts0

    def hypertan(self, utime, glat, glon, galt):
        """
        Temperature increases in altitude as a hyperbolic tangent.

        Parameters
        ----------
        maxtemp: float
            The temperature at infinity to asymptope to (K)
        scale_height: float
            Vertical scale height (m)
        """

        Ts = self.maxtemp*np.tanh(galt/self.scale_height)

        # s = (utime.shape[0],)+galt.shape
        # Ts0 = np.full(s, Ts)

        s = output_shape(utime, galt)
        if not s:
            Ts0 = Ts
        else:
            Ts0 = np.broadcast_to(Ts, s)

        return Ts0


    def gemini(self, utime, glat, glon, galt):
        """
        Temperature output from GEMINI model.

        Parameters
        ----------
        gemini_output_dir: string
            Path to directory of GEMINI output files
        species: string
            Which species (ion or electron) should be read from GEMINI output ('Te' or 'Ti')
        """

        xg = read.grid(self.gemini_output_dir)

        if not utime.shape:
            dat = read.frame(self.gemini_output_dir, dt.datetime.utcfromtimestamp(utime), var=self.species)
            Ts = model2pointsgeogcoords(xg, dat[self.species], galt, glon, glat)
            Ts0 = Ts.reshape(galt.shape)

        else:
            s = output_shape(utime, galt)
            Ts0 = np.empty(s)
            for i, ut in enumerate(utime):
                dat = read.frame(self.gemini_output_dir, dt.datetime.utcfromtimestamp(ut), var=self.species)
                Ts = model2pointsgeogcoords(xg, dat[self.species], galt, glon, glat)
                Ts0[i] = Ts.reshape(galt.shape)

        return Ts0



class Velocity(object):
    def __init__(self, type, params, utime0, apex=None):
        # set density function
        self.Vi_function = getattr(self, type)
        # set starttime
        self.utime0 = utime0

        # initialize Apex object for all mapping
        self.apex = apex

        # set parameters as class attributes
        for param, value in params.items():
            setattr(self, param, value)



    def __call__(self, utime, glat, glon, galt):
        return self.Vi_function(utime, glat, glon, galt)


    def uniform(self, utime, glat, glon, galt):
        """
        Uniform velocity at all points.  Uniform velocity is caluclated by translating the
        velocity vetor to all points in cartesain coordinates, then adjusting the vertical velocity
        to zero and renormalizing so they all have the same magnitude.  This creates the approximate
        effect of uniform velocity across the entire FoV regardless of the local coordinate system.

        Parameters
        ----------
        value: list
            The vector value to assign at all points [E, N, U] (m/s)
        """

        # alat, alon = self.apex.geo2apex(glat.ravel(), glon.ravel(), galt.ravel()/1000.)
        alat, alon = self.apex.geo2apex(glat, glon, galt/1000.)
        map_glat, map_glon, _ = self.apex.apex2geo(alat, alon, 300.)

        # Find ECEF velocity components for given geodetic velocity at center of points
        u, v, w = pm.enu2uvw(self.value[0], self.value[1], self.value[2], np.nanmean(map_glat), np.nanmean(map_glon))
        # Find ENU components for same velocity translated to all mapped locations
        e, n, u = pm.uvw2enu(u, v, w, map_glat, map_glon)
        u = np.zeros(u.shape)   # set up component to zero
        V_map = np.array([e,n,u])
        # rescale velocities so that they all have the same magnitude as the original vector
        V_scale = V_map*np.linalg.norm(self.value)/np.linalg.norm(V_map, axis=0)
        # map velocities along field lines back to the original heights
        if not np.isscalar(galt):
            # map_V_to_height doesn't handle multidimensional arrays, so must flatten and reform
            V0 = self.apex.map_V_to_height(alat.ravel(), alon.ravel(), 300., galt.ravel()/1000., V_scale.reshape(3,-1))
            # reform original array shape
            V0 = V0.T.reshape(galt.shape+(3,))
        else:
            V0 = self.apex.map_V_to_height(alat, alon, 300., galt/1000., V_scale)

        s = output_shape(utime, galt)
        if not s:
            VE0 = V0
        else:
            VE0 = np.broadcast_to(V0, s+(3,))

        # s = (utime.shape[0],)+galt.shape+(3,)
        # VE0 = np.broadcast_to(VE, s)

        return VE0

    def uniform_mlat_aligned(self, utime, glat, glon, galt):
        """
        Velocity will have the same Apex magnetic components at all points.  This is useful for flows
        expected to align with magnetic meridians (i.e., in the auroral zone).

        Parameters
        ----------
        value: list
            The Apex vector value to assign at all points [mag E, mag N, mag U] (m/s)
        """

        Ve1, Ve2, Ve3 = self.value

        # Find base vector at each location
        _, _, _, _, _, _, _, _, _, e1, e2, e3 = self.apex.basevectors_apex(glat.ravel(), glon.ravel(), galt.ravel()/1000.)
        # reshape basevector arrays to match the original input
        e1 = e1.T.reshape(glat.shape+(3,))
        e2 = e2.T.reshape(glat.shape+(3,))
        e3 = e3.T.reshape(glat.shape+(3,))

        # calculate V in geodetic coordinates
        VE = Ve1*e1 + Ve2*e2 + Ve3*e3

        # s = (utime.shape[0],)+galt.shape+(3,)
        # VE0 = np.broadcast_to(VE, s)

        s = output_shape(utime, galt)
        if not s:
            VE0 = VE
        else:
            VE0 = np.broadcast_to(VE, s+(3,))


        return VE0


    def uniform_glat_aligned(self, utime, glat, glon, galt):
        """
        Velocity will have the same Geodetic components at all points.  CAUTION: This is probably not
        the most appropriate function for creating a "uniform" field close to the poles.

        Parameters
        ----------
        value: list
            The vector value to assign at all points [E, N, U] (m/s)
        """

        # s = (utime.shape[0],)+galt.shape+(3,)
        # VE0 = np.broadcast_to(self.value, s)

        s = output_shape(utime, galt)
        if not s:
            VE0 = self.value
        else:
            VE0 = np.broadcast_to(self.value, s+(3,))


        return VE0

    def gemini(self, utime, glat, glon, galt):
        """
        Velocity output from GEMINI model.

        Parameters
        ----------
        gemini_output_dir: string
            Path to directory of GEMINI output files
        """

        xg = read.grid(self.gemini_output_dir)

        if not utime.shape:
            dat = read.frame(self.gemini_output_dir, dt.datetime.utcfromtimestamp(utime), var='v1')
            V1 = model2pointsgeogcoords(xg, dat['v1'], galt, glon, glat)
            dat = read.frame(self.gemini_output_dir, dt.datetime.utcfromtimestamp(utime), var='v2')
            V2 = model2pointsgeogcoords(xg, dat['v2'], galt, glon, glat)
            dat = read.frame(self.gemini_output_dir, dt.datetime.utcfromtimestamp(utime), var='v3')
            V3 = model2pointsgeogcoords(xg, dat['v3'], galt, glon, glat)

            V = np.array([V1, V2, V3]).T
            Vi0 = V.reshape(galt.shape+(3,))

        else:
            s = output_shape(utime, galt)+(3,)
            Vi0 = np.empty(s)
            for i, ut in enumerate(utime):
                dat = read.frame(self.gemini_output_dir, dt.datetime.utcfromtimestamp(utime[i]), var='v1')
                V1 = model2pointsgeogcoords(xg, dat['v1'], galt, glon, glat)
                dat = read.frame(self.gemini_output_dir, dt.datetime.utcfromtimestamp(utime[i]), var='v2')
                V2 = model2pointsgeogcoords(xg, dat['v2'], galt, glon, glat)
                dat = read.frame(self.gemini_output_dir, dt.datetime.utcfromtimestamp(utime[i]), var='v3')
                V3 = model2pointsgeogcoords(xg, dat['v3'], galt, glon, glat)

                V = np.array([V1, V2, V3]).T
                Vi0[i] = V.reshape(galt.shape+(3,))

        return Vi0
