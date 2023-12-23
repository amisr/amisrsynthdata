import numpy as np
import pymap3d as pm
import datetime as dt
import warnings
# from apexpy import Apex
from .utils import output_shape, gemini_helper


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
        Uniform velocity at all points.  Uniform velocity is caluclated by
        translating the velocity vetor to all points in cartesain coordinates,
        then adjusting the vertical velocity to zero and renormalizing so they
        all have the same magnitude.  This creates the approximate effect of
        uniform velocity across the entire FoV regardless of the local
        coordinate system.

        Parameters
        ----------
        value: list
            The vector value to assign at all points [E, N, U] (m/s)
        cent_lat: float
            Geodetic latitude of initial specified vector (deg)
        cent_lon: float
            Geodetic longitude of initial specified vector (deg)
        """

        # warning filters are needed to supress apexpy warnings with NaN input
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            alat, alon = self.apex.geo2apex(glat, glon, galt / 1000.)
            map_glat, map_glon, _ = self.apex.apex2geo(alat, alon, 300.)

        # Find ECEF velocity components for given geodetic velocity at center
        # of points
        u, v, w = pm.enu2uvw(
            self.value[0], self.value[1], self.value[2],
            self.cent_lat, self.cent_lon)
        # Find ENU components for same velocity translated to all mapped
        # locations
        e, n, u = pm.uvw2enu(u, v, w, map_glat, map_glon)
        u = np.zeros(u.shape)   # set up component to zero
        V_map = np.array([e, n, u])
        # rescale velocities so that they all have the same magnitude as the
        # original vector
        V_scale = V_map * \
            np.linalg.norm(self.value) / np.linalg.norm(V_map, axis=0)
        # map velocities along field lines back to the original heights
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            if not np.isscalar(galt):
                # map_V_to_height doesn't handle multidimensional arrays, so
                # must flatten and reform
                V0 = self.apex.map_V_to_height(alat.ravel(), alon.ravel(
                ), 300., galt.ravel() / 1000., V_scale.reshape(3, -1))
                # reform original array shape
                V0 = V0.T.reshape(galt.shape + (3,))
            else:
                V0 = self.apex.map_V_to_height(
                    alat, alon, 300., galt / 1000., V_scale)

        s = output_shape(utime, galt)
        if not s:
            VE0 = V0
        else:
            VE0 = np.broadcast_to(V0, s + (3,))

        return VE0

    def uniform_mlat_aligned(self, utime, glat, glon, galt):
        """
        Velocity will have the same Apex magnetic components at all points.
        This is useful for flows expected to align with magnetic meridians
        (i.e., in the auroral zone).

        Parameters
        ----------
        value: list
            The Apex vector value to assign at all points
            [mag E, mag N, mag U] (m/s)
        """

        Ve1, Ve2, Ve3 = self.value

        # Find base vector at each location
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            if np.isscalar(galt):
                (_, _, _, _, _, _,
                    _, _, _, e1, e2, e3) = self.apex.basevectors_apex(
                    glat, glon, galt / 1000.)
            else:
                (_, _, _, _, _, _,
                    _, _, _, e1, e2, e3) = self.apex.basevectors_apex(
                    glat.ravel(), glon.ravel(), galt.ravel() / 1000.)
                # reshape basevector arrays to match the original input
                e1 = e1.T.reshape(glat.shape + (3,))
                e2 = e2.T.reshape(glat.shape + (3,))
                e3 = e3.T.reshape(glat.shape + (3,))

        # calculate V in geodetic coordinates
        VE = Ve1 * e1 + Ve2 * e2 + Ve3 * e3

        s = output_shape(utime, galt)
        if not s:
            VE0 = VE
        else:
            VE0 = np.broadcast_to(VE, s + (3,))

        return VE0

    def uniform_glat_aligned(self, utime, glat, glon, galt):
        """
        Velocity will have the same Geodetic components at all points.
        CAUTION: This is probably not the most appropriate function for
        creating a "uniform" field close to the poles.

        Parameters
        ----------
        value: list
            The vector value to assign at all points [E, N, U] (m/s)
        """

        s = output_shape(utime, galt)
        if not s:
            if np.isnan(galt):
                VE0 = np.array([np.nan, np.nan, np.nan])
            else:
                VE0 = np.array(self.value)
        else:
            VE0 = np.broadcast_to(self.value, s + (3,)).copy()
            VE0[..., np.isnan(galt), 0] = np.nan
            VE0[..., np.isnan(galt), 1] = np.nan
            VE0[..., np.isnan(galt), 2] = np.nan

        return VE0

    def gemini(self, utime, glat, glon, galt):
        """
        Velocity output from GEMINI model.

        Parameters
        ----------
        gemini_output_dir: string
            Path to directory of GEMINI output files
        """

        gh = gemini_helper(self.gemini_output_dir, glat, glon, galt)

        if not utime.shape:
            V1 = gh.query_model(dt.datetime.utcfromtimestamp(utime), 'v1')
            V2 = gh.query_model(dt.datetime.utcfromtimestamp(utime), 'v2')
            V3 = gh.query_model(dt.datetime.utcfromtimestamp(utime), 'v3')
            Vi0 = np.moveaxis([V1, V2, V3], 0, -1)

        else:
            s = output_shape(utime, galt) + (3,)
            Vi0 = np.empty(s)
            for i, ut in enumerate(utime):
                V1 = gh.query_model(dt.datetime.utcfromtimestamp(ut), 'v1')
                V2 = gh.query_model(dt.datetime.utcfromtimestamp(ut), 'v2')
                V3 = gh.query_model(dt.datetime.utcfromtimestamp(ut), 'v3')
                Vi0[i] = np.moveaxis([V1, V2, V3], 0, -1)

        return Vi0
