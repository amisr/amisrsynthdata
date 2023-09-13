import numpy as np
import datetime as dt
from .utils import output_shape, gemini_helper

# Because ion and electron temperature functions tend to be similar, they are
#   both caputred within this function.  Ions and Electrons DO NOT have to use
#   the same background function or parameters, these are specified seperately
#   in the config file.


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

        s = output_shape(utime, galt)
        if not s:
            if np.isnan(galt):
                Ts0 = np.nan
            else:
                Ts0 = self.value
        else:
            Ts0 = np.full(s, self.value)
            Ts0[..., np.isnan(galt)] = np.nan

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

        Ts = self.maxtemp * np.tanh(galt / self.scale_height)

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
            Which species (ion or electron) should be read from GEMINI output
            ('Te' or 'Ti')
        """

        gh = gemini_helper(self.gemini_output_dir, glat, glon, galt)

        if not utime.shape:
            Ts0 = gh.query_model(
                dt.datetime.utcfromtimestamp(utime),
                self.species)

        else:
            s = output_shape(utime, galt)
            Ts0 = np.empty(s)
            for i, ut in enumerate(utime):
                Ts0[i] = gh.query_model(
                    dt.datetime.utcfromtimestamp(ut), self.species)

        return Ts0
