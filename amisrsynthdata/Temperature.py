# Temperature.py
import numpy as np

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

        s = (utime.shape[0],)+galt.shape
        Ts0 = np.full(s, self.value)

        return Ts0

    def hypertan(self, utime, glat, glon, galt):

        Ts = self.maxtemp*np.tanh(galt/self.scale_height)

        s = (utime.shape[0],)+galt.shape
        Ts0 = np.full(s, Ts)

        return Ts0
