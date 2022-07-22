# Temperature.py
import numpy as np

# Because ion and electron temperature functions tend to be similar, they are both caputred
#   within this function.  Ions and Electrons DO NOT have to use the same background function
#   or parameters, these are specified seperately in the config file.
class Temperature(object):
    def __init__(self, utime0, config_params):
        # set temperature function
        self.Ts_function = getattr(self, config_params['type'])
        # set starttime
        self.utime0 = utime0

        # assign remaining config options to parameters to be handled by each function
        config_params.pop('type')
        self.params = config_params


    def __call__(self, utime, glat, glon, galt):
        return self.Ts_function(utime, glat, glon, galt)


    def uniform(self, utime, glat, glon, galt):
        Ts = float(self.params['value'])

        s = (utime.shape[0],)+galt.shape
        Ts0 = np.full(s, Ts)

        return Ts0

    def hypertan(self, utime, glat, glon, galt):
        maxTs = float(self.params['maxtemp'])
        scale_height = float(self.params['scale_height'])

        Ts = maxTs*np.tanh(galt/scale_height)

        s = (utime.shape[0],)+galt.shape
        Ts0 = np.full(s, Ts)

        return Ts0
