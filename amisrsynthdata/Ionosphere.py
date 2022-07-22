# Ionosphere.py
import numpy as np
import datetime as dt
from apexpy import Apex
from .Density import Density
from .Velocity import Velocity
from .Temperature import Temperature


# try:
#     import ConfigParser as configparser
# except ImportError:
#     import configparser

# Figure out how to specify multiple species sensibly
# this can probably wait as a feature - unlikely to come up often
# may be nice to be able to have some kind of placeholder just to make the arrays the full shape though

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

    def density(self, utime, glat, glon, galt):
        # call each density instance and sum the results
        dens = np.zeros((utime.shape[0],)+galt.shape)
        for fun in self.density_functions:
            dens = dens + fun(utime, glat, glon, galt)
        return dens

    def velocity(self, utime, glat, glon, galt):
        vel = np.zeros((utime.shape[0],)+galt.shape+(3,))
        for fun in self.velocity_functions:
            vel = vel + fun(utime, glat, glon, galt)
        return vel

    def etemp(self, utime, glat, glon, galt):
        # call each electron temperature instance and sum the results
        etemp = np.zeros((utime.shape[0],)+galt.shape)
        for fun in self.etemp_functions:
            etemp = etemp + fun(utime, glat, glon, galt)
        return etemp

    def itemp(self, utime, glat, glon, galt):
        # call each ion temperature instance and sum the results
        itemp = np.zeros((utime.shape[0],)+galt.shape)
        for fun in self.itemp_functions:
            itemp = itemp + fun(utime, glat, glon, galt)
        return itemp
