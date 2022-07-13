# Ionosphere.py
import numpy as np
import datetime as dt
from apexpy import Apex
from .Density import Density
from .Velocity import Velocity
from .Temperature import Temperature


try:
    import ConfigParser as configparser
except ImportError:
    import configparser

# Figure out how to specify multiple species sensibly
# this can probably wait as a feature - unlikely to come up often
# may be nice to be able to have some kind of placeholder just to make the arrays the full shape though

class Ionosphere(object):

    def __init__(self, config_file):

        config = configparser.ConfigParser()
        config.read(config_file)

        starttime = dt.datetime.fromisoformat(config['GENERAL']['STARTTIME'])
        start_utime = (starttime-dt.datetime.utcfromtimestamp(0)).total_seconds()

        # initialize Apex object
        self.apex = Apex(date=starttime)

        # intialize parameter functions
        self.density = Density(start_utime, config['DENSITY'])
        self.velocity = Velocity(start_utime, config['VELOCITY'], apex=self.apex)
        self.etemp = Temperature(start_utime, config['ETEMP'])
        self.itemp = Temperature(start_utime, config['ITEMP'])

        self.ion_mass = np.array([float(im) for im in config['GENERAL']['ION_MASS'].split(',')])
