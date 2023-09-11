# ionosphere.py

import numpy as np
import datetime as dt
from apexpy import Apex
from .state_functions import Density, Temperature, Velocity
from .state_functions.utils import output_shape


class Ionosphere(object):

    def __init__(self, config):

        starttime = config['GENERAL']['starttime']
        start_utime = (
            starttime -
            dt.datetime.utcfromtimestamp(0)).total_seconds()

        # initialize Apex object
        self.apex = Apex(date=starttime)

        self.ion_mass = config['GENERAL']['ion_mass']

        # create lists of all the functions that will be used for each
        # ionospheric state
        self.density_functions = [
            Density(
                name,
                params,
                start_utime) for state_funct in config['DENSITY'] for name,
            params in state_funct.items()]
        self.velocity_functions = [
            Velocity(
                name,
                params,
                start_utime,
                apex=self.apex) for state_funct in config['VELOCITY'] for name,
            params in state_funct.items()]
        self.etemp_functions = [
            Temperature(
                name,
                params,
                start_utime) for state_funct in config['ETEMP'] for name,
            params in state_funct.items()]
        self.itemp_functions = [
            Temperature(
                name,
                params,
                start_utime) for state_funct in config['ITEMP'] for name,
            params in state_funct.items()]

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
        dens = self.zero_array(utime, galt)
        for fun in self.density_functions:
            dens = dens + fun(utime, glat, glon, galt)
        return dens

    def velocity(self, utime, glat, glon, galt):
        vel = self.zero_array(utime, galt, vec=True)
        for fun in self.velocity_functions:
            vel = vel + fun(utime, glat, glon, galt)
        return vel

    def etemp(self, utime, glat, glon, galt):
        etemp = self.zero_array(utime, galt)
        for fun in self.etemp_functions:
            etemp = etemp + fun(utime, glat, glon, galt)
        return etemp

    def itemp(self, utime, glat, glon, galt):
        itemp = self.zero_array(utime, galt)
        for fun in self.itemp_functions:
            itemp = itemp + fun(utime, glat, glon, galt)
        return itemp
