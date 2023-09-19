# ionosphere.py

import numpy as np
import datetime as dt
from apexpy import Apex
from .state_functions import Density, Temperature, Velocity
from .state_functions.utils import output_shape


class Ionosphere(object):
    """
    Class describing the ionosphere used to create synthetic data.

    Parameters
    ----------
    config : :obj:`yaml`
        Ionosphere configuration parameters
    """

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
        """
        Generate an array of zeros of the correct output shape based on input
        time an coordinate arrays.  This can then be indexed and populated
        with values.

        Parameters
        ----------
        ut : float or np.ndarray
            time or time array
        x : float or np.ndarray
            coordinate or coordinate array
        vec: bool
            Flag indicating whether or not the returned array should add a
            dimension for vector components

        Returns
        -------
        float or np.ndarray
            Either float 0.0 for scalar input or an array full of zeros
            of the correct shape
        """

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
        """
        Calculate density at the requested times and locations

        Parameters
        ----------
        utime : float or np.ndarray
            Unix Time
        glat : float or np.ndarray
            Geodetic Latitude
        glon : float or np.ndarray
            Geodetic Longitude
        galt : float or np.ndarray
            Geodetic Altitude

        Returns
        -------
        float or np.ndarray
            Density values at each of the input times and locations

        """
        dens = self.zero_array(utime, galt)
        for fun in self.density_functions:
            dens = dens + fun(utime, glat, glon, galt)
        return dens

    def itemp(self, utime, glat, glon, galt):
        """
        Calculate ion temperature at the requested times and locations

        Parameters
        ----------
        utime : float or np.ndarray
            Unix Time
        glat : float or np.ndarray
            Geodetic Latitude
        glon : float or np.ndarray
            Geodetic Longitude
        galt : float or np.ndarray
            Geodetic Altitude

        Returns
        -------
        float or np.ndarray
            Ion temperature values at each of the input times and locations

        """
        itemp = self.zero_array(utime, galt)
        for fun in self.itemp_functions:
            itemp = itemp + fun(utime, glat, glon, galt)
        return itemp

    def etemp(self, utime, glat, glon, galt):
        """
        Calculate electron temperature at the requested times and locations

        Parameters
        ----------
        utime : float or np.ndarray
            Unix Time
        glat : float or np.ndarray
            Geodetic Latitude
        glon : float or np.ndarray
            Geodetic Longitude
        galt : float or np.ndarray
            Geodetic Altitude

        Returns
        -------
        float or np.ndarray
            Electron temperature values at each of the input times and
            locations

        """
        etemp = self.zero_array(utime, galt)
        for fun in self.etemp_functions:
            etemp = etemp + fun(utime, glat, glon, galt)
        return etemp

    def velocity(self, utime, glat, glon, galt):
        """
        Calculate velocity at the requested times and locations

        Parameters
        ----------
        utime : float or np.ndarray
            Unix Time
        glat : float or np.ndarray
            Geodetic Latitude
        glon : float or np.ndarray
            Geodetic Longitude
        galt : float or np.ndarray
            Geodetic Altitude

        Returns
        -------
        np.ndarray
            Velocity components at each of the input times and locations

        Note
        ----
        The returned velocity array will be one dimension larger than the
        input time and coordinate arrays to account for the three velocity
        components (east, north, up)

        """
        vel = self.zero_array(utime, galt, vec=True)
        for fun in self.velocity_functions:
            vel = vel + fun(utime, glat, glon, galt)
        return vel
