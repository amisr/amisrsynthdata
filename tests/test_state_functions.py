"""
Tests to run on all state functions
These mostly ensure functions can all handle expected input and output.
"""

import pytest
import yaml
import numpy as np
import datetime as dt
from apexpy import Apex
import os
from amisrsynthdata.state_functions import Density, Temperature, Velocity
from amisrsynthdata.state_functions import utils

# Set up time arrays
sut = np.datetime64('2020-01-01T06:00:00').astype(int)
eut = np.datetime64('2020-01-01T07:00:00').astype(int)
time_0d = np.datetime64('2020-01-01T06:30:00').astype(int)
time_1d = np.arange(sut, eut, 60.)

# Set up position arrays
glat_0d = 66.0
glon_0d = -147.0
galt_0d = 300. * 1000.

glat_1d = np.linspace(65.0, 67.0, 10)
glon_1d = np.linspace(-145.0, -149.0, 10)
galt_1d = np.linspace(300.0, 400.0, 10) * 1000.

glat_2d, glon_2d = np.meshgrid(glat_1d, glon_1d)
galt_2d = np.full(glat_2d.shape, 300. * 1000.)

glat_3d, glon_3d, galt_3d = np.meshgrid(glat_1d, glon_1d, galt_1d)


# Set up state function inputs
config_file = os.path.join(
    os.path.dirname(__file__),
    'state_functions_config.yaml')
with open(config_file, 'r') as cf:
    config = yaml.load(cf, Loader=yaml.FullLoader)

dens_functions = [(name, params) for sf in config['DENSITY']
                  for name, params in sf.items()]
temp_functions = [(name, params) for sf in config['TEMPERATURE']
                  for name, params in sf.items()]
vel_functions = [(name, params) for sf in config['VELOCITY']
                 for name, params in sf.items()]


@pytest.mark.parametrize('utime', [time_0d, time_1d])
@pytest.mark.parametrize('glat,glon,galt',
                         [(glat_0d,
                           glon_0d,
                           galt_0d),
                          (glat_1d,
                           glon_1d,
                           galt_1d),
                             (glat_2d,
                              glon_2d,
                              galt_2d),
                             (glat_3d,
                              glon_3d,
                              galt_3d)])
class TestStateFunctionsShape:

    @pytest.mark.parametrize('name,params', dens_functions)
    def test_dens_array_shape(self, name, params, utime, glat, glon, galt):
        dens = Density(name, params, sut)
        Ne = dens(utime, glat, glon, galt)
        expected_shape = utils.output_shape(utime, galt)
        if not expected_shape:
            assert np.isscalar(Ne)
        else:
            assert Ne.shape == expected_shape

    @pytest.mark.parametrize('name,params', temp_functions)
    def test_temp_array_shape(self, name, params, utime, glat, glon, galt):
        temp = Temperature(name, params, sut)
        Ts = temp(utime, glat, glon, galt)
        expected_shape = utils.output_shape(utime, galt)
        if not expected_shape:
            assert np.isscalar(Ts)
        else:
            assert Ts.shape == expected_shape

    @pytest.mark.parametrize('name,params', vel_functions)
    def test_vel_array_shape(self, name, params, utime, glat, glon, galt):
        apex = Apex(sut.astype('datetime64[s]').astype(dt.datetime))
        # apex = Apex(2020)
        vel = Velocity(name, params, sut, apex=apex)
        V = vel(utime, glat, glon, galt)
        expected_shape = utils.output_shape(utime, galt)
        if not expected_shape:
            expected_shape = (3,)
        else:
            expected_shape = expected_shape + (3,)
        assert V.shape == expected_shape


glat_nan_0d = np.nan
glon_nan_0d = np.nan
galt_nan_0d = np.nan

glat_nan_1d = glat_1d.copy()
glat_nan_1d[5] = np.nan
glon_nan_1d = glon_1d.copy()
glon_nan_1d[5] = np.nan
galt_nan_1d = galt_1d.copy()
galt_nan_1d[5] = np.nan

glat_nan_2d = glat_2d.copy()
glat_nan_2d[4, 6] = np.nan
glon_nan_2d = glon_2d.copy()
glon_nan_2d[4, 6] = np.nan
galt_nan_2d = galt_2d.copy()
galt_nan_2d[4, 6] = np.nan

glat_nan_3d = glat_3d.copy()
glat_nan_3d[5, 3, 2] = np.nan
glon_nan_3d = glon_3d.copy()
glon_nan_3d[5, 3, 2] = np.nan
galt_nan_3d = galt_3d.copy()
galt_nan_3d[5, 3, 2] = np.nan


@pytest.mark.parametrize('utime', [time_0d, time_1d])
@pytest.mark.parametrize('glat,glon,galt',
                         [(glat_nan_0d,
                           glon_nan_0d,
                           galt_nan_0d),
                          (glat_nan_1d,
                           glon_nan_1d,
                           galt_nan_1d),
                             (glat_nan_2d,
                              glon_nan_2d,
                              galt_nan_2d),
                             (glat_nan_3d,
                              glon_nan_3d,
                              galt_nan_3d)])
class TestStateFunctionsNaN:

    @pytest.mark.parametrize('name,params', dens_functions)
    def test_dens_nan_input(self, name, params, utime, glat, glon, galt):
        dens = Density(name, params, sut)
        Ne = dens(utime, glat, glon, galt)
        if np.isscalar(galt) or np.isscalar(utime):
            truth_nan = np.isnan(galt)
        else:
            truth_nan = np.isnan(np.broadcast_to(galt, Ne.shape))
        np.testing.assert_array_equal(np.isnan(Ne), truth_nan)

    @pytest.mark.parametrize('name,params', temp_functions)
    def test_temp_nan_input(self, name, params, utime, glat, glon, galt):
        temp = Temperature(name, params, sut)
        Ts = temp(utime, glat, glon, galt)
        if np.isscalar(galt) or np.isscalar(utime):
            truth_nan = np.isnan(galt)
        else:
            truth_nan = np.isnan(np.broadcast_to(galt, Ts.shape))
        np.testing.assert_equal(np.isnan(Ts), truth_nan)

    @pytest.mark.parametrize('name,params', vel_functions)
    def test_vel_nan_input(self, name, params, utime, glat, glon, galt):
        apex = Apex(sut.astype('datetime64[s]').astype(dt.datetime))
        vel = Velocity(name, params, sut, apex=apex)
        V = vel(utime, glat, glon, galt)
        if np.isscalar(galt) or np.isscalar(utime):
            truth_nan = np.isnan(galt)
        else:
            truth_nan = np.isnan(np.broadcast_to(galt, V.shape[:-1]))
        np.testing.assert_equal(np.all(np.isnan(V), axis=-1), truth_nan)
