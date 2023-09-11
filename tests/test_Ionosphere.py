"""
Test of amisrsynthdata Ionosphere class
"""

import pytest
import yaml
import os
import numpy as np
from apexpy import Apex
from amisrsynthdata.ionosphere import Ionosphere
from amisrsynthdata.state_functions import utils


@pytest.fixture
def config():
    config_file = os.path.join(os.path.dirname(__file__), 'config.yaml')
    with open(config_file, 'r') as cf:
        config = yaml.load(cf, Loader=yaml.FullLoader)
    return config


@pytest.fixture
def ionosphere(config):
    iono = Ionosphere(config)
    return iono


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


def test_init(ionosphere, config):
    assert isinstance(ionosphere.apex, Apex)
    assert ionosphere.ion_mass == config['GENERAL']['ion_mass']


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
class TestIonosphereFunctions:

    def test_density(self, ionosphere, config, utime, glat, glon, galt):
        dens = ionosphere.density(utime, glat, glon, galt)
        truth_dens = config['DENSITY'][0]['uniform']['value']
        np.testing.assert_allclose(dens, truth_dens)
        expected_shape = utils.output_shape(utime, galt)
        if expected_shape:
            assert dens.shape == expected_shape

    def test_itemp(self, ionosphere, config, utime, glat, glon, galt):
        itemp = ionosphere.itemp(utime, glat, glon, galt)
        truth_itemp = config['ITEMP'][0]['uniform']['value']
        np.testing.assert_allclose(itemp, truth_itemp)
        expected_shape = utils.output_shape(utime, galt)
        if expected_shape:
            assert itemp.shape == expected_shape

    def test_etemp(self, ionosphere, config, utime, glat, glon, galt):
        etemp = ionosphere.etemp(utime, glat, glon, galt)
        truth_etemp = config['ETEMP'][0]['uniform']['value']
        np.testing.assert_allclose(etemp, truth_etemp)
        expected_shape = utils.output_shape(utime, galt)
        if expected_shape:
            assert etemp.shape == expected_shape

    def test_velocity(self, ionosphere, config, utime, glat, glon, galt):
        vel = ionosphere.velocity(utime, glat, glon, galt)
        truth_vel0 = config['VELOCITY'][0]['uniform_glat_aligned']['value']
        truth_vel = np.broadcast_to(truth_vel0, vel.shape)
        np.testing.assert_allclose(vel, truth_vel)
        expected_shape = utils.output_shape(utime, galt)
        if expected_shape:
            assert vel.shape == expected_shape + (3,)

    def test_sum_density(self, config, utime, glat, glon, galt):
        config['DENSITY'].append({'uniform': {'value': 5.e10}})
        iono = Ionosphere(config)
        dens = iono.density(utime, glat, glon, galt)
        truth_dens = config['DENSITY'][0]['uniform']['value'] + 5.e10
        np.testing.assert_allclose(dens, truth_dens)
        expected_shape = utils.output_shape(utime, galt)
        if expected_shape:
            assert dens.shape == expected_shape

    def test_sum_itemp(self, config, utime, glat, glon, galt):
        config['ITEMP'].append({'uniform': {'value': 2000.}})
        iono = Ionosphere(config)
        itemp = iono.itemp(utime, glat, glon, galt)
        truth_itemp = config['ITEMP'][0]['uniform']['value'] + 2000.
        np.testing.assert_allclose(itemp, truth_itemp)
        expected_shape = utils.output_shape(utime, galt)
        if expected_shape:
            assert itemp.shape == expected_shape

    def test_sum_etemp(self, config, utime, glat, glon, galt):
        config['ETEMP'].append({'uniform': {'value': 2000.}})
        iono = Ionosphere(config)
        etemp = iono.etemp(utime, glat, glon, galt)
        truth_etemp = config['ETEMP'][0]['uniform']['value'] + 2000.
        np.testing.assert_allclose(etemp, truth_etemp)
        expected_shape = utils.output_shape(utime, galt)
        if expected_shape:
            assert etemp.shape == expected_shape

    def test_sum_velocity(self, config, utime, glat, glon, galt):
        config['VELOCITY'].append(
            {'uniform_glat_aligned': {'value': [200., 0., 0]}})
        iono = Ionosphere(config)
        vel = iono.velocity(utime, glat, glon, galt)
        truth_vel0 = config['VELOCITY'][0]['uniform_glat_aligned']['value']
        truth_vel0 = truth_vel0 + np.array([200., 0., 0.])
        truth_vel = np.broadcast_to(truth_vel0, vel.shape)
        np.testing.assert_allclose(vel, truth_vel)
        expected_shape = utils.output_shape(utime, galt)
        if expected_shape:
            assert vel.shape == expected_shape + (3,)


@pytest.mark.parametrize('utime', [time_0d, time_1d])
@pytest.mark.parametrize('galt', [galt_0d, galt_1d, galt_2d, galt_3d])
class TestZeroArray:

    def test_zero_array(self, ionosphere, utime, galt):
        zout = ionosphere.zero_array(utime, galt)
        expected_shape = utils.output_shape(utime, galt)
        if expected_shape:
            np.testing.assert_array_equal(np.zeros(expected_shape), zout)
        else:
            assert zout == 0.0

    def test_zero_array_vector(self, ionosphere, utime, galt):
        zout = ionosphere.zero_array(utime, galt, vec=True)
        expected_shape = utils.output_shape(utime, galt)
        if expected_shape:
            np.testing.assert_array_equal(
                np.zeros(expected_shape + (3,)), zout)
        else:
            np.testing.assert_array_equal(np.zeros(3), zout)
