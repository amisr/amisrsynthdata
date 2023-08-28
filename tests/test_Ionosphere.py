"""
Test of amisrsynthdata Ionosphere class
"""

import pytest
#import amisrsynthdata
import yaml
import os
import numpy as np
from apexpy import Apex
from amisrsynthdata.ionosphere import Ionosphere

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

@pytest.fixture
def varin():
    sut = np.datetime64('2020-01-01T06:00:00').astype(int)
    eut = np.datetime64('2020-01-01T07:00:00').astype(int)
    utime = np.arange(sut, eut, 60.)
    glat = np.linspace(65.0, 67.0, 10)
    glon = np.linspace(-145.0, 149.0, 10)
    galt = np.linspace(300.0, 400.0, 10)*1000.
    return [utime, glat, glon, galt]


def test_init(ionosphere, config):
    assert isinstance(ionosphere.apex, Apex)
    assert ionosphere.ion_mass ==  config['GENERAL']['ion_mass']

def test_density(ionosphere, config, varin):

    dens = ionosphere.density(*varin)

    assert np.allclose(dens, config['DENSITY'][0]['uniform']['value'])

def test_velocity(ionosphere, config, varin):

    vel = ionosphere.velocity(*varin)

    assert np.allclose(vel, config['VELOCITY'][0]['uniform_glat_aligned']['value'])


def test_etemp(ionosphere, config, varin):

    etemp = ionosphere.etemp(*varin)

    assert np.allclose(etemp, config['ETEMP'][0]['uniform']['value'])


def test_itemp(ionosphere, config, varin):

    itemp = ionosphere.itemp(*varin)

    assert np.allclose(itemp, config['ITEMP'][0]['uniform']['value'])

def test_sum_density(config, varin):

    config['DENSITY'].append({'uniform': {'value': 5.e10}})
    iono = Ionosphere(config)
    dens = iono.density(*varin)
    truth_dens = config['DENSITY'][0]['uniform']['value'] + 5.e10

    assert np.allclose(dens, truth_dens)

def test_sum_velocity(config, varin):

    config['VELOCITY'].append({'uniform_glat_aligned': {'value': [200., 0., 0]}})
    iono = Ionosphere(config)
    vel = iono.velocity(*varin)
    truth_vel = np.array(config['VELOCITY'][0]['uniform_glat_aligned']['value']) + np.array([200., 0., 0.])
    
    assert np.allclose(vel, truth_vel)

def test_sum_etemp(config, varin):

    config['ETEMP'].append({'uniform': {'value': 2000.}})
    iono = Ionosphere(config)
    etemp = iono.etemp(*varin)
    truth_etemp = config['ETEMP'][0]['uniform']['value'] + 2000.

    assert np.allclose(etemp, truth_etemp)

def test_sum_itemp(config, varin):

    config['ITEMP'].append({'uniform': {'value': 2000.}})
    iono = Ionosphere(config)
    itemp = iono.itemp(*varin)
    truth_itemp = config['ITEMP'][0]['uniform']['value'] + 2000.

    assert np.allclose(itemp, truth_itemp)

def test_zero_array(ionosphere, varin):

    # scalar input
    zout = ionosphere.zero_array(varin[0][0], varin[3][0])
    assert zout == 0.0

    # vector input
    zout = ionosphere.zero_array(varin[0], varin[3])
    assert np.array_equal(zout, np.zeros((len(varin[0]), len(varin[3]))))

    # vector option
    zout = ionosphere.zero_array(varin[0], varin[3], vec=True)
    assert np.array_equal(zout, np.zeros((len(varin[0]), len(varin[3]), 3)))


