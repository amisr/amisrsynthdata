"""
Test of amisrsynthdata.Radar class
"""

import pytest
import amisrsynthdata
import yaml
import h5py
import numpy as np
import pymap3d as pm

@pytest.fixture
def radar():

    config_file = 'config.yaml'
    with open(config_file, 'r') as cf:
        config = yaml.load(cf, Loader=yaml.FullLoader)

    rad = amisrsynthdata.Radar(config)

    return rad

@pytest.fixture
def datafile():
    
    filename = 'synthetic_data.h5'
    h5file = h5py.File(filename, 'r')
    
    return h5file


def test_init(radar):
    assert radar.radar_name == 'Poker Flat'
    assert radar.radar_abbrev == 'PFISR'
    assert radar.site_lat == 65.13
    assert radar.site_lon == -147.47
    assert radar.site_alt == 213.
    assert radar.integration_period == 60.


def test_beamcodes(radar, datafile):
    bc_array = datafile['BeamCodes'][:]
    assert np.array_equal(radar.beam_codes, bc_array, equal_nan=True)

def test_slant_range_p(radar, datafile):
    sr_array = datafile['NeFromPower/Range'][:]
    assert np.array_equal(radar.slant_range_p, sr_array)

def test_slant_range_p_coords(radar):

    for bc, lat_p, lon_p, alt_p in zip(radar.beam_codes, radar.lat_p, radar.lon_p, radar.alt_p):
        az = bc[1]
        el = bc[2]
        lat, lon, alt = pm.aer2geodetic(az, el, radar.slant_range_p, radar.site_lat, radar.site_lon, radar.site_alt)
        assert np.allclose(lat_p, lat)
        assert np.allclose(lon_p, lon)
        assert np.allclose(alt_p, alt)

def test_slant_range(radar, datafile):
    sr_array = datafile['FittedParams/Range'][:]
    assert np.allclose(radar.slant_range, sr_array, equal_nan=True)

def test_slat_range_coords(radar, datafile):
    lat_array = datafile['Geomag/Latitude'][:]
    lon_array = datafile['Geomag/Longitude'][:]
    alt_array = datafile['Geomag/Altitude'][:]
    assert np.allclose(radar.lat, lat_array, equal_nan=True)
    assert np.allclose(radar.lon, lon_array, equal_nan=True)
    assert np.allclose(radar.alt, alt_array, equal_nan=True)

def test_kvec(radar):

    kvec = radar.kvec_all_gates()

    for sr, bc, kv in zip(radar.slant_range, radar.beam_codes, kvec):
        az = bc[1]
        el = bc[2]
        # aer -> uvw
        e, n, u = pm.aer2enu(az, el, 1.0)
        u, v, w = pm.enu2uvw(e, n, u, radar.site_lat, radar.site_lon, radar.site_alt)
        # aer(slant rng) -> new position
        lat, lon, alt = pm.aer2geodetic(az, el, sr, radar.site_lat, radar.site_lon, radar.site_alt)
        # uvw(new position) -> aer
        ke, kn, ku = pm.uvw2enu(u, v, w, lat, lon)

        assert np.allclose(np.array([ke, kn, ku]).T, kv, equal_nan=True)

