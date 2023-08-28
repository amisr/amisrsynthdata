"""
Test of amisrsynthdata.Radar class
"""

import pytest
import yaml
import h5py
import numpy as np
import pymap3d as pm
import os

from amisrsynthdata.radar import Radar

@pytest.fixture
def radar():

    config_file = os.path.join(os.path.dirname(__file__), 'config.yaml')
    with open(config_file, 'r') as cf:
        config = yaml.load(cf, Loader=yaml.FullLoader)

    rad = Radar(config)

    return rad

@pytest.fixture
def datafile():
    
    filename = os.path.join(os.path.dirname(__file__), 'synthetic_data.h5')
    h5file = h5py.File(filename, 'r')
    
    return h5file


def test_init(radar):
    assert radar.radar_name == 'Poker Flat'
    assert radar.radar_abbrev == 'PFISR'
    assert radar.site_lat == 65.13
    assert radar.site_lon == -147.47
    assert radar.site_alt == 213.
    assert radar.integration_period == 60.


#def test_beamcodes(radar, datafile):
#    bc_array = datafile['BeamCodes'][:]
#    assert np.array_equal(radar.beam_codes, bc_array, equal_nan=True)

def test_beams_from_beam_codes(radar):
    beamcodes = [64016, 64157, 64964, 65066]
    truth_az = [  14.04, -154.3 ,  -34.69,   75.03]
    truth_el = [90.  , 77.5 , 66.09, 65.56]
    truth_ksys = [np.nan, np.nan, np.nan, np.nan]
    azimuth, elevation, ksys = radar.beams_from_beam_codes(beamcodes)
    assert np.array_equal(azimuth, truth_az)
    assert np.array_equal(elevation, truth_el)
    assert np.array_equal(ksys, truth_ksys, equal_nan=True)

def test_beams_from_az_el(radar):
    azimuth = [  14.04, -154.3 ,  -34.69,   75.03]
    elevation = [90.  , 77.5 , 66.09, 65.56]
    truth_bc = [90000, 90001, 90002, 90003]
    truth_ksys = [np.nan, np.nan, np.nan, np.nan]
    beamcodes, ksys = radar.beams_from_az_el(azimuth, elevation)
    assert np.array_equal(beamcodes, truth_bc)
    assert np.array_equal(ksys, truth_ksys, equal_nan=True)
    
def test_caculate_acf_gates(radar):
    slant_range_config = [80000., 800000., 3000.]
    truth_sr = np.arange(*slant_range_config)
    truth_lat = list()
    truth_lon = list()
    truth_alt = list()
    for az, el in zip(radar.beam_azimuth, radar.beam_elevation):
        lat, lon, alt = pm.aer2geodetic(az, el, truth_sr, radar.site_lat, radar.site_lon, radar.site_alt)
        truth_lat.append(lat)
        truth_lon.append(lon)
        truth_alt.append(alt)

    slant_range, lat, lon, alt = radar.calculate_acf_gates(slant_range_config)

    assert np.array_equal(slant_range, truth_sr)
    assert np.array_equal(lat, truth_lat)
    assert np.array_equal(lon, truth_lon)
    assert np.array_equal(alt, truth_alt)

def test_calculate_gates(radar):

    altbins = np.arange(100., 800., 50.)*1000.
    truth_sr = np.empty((len(radar.beam_codes), len(altbins)-1))
    truth_lat = np.empty((len(radar.beam_codes), len(altbins)-1))
    truth_lon = np.empty((len(radar.beam_codes), len(altbins)-1))
    truth_alt = np.empty((len(radar.beam_codes), len(altbins)-1))
    for b in range(len(radar.beam_codes)):
        for i in range(len(altbins)-1):
            abidx = np.argwhere((radar.acf_alt[b]>=altbins[i]) & (radar.acf_alt[b]<altbins[i+1])).flatten()
            sr = np.mean(radar.acf_slant_range[abidx])
            lat, lon, alt = pm.aer2geodetic(radar.beam_azimuth[b], radar.beam_elevation[b], sr, radar.site_lat, radar.site_lon, radar.site_alt)
            truth_sr[b, i] = sr
            truth_lat[b, i] = lat
            truth_lon[b, i] = lon
            truth_alt[b, i] = alt

    slant_range, lat, lon, alt = radar.calculate_gates([[100000., 800000., 50000.]])

    assert np.array_equal(slant_range, truth_sr, equal_nan=True)
    assert np.array_equal(lat, truth_lat, equal_nan=True)
    assert np.array_equal(lon, truth_lon, equal_nan=True)
    assert np.array_equal(alt, truth_alt, equal_nan=True)

#def test_slant_range_p(radar, datafile):
#    sr_array = datafile['NeFromPower/Range'][:]
#    assert np.array_equal(radar.acf_slant_range, sr_array)
#
#def test_slant_range_p_coords(radar):
#
#    for bc, lat_p, lon_p, alt_p in zip(radar.beam_codes, radar.lat_p, radar.lon_p, radar.alt_p):
#        az = bc[1]
#        el = bc[2]
#        lat, lon, alt = pm.aer2geodetic(az, el, radar.slant_range_p, radar.site_lat, radar.site_lon, radar.site_alt)
#        assert np.allclose(lat_p, lat)
#        assert np.allclose(lon_p, lon)
#        assert np.allclose(alt_p, alt)
#
#def test_slant_range(radar, datafile):
#    sr_array = datafile['FittedParams/Range'][:]
#    assert np.allclose(radar.slant_range, sr_array, equal_nan=True)
#
#def test_slat_range_coords(radar, datafile):
#    lat_array = datafile['Geomag/Latitude'][:]
#    lon_array = datafile['Geomag/Longitude'][:]
#    alt_array = datafile['Geomag/Altitude'][:]
#    assert np.allclose(radar.lat, lat_array, equal_nan=True)
#    assert np.allclose(radar.lon, lon_array, equal_nan=True)
#    assert np.allclose(radar.alt, alt_array, equal_nan=True)

def test_kvec(radar):

    kvec = radar.kvec_all_gates()

    truth_kvec = list()
    for sr, az, el in zip(radar.slant_range, radar.beam_azimuth, radar.beam_elevation):
        #az = bc[1]
        #el = bc[2]
        # aer -> uvw
        e, n, u = pm.aer2enu(az, el, 1.0)
        u, v, w = pm.enu2uvw(e, n, u, radar.site_lat, radar.site_lon, radar.site_alt)
        # aer(slant rng) -> new position
        lat, lon, alt = pm.aer2geodetic(az, el, sr, radar.site_lat, radar.site_lon, radar.site_alt)
        # uvw(new position) -> aer
        ke, kn, ku = pm.uvw2enu(u, v, w, lat, lon)
        truth_kvec.append(np.array([ke, kn, ku]).T)


    assert np.array_equal(kvec, truth_kvec, equal_nan=True)

