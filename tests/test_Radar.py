"""
Test of amisrsynthdata.Radar class
"""

import pytest
import yaml
import numpy as np
import pymap3d as pm
import os
import warnings

from amisrsynthdata.radar import Radar


@pytest.fixture
def radar():

    config_file = os.path.join(os.path.dirname(__file__), 'config.yaml')
    with open(config_file, 'r') as cf:
        config = yaml.load(cf, Loader=yaml.FullLoader)

    rad = Radar(config)

    return rad


def test_init(radar):
    assert radar.radar_name == 'Poker Flat'
    assert radar.radar_abbrev == 'PFISR'
    assert radar.site_lat == 65.13
    assert radar.site_lon == -147.47
    assert radar.site_alt == 213.
    assert radar.integration_period == 60.


def test_beams_from_beam_codes(radar):
    beamcodes = [64016, 64157, 64964, 65066]
    truth_az = [14.04, -154.3, -34.69, 75.03]
    truth_el = [90., 77.5, 66.09, 65.56]
    truth_ksys = [np.nan, np.nan, np.nan, np.nan]
    azimuth, elevation, ksys = radar.beams_from_beam_codes(beamcodes)
    np.testing.assert_allclose(azimuth, truth_az)
    np.testing.assert_allclose(elevation, truth_el)
    np.testing.assert_allclose(ksys, truth_ksys)


def test_beams_from_az_el(radar):
    azimuth = [14.04, -154.3, -34.69, 75.03]
    elevation = [90., 77.5, 66.09, 65.56]
    truth_bc = [90001, 90002, 90003, 90004]
    truth_ksys = [np.nan, np.nan, np.nan, np.nan]
    beamcodes, ksys = radar.beams_from_az_el(azimuth, elevation)
    np.testing.assert_allclose(beamcodes, truth_bc)
    np.testing.assert_allclose(ksys, truth_ksys)


def test_caculate_acf_gates(radar):
    slant_range_config = [80000., 800000., 3000.]
    truth_sr = np.arange(*slant_range_config)
    truth_lat = list()
    truth_lon = list()
    truth_alt = list()
    for az, el in zip(radar.beam_azimuth, radar.beam_elevation):
        lat, lon, alt = pm.aer2geodetic(
            az, el, truth_sr, radar.site_lat, radar.site_lon, radar.site_alt)
        truth_lat.append(lat)
        truth_lon.append(lon)
        truth_alt.append(alt)

    slant_range, lat, lon, alt = radar.calculate_acf_gates(slant_range_config)

    np.testing.assert_allclose(slant_range, truth_sr)
    np.testing.assert_allclose(lat, truth_lat)
    np.testing.assert_allclose(lon, truth_lon)
    np.testing.assert_allclose(alt, truth_alt)


def test_calculate_gates(radar):

    altbins = np.arange(100., 800., 50.) * 1000.
    truth_sr = np.empty((len(radar.beam_codes), len(altbins) - 1))
    truth_lat = np.empty((len(radar.beam_codes), len(altbins) - 1))
    truth_lon = np.empty((len(radar.beam_codes), len(altbins) - 1))
    truth_alt = np.empty((len(radar.beam_codes), len(altbins) - 1))
    for b in range(len(radar.beam_codes)):
        for i in range(len(altbins) - 1):
            abidx = np.argwhere((radar.acf_alt[b] >= altbins[i]) & (
                radar.acf_alt[b] < altbins[i + 1])).flatten()
            # supress 'Mean of empty slice' warning so it doesn't clutter
            # output - we expect empty slices at high altitudes
            with warnings.catch_warnings():
                warnings.filterwarnings(action='ignore')
                sr = np.mean(radar.acf_slant_range[abidx])
            lat, lon, alt = pm.aer2geodetic(
                radar.beam_azimuth[b], radar.beam_elevation[b], sr,
                radar.site_lat, radar.site_lon, radar.site_alt)
            truth_sr[b, i] = sr
            truth_lat[b, i] = lat
            truth_lon[b, i] = lon
            truth_alt[b, i] = alt

    slant_range, lat, lon, alt = radar.calculate_gates(
        [[100000., 800000., 50000.]])

    np.testing.assert_allclose(slant_range, truth_sr)
    np.testing.assert_allclose(lat, truth_lat)
    np.testing.assert_allclose(lon, truth_lon)
    np.testing.assert_allclose(alt, truth_alt)


def test_kvec(radar):

    kvec = radar.kvec_all_gates()

    truth_kvec = list()
    for sr, az, el in zip(
            radar.slant_range, radar.beam_azimuth, radar.beam_elevation):
        # aer -> uvw
        e, n, u = pm.aer2enu(az, el, 1.0)
        u, v, w = pm.enu2uvw(e, n, u, radar.site_lat,
                             radar.site_lon, radar.site_alt)
        # aer(slant rng) -> new position
        lat, lon, alt = pm.aer2geodetic(
            az, el, sr, radar.site_lat, radar.site_lon, radar.site_alt)
        # uvw(new position) -> aer
        ke, kn, ku = pm.uvw2enu(u, v, w, lat, lon)
        truth_kvec.append(np.array([ke, kn, ku]).T)

    np.testing.assert_allclose(kvec, truth_kvec)
