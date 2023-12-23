"""
Test of amisrsynthdata SyntheticData class
"""

import pytest
import yaml
import h5py
import numpy as np
import datetime as dt
import filecmp
import os
import warnings

from amisrsynthdata.syntheticdata import SyntheticData
from amisrsynthdata.ionosphere import Ionosphere
from amisrsynthdata.radar import Radar


@pytest.fixture
def config():
    config_file = os.path.join(os.path.dirname(__file__), 'config.yaml')
    with open(config_file, 'r') as cf:
        config = yaml.load(cf, Loader=yaml.FullLoader)
    return config


@pytest.fixture
def synthdata(config):
    sd = SyntheticData(config)
    return sd


@pytest.fixture
def datafile():

    filename = os.path.join(os.path.dirname(__file__), 'synthetic_data.h5')
    h5file = h5py.File(filename, 'r')

    return h5file


def test_init(synthdata):
    assert isinstance(synthdata.radar, Radar)
    assert isinstance(synthdata.iono, Ionosphere)
    assert isinstance(synthdata.utime, np.ndarray)
    assert isinstance(synthdata.time, np.ndarray)
    assert isinstance(synthdata.ne, np.ndarray)
    assert isinstance(synthdata.ti, np.ndarray)
    assert isinstance(synthdata.te, np.ndarray)
    assert isinstance(synthdata.vlos, np.ndarray)
    assert isinstance(synthdata.ne_notr, np.ndarray)
    assert isinstance(synthdata.ne_err, np.ndarray)
    assert isinstance(synthdata.ti_err, np.ndarray)
    assert isinstance(synthdata.te_err, np.ndarray)
    assert isinstance(synthdata.vlos_err, np.ndarray)
    assert isinstance(synthdata.ne_notr_err, np.ndarray)


def test_generate_time_array(synthdata, config):
    starttime = '2020-01-01T12:00:00'
    endtime = '2020-01-02T02:00:00'

    st = np.datetime64(starttime)
    et = np.datetime64(endtime)
    step = np.timedelta64(int(config['RADAR']['integration_period']), 's')
    times = np.array(
        [np.arange(st, et, step), np.arange(st + step, et + step, step)]).T
    truth_utime = times.astype('datetime64[s]').astype('int')
    truth_time = times.astype(dt.datetime)

    dt_st = dt.datetime.strptime(starttime, '%Y-%m-%dT%H:%M:%S')
    dt_et = dt.datetime.strptime(endtime, '%Y-%m-%dT%H:%M:%S')
    utime, time = synthdata.generate_time_array(dt_st, dt_et)

    np.testing.assert_allclose(utime, truth_utime)
    np.testing.assert_equal(time, truth_time)


def test_generate_radar_measurements(synthdata, config):

    ne, ti, te, vlos, ne_notr = synthdata.generate_radar_measurements()

    kvec = synthdata.radar.kvec_all_gates()
    truth_vlos = np.dot(
        kvec, config['VELOCITY'][0]['uniform_glat_aligned']['value'])

    np.testing.assert_allclose(ne[np.isfinite(ne)],
                               config['DENSITY'][0]['uniform']['value'])
    np.testing.assert_allclose(ti[np.isfinite(ti)],
                               config['ITEMP'][0]['uniform']['value'])
    np.testing.assert_allclose(te[np.isfinite(te)],
                               config['ETEMP'][0]['uniform']['value'])
    np.testing.assert_allclose(vlos, np.broadcast_to(truth_vlos, vlos.shape))
    # Also check that these are the same shape as radar arrays


def test_generate_errors(synthdata, config, datafile):
    truth_ne_err = datafile['FittedParams/dNe'][:]
    truth_ti_err = datafile['FittedParams/Errors'][:, :, :, 0, 1]
    truth_te_err = datafile['FittedParams/Errors'][:, :, :, -1, 1]
    truth_vlos_err = datafile['FittedParams/Errors'][:, :, :, 0, 3]
    truth_ne_notr_err = datafile['NeFromPower/Ne_NoTr'][:] * \
        datafile['NeFromPower/dNeFrac'][:]

    ne_err, ti_err, te_err, vlos_err, ne_notr_err = synthdata.generate_errors(
        config['GENERAL']['rel_err'], config['GENERAL']['err_ref_rng'])

    np.testing.assert_allclose(
        np.broadcast_to(
            ne_err,
            truth_ne_err.shape),
        truth_ne_err)
    np.testing.assert_allclose(
        np.broadcast_to(
            ti_err,
            truth_ti_err.shape),
        truth_ti_err)
    np.testing.assert_allclose(
        np.broadcast_to(
            te_err,
            truth_te_err.shape),
        truth_te_err)
    np.testing.assert_allclose(
        np.broadcast_to(
            vlos_err,
            truth_vlos_err.shape),
        truth_vlos_err)
    np.testing.assert_allclose(
        np.broadcast_to(
            ne_notr_err,
            truth_ne_notr_err.shape),
        truth_ne_notr_err)


def test_noisy_measurements(synthdata):

    ne, ti, te, vlos, ne_notr = synthdata.noisy_measurements()

    # Noisy arrays should be the same shape but not same values
    assert ne.shape == synthdata.ne.shape
    assert not np.array_equal(ne, synthdata.ne)
    assert ti.shape == synthdata.ti.shape
    assert not np.array_equal(ti, synthdata.ti)
    assert te.shape == synthdata.te.shape
    assert not np.array_equal(te, synthdata.te)
    assert vlos.shape == synthdata.vlos.shape
    assert not np.array_equal(vlos, synthdata.vlos)
    assert ne_notr.shape == synthdata.ne_notr.shape
    assert not np.array_equal(ne_notr, synthdata.ne_notr)

    # 3-sigma test
    # 99.7% of data points should be within 3 standard deviations of the mean
    # Test at least 90% data is within this limit
    ne_sf = np.sum(np.abs(ne - synthdata.ne) < 3 * synthdata.ne_err) / ne.size
    ti_sf = np.sum(np.abs(ti - synthdata.ti) < 3 * synthdata.ti_err) / ti.size
    te_sf = np.sum(np.abs(te - synthdata.te) < 3 * synthdata.te_err) / te.size
    vlos_sf = np.sum(np.abs(vlos - synthdata.vlos) <
                     3 * synthdata.vlos_err) / vlos.size
    ne_notr_sf = np.sum(np.abs(ne_notr - synthdata.ne_notr)
                        < 3 * synthdata.ne_notr_err) / ne_notr.size

    assert ne_sf >= 0.9
    assert ti_sf >= 0.9
    assert te_sf >= 0.9
    assert vlos_sf >= 0.9
    assert ne_notr_sf >= 0.9


def hdf52dict(h5):
    out = dict()
    for key, value in h5.items():
        if isinstance(value, h5py._hl.dataset.Dataset):
            try:
                out[key] = value[:]
            except ValueError:
                out[key] = value[()]
        else:
            continue
    return out


def assert_dict_equal(dict1, dict2, rtol=1.e-2, atol=0.):
    np.testing.assert_equal(dict1.keys(), dict2.keys())
    for k in dict1.keys():
        np.testing.assert_allclose(
            dict1[k],
            dict2[k],
            rtol=rtol,
            atol=atol,
            err_msg=f'Error in {k}')


def test_generate_beamcodes(synthdata, datafile):
    beamcodes = synthdata.generate_beamcodes()
    truth_beamcodes = datafile['BeamCodes'][:]
    np.testing.assert_allclose(beamcodes, truth_beamcodes)


def test_generate_fitted_params(synthdata, datafile):

    fitted_params, fit_info, ne_from_power = synthdata.generate_fitted_params()

    truth_fitted_params = hdf52dict(datafile['FittedParams'])
    truth_fit_info = hdf52dict(datafile['FittedParams/FitInfo'])
    truth_ne_from_power = hdf52dict(datafile['NeFromPower'])

    assert_dict_equal(fitted_params, truth_fitted_params)
    assert_dict_equal(fit_info, truth_fit_info)
    assert_dict_equal(ne_from_power, truth_ne_from_power)


def test_generate_time(synthdata, datafile):
    time = synthdata.generate_time()
    truth_time = hdf52dict(datafile['Time'])
    assert_dict_equal(time, truth_time)


def test_generate_geomag(synthdata, datafile):
    geomag = synthdata.generate_geomag()
    truth_geomag = hdf52dict(datafile['Geomag'])
    assert_dict_equal(geomag, truth_geomag)


def test_generate_site(synthdata, datafile):
    site = synthdata.generate_site()
    truth_site = hdf52dict(datafile['Site'])
    # Remove 'Name' - not generated by function
    del truth_site['Name']
    assert_dict_equal(site, truth_site)


def test_create_hdf5_output(synthdata):
    file = 'temp_out.h5'
    synthdata.create_hdf5_output(file)
    truth_file = os.path.join(os.path.dirname(__file__), 'synthetic_data.h5')
    assert os.path.isfile(file)
    if not filecmp.cmp(file, truth_file):
        warnings.warn(
            Warning("Generated file is not identical to comparison file."))
    os.remove(file)

# Won't run on GHActions for CI
# Issue with installing caropy (specifically dependency on proj)


@pytest.mark.skipif(os.getenv("GITHUB_ACTIONS") == "true",
                    reason="Test doesn't work in Github Actions due to "
                           "challenges installing cartopy.")
def test_create_summary_plots(synthdata, config):
    synthdata.create_summary_plots(**config['SUMMARY_PLOT'])
    prefix = config['SUMMARY_PLOT']['output_prefix']

    for s in ['ne', 'ti', 'te', 'vlos']:
        plotname = f'{prefix}{s}.png'
        assert os.path.isfile(plotname)
        os.remove(plotname)
