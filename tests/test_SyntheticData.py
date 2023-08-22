"""
Test of amisrsynthdata SyntheticData class
"""

import pytest
import amisrsynthdata
import yaml
import h5py
import numpy as np
import datetime as dt
import filecmp
import os

@pytest.fixture
def config():
    config_file = 'config.yaml'
    with open(config_file, 'r') as cf:
        config = yaml.load(cf, Loader=yaml.FullLoader)
    return config


@pytest.fixture
def synthdata(config):
    sd = amisrsynthdata.SyntheticData(config)
    return sd

@pytest.fixture
def datafile():
    
    filename = 'synthetic_data.h5'
    h5file = h5py.File(filename, 'r')
    
    return h5file


def test_init(synthdata):
    assert isinstance(synthdata.radar, amisrsynthdata.Radar)
    assert isinstance(synthdata.iono, amisrsynthdata.Ionosphere)
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
    times = np.array([np.arange(st, et, step), np.arange(st+step, et+step, step)]).T
    truth_utime = times.astype('datetime64[s]').astype('int')
    truth_time = times.astype(dt.datetime)
    
    utime, time = synthdata.generate_time_array(dt.datetime.fromisoformat(starttime), dt.datetime.fromisoformat(endtime))

    assert np.array_equal(utime, truth_utime)
    assert np.array_equal(time, truth_time)

def test_generate_radar_measurements(synthdata, config):

    ne, ti, te, vlos, ne_notr = synthdata.generate_radar_measurements()

    kvec = synthdata.radar.kvec_all_gates()
    truth_vlos = np.dot(kvec, config['VELOCITY'][0]['uniform_glat_aligned']['value'])

    print(vlos.shape, truth_vlos.shape)
   
    assert np.allclose(ne, config['DENSITY'][0]['uniform']['value'])
    assert np.allclose(ti, config['ITEMP'][0]['uniform']['value'])
    assert np.allclose(te, config['ETEMP'][0]['uniform']['value'])
    assert np.allclose(vlos, truth_vlos, equal_nan=True)

def test_generate_errors(synthdata, config, datafile):
    truth_ne_err = datafile['FittedParams/dNe'][:]
    truth_ti_err = datafile['FittedParams/Errors'][:,:,:,0,1]
    truth_te_err = datafile['FittedParams/Errors'][:,:,:,-1,1]
    truth_vlos_err = datafile['FittedParams/Errors'][:,:,:,0,3]
    truth_ne_notr_err = datafile['NeFromPower/Ne_NoTr'][:]*datafile['NeFromPower/dNeFrac'][:]

    ne_err, ti_err, te_err, vlos_err, ne_notr_err = synthdata.generate_errors(config['GENERAL']['err_coef'])

    assert np.allclose(ne_err, truth_ne_err, equal_nan=True)
    assert np.allclose(ti_err, truth_ti_err, equal_nan=True)
    assert np.allclose(te_err, truth_te_err, equal_nan=True)
    assert np.allclose(vlos_err, truth_vlos_err, equal_nan=True)
    #assert np.allclose(ne_notr_err, truth_ne_notr_err)

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
    ne_sf = np.sum(np.abs(ne-synthdata.ne)<3*synthdata.ne_err)/ne.size
    ti_sf = np.sum(np.abs(ti-synthdata.ti)<3*synthdata.ti_err)/ti.size
    te_sf = np.sum(np.abs(te-synthdata.te)<3*synthdata.te_err)/te.size
    vlos_sf = np.sum(np.abs(vlos-synthdata.vlos)<3*synthdata.vlos_err)/vlos.size
    ne_notr_sf = np.sum(np.abs(ne_notr-synthdata.ne_notr)<3*synthdata.ne_notr_err)/ne_notr.size

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


def test_generate_beamcodes(synthdata, datafile):

    beamcodes = synthdata.generate_beamcodes()
    truth_beamcodes = datafile['BeamCodes'][:]

    assert np.array_equal(beamcodes, truth_beamcodes, equal_nan=True)



def test_generate_fitted_params(synthdata, datafile):

    fitted_params, fit_info, ne_from_power = synthdata.generate_fitted_params()

    truth_fitted_params = hdf52dict(datafile['FittedParams'])
    truth_fit_info = hdf52dict(datafile['FittedParams/FitInfo'])
    truth_ne_from_power = hdf52dict(datafile['NeFromPower'])

    np.testing.assert_equal(fitted_params, truth_fitted_params)
    np.testing.assert_equal(fit_info, truth_fit_info)
    np.testing.assert_equal(ne_from_power, truth_ne_from_power)

def test_generate_time(synthdata, datafile):

    time = synthdata.generate_time()
    truth_time = hdf52dict(datafile['Time'])
    np.testing.assert_equal(time, truth_time)

def test_generate_geomag(synthdata, datafile):
    geomag = synthdata.generate_geomag()
    truth_geomag = hdf52dict(datafile['Geomag'])
    np.testing.assert_equal(geomag, truth_geomag)

def test_generate_site(synthdata, datafile):
    site = synthdata.generate_site()
    truth_site = hdf52dict(datafile['Site'])
    # Remove 'Name' - not generated by function
    del truth_site['Name']
    np.testing.assert_equal(site, truth_site)

def test_create_hdf5_output(synthdata):
    synthdata.create_hdf5_output('temp_out.h5')
    assert filecmp.cmp('temp_out.h5', 'synthetic_data.h5')
    os.remove('temp_out.h5')

def test_create_summary_plots(synthdata, config):
    synthdata.create_summary_plots(**config['SUMMARY_PLOT'])
    prefix = config['SUMMARY_PLOT']['output_prefix']

    for s in ['ne', 'ti', 'te', 'vlos']:
        plotname = f'{prefix}{s}.png'
        assert os.path.isfile(plotname)
        os.remove(plotname)

