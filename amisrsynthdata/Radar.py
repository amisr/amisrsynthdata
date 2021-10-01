# Radar.py
import numpy as np
import pymap3d as pm

try:
    import ConfigParser as configparser
except ImportError:
    import configparser

# NEEDS MAJOR REFACTORING FOR EFFICIENCY/MODULARIZING

class Radar(object):

    def __init__(self, config):

        self.read_config(config)

        # beams defined by standard beam code
        bc_data = np.loadtxt(self.beamcode_filename)
        idx = np.where(np.in1d(bc_data[:,0],self.beamcodes))[0]
        beams_bc = np.array([bc_data[idx,0], bc_data[idx,1], bc_data[idx,2], np.full(len(idx), np.nan)]).T

        # beams defined by azimuth and elevation
        beams_ae = np.array([np.arange(len(self.beam_azimuth))+90000, self.beam_azimuth, self.beam_elevation, np.full(len(self.beam_azimuth), np.nan)]).T

        # combine both beam arrays
        self.beam_codes = np.concatenate((beams_bc, beams_ae), axis=0)

        az = self.beam_codes[:,1]
        el = self.beam_codes[:,2]


        if len(self.altbins) == 3:
            self.altbins = np.arange(self.altbins[0], self.altbins[1], self.altbins[2])


        self.slant_range = np.arange(self.range_start,self.range_end, self.range_step)  # move start/end range to config file
        self.lat_nb, self.lon_nb, self.alt_nb = pm.aer2geodetic(az[:,None], el[:,None], self.slant_range[None,:], self.site_lat, self.site_lon, self.site_alt)

        self.fit_slant_range = np.array([[np.nanmean(np.where((beam>=self.altbins[i]) & (beam<self.altbins[i+1]), self.slant_range, np.nan)) for i in range(len(self.altbins)-1)] for beam in self.alt_nb])
        self.lat, self.lon, self.alt = pm.aer2geodetic(az[:,None], el[:,None], self.fit_slant_range, self.site_lat, self.site_lon, self.site_alt)

        ke, kn, ku = pm.aer2enu(az, el, 1.0)
        self.kvec = np.array([ke, kn, ku]).T


    def kvec_all_gates(self):

        kx, ky, kz = pm.enu2uvw(self.kvec[:,0], self.kvec[:,1], self.kvec[:,2], self.site_lat, self.site_lon)
        ke, kn, ku = pm.uvw2enu(kx[:,None], ky[:,None], kz[:,None], self.lat, self.lon)
        kvec = np.array([ke, kn, ku]).transpose(1,2,0)

        return kvec


    def read_config(self, config_file):

        config = configparser.ConfigParser()
        config.read(config_file)

        self.site_lat, self.site_lon, self.site_alt = [float(i) for i in config['RADAR'].get('SITE_COORDS').split(',')]
        self.beamcode_filename = config['RADAR'].get('BEAMCODE_FILENAME', None)
        self.beamcodes = [float(i) for i in config['RADAR'].get('BEAMCODES','').split(',') if i]
        self.beam_azimuth = [float(i) for i in config['RADAR'].get('BEAM_AZIMUTH','').split(',') if i]
        self.beam_elevation = [float(i) for i in config['RADAR'].get('BEAM_ELEVATION','').split(',') if i]
        self.altbins = np.array([float(i) for i in config['RADAR'].get('ALTBINS').split(',')])
        self.range_step = config.getfloat('RADAR','RANGE_STEP')
        self.range_start = config.getfloat('RADAR','RANGE_START')
        self.range_end = config.getfloat('RADAR','RANGE_END')
        self.integration_period = config.getfloat('RADAR','INTEGRATION_PERIOD')
        self.vel_error = config.getfloat('RADAR','VEL_ERROR')
        self.radar_name = config.get('RADAR', 'NAME')
