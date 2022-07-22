# Radar.py
import numpy as np
import pymap3d as pm

# try:
#     import ConfigParser as configparser
# except ImportError:
#     import configparser

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


        self.slant_range_p = np.arange(self.range_start,self.range_end, self.range_step)  # move start/end range to config file
        self.lat_p, self.lon_p, self.alt_p = pm.aer2geodetic(az[:,None], el[:,None], self.slant_range_p[None,:], self.site_lat, self.site_lon, self.site_alt)

        self.slant_range = np.array([[np.nanmean(np.where((beam>=self.altbins[i]) & (beam<self.altbins[i+1]), self.slant_range_p, np.nan)) for i in range(len(self.altbins)-1)] for beam in self.alt_p])
        self.lat, self.lon, self.alt = pm.aer2geodetic(az[:,None], el[:,None], self.slant_range, self.site_lat, self.site_lon, self.site_alt)

        ke, kn, ku = pm.aer2enu(az, el, 1.0)
        self.kvec = np.array([ke, kn, ku]).T


    def kvec_all_gates(self):

        kx, ky, kz = pm.enu2uvw(self.kvec[:,0], self.kvec[:,1], self.kvec[:,2], self.site_lat, self.site_lon)
        ke, kn, ku = pm.uvw2enu(kx[:,None], ky[:,None], kz[:,None], self.lat, self.lon)
        kvec = np.array([ke, kn, ku]).transpose(1,2,0)

        return kvec


    def read_config(self, config):

        # config = configparser.ConfigParser()
        # config.read(config_file)
        print(config['RADAR'])

        self.site_lat, self.site_lon, self.site_alt = config['RADAR']['site_coords']
        self.beamcode_filename = config['RADAR']['beamcode_filename']
        self.beamcodes = config['RADAR']['beamcodes']
        self.beam_azimuth = config['RADAR'].get('beam_azimuth', [])
        self.beam_elevation = config['RADAR'].get('beam_elevation', [])
        self.altbins = config['RADAR']['altbins']
        self.range_step = config['RADAR']['range_step']
        self.range_start = config['RADAR']['range_start']
        self.range_end = config['RADAR']['range_end']
        self.integration_period = config['RADAR']['integration_period']
        self.radar_name = config['RADAR']['name']
