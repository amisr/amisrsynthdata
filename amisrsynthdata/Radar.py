# Radar.py
import numpy as np
import pymap3d as pm
from importlib_resources import files

# NEEDS MAJOR REFACTORING FOR EFFICIENCY/MODULARIZING

class Radar(object):

    def __init__(self, config):

        self.read_config(config)

        # beams defined by standard beam code (beamcode files in package data)
        bc_file = files('amisrsynthdata.beamcodes').joinpath('bcotable_{}.txt'.format(self.radar_abbrev.lower().replace('-','')))
        bc_data = np.loadtxt(bc_file)


        idx = np.where(np.in1d(bc_data[:,0],self.beamcodes))[0]
        beams_bc = np.array([bc_data[idx,0], bc_data[idx,1], bc_data[idx,2], np.full(len(idx), np.nan)]).T

        # beams defined by azimuth and elevation
        beams_ae = np.array([np.arange(len(self.beam_azimuth))+90000, self.beam_azimuth, self.beam_elevation, np.full(len(self.beam_azimuth), np.nan)]).T

        # combine both beam arrays
        self.beam_codes = np.concatenate((beams_bc, beams_ae), axis=0)

        az = self.beam_codes[:,1]
        el = self.beam_codes[:,2]

        # form list of altitude bins
        altbins = list()
        for segment in self.altitude_bins:
            altbins.extend(np.arange(*segment))

        # form slant range bins
        self.slant_range_p = np.arange(*self.acf_slant_range)
        self.lat_p, self.lon_p, self.alt_p = pm.aer2geodetic(az[:,None], el[:,None], self.slant_range_p[None,:], self.site_lat, self.site_lon, self.site_alt)

        self.slant_range = np.array([[np.nanmean(np.where((beam>=altbins[i]) & (beam<altbins[i+1]), self.slant_range_p, np.nan)) for i in range(len(altbins)-1)] for beam in self.alt_p])
        self.lat, self.lon, self.alt = pm.aer2geodetic(az[:,None], el[:,None], self.slant_range, self.site_lat, self.site_lon, self.site_alt)

        ke, kn, ku = pm.aer2enu(az, el, 1.0)
        self.kvec = np.array([ke, kn, ku]).T


    def kvec_all_gates(self):

        kx, ky, kz = pm.enu2uvw(self.kvec[:,0], self.kvec[:,1], self.kvec[:,2], self.site_lat, self.site_lon)
        ke, kn, ku = pm.uvw2enu(kx[:,None], ky[:,None], kz[:,None], self.lat, self.lon)
        kvec = np.array([ke, kn, ku]).transpose(1,2,0)

        return kvec


    def read_config(self, config):

        self.radar_name = config['RADAR']['full_name']
        self.radar_abbrev = config['RADAR']['abbreviation']
        self.site_lat, self.site_lon, self.site_alt = config['RADAR']['site_coords']
        self.beamcodes = config['RADAR'].get('beamcodes', [])
        self.beam_azimuth = config['RADAR'].get('beam_azimuth', [])
        self.beam_elevation = config['RADAR'].get('beam_elevation', [])
        self.altitude_bins = config['RADAR']['altitude_bins']
        self.acf_slant_range = config['RADAR']['acf_slant_range']
        self.integration_period = config['RADAR']['integration_period']
