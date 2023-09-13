# radar.py
import numpy as np
import pymap3d as pm
from importlib_resources import files
import warnings


class Radar(object):

    def __init__(self, config):

        # self.read_config(config)

        # These can be initalization parameters
        self.radar_name = config['RADAR']['full_name']
        self.radar_abbrev = config['RADAR']['abbreviation']
        (self.site_lat,
         self.site_lon,
         self.site_alt) = config['RADAR']['site_coords']

        beam_bc1 = config['RADAR'].get('beamcodes', [])
        beam_az1, beam_el1, beam_ksys1 = self.beams_from_beam_codes(beam_bc1)

        beam_az2 = config['RADAR'].get('beam_azimuth', [])
        beam_el2 = config['RADAR'].get('beam_elevation', [])
        beam_bc2, beam_ksys2 = self.beams_from_az_el(beam_az2, beam_el2)

        # combine both beam arrays
        self.beam_codes = np.concatenate((beam_bc1, beam_bc2))
        self.beam_azimuth = np.concatenate((beam_az1, beam_az2))
        self.beam_elevation = np.concatenate((beam_el1, beam_el2))
        self.beam_ksys = np.concatenate((beam_ksys1, beam_ksys2))

        alt_bin_config = config['RADAR']['altitude_bins']
        slant_range_config = config['RADAR']['acf_slant_range']

        (self.acf_slant_range,
         self.acf_lat,
         self.acf_lon,
         self.acf_alt) = self.calculate_acf_gates(slant_range_config)
        (self.slant_range,
         self.lat,
         self.lon,
         self.alt) = self.calculate_gates(alt_bin_config)

        self.integration_period = config['RADAR']['integration_period']

    def beams_from_beam_codes(self, beamcodes):
        # beams defined by standard beam code (beamcode files in package data)
        bc_file = files('amisrsynthdata.beamcodes').joinpath(
            'bcotable_{}.txt'.format(
                self.radar_abbrev.lower().replace('-', '')))
        bc_data = np.loadtxt(bc_file)

        idx = np.where(np.in1d(bc_data[:, 0], beamcodes))[0]
        beam_azimuth = bc_data[idx, 1]
        beam_elevation = bc_data[idx, 2]
        beam_ksys = np.full(len(idx), np.nan)

        return beam_azimuth, beam_elevation, beam_ksys

    def beams_from_az_el(self, beam_azimuth, beam_elevation):
        beamcodes = np.arange(len(beam_azimuth)) + 90000
        beam_ksys = np.full(len(beam_azimuth), np.nan)

        return beamcodes, beam_ksys

    def calculate_acf_gates(self, slant_range_config):

        # form slant range bins
        slant_range_p = np.arange(*slant_range_config)
        lat_p, lon_p, alt_p = pm.aer2geodetic(
            self.beam_azimuth[:, None], self.beam_elevation[:, None],
            slant_range_p[None, :],
            self.site_lat, self.site_lon, self.site_alt)
        return slant_range_p, lat_p, lon_p, alt_p

    def calculate_gates(self, alt_bins_config):
        # form list of altitude bins
        altbins = list()
        for segment in alt_bins_config:
            altbins.extend(np.arange(*segment))

        # supress 'Mean of empty slice' warning so it doesn't clutter output -
        # we expect empty slices at high altitudes
        with warnings.catch_warnings():
            warnings.filterwarnings(
                action='ignore', message='Mean of empty slice')
            slant_range = np.array([[np.nanmean(
                np.where((beam >= altbins[i]) & (beam < altbins[i + 1]),
                         self.acf_slant_range, np.nan))
                for i in range(len(altbins) - 1)] for beam in self.acf_alt])
        lat, lon, alt = pm.aer2geodetic(
            self.beam_azimuth[:, None], self.beam_elevation[:, None],
            slant_range, self.site_lat, self.site_lon, self.site_alt)
        return slant_range, lat, lon, alt

    def kvec_all_gates(self):

        ke, kn, ku = pm.aer2enu(self.beam_azimuth, self.beam_elevation, 1.)
        kx, ky, kz = pm.enu2uvw(ke, kn, ku, self.site_lat, self.site_lon)
        ke_all, kn_all, ku_all = pm.uvw2enu(
            kx[:, None], ky[:, None], kz[:, None], self.lat, self.lon)
        kvec = np.array([ke_all, kn_all, ku_all]).transpose(1, 2, 0)

        return kvec
