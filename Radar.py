# Radar.py
import numpy as np
# import coord_convert as cc
import pymap3d as pm
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from scipy.stats import binned_statistic

try:
    import ConfigParser as configparser
except ImportError:
    import configparser

# NEEDS MAJOR REFACTORING FOR EFFICIENCY/MODULARIZING

class Radar(object):

    def __init__(self, config):

        self.read_config(config)

        bc_data = np.loadtxt(self.beamcode_filename)
        idx = np.where(np.in1d(bc_data[:,0],self.beamcodes))[0]
        bc = bc_data[idx,0]
        az = bc_data[idx,1]
        el = bc_data[idx,2]
        self.beam_codes = np.array([bc, az, el, np.full(bc.shape, np.nan)]).T
        print('BEAM CODES:', self.beam_codes.shape)


        self.slant_range = np.arange(80.,800., self.range_step)*1000.   # move start/end range to config file
        self.lat_nb, self.lon_nb, self.alt_nb = pm.aer2geodetic(az[:,None], el[:,None], self.slant_range[None,:], self.site_coords[0], self.site_coords[1], self.site_coords[2]*1000.)
        self.x_nb, self.y_nb, self.z_nb = pm.geodetic2ecef(self.lat_nb, self.lon_nb, self.alt_nb)

        print('SLANT RANGE:', self.slant_range.shape, self.alt_nb.shape)

        # srange_fit = np.digitize(self.alt_nb, self.altbins*1000.)
        # srange_fit, _, _ = binned_statistic(self.alt_nb, self.slant_range, statistic='mean', bins=self.altbins*1000.)
        # print(srange_fit)

        self.fit_slant_range = np.array([[np.nanmean(np.where((beam>=self.altbins[i]*1000.) & (beam<self.altbins[i+1]*1000.), self.slant_range, np.nan)) for i in range(len(self.altbins)-1)] for beam in self.alt_nb])
        print('FIT SLANT RANGE:', self.fit_slant_range.shape)

        self.x, self.y, self.z = pm.aer2ecef(az[:,None], el[:,None], self.fit_slant_range, self.site_coords[0], self.site_coords[1], self.site_coords[2]*1000.)

        ke, kn, ku = pm.aer2enu(az, el, 1.0)
        kx, ky, kz = pm.enu2uvw(ke, kn, ku, self.site_coords[0], self.site_coords[1])

        # for a in self.alt_nb:
        #     for i in range(len(self.altbins)-1):
        #         sr_in_bin = np.where((a>=self.altbins[i]*1000.) & (a<self.altbins[i+1]*1000.), self.slant_range, np.nan)
        #         print(self.altbins[i], np.nanmean(sr_in_bin))




        self.X0, self.Y0, self.Z0 = pm.geodetic2ecef(self.site_coords[0], self.site_coords[1], self.site_coords[2]*1000.)
        # self.x, self.y, self.z, kx, ky, kz = self.get_gate_locations(self.beam_codes[:,1], self.beam_codes[:,2], self.range_step)
        print(self.x.shape)
        # form array of k vector at each range gate
        # move this reordering of arrays into the saving synthetic file script?
        self.kvec = np.array([np.tile(np.array([x,y,z]), (self.x.shape[1],1)) for x, y, z in zip(kx,ky,kz)])
        self.kx = np.tile(kx, (self.x.shape[1],1)).T
        self.ky = np.tile(ky, (self.x.shape[1],1)).T
        self.kz = np.tile(kz, (self.x.shape[1],1)).T

        self.lat, self.lon, self.alt, self.ke, self.kn, self.ku = self.geodetic_locations(self.x, self.y, self.z, self.kx, self.ky, self.kz)
        # self.plot_radar()

        # self.x_nb, self.y_nb, self.z_nb, kx, ky, kz = self.get_gate_locations(self.beam_codes[:,1], self.beam_codes[:,2], self.bin_step)
        # _, _, self.alt_nb, _, _, _ = self.geodetic_locations(self.x_nb, self.y_nb, self.z_nb, kx[:,None], ky[:,None], kz[:,None])


    def read_config(self, config_file):

        config = configparser.ConfigParser()
        config.read(config_file)

        getfloat = lambda x: config.getfloat('RADAR', x)
        getarray = lambda x: np.array(eval(config.get('RADAR', x)))
        getstring = lambda x: config.get('RADAR', x)

        self.site_coords = [float(i) for i in config.get('RADAR','SITE_COORDS').split(',')]
        self.beamcode_filename = config.get('RADAR', 'BEAMCODE_FILENAME')
        self.beamcodes = [float(i) for i in config.get('RADAR','BEAMCODES').split(',')]
        self.altbins = np.array([float(i) for i in config.get('RADAR','ALTBINS').split(',')])
        self.range_step = config.getfloat('RADAR','RANGE_STEP')
        self.bin_step = config.getfloat('RADAR','BIN_STEP')
        self.output_filename = config.get('RADAR','OUTPUT_FILENAME')
        self.integration_period = config.getfloat('RADAR','INTEGRATION_PERIOD')
        self.vel_error = config.getfloat('RADAR','VEL_ERROR')
        self.radar_name = config.get('RADAR', 'NAME')


    def get_gate_locations(self, az, el, rs):

        # create array of ranges
        ranges = np.arange(80.,800., rs)*1000.

        # find E, N, U components of k (find k vector in geodetic coordinates)
        ke, kn, ku = pm.aer2enu(az, el, 1.0)
        # az = np.array(az)*np.pi/180.
        # el = np.array(el)*np.pi/180.
        # ke = np.cos(el)*np.sin(az)
        # kn = np.cos(el)*np.cos(az)
        # ku = np.sin(el)

        # convert geodetic k vector to ECEF
        kx, ky, kz = pm.enu2uvw(ke, kn, ku, self.site_coords[0], self.site_coords[1])

        # calculate position of each range gate in ECEF
        x = kx[:,None]*ranges + self.X0
        y = ky[:,None]*ranges + self.Y0
        z = kz[:,None]*ranges + self.Z0

        return x, y, z, kx, ky, kz



    def geodetic_locations(self, x, y, z, kx, ky, kz):
        print(x.shape, kx.shape)
        # calculate gate position and k vectors in geodetic coordinates
        lat, lon, alt = pm.ecef2geodetic(x, y, z)
        ke, kn, ku = pm.uvw2enu(kx, ky, kz, lat, lon)
        return lat, lon, alt, ke, kn, ku

    def plot_radar(self):

        fig = plt.figure(figsize=(10,10))
        ax = fig.add_subplot(111,projection='3d')

        for x, y, z, kx, ky, kz in zip(self.X, self.Y, self.Z, self.kx, self.ky, self.kz):
            ax.quiver(x, y, z, kx, ky, kz, length=100000, color='orange')


        plt.show()
