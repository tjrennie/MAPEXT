#!/usr/bin/env python3
#===============================================================================
# SciptName     : mapClass.py
# ScriptDir     : /src/core/mapClass.py
# Author(s)     : T.J.Rennie
# Description   : astroMap class for holding astronomical maps in processing
#===============================================================================

# P-I: IMPORTS
from astropy.io import fits
from astropy_healpix import HEALPix
import astropy.units as u
from astropy.wcs import WCS
from astropy.cosmology import Planck15
import healpy
import numpy as np
import matplotlib.pyplot as plt
from astropy.convolution import convolve, Gaussian2DKernel
from skimage.measure import block_reduce
from copy import deepcopy
import reproject
import os

# I: PROCESS CLASS
class astroMap():

    def __init__(self,
                 survey=None,name=None,reference=None,
                 filename=None,
                 frequency=None,wavelength=None,beamwidth=None,
                 resolution=None,projection=None,unit=None,error={},
                 note=[],
                 test=False):
        if test: # TEST VARIABLES FOR PIPELINE CHECKS
            # NAMES
            self.SURV = 'TEST'
            self.NAME = 'TEST'
            self.REFS = 'Rennie 2022'
            # FILE REFERENCES
            self.FILE = 'PLACEHOLDER'
            # OBSERVATION
            self.FREQ = 30.5e9*u.Hz
            self.WLEN = 3e8*(u.m*u.Hz)/self.FREQ
            self.BEAM = [4.5*u.arcmin,4.5*u.arcmin]
            self.BMEA = (np.pi*self.BEAM[0]*self.BEAM[1]).to(u.sr)
            self.BAND = 1e9*u.Hz
            # MAP INFO
            self.RESO = 4.5*u.arcmin
            self.PROJ = 'CART' #CART or HPX
            self.ERRS = {'CALIBRATION':1*u.pc}
            # AOB
            self.NOTE = note
            # MAP
            self.MAP = np.ones((101,201))*(u.MJy/u.sr)
            self.WCS = WCS(naxis=2)
            self.WCS.wcs.cdelt = [1/60,1/60]
            self.WCS.wcs.crpix = [100,50]
            self.WCS.wcs.crval = [30,0]
            self.WCS.wcs.ctype = ['GLON-CYP','GLAT-CYP']
        else:
            # NAMES
            self.SURV = survey
            self.NAME = name
            self.REFS = reference
            # FILE REFERENCES
            self.FILE = filename
            # OBSERVATION
            if frequency != None:
                self.FREQ = frequency
                self.WLEN = 3e8*(u.m*u.Hz)/self.FREQ
            else:
                self.WLEN = wavelength
                self.FREQ = 3e8*(u.m*u.Hz)/self.WLEN
            self.BEAM = beamwidth
            self.BMEA = (np.pi*self.BEAM[0]*self.BEAM[1]).to(u.sr)
            # MAP INFO
            self.RESO = resolution
            self.PROJ = projection
            self.ERRS = error
            # AOB
            self.NOTE = note
            self.WCS, self.MAP = self.load_map()
            if unit=='KCMB':
                self.MAP *= u.K
                self.UNIT = 'KCMB'
            else:
                self.MAP *= unit.value
                self.MAP *= unit.unit
                self.UNIT = self.MAP.unit
        self.ID = self.SURV+'_'+self.NAME

    def load_map(self,
                 n_wcs=0,
                 n_map=0):
        if self.PROJ[:3] == 'WCS':
            f = fits.open(self.FILE)
            map = f[n_map].data[:]
            map[map==0] = np.nan
            wcs = WCS(f[n_wcs].header)[:]
        elif self.PROJ[:3] == 'HPX':
            f = fits.open(self.FILE)
            wcs = {}
            wcs['NSIDE'] = f[1].header['NSIDE']
            wcs['ORDER'] = f[1].header['ORDERING'].lower()
            if self.SURV.lower() == 'parkes':
                # print('PARKES: FORCE ORDERING:RING')
                wcs['ORDER'] = 'ring'
            map = np.array(healpy.read_map(self.FILE))
            map[map<=-1.5e30] = np.nan
            map[map==0] = np.nan
        return wcs, map

    def convert_unit(self,
                     units=u.Jy/u.sr,
                     update=True):
        if units==self.UNIT:
            map_final = self.MAP[:]
        else:
            # first convert to Jy.sr
            if self.UNIT in [u.W/(u.m**2 * u.sr)]:
                map_int = self.MAP.to(u.Jy/u.sr, equivalencies=u.spectral_density(self.FREQ))[:]
            elif self.UNIT in [u.K]:
                map_int = self.MAP.to(u.Jy/u.sr, equivalencies=u.brightness_temperature(self.FREQ))[:]
            elif self.UNIT in [u.Jy/u.sr]:
                map_int = self.MAP[:]
            elif self.UNIT in [u.Jy/u.beam]:
                map_int = self.MAP.to(u.Jy/u.sr, equivalencies=u.beam_angular_area(self.BMEA))[:]
            elif self.UNIT == 'KCMB':
                print('KCMB')
                map_int = self.MAP.to(u.Jy/u.sr, equivalencies=u.thermodynamic_temperature(self.FREQ, Planck15.Tcmb0))
            else:
                print('UNIT NOT RECOGNISED {}'.format(self.UNIT))

            # convert to final units
            # first convert to Jy.sr
            if units in [u.W/(u.m**2 * u.sr)]:
                map_final = map_int.to(units, equivalencies=u.spectral_density(self.FREQ))[:]
            elif units in [u.K]:
                map_final = map_int.to(units, equivalencies=u.brightness_temperature(self.FREQ))[:]
            elif units in [u.Jy/u.sr]:
                map_final = map_int[:]
            elif units in [u.Jy/u.beam]:
                map_final = map_int.to(units, equivalencies=u.beam_angular_area(self.BMEA))[:]
            elif units in [u.Jy/(u.pixel**2)]:
                map_final = map_int.value * u.Jy *abs(np.radians(self.WCS.wcs.cdelt[0])*np.radians(self.WCS.wcs.cdelt[1]))/(u.pixel**2)
            else:
                print('UNIT NOT RECOGNISED {}'.format(self.UNIT))

            # remove intermediate map and output
            del map_int

        if update:
            self.MAP = map_final
            self.NOTE.append('UNIT CONVERTED FROM {} TO {}'.format(self.UNIT, units))
            self.UNIT = units
            return
        else:
            return map_final, self.WCS

    def rtn_coords(self):
        if type(self.WCS) is dict:
            hp = HEALPix(nside=self.WCS['NSIDE'],
                         order=self.WCS['ORDER'],
                         frame=Galactic())
            coords = hp.healpix_to_skycoord(np.arange(self.map.shape[0]))
            return coords
        else:
            x,y = np.meshgrid(np.arange(self.MAP.shape[1]),
                              np.arange(self.MAP.shape[0]),
                              sparse=True)
            coords = self.WCS.pixel_to_world(x,y)
            return coords

    def downsample(self,
                   target_pix=1*u.arcmin,
                   minScale=1,
                   update_map=True):
        current_pix = abs(self.WCS.wcs.cdelt[0])*u.degree.to(u.arcmin)
        sf = int((target_pix/(self.WCS.wcs.cdelt[0]*u.degree)).decompose().value)
        if sf>minScale:
            print('Downscaling by a factor of {}'.format(sf))
            new_map = block_reduce(self.MAP.value,
                                   block_size=(sf,sf),
                                   func=np.nanmean,
                                   cval=np.nan)*self.UNIT
            new_wcs = WCS(naxis=2)
            new_wcs.wcs.crpix = self.WCS.wcs.crpix/sf
            new_wcs.wcs.cdelt = self.WCS.wcs.cdelt*sf
            new_wcs.wcs.crval = self.WCS.wcs.crval
            new_wcs.wcs.ctype = self.WCS.wcs.ctype
        else:
            if self.genParams['verbose']:
                print('NO RESCALING REQUIRED')

        if update_map:
            self.MAP = new_map
            self.NOTE.append('MAP DOWNSAMPLED FROM {} TO {}'.format(current_pix,target_pix))
            self.WCS = new_wcs
            return
        else:
            return new_map, new_wcs

    def smooth_to(self,
                  target_res=5*u.arcmin,
                  update_map=True):
        if self.RESO < target_res:
            map_in = deepcopy(self.MAP[:])
            map_in[np.isfinite(map_in)==False] = 0.
            std = np.sqrt((target_res**2-self.RESO**2)/(8*np.log(2)))
            sd = abs((std.to(u.degree).value)/self.WCS.wcs.cdelt[0])
            kernel = Gaussian2DKernel(x_stddev=sd)
            map_final = convolve(map_in,kernel)
            map_final[map_in==0.] = np.nan
        else:
            map_final = np.zeros(map_in.shape)
            map_final[:] = np.nan
            print('aperture error')
        if update_map:
            self.MAP = map_final
            self.NOTE.append('MAP SMOOTHED FROM {} TO {}'.format(self.RESO, target_res))
            self.RESO = target_res
            return
        else:
            return map_final, target_resolution

    def quickplot(self,map_kwargs={}):
        ax = plt.subplot(projection=self.WCS)
        ax.imshow(self.MAP.value,**map_kwargs)
        ax.coords['glon'].set_major_formatter('d.dd')
        ax.coords['glat'].set_major_formatter('d.dd')
        plt.show()

    def reproject(self,
                  cdelt=None,ctype=None,crpix=None,crval=None,shape_out=None,
                  healpix=None,nside=2048,order='ring',
                  update=True):
        # define initial projection
        if type(self.WCS) is WCS:
            initial_type = 'wcs'
        else:
            initial_type = 'hpx'
        # define final projection
        if healpix == None:
            final_type = 'wcs'
            # create new wcs to go to
            new_wcs = WCS(naxis=2)
            if cdelt!=None:
                new_wcs.wcs.cdelt = cdelt
            else:
                new_wcs.wcs.cdelt = self.WCS.wcs.cdelt

            if ctype!=None:
                new_wcs.wcs.ctype = ctype
            else:
                new_wcs.wcs.ctype = self.WCS.wcs.ctype

            if crpix!=None:
                new_wcs.wcs.crpix = crpix
            else:
                new_wcs.wcs.crpix = self.WCS.wcs.crpix

            if crval!=None:
                new_wcs.wcs.crval = crval
            else:
                new_wcs.wcs.crval = self.WCS.wcs.crval

            if shape_out==None:
                shape_out=[int(new_wcs.wcs.crpix[1]*2-1),
                           int(new_wcs.wcs.crpix[0]*2-1)]
        else:
            final_type = 'hpx'
            new_wcs = {}
            new_wcs['NSIDE'] = healpix
            new_wcs['ORDERING'] = order

        # case-by-case sampling
        if (initial_type=='wcs') & (final_type=='wcs'):
            new_map,footprint=reproject.reproject_interp((self.MAP.value,
                                                          self.WCS),
                                                         new_wcs,
                                                         shape_out=shape_out)

        elif (initial_type=='wcs') & (final_type=='hpx'):
            new_map,footprint=reproject.reproject_to_healpix((self.MAP.value,
                                                              self.WCS),
                                                             coord_system_out='galactic',
                                                             nside=new_wcs['NSIDE'],
                                                             nested=new_wcs['ORDER']=='nested')

        elif (initial_type=='hpx') & (final_type=='wcs'):
            new_x,new_y = np.meshgrid(np.arange(shape_out[0]),
                                  np.arange(shape_out[1]),
                                  sparse=True)
            new_coords = new_wcs.pixel_to_world(new_x,new_y)
            #HEALPIX projection
            hp = HEALPix(nside=self.WCS['NSIDE'],
                         order='RING',#self.WCS['ORDER'],
                         frame='galactic')
            new_coords_hpx = hp.lonlat_to_healpix(new_coords.l,
                                                  new_coords.b)
            new_map = np.zeros(shape_out)
            new_map[:] = np.nan
            m = self.MAP.value
            m[m==0] = np.nan
            m[np.isfinite(m) == False] = np.nan
            if self.MAP[0] != False:
                new_map = self.MAP.value[new_coords_hpx]
            else:
                hpx1 = fits.open(self.FILE)[1].data
                idx = hpx1['PIXEL']
                data = hpx1['SIGNAL']
                new_map = np.zeros(shape_out)
                i = np.searchsorted(idx,new_coords_hpx).astype(int)
                new_map = data[i]
            footprint = np.isfinite(new_map)

            # new_map,footprint=reproject.reproject_from_healpix(self.FILE,
            #                                                    new_wcs,
            #                                                    shape_out=shape_out)
        else:
            print('TRANSFORM NOT CURRENTLY SUPPORTED')
            return

        new_map[footprint==0] = np.nan
        new_map *= self.MAP.unit

        if update:
            self.MAP = new_map
            self.WCS = new_wcs
            self.NOTE.append('MAP REPROJECTED TO CURRENT WCS')
            return
        else:
            return new_map, new_wcs

    def return_dictionary(self):
        out_dict =  {'general':    {'SURV':    self.SURV,
                                    'NAME':    self.NAME,
                                    'REFS':    self.REFS,
                                    'FILE':    self.FILE},
                     'spectral':   {'FREQ':    self.FREQ.to(u.Hz).value,
                                    'WLEN':    self.WLEN.to(u.m).value},
                     'beam':       {'BEAM':    [self.BEAM[0].to(u.arcmin).value,
                                                self.BEAM[1].to(u.arcmin).value],
                                    'BMEA':    self.BMEA.to(u.sr).value},
                     'mapspace':   {'RESO':    self.RESO.to(u.arcmin).value,
                                    'PROJ':    self.PROJ,
                                    'ERRS':    self.ERRS},
                     'note':        np.array(self.NOTE,dtype='<S200')}
        return out_dict

    def save_out(self,ow=True,suffix='proc20220308'):
        hdu = fits.PrimaryHDU(self.MAP.value,header=self.WCS.to_header())
        hdul = fits.HDUList(hdu)
        # Add keywords
        kw = self.return_dictionary()
        for _ in ['general','spectral','beam','mapspace']:
            for __ in kw[_]:
                hdul[0].header[__] = str(kw[_][__])
        for _l,l in enumerate(self.NOTE):
            hdul[0].header['COMMENT'] = self.NOTE[_l]
        # Save out
        outfilename = '{}_{}.fits'.format(self.ID,suffix)
        try:
            hdul.writeto(outfilename,overwrite=ow)
        except:
            print('File already exists - not overwritten')
        del hdul
        return
