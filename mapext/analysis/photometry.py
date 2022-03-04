
import astropy.wcs as wcs
import astropy.units as u
import numpy as np
from astropy.coordinates import SkyCoord

import matplotlib.pyplot as plt

from mapext.core.processClass import Process
from mapext.core.mapClass import astroMap
from mapext.core.srcClass import astroSrc

class aperture(Process):

    def run(self,map,src_lst,apertures=10*u.arcmin):
        # load map in Jy/pix units
        MAP, WCS = map.convert_unit(units=u.Jy/(u.pixel**2),update=False)

        # ensure src_lst is a list of sources and retrieve coordinates of each one
        if type(src_lst) is astroSrc:
            src_lst = [src_lst]
        src_skyCoord = [src.COORD for src in src_lst]
        map_skyCoord = map.rtn_coords()

        coords_r = np.array([map_skyCoord.separation(sc).degree for sc in src_skyCoord])
        mask = np.zeros(coords_r.shape).astype(int)
        apertures_deg = apertures.to(u.degree).value
        mask[coords_r<=apertures_deg] = 1
        mask[np.all([coords_r>=apertures_deg*(4/3),coords_r<=apertures_deg*(5/3)],axis=0)] = 2

        ap_count = np.sum(mask==1,axis=(1,2))*u.pixel**2
        an_count = np.sum(mask==2,axis=(1,2))*u.pixel**2

        ap = (mask==1).astype(float)*MAP[np.newaxis,:,:]
        ap[ap==0] = np.nan

        an = (mask==2).astype(float)*MAP[np.newaxis,:,:]
        an[an==0] = np.nan

        ap_sum   = np.nansum(ap,axis=(1,2))*u.pixel**2
        an_mean  = np.nanmean(an,axis=(1,2))
        an_med   = np.nanmedian(an,axis=(1,2))
        an_std   = np.nanstd(an,axis=(1,2))

        Sv = ap_sum - an_med*ap_count
        Sv_e = (an_std*np.sqrt(an_count + ap_count))

        for _ in range(Sv.shape[0]):
            print(_,' : Sv = ',Sv[_].value,'±',Sv_e[_].value)
            # OUTPUT FLUX INFO
            flux_info = [('ApPhoto',map.SURV,map.NAME,
                          (map.FREQ/u.Hz).value,Sv[_].value,Sv_e[_].value)]
            src_lst[_].add_flux(flux_info)

            with open('{}_{}_apho.txt'.format(src_lst[_].NAME,map.NAME), 'a') as the_file:
                the_file.write('Sv = {} ± {}'.format(Sv[_].value,Sv_e[_].value))
