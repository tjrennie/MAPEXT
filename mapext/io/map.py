# from skimage.features import hessian_matrix, hessian_matrix_eigvals
import numpy as np
import astropy.units as u
import matplotlib.pyplot as plt
from reproject import reproject_interp
from astropy.wcs import wcs

# def detect_ridges(gray, sigma=1.0):
#     H_elems = hessian_matrix(gray, sigma=sigma, order='rc')
#     maxima_ridges, minima_ridges = hessian_matrix_eigvals(H_elems)
#     return maxima_ridges, minima_ridges


from mapext.core.processClass import Process
from mapext.core.mapClass import astroMap
from mapext.core.srcClass import astroSrc
from mapext.io.outFile import outFileHandler

class cutout(Process):

    def run(self,map,src_lst,dist=15*u.arcmin):
        # load map in Jy/pix units
        MAP, WCS = map.convert_unit(units=u.Jy/(u.pixel**2),update=False)
        shape_out = [int(dist.to(u.arcmin).value*2 + 1),
                     int(dist.to(u.arcmin).value*2 + 1)]
        crpix_out = [int(dist.value),
                     int(dist.value)]

        for s in src_lst:
            wcs_out = wcs.WCS(naxis=2)
            wcs_out.wcs.cdelt = [-1/60,1/60]
            wcs_out.wcs.crval = [s.COORD.l.degree,s.COORD.b.degree]
            wcs_out.wcs.crpix = crpix_out
            wcs_out.wcs.ctype = ['GLON-CYP','GLAT-CYP']
            map_out, footprint = reproject_interp((MAP.value,WCS),wcs_out,shape_out=shape_out)
            map_out *= MAP.unit

            if 'OUTFILE' in self.RUNDICT:
                self.RUNDICT['OUTFILE'].save_map(s,map_out,wcs_out,map.ID)

        # i1, i2 = detect_ridges(MAP.value, sigma = 3.0)

        # fig, axs = plt.subplots(2,2)
        # axs[0,0].imshow(MAP.value)
        # axs[1,0].imshow(i1)
        # axs[1,1].imshow(i2)
        #
        # plt.show()
