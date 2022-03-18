from skimage.feature import hessian_matrix, hessian_matrix_eigvals
import numpy as np
import astropy.units as u
import matplotlib.pyplot as plt
from reproject import reproject_interp
from astropy.wcs import wcs
from matplotlib.colors import SymLogNorm

# def detect_ridges(gray, sigma=1.0):
#     H_elems = hessian_matrix(gray, sigma=sigma, order='rc')
#     maxima_ridges, minima_ridges = hessian_matrix_eigvals(H_elems)
#     return maxima_ridges, minima_ridges


from mapext.core.processClass import Process
from mapext.core.mapClass import astroMap
from mapext.core.srcClass import astroSrc
from mapext.io.outFile import outFileHandler
from mapext.core.usefulFunctions import set_lims

tab20_clrs = np.loadtxt('mapext/core/tab20.txt', delimiter=' ')/255
tab10_clrs = np.loadtxt('mapext/core/tab10.txt', delimiter=' ')/255

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

class map_with_srcs(Process):

    def run(self,map,srcs):
        # load map in Jy/pix units
        MAP, WCS = map.convert_unit(units=u.Jy/(u.pixel**2),update=False)

        plt.figure()
        ax = plt.subplot(projection=WCS)

        vmin,vmax = set_lims(MAP.value,lower=0.05)

        plt.imshow(MAP.value,cmap='Greys',vmin=vmin,vmax=vmax)
        # if type(srcs) is dict:
        for __,_ in enumerate(srcs.keys()):
            slst = srcs[_]
            for _s,s in enumerate(slst):
                if _s ==0:
                    plt.scatter(s.COORD.l.degree,s.COORD.b.degree,
                                c=[tab10_clrs[__]],s=10, transform=ax.get_transform('world'),label=_)
                else:
                    plt.scatter(s.COORD.l.degree,s.COORD.b.degree,
                                c=[tab10_clrs[__]],s=10, transform=ax.get_transform('world'))
        plt.legend()
        plt.show()



class interest_map(Process):

    def run(self,map,sigma=3.0):

        MAP, WCS = map.convert_unit(units=u.Jy/(u.pixel**2),update=False)

        H_elems = hessian_matrix(MAP.value+1e3, sigma=sigma, order='rc')
        srcs, ridges = hessian_matrix_eigvals(H_elems)

        fig, ax = plt.subplots(1,3,sharex=True,sharey=True,subplot_kw={'projection':WCS})
        ax[0].imshow(-ridges,norm=SymLogNorm(1e-3,vmin=-1e-1,vmax=1e-1),cmap='bwr')
        ax[1].imshow(-srcs,norm=SymLogNorm(1e-3,vmin=-1e-2,vmax=1e-2),cmap='bwr')

        VAL = 10**(np.log10(-ridges) - np.log10(-srcs))

        ax[2].imshow(VAL,cmap='Greens')
        # ax.imshow(srcs,norm=LogNorm(vmin=1e-4,vmax=1e-2),cmap='Greens')
        ax[0].contour(MAP.value,vmin=-0.03,vmax=0.05,levels=60,colors='black',lw=0.5)
        ax[1].contour(MAP.value,vmin=-0.03,vmax=0.05,levels=60,colors='black',lw=0.5)
        ax[2].contour(MAP.value,vmin=-0.03,vmax=0.05,levels=60,colors='black',lw=0.5)

        plt.show()
