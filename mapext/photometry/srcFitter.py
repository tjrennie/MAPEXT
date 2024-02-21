import numpy as np
from astropy.modeling import models, fitting
from astropy.units import Quantity
from tqdm import tqdm
from astropy import units as u
from mapext.core import astroMap, astroSrc
import matplotlib.pyplot as plt

__all__ = ['gaussian_fitter']

def gaussian_fitter(mList, sList, verbose=False, msiz='10 arcmin', progress=False):
    """ Gaussian source fitter.


    :param mList: List of maps to reduce
    :type mList: list of mapext.astroMap
    :param sList: list of sources to interrogate
    :type sList: list of mapext.astroSrc
    :param verbose: Set verbosity of function, defaults to False
    :type verbose: bool, optional
    :param progress: Use tqdm to see progress through loops, defaults to False
    :type progress: bool, optional
    :param msiz: Map size to use in fitting, defaults to '10 arcmin'
    :type msiz: str, optional
    
    :raises ValueError: Review source list format
    """
    # sort list of maps and iterations required
    if type(mList) != list:
        mList = [mList]
    maps_to_iterate = {}
    nmaps = 0
    for _map,map in enumerate(mList):
        maps_to_iterate[_map] = map.stokes
        nmaps += len(map.stokes)
    # sort list of sources into one list of center coordinates
    if type(sList) != list:
        sList = [sList]
    if type(sList[0]) == astroSrc:
        sCenters = [[_.coord.l.degree,_.coord.b.degree] for _ in sList]
    elif type(sList[0]) == list:
        sCenters = [None]*len(sList)
        for _src,src in enumerate(sList):
            sCenters[_src] = [src[0], src[1]]
    else:
        raise ValueError('Review source list format')
    sCenters = np.array(sCenters)
    msiz = Quantity(msiz)
    # iter through maps:
    if progress:
        maps_pbar = tqdm(total=nmaps, desc="Maps")
    for aMap in mList:
        stokes = aMap.stokes
        for stokes_comp in stokes:
            stokes_map = getattr(aMap, stokes_comp, None)
            stokes_map.convert_units(new_units='Jy/pixel^2')
            # Convert longitudes and latitudes to pixel coordinates
            pixel_coords = stokes_map.proj.all_world2pix(np.array(sCenters[:,0]), np.array(sCenters[:,1]), 0)
            pixel_coords = np.round(pixel_coords).astype(int).T
            #
            pixscale = np.abs(stokes_map.proj.wcs.cdelt[0])
            delt = msiz.to(u.degree).value / pixscale
            X, Y = np.meshgrid(np.arange(2*delt + 1) - delt,
                               np.arange(2*delt + 1) - delt)
            R = np.sqrt(X*X + Y*Y)
            mask = R<delt
            std_init = (stokes_map.reso.to(u.degree).value / pixscale / 2.355)
            s = np.zeros(sCenters.shape[0])
            s_e = np.zeros(sCenters.shape[0])
            for _s,(l,b) in enumerate(pixel_coords):
                dat = stokes_map.data[int(b-delt):int(b+delt+2), int(l-delt):int(l+delt+2)]
                model = models.Gaussian2D() + models.Polynomial2D(degree=1)
                model.x_mean_0 = 0
                model.y_mean_0 = 0
                model.amplitude_0 = np.nanmax(dat) - np.nanmedian(dat)
                model.x_stddev_0 = std_init
                model.y_stddev_0 = std_init
                # Perform fit
                fitter = fitting.LevMarLSQFitter(calc_uncertainties=True)
                fit = fitter(model, X[mask], Y[mask], dat[mask])
                try:
                    std = np.sqrt(np.diag(np.array(fit.cov_matrix.cov_matrix)))
                except:
                    std = np.zeros(6)
                s[_s] = 2 * np.pi * fit.amplitude_0.value * fit.x_stddev_0.value * fit.y_stddev_0.value
                s_var = (s[_s]**2 * ( (std[0]/fit.amplitude_0.value)**2 + (std[3]/fit.x_stddev_0.value)**2 + (std[4]/fit.y_stddev_0.value)**2 )) + (s[_s] * stokes_map.cal * 1e-2)**2
                s_e[_s] = np.sqrt(s_var)
            # output = []
            if type(sList[0])!=astroSrc:
                verbose=True
            if verbose:
                print(f'{stokes_map.name:20.20s} : {stokes_map.freq}')
                print(f'GLON    GLAT    S      Se')
            for _s, (S,Se,src,cen) in enumerate(zip(s,s_e,sList,sCenters)):
                if type(src) == astroSrc:
                    src.append_flux_measure([(
                        stokes_map.name,
                        stokes_comp,
                        stokes_map.freq.to(u.Hz).value,
                        stokes_map.reso.to(u.arcmin).value,
                        S, Se, 'gaussianFit')])
                    src.append_model([(
                        stokes_map.name,
                        stokes_comp,
                        stokes_map.freq.to(u.Hz).value,
                        stokes_map.reso.to(u.arcmin).value,
                        ('gaussianModel',
                         fit.amplitude_0.value,
                         cen[0]+(pixscale*fit.x_mean_0.value),
                         cen[1]+(pixscale*fit.y_mean_0.value),
                         fit.x_stddev_0.value,
                         fit.y_stddev_0.value,
                         fit.theta_0.value),
                        ('gaussianModel',
                         std[0],
                         std[1],
                         std[2],
                         std[3],
                         std[4],
                         std[5]),
                        'gaussianFit' )])
                if verbose:
                    print(f'{cen[0]:7.3f} {cen[1]:7.3f} {S:6.2f} {Se:5.2f}')
            if progress:
                maps_pbar.update(1)
    if progress:
        maps_pbar.close()
    return