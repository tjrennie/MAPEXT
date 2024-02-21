import numpy as np
from astropy.coordinates import SkyCoord
from tqdm import tqdm
from astropy import units as u
from mapext.core import astroMap, astroSrc

__all__ = ['aperture_photometry']

def aperture_photometry(mList, sList, verbose=False, progress=False, aperture=9, circles=[1,1.333,1.667]):
    """ Aperture photometry function for circular apertures and a ring-shaped annulus.

    :param mList: List of maps to reduce
    :type mList: list of mapext.astroMap
    :param sList: list of sources to interrogate
    :type sList: list of mapext.astroSrc
    :param verbose: Set verbosity of function, defaults to False
    :type verbose: bool, optional
    :param progress: Use tqdm to see progress through loops, defaults to False
    :type progress: bool, optional
    :param aperture: Set aperture size, defaults to 9
    :type aperture: float, optional
    :param circles: Set ratio of aperture to inner annulus to outer annulus, defaults to [1,4/3,5/3] (so the area of the aperture and annulus are the same)
    :type circles: list, optional
    
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
        sCenters = [_.coord for _ in sList]
    elif type(sList[0]) == list:
        sCenters = [None]*len(sList)
        for _src,src in enumerate(sList):
            sCenters[_src] = SkyCoord(src[0], src[1], unit='degree', frame='galactic')
    else:
        raise ValueError('Review source list format')
    # iter through maps:
    if progress:
        maps_pbar = tqdm(total=nmaps, desc="Maps")
    for aMap in mList:
        stokes = aMap.stokes
        for stokes_comp in stokes:
            stokes_map = getattr(aMap, stokes_comp, None)
            stokes_map.convert_units(new_units='Jy/pixel^2')
            coords = stokes_map.get_pixel_coords()
            r = np.zeros([len(sCenters),stokes_map.data.shape[0], stokes_map.data.shape[1]])
            for _c,c in enumerate(sCenters):
                r[_c] = c.separation(coords).arcmin
            radii = aperture*np.array(circles)
            mask = (r<radii[0]).astype(int)
            mask[np.all([r>=radii[1],r<radii[2]],axis=0)] = 2
            # Aperture statistics
            ap_sum = np.nansum((mask==1).astype(float)*stokes_map.data, axis=(-1,-2))
            ap_hit = np.nansum((mask==1).astype(float), axis=(-1,-2))
            # Annulus statistics
            s_temp = (mask==2).astype(float)*stokes_map.data
            s_temp[s_temp==0] = np.nan
            an_med = np.nanmedian(s_temp, axis=(-1,-2))
            an_std = np.nanstd(s_temp, axis=(-1,-2))
            an_hit = np.nansum((mask==2).astype(float), axis=(-1,-2))
            # final calculations
            s = ap_sum - (ap_hit*an_med)
            s_var = (an_std**2 * ap_hit * (1 + (np.pi/2)*(ap_hit/an_hit))) + (s * stokes_map.cal * 1e-2)**2
            s_e = np.sqrt(s_var)
            output = []
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
                        S, Se, 'apPhoto')])
                if verbose:
                    print(f'{cen.l.deg:7.3f} {cen.b.deg:7.3f} {S:6.2f} {Se:5.2f}')
            if progress:
                maps_pbar.update(1)
    if progress:
        maps_pbar.close()
    return