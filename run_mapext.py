from mapext.core.mapClass import astroMap
from mapext.core.srcClass import astroSrc
from mapext.analysis import photometry,sourcefit
from mapext.io import map
import numpy as np
from astropy.coordinates import SkyCoord
import astropy.units as u
from mapext.io.outFile import outFileHandler
from tqdm import tqdm

mapDict = {
    'MAP1' : {
          'survey' : 'SURVEYNAME',
            'name' : 'MAPNAME',
       'reference' : 'REFERENCE IN BIB FILE',
        'filename' : 'PATHTOFILE',
       'frequency' : ?*u.Hz, # SUPPLY EITHER FREQUENCY OR WAVELENGTH - MAPEXT WILL CALCULATE THE OTHER
      'wavelength' : ?*u.m,
       'beamwidth' : [?*u.arcmin,?*u.arcmin], # GIVE 2D BEAM PROFILE
      'resolution' : ?*u.arcmin, # RESOLUTION OF MAP AS FWHM OF GAUSSIAN BEAM PROFILE
      'projection' : 'WCS', # 'WCS' or 'HPX<nside>'
            'unit' : 1*u.K, # UNITS OF THE MAP (for CMB TEMP. USE 'KCMB')
           'error' : {'SOURCE 1':1e-3*u.K}, # DICTIONARY OF ERRORS AND SOURCES
           }
}

srcDict = [
            {'coord':SkyCoord(039.22407*u.degree,-00.31469*u.degree,frame='galactic')}
          ]

# Initialise outFile to output results in
outFile = outFileHandler(update=True)#,use_date=False,name='mapext_run_file_2022-03-01')
rundict = {'OUTFILE':outFile}

# initialise source object list to query and add sources to outFile
s = []
for _ in range(len(srcDict)):
    src = astroSrc(**srcDict[_])
    s.append(src)
    outFile.add_source(src)

# Loop over maps
for _ in tqdm(mapDict):
    m = astroMap(**mapDict[_])
    # convert healpix to cartesian if needbe
    if m.PROJ[:3] == 'HPX':
        m.reproject(crval=[25,0.0],crpix=[360,240],cdelt=[1/60,1/60],ctype=['GLON-CYP','GLAT-CYP'],shape_out=[3001,481])
    m.smooth_to(target_res=5*u.arcmin) # smooth to 5 arcmin
    m.convert_unit(units=u.K) # convert units to K
    # m.quickplot(map_kwargs={'cmap':'jet'}) # Quickplot of the survey to check the map looks good
    photometry.aperture(m,s,RUNDICT=rundict) # Aperture photometry on the source list
    map.cutout(m,s,RUNDICT=rundict) # Retrieve cutouts for sources and save them
    sourcefit.single_source(m,s) # Fit gaussian model to sources

    # Update map section of outFile
    rtnParam = m.return_dictionary()
    outFile.update_file({'maps/'+m.ID:rtnParam})
    # Delete map once used
    del m

# update sources in outfile
for _ in s:
    outFile.add_source(_)
    print(_.flux)
