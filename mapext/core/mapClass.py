import numpy as np
import astropy.units as u
from astropy.units import Quantity, Unit
from astropy.convolution import convolve, Gaussian2DKernel
from astropy.io import fits
from astropy.wcs import WCS
from astropy.wcs.utils import pixel_to_skycoord
from astropy_healpix import HEALPix
from reproject import reproject_interp,reproject_from_healpix
from mapext import _version

__all__ = ['astroMap', 'mapObj']

def get_astropy_quantity(quantity_string):
    if isinstance(quantity_string, str):
        return Quantity(quantity_string)
    elif isinstance(quantity_string, Quantity):
        return quantity_string
    else:
        raise ValueError('Cannot convert input to astropy Quantity object.')
        
def get_astropy_unit(unit_string):
    if isinstance(unit_string, str):
        return Unit(unit_string)
    elif isinstance(unit_string, Unit):
        return unit_string
    else:
        raise ValueError('Cannot convert input to astropy Unit object.')

coorddict = {'g' : 'galactic'}

class astroMap():
    """Class to hold one or many astronomical maps alongside metadata describing them. Each map held inside an `astroMap` object is contained in a `mapObj` object inside the `astroMap` instance.
    
    Metadata is split into overarching/common metadata (such as origin, reference, frequency) which is held in the `astroMap`, and map-specific metadata (such as Stokes component, projection and resolution) which is held in the `mapObj`.
    
    :param load_maps: Stokes parameters of maps to preload from supplied map details. Defaults to 'IQU'.
    :type load_maps: str, kwarg, optional
    :param Metadata: Dictionary containing general metadata about the maps to be initialised. Must contain `name` and `freq` keys.
    :type Metadata: dict, kwarg    
    :param I: Dictionary containing metadata for the Stokes I map to be initialised. Must contain `fname`, `unit` and `resolution` keys, can contain `header`, `tfield` and `temp` keys. All keys are passed onto a `mapObj` instance without interpretation.
    :type I: dict, kwarg, optional    
    :param Q: Dictionary containing metadata for the Stokes Q map to be initialised. Must contain `fname`, `unit` and `resolution` keys, can contain `header`, `tfield` and `temp` keys. All keys are passed onto a `mapObj` instance without interpretation.
    :type Q: dict, kwarg, optional
    :param U: Dictionary containing metadata for the Stokes U map to be initialised. Must contain `fname`, `unit` and `resolution` keys, can contain `header`, `tfield` and `temp` keys. All keys are passed onto a `mapObj` instance without interpretation.
    :type U: dict, kwarg, optional
    
    :raises ValueError: Key `freq` must be specified in Metadata kwarg dictionary.
    
    :warnings: All functions are defined for cartesian map projections, but may not be for native HEALPix map operation. Please check that your maps are suitable to be processed in cartesian before reprojecting for compatibility.
    """

    def __init__(self, load_maps='IQU', **kwargs):
        """ Initialise astroMap object with parameters given in keyword arguments. Map objects are initialised into mapObj objects inside this class ready for processing.
        """
        # INITIALISE MAP METADATA
        metadata = kwargs['Metadata']
        self.name = metadata.get('name', 'mapobject')
        self.ID = metadata.get('id', self.name)
        self.freq = get_astropy_quantity(metadata.get('freq', None))
        if self.freq == None:
            raise ValueError('Key `freq` must be specified in Metadata kwarg dictionary.')
        self.stokes = list(load_maps)

        # run through IQU maps
        for stokes_parameter in ['I', 'Q', 'U']:
            map_params = kwargs.get(stokes_parameter, None)
            if (map_params is None) or (stokes_parameter not in load_maps):
                setattr(self, stokes_parameter, None)
            else:
                setattr(self, stokes_parameter, mapObj(name=f'{self.ID}_{stokes_parameter}', freq=self.freq, **map_params))
        return
    
    def set_map_projection(self, **kwargs):
        """ Method to set map projection for all available maps in this object. For details see :func:`mapext.core.mapClass.mapObj.reproject`.
        """
        for comp_stokes in self.stokes:
            comp_obj = getattr(self, comp_stokes)
            if comp_obj != None:
                comp_obj.reproject(**kwargs)
        return
    
    def set_map_resolution(self, **kwargs):
        """ Method to set map resolution for all available maps in this object. For details see :func:`mapext.core.mapClass.mapObj.smooth_to`.
        """
        for comp_stokes in self.stokes:
            comp_obj = getattr(self, comp_stokes)
            if comp_obj != None:
                comp_obj.smooth_to(**kwargs)
        return
    
    def set_map_units(self, **kwargs):
        """ Method to set map units for all available maps in this object. For details see :func:`mapext.core.mapClass.mapObj.convert_units`.
        """
        for comp_stokes in self.stokes:
            comp_obj = getattr(self, comp_stokes)
            if comp_obj != None:
                comp_obj.smooth_to(**kwargs)
        return


class mapObj():
    """
    Class to hold a single astronomical map as part of an astroMap object. Holds map-specific metadata that can change between Stokes parameter maps from teh same source. 
    
    :param name: Map name. If initialised as part of `astroMap` object, will be set to `{astromapName}_{stokesParameter}`.
    :type name: string, kwarg
    :param freq: Map central frequency.
    :type freq: string or astropy.units.Quantity, kwarg
    :param fname: Filename for data to be read from.
    :type fname: str, kwarg
    :param resolution: Resolution of the map as read in.
    :type resolution: str, kwarg
    :param beamwidth: Beamwidth of original telescope beam before any postprocessing or smoothing. Defaults to resolution value.
    :type beamwidth: string or astropy.units.Quantity, kwarg, optional
    :param cal: Calibration uncertainty as percentage (e.g. 5% should be inputted at 5). Defaults to 8%.
    :type cal: float, kwarg, optional
    :param unit: Units of teh input map. If units are Kelvin (K), use temp keyword to specify temperature type.
    :type unit: string or astropy.units.Unit, kwarg, optional
    :param temp: Type of temperature specified in units: `cmb` for CMB temperature (as defined by Planck) or `rj` for Rayleigh-Jeans temperature.
    :type temp: str, kwarg, optional
    :param zaxis: index of map to read in if specified in an image cube 9greater than 2 dimensions).
    :type zaxis: int, kwarg, optional
    """

    def __init__(self, name=None, freq=None, **kwargs):
        """Initialises mapObj object.
        """
        self.name = name
        self.freq = get_astropy_quantity(freq)
        self.reso = get_astropy_quantity(kwargs.get('resolution', None))
        self.bwid = get_astropy_quantity(kwargs.get('beamwidth', str(self.reso)))
        self.unit = get_astropy_unit(kwargs.get('unit', None))
        if self.unit in [u.K, u.mK, u.uK, u.nK]:
            self.temp = kwargs.get('temp', 'rj')
        else:
            self.temp = None
        self.cal = kwargs.get('calibration', 8)
        self.tfield = kwargs.get('tfield', None)
        self.zaxis = kwargs.get('zaxis', None)
        self.fname = str(kwargs.get('fname', None))
        filetype = self.determine_map_type()
        if filetype == 'wcs':
            self.hdrno = int(kwargs.get('header', 0))
            self.data, self.proj = self.read_wcs_map()
        elif filetype == 'healpix':
            self.hdrno = int(kwargs.get('header', 1))
            self.data, self.proj = self.read_hpx_map()
    
    def determine_map_type(self):
        """ Method to determine if a fits file is a HEALPix file or a WCS (cartesian) file

        :return: Pixelisation scheme, `wcs` if a cartesian WCS-defined object or `healpix` if HEALPix pixelisation.
        :rtype: string
        """
        try:
            with fits.open(self.fname) as hdul:
                header = hdul[1].header
                test = header["ORDERING"]
                return 'healpix'
        except Exception:
            return 'wcs'
        
    def read_wcs_map(self):
        """ Method to read a fits file containing a map with WCS (cartesian) projection to memory

        :return: Map data, Map projection.
        :rtype: numpy.array, astropy.wcs.WCSObject
        """
        hdu = fits.open(self.fname)[self.hdrno]
        proj = WCS(hdu.header)
        if self.tfield != None:
            data = np.array(hdu.data[self.tfield])
        else:
            data = np.array(hdu.data)
        if self.zaxis != None:
            for _ in self.zaxis:
                data = data[_]
            proj = proj.celestial
        return data, proj
        
    def read_hpx_map(self):
        """ Method to read a fits file containing a map with HEALPix projection to memory

        :return: Map data, Map projection.
        :rtype: numpy.array, astropy_healpix.HEALPix
        """
        hdu = fits.open(self.fname)[self.hdrno]
        coordsys = None
        if "COORDSYS" in hdu.header.keys():
            if hdu.header["COORDSYS"].lower() in coorddict:
                coordsys = coorddict[hdu.header["COORDSYS"].lower()]
            else:
                coordsys = hdu.header["COORDSYS"].lower()
    
        else:
            coordsys = 'galactic'
        proj = HEALPix(nside=hdu.header['NSIDE'], order=hdu.header['ORDERING'], frame=coordsys)
        if self.tfield != None:
            data = np.array(hdu.data[self.tfield])
        else:
            data = np.array(hdu.data)
        return data, proj
    
    def reproject(self, wcs=None, shape_out=None, order='nearest_neighbour'):
        """ Method to reproject current map data into another projection specified by the keyword arguments supplied

        :param wcs: WCS object describing final projection. Defaults to None.
        :type wcs: astropy.wcs.WCSobject, optional
        :param shape_out: Desired final map shape supplied as a length-2 list. Defaults to None.
        :type shape_out: numpy.ndarray, optional)
        :param order: Interpolation scheme to use. Defaults to 'nearest_neighbour'.
        :type order: str, optional

        :todo: Add support for HEALPix -> HEALPix transformations
        :todo: Add support for WCS -> HEALPix transformations
        """
        if type(self.proj) == WCS:
            new_map, new_footprint = reproject_interp((self.data, self.proj), wcs, shape_out=shape_out)
            new_map[new_footprint==0] = np.nan
            self.data, self.proj = new_map, wcs
        elif type(self.proj) == HEALPix:
            new_map, new_footprint = reproject_from_healpix((self.data, str(self.proj.frame.name).lower()), wcs, shape_out=shape_out, nested=True)
            # new_map[new_footprint==0] = np.nan
            self.data, self.proj = new_map, wcs
        return
    
    def smooth_to(self, reso=None):
        """ Method to smooth map to desired resolution based on arguments supplied

        Args:
            reso (astropy quantity, optional): Resolution of the final map. Need not be supplied if fwhm given, but priorotised if both given. Defaults to None.

        :param reso: Resolution of the final map. Need not be supplied if fwhm given, but priorotised if both given. Defaults to None.
        :type reso: astropy.units.Quantity
        :param fwhm: Full-width half maximum of a point source is the (desired) final map.
        :type fwhm: astropy.units.Quantity

        :raises ValueError: Map resolution lower than target resolution ({reso} > {self.reso}). Please re-evaluate use of this map in this run.
        :raises ValueError: Operation not yet supported.
        
        :todo: Add support for smoothing HEALPix maps
        """
        new_reso = reso
        if new_reso < self.reso:
            raise ValueError(f"Map resolution lower than target resolution ({reso} > {self.reso}). Please re-evaluate use of this map in this run.")
        if new_reso == self.reso:
            return
        else:
            if type(self.proj) == WCS:
                dr = np.sqrt(new_reso**2 - self.reso**2).to(u.degree).value/2.355
                kernel = Gaussian2DKernel(
                    x_stddev=dr/np.abs(self.proj.wcs.cdelt[0]),
                    y_stddev=dr/np.abs(self.proj.wcs.cdelt[1]))
                new_map = convolve(self.data, kernel)
                self.data, self.reso = new_map, new_reso
            elif type(self.proj) == HEALPix:
                raise ValueError("Operation not yet supported")
        return
    
    def convert_units(self, new_units=u.Jy/(u.pixel**2)):
        """ Method to convert map to to desired units based on arguments supplied

        :param new_units: (_type_, optional): _description_. Defaults to u.Jy/(u.pixel**2).

        :raises ValueError: Unit not recognised as flux density units
        """
        if new_units != self.unit:
            # first convert to Jy.sr
            if self.unit in [u.W/(u.m**2 * u.sr)]:
                map_int = (self.data*self.unit).to(u.Jy/u.sr, equivalencies=u.spectral_density(self.freq))[:]
            elif (self.unit in [u.K, u.mK, u.uK, u.nK]) and (self.temp == 'rj'):
                map_int = (self.data*self.unit).to(u.Jy/u.sr, equivalencies=u.brightness_temperature(self.freq))[:]
            elif (self.unit in [u.K, u.mK, u.uK, u.nK]) and (self.temp == 'cmb'):
                from astropy.cosmology import Planck15
                map_int = (self.data*self.unit).to(
                    u.Jy/u.sr, equivalencies=u.thermodynamic_temperature(self.freq, Planck15.Tcmb0))
            elif (self.unit in [u.Jy/u.sr, get_astropy_unit('MJy/sr')]):
                map_int = (self.data * self.unit).to(u.Jy/u.sr)[:]
            elif self.unit in [u.Jy/u.beam, u.mJy/u.beam]:
                map_int = (self.data*self.unit).to(u.Jy/u.sr, equivalencies=u.beam_angular_area(self.get_beam_area()))[:]
            elif self.unit in [u.Jy/(u.pixel**2)]:
                pixel = abs(self.proj.wcs.cdelt[0]*self.proj.wcs.cdelt[1]) * (u.degree/u.pixel)**2
                map_int = ((self.data * self.unit) / pixel).to(u.Jy/u.sr)
            else:
                raise ValueError(f'Unit {self.unit} not recognised as flux density units')
            # convert to final units
            # first convert to Jy.sr
            if new_units in [u.W/(u.m**2 * u.sr)]:
                map_final = map_int.to(new_units, equivalencies=u.spectral_density(self.freq))[:]
            elif new_units in [u.K]:
                map_final = map_int.to(new_units, equivalencies=u.brightness_temperature(self.freq))[:]
            elif new_units in [u.Jy/u.sr]:
                map_final = map_int[:]
            elif new_units in [u.Jy/u.beam]:
                map_final = map_int.to(new_units, equivalencies=u.beam_angular_area(self.get_beam_area()))[:]
            elif new_units in [u.Jy/(u.pixel**2)]:
                pixel = abs(self.proj.wcs.cdelt[0]*self.proj.wcs.cdelt[1]) * (u.degree/u.pixel)**2
                map_final = (map_int * pixel).to(u.Jy/(u.pixel**2))
            else:
                raise ValueError(f'Unit {new_units} not recognised as flux density units')
            self.data, self.unit = map_final.value, map_final.unit
        return

    def get_pixel_coords(self):
        """Method to return an array of center coordinates for each pixel

        :return: Array of SkyCoord objects representing sky coordinates of the center of each map pixel.
        :rtype: numpy.ndarray of astropy.coords.SkyCoord
        """
        Xp,Yp = np.meshgrid(np.arange(self.data.shape[1]),
                            np.arange(self.data.shape[0]))
        sc = pixel_to_skycoord(Xp, Yp, self.proj)
        return sc
    
    def get_beam_area(self):
        """Method to return beam area

        :return: Beam area.
        :rtype: astropy.unit.Quantity
        """
        return np.pi * (self.bwid/2.355)**2
    
    def meta_to_hdr(self):
        unit = ''
        if self.unit == u.K:
            unit = f'{self.unit}_{self.temp.upper()}'
        else:
            unit = f'{self.unit}'
        include = [
            f'Resolution : {self.reso}',
            f'Beamwidth : {self.bwid}',
            f'Frequency : {self.freq}',
            f'Units : {unit}',
            f'calUncertainty : {self.uncert_cal}']
        return include
    
    def save_out(self, filename=None, dir = '.'):
        """Method to save out the current state of a map.
        """
        hdu = fits.PrimaryHDU(self.data, header=self.proj.to_header())
        hdul = fits.HDUList(hdu)
        
        hdul[0].header['COMMENT'] = f'Processed with MAPEXT {_version.__version__}'
        hdul[0].header['COMMENT'] = '-'*80
        for _ in self.meta_to_hdr():
            hdul[0].header['COMMENT'] = _
            
        if filename is not None:
            outfilename = filename
        else:
            outfilename = '{}-{}.fits'.format(self.name, self.reso)
        try:
            hdul.writeto(f'{dir}/{outfilename}', overwrite=True)
        except:
            print('File already exists - not overwritten')
        del hdul
        return
        