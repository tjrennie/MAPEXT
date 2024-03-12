import numpy as np
from astropy.coordinates import SkyCoord

__all__ = ['astroSrc']

class astroSrc(baseSrc):
    """Class to hold information pertaining to a specific astronomical source. this class is currently aimed at point sources, although will be expanded to be more flexible.
    
    :param name: Source name
    :type name: string, optional
    :param coords: Center coordinates for the source.
    :type coords: list of astropy.coords.SkyCoord or list of floats
    :param frame: Coordinate frame if list of floats supplied for coord.
    :type frame: string, optional
    :param flux: Array of flux values (represented by tuples). Each tuple should be of the form (<name, str>, <stokes_parameter, str>, <frequency, Hz>, <resolution, arcmin>, <flux, Jy>, <std_flux, Jy>, <origin/reference, str>).
    :type flux: list of tuples
    :param model: Array of source models (represented by tuples). Each tuple should be of the form (<name, str>, <stokes_parameter, str>, <frequency, Hz>, <resolution, arcmin>, <model, obj>, <std_model, obj>, <origin/reference, str>).
    :type model: list of tuples
    :param flux_model: Array of source models (represented by tuples). Each tuple should be of the form (<name, str>, <model, obj>, <origin/reference, str>).
    :type flux_model: list of tuples
    """

    def __init__(self, **kwargs):
        """Initialisation function.
        """
        name = kwargs.get('coords', None)
        coords = kwargs.get('coords', None)
        frame = kwargs.get('frame', 'galactic')
        if isinstance(coords, SkyCoord) == False:
            self.coord = SkyCoord(*coords, frame=frame, unit='degree')
        self.flux = np.array([],
            dtype=[('name','<U40'),('stokes','<U40'),('freq',float),('resolution',float),('Sv',float),('Sv_e',float),('origin','<U40')])
        self.model = np.array([],
            dtype=[('name','<U40'),('stokes','<U40'),('freq',float),('resolution',float),('model', object),('model_e', object),('origin','<U40')])
        self.flux_model = np.array([],
            dtype=[('name','<U40'),('model', object),('origin','<U40')])
        return
    
    def append_flux_measure(self, fluxentry):
        """ Function to append new flux measure entry to table.
        
        :param fluxentry: Array of flux values (represented by tuples). Each tuple should be of the form (<name, str>, <stokes_parameter, str>, <frequency, Hz>, <resolution, arcmin>, <flux, Jy>, <std_flux, Jy>, <origin/reference, str>).
        :type fluxentry: list of tuples
        """
        new_flux = np.array(fluxentry, dtype=[('name','<U40'),('stokes','<U40'),('freq',float),('resolution',float),('Sv',float),('Sv_e',float),('origin','<U40')])
        self.flux = np.concatenate([self.flux,new_flux],axis=0)
        return
    
    def append_model(self, modelentry):
        """ Function to append new source model to table.
        
        :param modelentry: Array of source models (represented by tuples). Each tuple should be of the form (<name, str>, <stokes_parameter, str>, <frequency, Hz>, <resolution, arcmin>, <model, obj>, <std_model, obj>, <origin/reference, str>).
        :type modelentry: list of tuples"""
        new_model = np.array(modelentry, dtype=[('name','<U40'),('stokes','<U40'),('freq',float),('resolution',float),('model', object),('model_e', object),('origin','<U40')])
        self.model = np.concatenate([self.model,new_model],axis=0)
        return
    
    def append_flux_model(self, fluxmodelentry):
        """ Function to append new flux model to table.

        :param flux_model: Array of source models (represented by tuples). Each tuple should be of the form (<name, str>, <model, obj>, <origin/reference, str>).
        :type flux_model: list of tuples
        """
        new_model = np.array(fluxmodelentry, dtype=[('name','<U40'),('model', object),('origin','<U40')])
        self.flux_model = np.concatenate([self.flux_model,new_model],axis=0)
        return