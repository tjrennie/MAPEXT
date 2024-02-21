import numpy as np
import h5py as h5
from datetime import datetime as dt

import mapext.settings as settings

__all__ = ['mapext_runfile', 'initialise_output']

def initialise_output(*args):
    # If method is passed a path, initialise dir as 
    if len(args)==1:
        settings.saving.path = args[0]
        print(settings.saving.path)
    print(settings.saving.path)
    # Check a path exists
    if settings.saving.path == None:
        raise ValueError("A path is required to initialise the output file.\nPlease specify this either in the settings file, or by setting 'mapext.settings.saving.path' directly.")
    # Initialise file
    runfile = mapext_runfile()
    # Set outfile in settings
    settings.saving.outfile = runfile
    return True

class mapext_runfile():
    """
    Class to manage exportinga and importing astroMap and astroSrc objects to and from one (or many) hdf5 file(s).
    
    :param path: Path for runfile to be saved to. If None, reverts to setting given in mapext.settings. Defaults to None.
    :type path: string
    :param naming: Naming scheme for runfile and ancillary files. If None, reverts to setting given in mapext.settings. Defaults to None.
    :type naming: string
    """

    def __init__(self,
                 path=None,
                 naming=None):
        """ Initialisation method for mapext_runfile.
        """
        if path == None:
            self.path = settings.saving.path
        else:
            self.path = path
        
        if naming == None:
            self.naming_scheme = settings.saving.default_naming_scheme
        else:
            self.naming_scheme = naming

        if self.path[-1] != '/':
            self.path += '/'
        self.mainname = self.get_fname('main')
        # Create root file
        rootfile = h5.File(self.path+self.mainname+'.hdf5','a')
        rundets = rootfile.create_group("RunDetails")
        rundets.attrs['DateTime'] = dt.now().strftime('%Y%m%dT%H:%M:%S')
        rundets.attrs['Author'] = settings.saving.author
        rundets.attrs['Description'] = settings.saving.description
        rootfile.close()
        return 

    def get_fname(self, ftype, **kwargs):
        """ Method to retrieve names for runfile and ancillary data files

        :param ftype: Filetype to be named
        :type ftype: string
        
        :raise ValueError: Do not recognise naming scheme
        
        :return: Filename
        :rtype: string
        """
        name_format = settings.saving.naming_schemes[self.naming_scheme][ftype]
        fname = ''
        for _ in name_format:
            if isinstance(_, list):
                # instruction set
                if _[0] == 'DateTime':
                    fname += dt.now().strftime(_[1])
                else:
                    ValueError(f'Do not recognise naming scheme "{_[0]}"')
            else:
                fname += str(_)
        return fname