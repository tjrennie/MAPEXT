import numpy as np
import h5py
from mapext.core.usefulFunctions import h5PushDict, ensure_dir

class outFileHandler():

    def __init__(self,name='mapext_run_file',dir='./',use_date=True,update=False):
        if use_date:
            name += '_'+str(np.datetime64('today','D'))
        self.NAME = name+'.hdf5'
        self.DIR  = dir
        if self.DIR[-1] == '/':
            self.DIR = self.DIR[:-1]
        ensure_dir(self.DIR)
        f = h5py.File(self.DIR+'/'+self.NAME,'a')
        if ('run_params' in f) and (update):
                del f['run_params']
                f['run_params/time'] = np.datetime64('now','s').astype('<S99')
        else:
            f['run_params/time'] = np.datetime64('now','s').astype('<S99')
        f.close()

    def update_file(self,changes):
        f = h5py.File(self.DIR+'/'+self.NAME,'a')
        h5PushDict(f,changes,ow=True)
        f.close()

    def add_source(self,source):
        ensure_dir(self.DIR+'/SRCHDF5/')
        # create and populate source file
        sf = h5py.File(self.DIR+'/SRCHDF5/'+source.NAME,'a')
        h5PushDict(sf,source.return_dictionary(),ow=True)
        sf.close()
        # create link to file
        f = h5py.File(self.DIR+'/'+self.NAME,'a')
        if '/sources' not in f:
            f.create_group('sources')
        if source.NAME not in f['/sources/']:
            f['/sources/'+source.NAME] = h5py.ExternalLink('SRCHDF5/'+source.NAME, '/')
        f.close()

    def save_map(self, source, map, wcs, name):
        sf = h5py.File(self.DIR+'/SRCHDF5/'+source.NAME,'a')
        dir = '/maps/'+name
        if dir in sf:
            del sf[dir]
        sf[dir] = map.value
        sf[dir].attrs['cdelt'] = wcs.wcs.cdelt
        sf[dir].attrs['crpix'] = wcs.wcs.crpix
        sf[dir].attrs['crval'] = wcs.wcs.crval
        sf[dir].attrs['ctype'] = np.array(wcs.wcs.ctype,dtype='<S30')
        sf[dir].attrs['unit'] = np.array(str(map.unit),dtype='<S30')
        sf.close()
