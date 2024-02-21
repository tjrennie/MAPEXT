import unittest
import mapext
import toml
import numpy as np
from astropy import units as u

fname_map_dict = 'mapext/tests/mapTests_config.toml'
list_maps = toml.load(fname_map_dict)
map_param = list_maps['Test353']

class Test_TestSEDAnalysis(unittest.TestCase):

    def test_smoothing(self):
        resos = [6,7,8,9,10,11,12,13,14]
        src = mapext.core.astroSrc(coords=[0,0])
        mapObj = mapext.core.astroMap(load_maps='I', **map_param)
        mapext.photometry.gaussian_fitter([mapObj], [src])
        init_flux = src.flux[0]
        for r in resos:
            mapObj.set_map_resolution(reso = r*u.arcmin)
            mapext.photometry.gaussian_fitter([mapObj], [src])
            current_flux = src.flux[-1]
            delta_s = np.abs(current_flux['Sv'] - init_flux['Sv'])
            combined_uncert = np.sqrt(current_flux['Sv_e']**2 + init_flux['Sv_e']**2)
            reduced_diff = delta_s/combined_uncert
            self.assertLessEqual(reduced_diff, 1)

if __name__ == '__main__':
    unittest.main()