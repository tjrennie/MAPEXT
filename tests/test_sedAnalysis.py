import unittest
import mapext
from mapext.emission import models
import toml
import numpy as np

fname_map_dict = '/Users/tjrennie/MAPEXT/mapext/tests/mapTests_config.toml'
real_Sv = np.array([21.39, 19.44, 17.96, 17.52, 22.71, 22.58, 22.32, 21.99, 24.74, 55.08, 214.01, 572.24, 1686.33, 4537.50, 12205.65, 12311.73, 9005.56, 7230.43])

class Test_TestSEDAnalysis(unittest.TestCase):

    def test_aperture_photometry(self):
        list_maps = toml.load(fname_map_dict)
        mapobjs = []
        src = mapext.core.astroSrc(coords=[0,0])
        for _ in list_maps.keys():
            mapObj = mapext.core.astroMap(load_maps='I', **list_maps[_])
            mapobjs.append(mapObj)
            mapext.photometry.aperture_photometry([mapObj], [src], verbose=True)
        deltas = src.flux['Sv'] - real_Sv
        deltas_scaled = deltas/src.flux['Sv_e']
        self.assertTrue(np.all(np.abs(deltas_scaled)<1,axis=0))

    def test_gaussian_fitting(self):
        list_maps = toml.load(fname_map_dict)
        mapobjs = []
        src = mapext.core.astroSrc(coords=[0,0])
        for _ in list_maps.keys():
            mapObj = mapext.core.astroMap(load_maps='I', **list_maps[_])
            mapobjs.append(mapObj)
            mapext.photometry.gaussian_fitter([mapObj], [src])
        deltas = src.flux['Sv'] - real_Sv
        deltas_scaled = deltas/src.flux['Sv_e']
        self.assertTrue(np.all(np.abs(deltas_scaled)<1,axis=0))
        
    def test_sedfit(self):
        list_maps = toml.load(fname_map_dict)
        mapobjs = []
        src = mapext.core.astroSrc(coords=[0,0])
        for _ in list_maps.keys():
            mapObj = mapext.core.astroMap(load_maps='I', **list_maps[_])
            mapobjs.append(mapObj)
            mapext.photometry.gaussian_fitter([mapObj], [src])
        model = models.freeFree_7000k(10) + models.ame_lognormal(10, 30, 0.5) + models.thermalDust(40, 1, 1)
        fitter = mapext.emission.sedFitter_LSQ(model, src, np.pi*np.radians(2/60)**2)
        fitter.fitSED_quickFit()
        
        J = fitter.fit_info.fit_jacobian(src.flux['freq']/1e9)
        C = fitter.fit_info.cov_matrix.cov_matrix
        S_cov = np.dot(J, np.dot(C, J.T))
        S_std = np.sqrt(np.diag(S_cov))

        S = fitter.fit_info(src.flux['freq']/1e9, fitter.beam_area)
        
        deltas = src.flux['Sv'] - S
        sigmas = np.sqrt(src.flux['Sv_e']**2 + S_std**2)
        deltas_scaled = deltas/sigmas
        
        print(deltas_scaled)
        
        self.assertTrue(np.all(np.abs(deltas_scaled)<1,axis=0))
        
if __name__ == '__main__':
    unittest.main()