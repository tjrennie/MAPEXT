import numpy as np
import astropy.units as u
from astropy.modeling.fitting import LevMarLSQFitter

from mapext.core import astroSrc

__all__ = ['sedFitter_LSQ']

class sedFitter_LSQ():
    """ Least Squares SED Fitter class.
    
    :param model: SED model to fit
    :type model: mapext.FittableEmissionModel or astropy.CompositeModel
    :param fluxes: Source flux measurements to fit to. Can be entered as a mapext.astroSrc object, the flux numpy.ndarray from a source, or as three lists for frequency [Hz] flux density and its uncertainty [Jy]
    :type fluxes: mapext.astroSrc object, numpy.ndarray or 3 arrays/lists
    :param beam_area: Beam area in steradians
    :type beam_area: float
    """
    
    def __init__(self, *args, **kwargs):
        # Check model
        self.model = args[0]
        # Check beam area
        self.beam_area = args[-1]
        if type(self.beam_area) == u.Quantity:
            self.beam_area = self.beam_area.to(u.sr).value
        # Check flux parameters
        self.astroSrc_obj = None
        nu, flux, fluxerr = None, None, None
        if len(args) == 3:
            if type(args[1]) is astroSrc:
                print('AstroSrc')
                self.astroSrc_obj = args[1]
                stokes = kwargs.get('fit_stokes', 'I')
                mask = np.array([_ in stokes for _ in args[1].flux['stokes']])
                nu, flux, fluxerr = args[1].flux['freq'][mask], args[1].flux['Sv'][mask], args[1].flux['Sv_e'][mask]
            elif type(args[1]) is np.ndarray:
                stokes = kwargs.get('fit_stokes', 'I')
                mask = np.array([_ in stokes for _ in args[1]['stokes']])
                nu, flux, fluxerr = args[1]['freq'][mask], args[1]['Sv'][mask], args[1]['Sv_e'][mask]
            else:
                ValueError('Input not recognised')
        elif len(args) == 5: 
            nu, flux, fluxerr = args[1:-1]
        else:
            ValueError('Input not recognised')
        self.flux = np.array([tuple(_) for _ in zip(nu/1e9, flux, fluxerr)], dtype = [('freq',float),('Sv',float),('Sv_e',float)])
        # Check input
        mask_nonfinite = np.all([
            np.isfinite(self.flux['freq']),
            np.isfinite(self.flux['Sv']),
            np.isfinite(self.flux['Sv_e'])],
            axis=0)
        self.flux = self.flux[mask_nonfinite]
        
    def write_model_to_astroSrc(self):
        """ Write the model to the astroSrc object.
        """
        self.astroSrc_obj.append_flux_model([('', self.fit_info, 'sedFitterLSQ')])
        return None
        
    def calc_jacobian(self, model, nu):
        """ Calculate the Jacobian matrix of derivitives for a compund model.

        :params model: Model to evalue the jacobian for
        :type model: astropy.modelling.CompositeModel
        :params nu: Frequency [Hz]
        :type nu: numpy.ndarray
        
        :return: Jacobian matrix of the model
        :rtype: np.ndarray
        """
        J = np.zeros([len(list(nu)), len(model.parameters)])
        idx = 0
        for _ in model:
            n_param = len(_.parameters)
            J[:,idx:idx+n_param] = np.array(_.deriv(nu/1e9, self.beam_area, *_.parameters)).T
            idx += n_param
        return J
    
    def fitSED(self):
        """ Run fitter.

        :return: Fit information
        :rtype: Same as model input
        """
        mainFit = LevMarLSQFitter(calc_uncertainties=True)
        self.fit_info = mainFit(self.model,
                                self.flux['freq'],
                                self.beam_area*np.ones(self.flux['freq'].shape),
                                self.flux['Sv'],
                                weights = 1/self.flux['Sv_e'])
        self.fit_info.fit_jacobian = lambda nu:self.calc_jacobian(self.fit_info,nu)
        
        if self.astroSrc_obj != None:
            self.write_model_to_astroSrc()
            
        return self.fit_info