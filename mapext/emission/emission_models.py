import numpy as np
from astropy.modeling import FittableModel, Parameter

__all__ = [
    'synchrotron_1comp',
    'freeFree','freeFree_7000k',
    'ame_lognormal',
    'thermalDust',
]

phys_const = {
    'c' : 299792458.,
    'k' : 1.3806488e-23,
    'h' : 6.62606957e-34,
}

class FittableEmissionModel(FittableModel):
    n_inputs = 2
    n_outputs = 1
    
    @staticmethod
    def fit_deriv(*args):
        return FittableEmissionModel.deriv(self, *args)

# Synchrotron emission model
class synchrotron_1comp(FittableEmissionModel):
    """Emission model for 1-component synchrotron emission without spectral break or curvature (power law).
    
    :param synch_S1: Synchrotron flux at 1GHz
    :type name: float
    :param synch_alp: Synchrotron spectral index
    :type name: float
    
    :Formulism:
        .. math::
            S^\\mathrm{\\,sync}_\\nu = S1 * \\nu^\\alpha
        For details see `Carol and Oslie (2007)`_ .

    :note: Either all or none of input ``nu``, ``area``, ``synch_S1`` and ``synch_alp`` must be provided consistently with compatible units or as unitless numbers.
    
    .. _Carol and Oslie (2007): https://ui.adsabs.harvard.edu/abs/1996ima..book.....C/abstract
    """
    synch_S1 = Parameter(default=1,
                         min=0,
                         description="Synchrotron flux at 1 GHz")
    synch_alp = Parameter(default=-0.7,
                          description='Synchrotron spectral index')

    @staticmethod
    def evaluate(nu, area, synch_S1, synch_alp):
        """Evaluate the emission model.

        :param nu: Frequency in GHz
        :type nu: float or numpy.ndarray
        :param area: Beam area in steradians
        :type name: float or numpy.ndarray
        :param synch_S1: Synchrotron flux at 1GHz
        :type name: float or numpy.ndarray
        :param synch_alp: Synchrotron spectral index
        :type name: float or numpy.ndarray

        :return: Evaluated function
        :rtype: float or numpy.ndarray
        """
        return synch_S1 * (nu**synch_alp)
    
    def deriv(self,nu, area, synch_S1, synch_alp):
        """Evaluate the first derivitives of emission model with respect to input parameters.

        :param nu: Frequency in GHz
        :type nu: float or numpy.ndarray
        :param area: Beam area in steradians
        :type name: float or numpy.ndarray
        :param synch_S1: Synchrotron flux at 1GHz
        :type name: float or numpy.ndarray
        :param synch_alp: Synchrotron spectral index
        :type name: float or numpy.ndarray

        :return: First derivitives with repect to input parameters in order.
        :rtype: list
        """
        d_synch_S1 = nu**synch_alp
        d_synch_alp = synch_S1 * (nu**synch_alp) * np.log(nu)
        return [d_synch_S1, d_synch_alp]

# Free-free emission model (Draine 2011)
class freeFree_7000k(FittableEmissionModel):
    """Emission model for free-free emission with an electron temperature of 7000K.
    
    :param ff_em: Free-free emission measure
    :type ff_em: float

    :Formulism:
        .. math::
            \\tau^\\mathrm{ff}_\\nu = 5.468\\times 10^{-2} \\cdot T_e^{-\\frac{3}{2}} \\left[ \\frac{\\nu}{\\mathrm{GHz}} \\right]^{-2} \\left[\\frac{EM}{\\mathrm{pc\\,cm}^-6}\\right]  g^\\mathrm{ff}_\\nu
        .. math::
            g^\\mathrm{ff}_\\nu = \\ln\\left(\\exp\\left\\{5.90 - \\frac{\\sqrt{3}}{\\pi}\\ln\\left(\\left[ \\frac{\\nu}{\\mathrm{GHz}} \\right] \\left[\\frac{T_e}{10^4\\,\\mathrm{K}}\\right] ^\\frac{3}{2}\\right)\\right\\} + 2.71828\\right)
        .. math::
            T^\\mathrm{ff}_\\nu = T_e \\left(1-e^{-\\tau^\\mathrm{ff}_\\nu}\\right)
        .. math::
            S^\\mathrm{ff}_\\nu = \\frac{2k_B\\Omega\\nu^2}{c^2} T^\\mathrm{ff}_\\nu
        For details see `Draine (2011)`_ .
    
    :note: Either all or none of input ``nu``, ``area`` and ``ff_em`` must be provided consistently with compatible units or as unitless numbers.
    
    .. _Draine (2011): https://ui.adsabs.harvard.edu/abs/2011piim.book.....D/abstract
    """
    ff_em = Parameter(default=100,
                      min=0,
                      description="Free-free emission measure")

    @staticmethod
    def evaluate(nu, area, ff_em):
        """Evaluate the emission model.

        :param nu: Frequency in GHz
        :type nu: float or numpy.ndarray
        :param area: Beam area in steradians
        :type name: float or numpy.ndarray
        :param ff_em: Free-free emission measure
        :type ff_em: float

        :return: Evaluated function
        :rtype: float or numpy.ndarray
        """
        T_e = 7000
        a = 0.366 * np.power(nu,0.1)* np.power(T_e,-0.15) * (np.log(np.divide(4.995e-2, nu)) + 1.5 * np.log(T_e))
        T_ff = 8.235e-2 * a * np.power(T_e,-0.35) * np.power(nu,-2.1) * (1. + 0.08) * ff_em
        S = 2. * phys_const['k'] * area * np.power(np.multiply(nu,1e9),2)  / phys_const['c']**2 * T_ff * 1e26
        return S
    
    def deriv(self,nu, area, ff_em):
        """Evaluate the first derivitives of emission model with respect to input parameters.

        :param nu: Frequency in GHz
        :type nu: float or numpy.ndarray
        :param area: Beam area in steradians
        :type name: float or numpy.ndarray
        :param ff_em: Free-free emission measure
        :type ff_em: float

        :return: First derivitives with repect to input parameters in order.
        :rtype: list
        """
        T_e = 7000
        a = 0.366 * np.power(nu, 0.1) * np.power(T_e, -0.15) * (np.log(np.divide(4.995e-2, nu)) + 1.5 * np.log(T_e))
        T_ff_derivative_ff_em = 8.235e-2 * a * np.power(T_e, -0.35) * np.power(nu, -2.1)
        d_ff_em = 2.0 * phys_const['k'] * area * np.power(np.multiply(nu, 1e9) / phys_const['c'], 2) * T_ff_derivative_ff_em * 1e26
        return [d_ff_em]

# Free-free emission model (Draine 2011)
class freeFree(FittableEmissionModel):
    """Emission model for free-free UCHii emission 
    
    :param ff_em: Free-free emission measure
    :type ff_em: float
    :param ff_Te: Free-free electron temperature
    :type ff_Te: float

    :Formulism:
        .. math::
            g^\\mathrm{ff}_\\nu = \\ln\\left(\\exp\\left\\{5.90 - \\frac{\\sqrt{3}}{\\pi}\\ln\\left(\\left[ \\frac{\\nu}{\\mathrm{GHz}} \\right] \\left[\\frac{T_e}{10^4\\,\\mathrm{K}}\\right] ^\\frac{3}{2}\\right)\\right\\} + 2.71828\\right)
        .. math::
            \\tau^\\mathrm{ff}_\\nu = 5.468\\times 10^{-2} \\cdot T_e^{-\\frac{3}{2}} \\left[ \\frac{\\nu}{\\mathrm{GHz}} \\right]^{-2} \\left[\\frac{EM}{\\mathrm{pc\\,cm}^-6}\\right]  g^\\mathrm{ff}_\\nu
        .. math::
            T^\\mathrm{ff}_\\nu = T_e \\left(1-e^{-\\tau^\\mathrm{ff}_\\nu}\\right)
        .. math::
            S^\\mathrm{ff}_\\nu = \\frac{2k_B\\Omega\\nu^2}{c^2} T^\\mathrm{ff}_\\nu
        For details see `Draine (2011)`_ .
    
    :note: Either all or none of input ``nu``, ``area`` and ``ff_em`` must be provided consistently with compatible units or as unitless numbers.
    
    .. _Draine (2011): https://ui.adsabs.harvard.edu/abs/2011piim.book.....D/abstract
    """
    ff_em = Parameter(default=100,
                      min=0,
                      description="Free-free emission measure")
    ff_Te = Parameter(default=7000,
                      min=0,
                      description="Free-free electron temperature")

    @staticmethod
    def evaluate(nu, area, ff_em, ff_Te):
        """Evaluate the emission model.

        :param nu: Frequency in GHz
        :type nu: float or numpy.ndarray
        :param area: Beam area in steradians
        :type name: float or numpy.ndarray
        :param ff_em: Free-free emission measure
        :type ff_em: float
        :param ff_Te: Free-free electron temperature
        :type ff_Te: float

        :return: Evaluated function
        :rtype: float or numpy.ndarray
        """
        g = np.log( np.exp( 5.90 - ( np.sqrt(3)/np.pi * np.log( nu * ((ff_Te/1e4)**1.5) ) ) ) + 2.71828 )
        tau = 5.468e-2 * (ff_Te**-1.5) * (nu**-2) * ff_em * g
        T_ff = ff_Te * (1 - np.exp(-1 * tau))
        S = 2. * phys_const['k'] * area * np.power(np.multiply(nu,1e9),2)  / phys_const['c']**2 * T_ff * 1e26
        return S
    
    def deriv(self, nu, area, ff_em, ff_Te):
        """Evaluate the first derivitives of emission model with respect to input parameters.

        :param nu: Frequency in GHz
        :type nu: float or numpy.ndarray
        :param area: Beam area in steradians
        :type name: float or numpy.ndarray
        :param ff_em: Free-free emission measure
        :type ff_em: float
        :param ff_Te: Free-free electron temperature
        :type ff_Te: float

        :return: Evaluated function
        :rtype: float or numpy.ndarray
        """
        g = np.log( np.exp( 5.90 - ( np.sqrt(3)/np.pi * np.log( nu * ((ff_Te/1e4)**1.5) ) ) ) + 2.71828 )
        tau = 5.468e-2 * (ff_Te**-1.5) * (nu**-2) * ff_em * g
        T_ff = ff_Te * (1 - np.exp(-1 * tau))
        S = 2. * phys_const['k'] * area * np.power(np.multiply(nu,1e9),2)  / phys_const['c']**2 * T_ff * 1e26
        
        dtau_dem = tau/ff_em
        dTff_dem = T_ff*dtau_dem/(np.exp(tau) - 1)
        dS_dem = S*dTff_dem/T_ff
        
        dg_dTe = 225693 / ((ff_Te * ((nu * (ff_Te**1.5))**(np.sqrt(3)/np.pi))) + (272908*ff_Te))
        dtau_dTe = (tau*dg_dTe/g) - (tau*1.5/ff_Te)
        dTff_dTe = (T_ff*dtau_dTe/(np.exp(tau) - 1)) + T_ff/ff_Te
        dS_dTe = S*dTff_dTe/T_ff
        
        return [dS_dem, dS_dTe]


# AME lognormal
class ame_lognormal(FittableEmissionModel):
    """Emission model for an AME lognormal source.
    
    :param ame_ampl: AME peak flux density
    :type ame_ampl: float
    :param ame_peak: AME peak frequency
    :type ame_peak: float
    :param ame_width: AME lognormal width
    :type ame_width: float

    :Formulaism:
        .. math::
            S^{\\mathrm{AME}}_\\nu = A_\\mathrm{AME} \\cdot \\exp\\left\\{ -\\frac{1}{2}\\left( \\frac{\\ln(\\nu/\\nu_\\mathrm{AME})}{W_\\mathrm{AME}} \\right)^2  \\right\\}
        For details see `Stevenson (2014)`_ .

    :notes: Either all or none of input ``nu``, ``area``, ``ame_ampl``, ``ame_peak`` and ``ame_width`` must be provided consistently with compatible units or as unitless numbers.
        
    .. _Stevenson (2014): https://ui.adsabs.harvard.edu/abs/2014ApJ...781..113S/abstract
    """
    ame_ampl = Parameter(default=10,
                         min=0,
                         description="AME peak flux density")
    ame_peak = Parameter(default=27,
                         min=5, max=60,
                         description='AME peak frequency')
    ame_width = Parameter(default=0.5,
                          min=0, max=1,
                          description="AME lognormal width")
    
    @staticmethod
    def evaluate(nu, area, ame_ampl, ame_peak, ame_width):
        """Evaluate the emission model.

        :param nu: Frequency in GHz
        :type nu: float or numpy.ndarray
        :param area: Beam area in steradians
        :type name: float or numpy.ndarray
        :param ame_ampl: AME peak flux density
        :type ame_ampl: float
        :param ame_peak: AME peak frequency
        :type ame_peak: float
        :param ame_width: AME lognormal width
        :type ame_width: float

        :return: Evaluated function
        :rtype: float or numpy.ndarray
        """
        nlog = np.log(nu)
        nmaxlog = np.log(ame_peak)
        return ame_ampl*np.exp(-0.5 * ((nlog-nmaxlog)/ame_width)**2)
    
    def deriv(self,nu, area, ame_ampl, ame_peak, ame_width):
        """Evaluate the first derivitives of emission model with respect to input parameters.

        :param nu: Frequency in GHz
        :type nu: float or numpy.ndarray
        :param area: Beam area in steradians
        :type name: float or numpy.ndarray
        :param ame_ampl: AME peak flux density
        :type ame_ampl: float
        :param ame_peak: AME peak frequency
        :type ame_peak: float
        :param ame_width: AME lognormal width
        :type ame_width: float

        :return: First derivitives with repect to input parameters in order.
        :rtype: list
        """
        def evaluate2(nu, area, ame_ampl, ame_peak, ame_width):
            nlog = np.log(nu)
            nmaxlog = np.log(ame_peak)
            return ame_ampl*np.exp(-0.5 * ((nlog-nmaxlog)/ame_width)**2)
        S = evaluate2(nu, area, ame_ampl, ame_peak, ame_width)
        d_ame_ampl = S / ame_ampl
        d_ame_peak = S * -1 * ((np.log(ame_peak) - np.log(nu))/(ame_peak*(ame_width**2)))
        d_ame_width = S*((np.log(ame_peak) - np.log(nu))**2 / (ame_width**3))
        return [d_ame_ampl, d_ame_peak, d_ame_width]

# Thermal dust
class thermalDust(FittableEmissionModel):
    """Emission model for the Planck modified thermal dust curve - a modified blackbody with opacity varying as frequency to some dust spectral index.

    :params tdust_Td: Thermal dust temperature
    :type tdust_Td: float
    :params tdust_tau: Thermal dust opacity (given as log_10(tau))
    :type tdust_tau: float
    :params tdust_beta: Thermal dust spectral index
    :type tdust_beta: float

    :Formulism:
        .. math:: S^\\mathrm{td}_\\nu = \\frac{2k_B\\Omega\\nu^3}{c^2} \\frac{1}{e^{h\\nu/k_BT_b}-1} \\cdot \\tau_{\\nu_0} \\cdot \\left(\\frac{\\nu}{\\nu_0}\\right)^\\beta
        For details see `Draine and Li (2001)`_ .
    
    :notes: Either all or none of input ``nu``, ``area``,  ``tdust_Td``, ``tdust_tau`` and ``tdust_beta`` must be provided consistently with compatible units or as unitless numbers.
        
    .. _Draine and Li (2001): https://ui.adsabs.harvard.edu/abs/2001ApJ...551..807D/abstract
    """
    tdust_Td = Parameter(default=20,
                         min=0,
                         description='Thermal dust temperature')
    tdust_tau = Parameter(default=-4,
                          description="Thermal dust opacity (given as log_10(tau))")
    tdust_beta = Parameter(default=1.5,
                           description="thermal dust spectral index")

    @staticmethod
    def evaluate(nu, area, tdust_Td, tdust_tau, tdust_beta):
        """Evaluate the emission model.

        :param nu: Frequency in GHz
        :type nu: float or numpy.ndarray
        :param area: Beam area in steradians
        :type name: float or numpy.ndarray
        :params tdust_Td: Thermal dust temperature
        :type tdust_Td: float
        :params tdust_tau: Thermal dust opacity (given as log_10(tau))
        :type tdust_tau: float
        :params tdust_beta: Thermal dust spectral index
        :type tdust_beta: float

        :return: Evaluated function
        :rtype: float or numpy.ndarray
        """
        nu9 = np.multiply(nu,1e9)
        planck = np.exp(phys_const['h']*nu9/phys_const['k']/tdust_Td) - 1.
        modify = 10**tdust_tau * (nu9/1.2e12)**tdust_beta
        return 2 * phys_const['h'] * nu9**3/phys_const['c']**2 /planck * modify * area * 1e26
    
    def deriv(self,nu, area, tdust_Td, tdust_tau, tdust_beta):
        """Evaluate the first derivitives of emission model with respect to input parameters.

        :param nu: Frequency in GHz
        :type nu: float or numpy.ndarray
        :param area: Beam area in steradians
        :type name: float or numpy.ndarray
        :params area: Thermal dust temperature
        :type tdust_Td: float
        :params tdust_tau: Thermal dust opacity (given as log_10(tau))
        :type tdust_tau: float
        :params tdust_beta: Thermal dust spectral index
        :type tdust_beta: float

        :return: First derivitives with repect to input parameters in order.
        :rtype: list
        """
        def evaluate2(nu, area, tdust_Td, tdust_tau, tdust_beta):
            nu9 = np.multiply(nu,1e9)
            planck = np.exp(phys_const['h']*nu9/phys_const['k']/tdust_Td) - 1.
            modify = 10**tdust_tau * (nu9/1.2e12)**tdust_beta
            return 2 * phys_const['h'] * nu9**3/phys_const['c']**2 /planck * modify * area * 1e26
        S = evaluate2(nu, area, tdust_Td, tdust_tau, tdust_beta)
        nu9 = np.multiply(nu,1e9)
        hvkT = phys_const['h']*nu9 / (phys_const['k']*tdust_Td)
        d_tdust_Td = S * (hvkT) * (1/tdust_Td) * (np.exp(hvkT) / (np.exp(hvkT) - 1))
        d_tdust_tau = S * np.log(10)
        d_tdust_beta = S * np.log(nu9/353e9)
        return [d_tdust_Td, d_tdust_tau, d_tdust_beta]