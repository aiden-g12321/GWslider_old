import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import trapezoid as integrate
from aLIGO_sensitivity.aLIGO_sensitivity import *
from constants import *
from wave_gen import times, waveform_inj, waveform_inj_FD


class GW_signal:

    def __init__(self, times, freqs, waveformTD, waveformFD, mass1, mass2, chiPlus, chiMinus):

        self.times= times
        self.freqs= freqs
        self.waveformTD= waveformTD
        self.waveformFD= waveformFD
        self.mass1= mass1
        self.mass2= mass2
        self.ratio= self.mass1/self.mass2
        self.chirp= mchirp_from_mass1_mass2(self.mass1, self.mass2)
        self.chi1= spin1z_from_mass1_mass2_chi_eff_chi_a(self.mass1, self.mass2, self.chiPlus, self.chiMinus)
        self.chi2= spin2z_from_mass1_mass2_chi_eff_chi_a(self.mass1, self.mass2, self.chiPlus, self.chiMinus)
        self.chiPlus= chiPlus
        self.chiMinus= chiMinus

        self.sqrtS = np.array([sqrtS(f) for f in self.freqs])
        self.df= self.freqs[1]-self.freqs[0]

    def plotTD(self):
        plt.plot(self.times, self.waveformTD)
        plt.show()
    
    def plotFD(self):
        plt.plot(self.freqs, self.waveformFD)
        plt.show()
    
    def PSD_LIGO(self):
         plt.loglog(self.freqs, self.sqrtS)
         plt.show()
    
    #innr product in the frequency domain
    def innerProduct(self,a):
        integrand = ((np.real(a)*np.real(self.waveformFD) + np.imag(a)*np.imag(self.waveformFD)))
        inner_prod = integrate(integrand, dx=self.df)
        return inner_prod * 4.
    



GW150914= GW_signal(data_times, freqs, data_strains, intrp_waveform, mass1_150914, mass2_150914, chiPlus_150914, chiMinus_150914)
simulatedGW= GW_signal(times, freqs_cut, waveform_inj, waveform_inj_FD, m1_inj, m2_inj,spin_plus_inj, spin_minus_inj)
 


simulatedGW.plotTD()
simulatedGW.plotFD()

