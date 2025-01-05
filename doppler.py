import numpy as np
from scipy import signal
import matplotlib.pyplot as plt
from my_plot import TEX_FONTS, LATEX_SUBPLOT_2, PlotStyle 

def PlotDoppler(wave: np.ndarray, flux: np.ndarray, flux_shift: np.ndarray, offsets: np.ndarray, cross: np.ndarray, name: str='doppler_0.png') -> None:
    '''
    Plots the original spectrum and the Doppler shifted spectrum along with the cross-correlation plot.
    '''
    with PlotStyle('science', TEX_FONTS):
        fig, ax = plt.subplots(1,2,figsize=LATEX_SUBPLOT_2)

        ax[0].plot(wave,flux,label='Original Spectrum',color='green')
        ax[0].plot(wave,flux_shift,label='Doppler Shifted Spectrum',color='C3')
        ax[0].legend(prop={'size': 6},loc=2)
        ax[0].set_ylim(3e-15,8.5e-15)
        ax[0].set_xlim(1.05,2.4)
        ax[0].set_xlabel(r'Wavelength ($\mu m$)')
        ax[0].set_ylabel(r'Flux ($W \: m^2 \: \mu m^{-1}$)')
        ax[0].grid()

        ax[1].plot(offsets/10000,cross,color='olivedrab')
        ax[1].set_xlim(-38,38)
        ax[1].grid()
        ax[1].set_xlabel('Velocity ($10^7$ m/s)')
        ax[1].axvline(x=0,color='r')
        ax[1].set_xticks(np.linspace(-30,30,7))
        
        pos = np.argmax(cross)  # Position of the maximum
        axins = ax[1].inset_axes([0.60, 0.60, 0.4, 0.4])
        # sub region of the original image
        axins.plot(offsets/10000,cross,color='olivedrab')
        axins.axvline(x=0,color='r')
        axins.set_xlim(offsets[pos]/10000-1, offsets[pos]/10000+1)
        axins.set_ylim(cross[pos]-.15e-24, cross[pos]+.1e-24)
        axins.grid()
        axins.set_yticks([])
        axins.tick_params(axis='both', which='major', labelsize=3.5)
        ax[1].indicate_inset_zoom(axins, edgecolor="black")

        plt.savefig(f'plots/{name}',bbox_inches='tight',dpi=300)

def CrossCorrelation(wave: np.ndarray, flux: np.ndarray, flux_shifted: np.ndarray) -> tuple:
    '''
    Cross-correlates the original spectrum with the Doppler shifted spectrum to recover
    the velocity of the shift.
    '''
    cross=signal.correlate(flux,flux_shifted)

    xcorr = np.arange(cross.size)

    # Convert this into lag units
    lags = xcorr - (flux.size-1)
    dist_per_lag = (wave[-1] - wave[0])/float(wave.size)/1.185  # This is just the x-spacing

    # Convert your lags into physical units
    offsets = -lags*dist_per_lag*C/1000   # C is the speed of light in m/s

    # Recover the velocity from the offset with an error
    pos = np.argmax(cross)
    print(f'Recovered Velocity = {offsets[pos]:.5}(+-{offsets[pos]-offsets[pos+1]:.2}) km/s')

    return offsets, cross

        
if __name__ == '__main__':
    from PyAstronomy import pyasl
    from constants import *

    wave1=np.load('spectral_data/wave1_convol.npy')
    flux1=np.load('spectral_data/flux1_convol.npy')

    wave2=np.load('spectral_data/wave2_convol.npy')
    flux2=np.load('spectral_data/flux2_convol.npy')

    velocity = 3000.0 # km/s
    print(f'Real Velocity = {velocity} km/s')

    print('Cross Correlation for Spectrum 1:')
    flux1_shifted, wave1_shifted = pyasl.dopplerShift(wave1, flux1, velocity, edgeHandling="firstlast")
    offsets1, cross1 = CrossCorrelation(wave1, flux1, flux1_shifted)
    PlotDoppler(wave1, flux1, flux1_shifted, offsets1, cross1, name='doppler_1.png')

    print('Cross Correlation for Spectrum 2:')
    flux2_shifted, wave2_shifted = pyasl.dopplerShift(wave2, flux2, velocity, edgeHandling="firstlast")
    offsets2, cross2 = CrossCorrelation(wave2, flux2, flux2_shifted)
    PlotDoppler(wave2, flux2, flux2_shifted, offsets2, cross2, name='doppler_2.png')