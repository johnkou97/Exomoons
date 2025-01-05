from tqdm import tqdm
import numpy as np
from scipy import signal, stats

def Rebin(wave: np.ndarray, flux: np.ndarray, start: float = 1, end: float = 2.5, dw: float = 1.5/100000) -> tuple:
    '''
    Rebin the spectrum to a new wavelength grid
    '''
    wave_rebin = np.arange(start, end, dw)
    wave_low = wave_rebin - dw / 2
    wave_up = wave_rebin + dw / 2
    flux_rebin = np.zeros_like(wave_rebin)
    
    for i,(l,u) in tqdm(enumerate(zip(wave_low,wave_up)),total=len(wave_rebin)):
        mask = (wave>l)*(wave<u)
        for j in range(len(mask)):
            if mask[j]:
                flux_rebin[i]=np.mean(flux[mask])
                break

    return wave_rebin, flux_rebin

def TukeyWindow(flux_rebin: np.ndarray, alpha: float = 0.03) -> np.ndarray:
    '''
    Apply Tukey window to the rebinned flux
    '''
    window = signal.windows.tukey(len(flux_rebin), alpha=alpha)
    flux_window = flux_rebin * window

    return flux_window

def GaussianKernel(mean: float = 0, standard_deviation: float = 1, n: int = 1080) -> np.ndarray:
    '''
    Create a Gaussian kernel
    '''
    x_values = np.linspace(-5 * standard_deviation, 5 * standard_deviation, n)
    y_values = stats.norm(mean, standard_deviation).pdf(x_values)
    kernel = y_values / sum(y_values)

    return kernel

def Convolve(wave_rebin: np.ndarray, flux_rebin: np.ndarray, kernel: np.ndarray) -> tuple:
    '''
    Convolve the flux with the kernel
    '''
    flux_convloved = signal.convolve(flux_rebin, kernel, mode='same')
    wave_convolved = np.linspace(wave_rebin[0], wave_rebin[-1], len(flux_convloved))

    return wave_convolved, flux_convloved

if __name__ == '__main__':
    import matplotlib.pyplot as plt
    from my_plot import PlotStyle

    highres=np.loadtxt('spectral_data/betapicb_highres.dat')
    wave1=highres[:,0]
    flux1=highres[:,1]

    spectrum=np.loadtxt('spectral_data/beta_pic_high_res_1.0_2.5_micron.dat')
    wave2=spectrum[:,0]
    flux2=spectrum[:,1]

    dw = 1.5/100000

    wave_rebin_1, flux_rebin_1 = Rebin(wave1, flux1, 1, 2.5, dw)
    np.save('spectral_data/wave1_rebin', wave_rebin_1)
    np.save('spectral_data/flux1_rebin', flux_rebin_1)

    wave_rebin_2, flux_rebin_2 = Rebin(wave2, flux2, 1, 2.5, dw)
    np.save('spectral_data/wave2_rebin', wave_rebin_2)
    np.save('spectral_data/flux2_rebin', flux_rebin_2)

    wave_rebin_1=np.load('spectral_data/wave1_rebin.npy')
    flux_rebin_1=np.load('spectral_data/flux1_rebin.npy')

    flux_window_1 = TukeyWindow(flux_rebin_1)
    flux_window_2 = TukeyWindow(flux_rebin_2)

    kernel = GaussianKernel(n=1080)

    wave_convolved_1, flux_convolved_1 = Convolve(wave_rebin_1, flux_window_1, kernel)
    np.save('spectral_data/wave1_convol', wave_convolved_1)
    np.save('spectral_data/flux1_convol', flux_convolved_1)

    wave_convolved_2, flux_convolved_2 = Convolve(wave_rebin_2, flux_window_2, kernel)
    np.save('spectral_data/wave2_convol', wave_convolved_2)
    np.save('spectral_data/flux2_convol', flux_convolved_2)

    # Plotting

    tex_fonts = {
        "axes.labelsize": 30,
        "font.size": 30,
        "legend.fontsize": 25,
        "xtick.labelsize": 25,
        "ytick.labelsize": 25,
    }

    # Load the data
    wave_convolved_1=np.load('spectral_data/wave1_convol.npy')
    flux_convolved_1=np.load('spectral_data/flux1_convol.npy')
    wave_rebin_1=np.load('spectral_data/wave1_rebin.npy')
    flux_rebin_1=np.load('spectral_data/flux1_rebin.npy')

    with PlotStyle('science', tex_fonts):
        plt.figure(figsize=(16,8))
        plt.scatter(wave1,flux1,marker='.',s=.3,label='original spectrum')
        plt.scatter(wave_rebin_1,flux_rebin_1,marker='.',s=.6,label='rebinned spectrum')
        plt.plot(wave_convolved_1,flux_convolved_1,color='C3',label='convolved spectrum')

        lgnd=plt.legend(scatterpoints=3,loc='upper right')
        lgnd.legendHandles[0]._sizes = [30]
        lgnd.legendHandles[1]._sizes = [30]

        plt.xlabel(r'Wavelength ($\mu m$)')
        plt.ylabel(r'Flux ($W \: m^2 \: \mu m^{-1}$)')
        plt.ylim(0,9.5e-15)
        plt.yticks(np.arange(0,10e-15,1e-15))
        plt.xlim(0.99,2.51)
        plt.savefig('plots/spectrum_1.png',dpi=300)

    # Load the data
    wave_convolved_2=np.load('spectral_data/wave2_convol.npy')
    flux_convolved_2=np.load('spectral_data/flux2_convol.npy')
    wave_rebin_2=np.load('spectral_data/wave2_rebin.npy')
    flux_rebin_2=np.load('spectral_data/flux2_rebin.npy')

    with PlotStyle('science', tex_fonts):
        plt.figure(figsize=(16,8))
        plt.scatter(wave2,flux2,marker='.',s=.3,label='original spectrum')
        plt.scatter(wave_rebin_2,flux_rebin_2,marker='.',s=.6,label='rebinned spectrum')
        plt.plot(wave_convolved_2,flux_convolved_2,color='C3',label='convolved spectrum')

        lgnd=plt.legend(scatterpoints=3,loc='upper right')
        lgnd.legendHandles[0]._sizes = [30]
        lgnd.legendHandles[1]._sizes = [30]

        plt.xlabel(r'Wavelength ($\mu m$)')
        plt.ylabel(r'Flux ($W \: m^2 \: \mu m^{-1}$)')
        plt.ylim(0,9.5e-15)
        plt.yticks(np.arange(0,10e-15,1e-15))
        plt.xlim(0.99,2.51)
        plt.savefig('plots/spectrum_2.png',dpi=300)