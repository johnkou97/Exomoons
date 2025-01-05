import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.optimize import curve_fit
from scipy.signal import lombscargle
from my_plot import TEX_FONTS, LATEX_PLOT, LATEX_SUBPLOT_2, PlotStyle
from constants import *

def FitFunction(t: np.ndarray, a: float, b: float, c: float, d: float) -> np.ndarray:
    '''
    Simple function to fit the data
    '''
    return   a*t + b* np.sin(c * t + d)

def Analyze(data: str = '11') -> tuple:
    '''
    Analyze the data from the csv files
    '''
    print(f'Analyzing data {data}')

    columns=['# Time (hours)',' RV',' -uncertainty',' +uncertainty']
    rvs_table = pd.read_csv(f'rv_data/preliminary_rvs_2022_11_{data}.csv', names=columns, header=None)
    rvs_fit = rvs_table.sort_values('# Time (hours)', ignore_index=True)

    time_beta = rvs_fit['# Time (hours)'].to_numpy()
    rv_beta = rvs_fit[' RV'].to_numpy()

    t = time_beta - (time_beta[0] + time_beta[-1]) / 2  # Centering the time
    y = rv_beta - np.mean(rv_beta)  # Centering the radial velocity
    errors = np.vstack((rvs_fit[' +uncertainty'], rvs_fit[' -uncertainty']))

    with PlotStyle('science', fonts=TEX_FONTS, figsize=LATEX_PLOT):
        plt.errorbar(t, y, errors, fmt='.', color='black')
        plt.ylabel('Radial Velocity (km/s)')
        plt.xlabel('Time (hours)')
        plt.savefig(f'plots_rossiter/rv_data_{data}.png', dpi=300)

    freqs = np.linspace(0.01, 30, 1000)
    amp = lombscargle(t, y, freqs)

    max_freq = freqs[np.argmax(amp)]
    period = 2 * np.pi / max_freq

    with PlotStyle('science', fonts=TEX_FONTS, figsize=LATEX_PLOT):
        plt.plot(freqs, amp, color='black')
        plt.ylabel('Amplitude')
        plt.xlabel('Period (hours)')
        plt.vlines(max_freq, 0, 100, color='r', linestyles='--')
        plt.ylim(0, np.amax(amp) + 1)
        plt.xlim(0, np.amax(freqs))
        plt.savefig(f'plots_rossiter/fft_{data}.png', dpi=300)

    print(f'Main Frequency: {max_freq}')
    print(f'Period: {period}')

    amp_est = (np.amax(y) - np.amin(y)) / 2
    print(f'Amplitude Estimation: {amp_est}')

    params, params_covariance = curve_fit(FitFunction, t, y, p0=[2, amp_est, max_freq, 0], absolute_sigma=True)

    with PlotStyle('science', fonts=TEX_FONTS, figsize=LATEX_PLOT):
        plt.errorbar(t, y, errors, fmt='.', color='black', label='Data')
        x = np.linspace(t[0], t[-1], 1000)
        plt.plot(x, FitFunction(x, params[0], params[1], params[2], params[3]), color='r', label='Fitted function')
        # also have the error bands
        plt.fill_between(x, FitFunction(x, params[0], params[1], params[2], params[3]), 
            FitFunction(x, params[0] + np.sqrt(params_covariance[0, 0]), params[1] + np.sqrt(params_covariance[1, 1]), params[2] + np.sqrt(params_covariance[2, 2]), params[3] + np.sqrt(params_covariance[3, 3]),), 
            color='r', alpha=0.3)
        plt.fill_between(x, FitFunction(x, params[0], params[1], params[2], params[3]), 
            FitFunction(x, params[0] - np.sqrt(params_covariance[0, 0]), params[1] - np.sqrt(params_covariance[1, 1]), params[2] - np.sqrt(params_covariance[2, 2]), params[3] - np.sqrt(params_covariance[3, 3]),), 
            color='r', alpha=0.3)
        plt.legend()
        plt.ylabel('Radial Velocity (km/s)')
        plt.xlabel('Time (hours)')
        plt.savefig(f'plots_rossiter/rv_fit_{data}.png', dpi=300)

    print(f'Main Frequency: {params[2] * 1000 / (60 * 60)}')
    print(f'Error of Main Frequency: {np.sqrt(params_covariance[2, 2] * 1000 / (60 * 60))}')
    print(f'Period: {2 * np.pi / params[2]}')

    print(f'Amplitude Estimation: {abs(params[1])}')
    print(f'Error of Amplitude: {np.sqrt(params_covariance[1, 1])}')

    return t, y, errors, freqs, amp, max_freq, params, params_covariance

if __name__ == '__main__':
    import os

    if not os.path.exists('plots_rossiter'):
        os.makedirs('plots_rossiter')

    t_11, y_11, errors_11, freqs_11, amp_11, max_freq_11, _, _ = Analyze('11')
    t_13, y_13, errors_13, freqs_13, amp_13, max_freq_13, _, _ = Analyze('13')

    # Plots for Thesis

    with PlotStyle('science', fonts=TEX_FONTS, figsize=None):
        fig, ax = plt.subplots(1, 2, figsize=LATEX_SUBPLOT_2)
        ax[0].errorbar(t_11, y_11, errors_11, fmt='.', color='black')
        ax[1].errorbar(t_13, y_13, errors_13, fmt='.', color='black')
        fig.add_subplot(111, frameon=False)
        plt.tick_params(labelcolor='none', which='both', top=False, bottom=False, left=False, right=False)
        plt.ylabel('Radial Velocity (km/s)')
        plt.xlabel('Time (hours)')
        plt.savefig('plots_rossiter/rv_full.png', bbox_inches='tight', dpi=300)

    with PlotStyle('science', fonts=TEX_FONTS, figsize=None):
        fig, ax = plt.subplots(1, 2, figsize=LATEX_SUBPLOT_2)
        ax[0].plot(freqs_11 * 1000 / (60 * 60), amp_11, color='black')
        ax[0].vlines(max_freq_11 * 1000 / (60 * 60), 0, 100, color='r', linestyles='--')
        ax[0].set_ylim(0, np.amax(amp_11) + 1)
        ax[1].plot(freqs_13 * 1000 / (60 * 60), amp_13, color='black')
        ax[1].vlines(max_freq_13 * 1000 / (60 * 60), 0, 100, color='r', linestyles='--')
        ax[1].set_ylim(0, np.amax(amp_13) + 1)
        fig.add_subplot(111, frameon=False)
        plt.tick_params(labelcolor='none', which='both', top=False, bottom=False, left=False, right=False)
        plt.ylabel('Amplitude')
        plt.xlabel('Angular Frequency (mHz)')
        plt.savefig('plots_rossiter/lombscargle_both.png', bbox_inches='tight', dpi=300)