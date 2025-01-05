import scienceplots						# Enables the 'science' style
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import corner
from contextlib import contextmanager	# Used to create custom plot context manager

'''
This script contains the custom plot settings for the project. It includes the following:
- TEX_FONTS: A dictionary containing the LaTeX font settings for the plots.
- CORNER_KWARGS: A dictionary containing the settings for the corner plot.
- SetSize: A function to set the size of the plot in inches.
- LATEX_PLOT: The size of the plot for LaTeX.
- LATEX_SUBPLOT_2: The size of the subplot for LaTeX.
- PlotStyle: A context manager to set the plot style, fonts, and figure size.
'''

TEX_FONTS = {
    "axes.labelsize": 12,
	"font.size": 12,
	"legend.fontsize": 10,
	"xtick.labelsize": 10,
	"ytick.labelsize": 10}

CORNER_KWARGS = dict(
	# smooth=0.9,
	# label_kwargs=dict(fontsize=16),
	# title_kwargs=dict(fontsize=16),
	quantiles=[0.16, 0.50, 0.84],
	levels=[0.20, 0.40, 0.85, 0.99],#(1 - np.exp(-0.5), 1 - np.exp(-2), 1 - np.exp(-9 / 2.)),
	plot_density=False,
	plot_datapoints=False,
	fill_contours=True,
	show_titles=False,
	# max_n_ticks=3,
	title_fmt=".2E")

def SetSize(width: float, fraction: float = 1, subplots: tuple = (1, 1)) -> tuple:
	"""
    Set figure dimensions to avoid scaling in LaTeX.
    Taken from https://jwalton.info/Embed-Publication-Matplotlib-Latex/

	Parameters
	----------
	width: float or string
			Document width in points, or string of predined document type
	fraction: float, optional
			Fraction of the width which you wish the figure to occupy
	subplots: array-like, optional
			The number of rows and columns of subplots.
	Returns
	-------
	fig_dim: tuple
			Dimensions of figure in inches
	"""
	# if width == 'thesis':
	# 	width_pt = 426.79135
	# elif width == 'beamer':
	# 	width_pt = 307.28987
	# else:
	width_pt = width

	# Width of figure (in pts)
	fig_width_pt = width_pt * fraction
	# Convert from pt to inches
	inches_per_pt = 1 / 72.27

	# Golden ratio to set aesthetic figure height
	# https://disq.us/p/2940ij3
	golden_ratio = (5**.5 - 1) / 2

	# Figure width in inches
	fig_width_in = fig_width_pt * inches_per_pt
	# Figure height in inches
	fig_height_in = fig_width_in * golden_ratio * (subplots[0] / subplots[1])

	return (fig_width_in, fig_height_in)

LATEX_PLOT = SetSize(390.0)
LATEX_SUBPLOT_2 = SetSize(390.0, subplots=(1, 2))

# Custom context manager
@contextmanager
def PlotStyle(style: str = 'science', fonts: dict = None, figsize: tuple = None):
    if style:
        mpl.style.use(style)
    if fonts:
        plt.rcParams.update(fonts)
    if figsize:
        plt.figure(figsize=figsize)
    yield
    plt.close()
    plt.style.use('default')  # Reset to default after use

if __name__ == '__main__':

	print(f'Latex size for a single plot: {LATEX_PLOT}')
	print(f'Latex size for a subplot of shape (1, 2): {LATEX_SUBPLOT_2}')