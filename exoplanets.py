if __name__ == '__main__':
    import os
    import pandas as pd
    import matplotlib.pyplot as plt
    from my_plot import TEX_FONTS, LATEX_PLOT, PlotStyle

    if not os.path.exists('plots'):
        os.makedirs('plots')

    df=pd.read_csv('exoplanet_catalog/exoplanet.eu_catalog.csv')
    rv=pd.read_csv('exoplanet_catalog/radial_velocity_exoplanets.csv')
    im=pd.read_csv('exoplanet_catalog/imaging_exoplanets.csv')
    tr=pd.read_csv('exoplanet_catalog/transit_exoplanets.csv')

    with PlotStyle('science', TEX_FONTS, LATEX_PLOT):
        plt.scatter(df['semi_major_axis'],df['mass'],s=1,color='grey',label='Other techniques')
        plt.scatter(tr['semi_major_axis'],tr['mass'],s=1,color='yellow',label='Transit')
        plt.scatter(rv['semi_major_axis'],rv['mass'],s=1,color='green',label='Radial Velocity')
        plt.scatter(im['semi_major_axis'],im['mass'],s=1,color='purple',label='Direct Imaging')
        plt.xscale('log')
        plt.yscale('log')
        plt.ylabel(r'Mass $(M_{jup})$')
        plt.xlabel(r'Semi-Major axis $(AU)$')
        plt.xticks([0.01,0.1,1,10,100,1000,10000])
        plt.yticks([0.0001,0.001,0.01,0.1,1,10,100])
        leg=plt.legend(fontsize=8,scatterpoints=3,frameon=True)
        leg.get_frame().set_edgecolor('black')
        plt.savefig('plots/exoplanets.png',dpi=600,bbox_inches='tight')
