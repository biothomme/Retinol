#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 11 11:07:14 2021

@author: Thomsn
"""
import matplotlib.pyplot as plt
import seaborn as sns

def default():
    import matplotlib as mpl
    mpl.rcParams['pdf.fonttype'] = 42
    return

def single_plot(color_df, p_threshold):
    COLORS = ['darkorchid',
              '#FCC200',
              'mediumblue',
              'darkgreen',
              'coral']
    
    species_list = set(color_df.keys().get_level_values(0))

    x = color_df.index
    fig = plt.figure(figsize=(8, 6))
    ax = plt.subplot(111)
        
    # wavelengths with significant differences - threshold .001!
    try:
        anova_results_wls = wlanova(color_df)
    except IndexError:
        significant_wls = []
    else:
        significant_wls = [yv for sw, yv in zip(anova_results_wls, x)
                           if sw < p_threshold]
        
    for significant_wl in significant_wls:
        if significant_wl > 380:
            color = wavelength_to_rgb(significant_wl, gamma=0.8)
        else:
            color = 'gray'
        plt.axvline(significant_wl, c=color, alpha=.05, linewidth=2)
        
    # lines and mean lines
    for i, species in enumerate(species_list):
        y = color_df[species]
        individuals = set(y.columns)
        for individual in individuals:
            ax.plot(x, y[individual], color = COLORS[i], alpha = .2)
        ax.plot(x, y.mean(axis=1).values, color = COLORS[i], label = species, linewidth = 2)

    ax.legend(loc='upper left', title='Species')
    ax.set_xlabel('Wavelength / nm')
    ax.set_ylabel('Diffuse reflexion')
    return fig


def wlanova(color_df, bonferroni=True):
    from scipy.stats import f_oneway
    from statsmodels.sandbox.stats.multicomp import multipletests

    species_list = list(set(color_df.columns.get_level_values(0)))
    p_list = []
    f_list = []
    for _, row in color_df.iterrows():
        f, p = f_oneway(row[species_list[0]].values,
                        row[species_list[1]].values)
        f_list.append(f)
        p_list.append(p)
    if bonferroni:
        p_list = multipletests(p_list, method='bonferroni')[1]
    return p_list


def wavelength_to_rgb(wavelength, gamma=0.8):
    '''This converts a given wavelength of light to an 
    approximate RGB color value. The wavelength must be given
    in nanometers in the range from 380 nm through 750 nm
    (789 THz through 400 THz).

    Based on code by Dan Bruton
    http://www.physics.sfasu.edu/astro/color/spectra.html
    '''
    wavelength = float(wavelength)
    if wavelength >= 380 and wavelength <= 440:
        attenuation = 0.3 + 0.7 * (wavelength - 380) / (440 - 380)
        R = ((-(wavelength - 440) / (440 - 380)) * attenuation) ** gamma
        G = 0.0
        B = (1.0 * attenuation) ** gamma
    elif wavelength >= 440 and wavelength <= 490:
        R = 0.0
        G = ((wavelength - 440) / (490 - 440)) ** gamma
        B = 1.0
    elif wavelength >= 490 and wavelength <= 510:
        R = 0.0
        G = 1.0
        B = (-(wavelength - 510) / (510 - 490)) ** gamma
    elif wavelength >= 510 and wavelength <= 580:
        R = ((wavelength - 510) / (580 - 510)) ** gamma
        G = 1.0
        B = 0.0
    elif wavelength >= 580 and wavelength <= 645:
        R = 1.0
        G = (-(wavelength - 645) / (645 - 580)) ** gamma
        B = 0.0
    elif wavelength >= 645 and wavelength <= 750:
        attenuation = 0.3 + 0.7 * (750 - wavelength) / (750 - 645)
        R = (1.0 * attenuation) ** gamma
        G = 0.0
        B = 0.0
    else:
        R = 0.0
        G = 0.0
        B = 0.0
    return (R, G, B)

def polygon_plot(received_signals_df, spectrum_loci_df, plot_type="hexagon",
                 axis_label=True, spectrum_loci_annotations=True):
    from convwale import Perceived_Signals
    import pandas as pd

    signals_received = Perceived_Signals(received_signals_df)
    signals_sl = Perceived_Signals(spectrum_loci_df)
    signals_received.get_x()
    signals_received.get_y()

    fig = plt.figure(figsize=(12, 9))
    ax = plt.subplot(111)

    ax.fill(signals_sl.get_x(), signals_sl.get_y(), alpha=.2, color="grey")
    plotting_frame = pd.DataFrame({
        "x": signals_received.data["x"],
        "y": signals_received.data["y"],
        "taxon": signals_received.get_taxa(),
        "leaf area":signals_received.data.index.get_level_values(1)})
    if len(pd.unique(plotting_frame["leaf area"])) == 1:
        sns.scatterplot(data=plotting_frame, x="x", y="y", ax=ax, s=60,
                        hue="taxon", style="leaf area", alpha=.9)
    else:
        sns.scatterplot(data=plotting_frame, x="x", y="y", ax=ax, s=60,
                        style="taxon", hue="leaf area", alpha=.9)
    if spectrum_loci_annotations:
        for i, row in signals_sl.data.iterrows():
            if i%50 == 0:
                ax.text(row["x"], row["y"], f"{i}", alpha=.5, 
                        verticalalignment="center", horizontalalignment="center")

    # ax.axis("off")
    ax.spines['left'].set_position('zero')
    ax.spines['left'].set_color('gray')
    ax.spines['bottom'].set_position('zero')
    ax.spines['bottom'].set_color('gray')
    ax.spines['right'].set_color('none')
    ax.spines['top'].set_color('none')
    ax.set_xlabel("")
    ax.set_ylabel("")
    ax.legend(loc="upper left", bbox_to_anchor=(.9, 1))
    if axis_label:
        ax.tick_params(color="gray", labelcolor="gray")
    else:
        ax.tick_params(color="gray", labelcolor="none")
    plt.gca().set_aspect('equal', adjustable='box')

    for i, receptor in enumerate(signals_received.data.columns.values[0:3]):
        receptor2 = signals_received.data.columns.values[0:3][(i+1)//3]
        va = "bottom" if i == 1 else "top"
        ha = "left" if i == 2 else ("center" if i == 1 else "right")
        ax.text(signals_received.TRIANGLE_COORDINATES[0][i],
                signals_received.TRIANGLE_COORDINATES[1][i],
                receptor, fontsize=10,
                verticalalignment=va, horizontalalignment=ha)
    if plot_type.lower() == "hexagon":
        for i, receptor in enumerate(signals_received.data.columns.values[0:3]):
            receptor2 = signals_received.data.columns.values[0:3][(i+1)%3]
            va = "top" if i == 2 else "bottom"
            ha = "center" if i == 2 else ("left" if i == 1 else "right")
            ax.text(signals_received.HEXAGON_COORDINATES[0][2*i+1],
                    signals_received.HEXAGON_COORDINATES[1][2*i+1],
                    f'{receptor}/\n{receptor2}', fontsize=10,
                    verticalalignment=va, horizontalalignment=ha)
        ax.plot(signals_received.HEXAGON_COORDINATES[0],
            signals_received.HEXAGON_COORDINATES[1], color="grey")
        ax.set_xlim([min(signals_received.HEXAGON_COORDINATES[0]),
                     max(signals_received.HEXAGON_COORDINATES[0])])
        ax.set_ylim([min(signals_received.HEXAGON_COORDINATES[1]),
                     max(signals_received.HEXAGON_COORDINATES[1])])
    elif plot_type.lower() == "triangle":
        ax.plot(signals_received.TRIANGLE_COORDINATES[0],
            signals_received.TRIANGLE_COORDINATES[1], color="grey")
        ax.set_xlim([min(signals_received.TRIANGLE_COORDINATES[0]),
                     max(signals_received.TRIANGLE_COORDINATES[0])])
        ax.set_ylim([min(signals_received.TRIANGLE_COORDINATES[1]),
                     max(signals_received.TRIANGLE_COORDINATES[1])])
    else:
        print(f"""Error: Plot type {plot_type} was not recognized. Please try 
              either 'hexagon' or 'triangle'.""")
        return
    return fig

default()