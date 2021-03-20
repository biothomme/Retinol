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
              'coral',
              'fuchsia',
              'y',
              'tomato',
              'teal',
              'papayawhip',
              'darkslategrey',
              'black',
              'dodgerblue',
              'cyan']
    
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

    ax.spines['left'].set_position('zero')
    ax.spines['left'].set_color('gray')
    ax.spines['bottom'].set_position('zero')
    ax.spines['bottom'].set_color('gray')
    ax.spines['right'].set_color('none')
    ax.spines['top'].set_color('none')
    ax.set_xlabel("")
    ax.set_ylabel("")
    ax.legend(loc="upper left", bbox_to_anchor=(1, 1))
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


def distance_heatmap(pairwise_color_df):
    import numpy as np
    pairwise_color_dist = pairwise_color_df.copy().sort_index().sort_index(axis=1)
    fig = plt.figure(figsize=(8, 6))
    ax = plt.subplot(111)

    im = ax.imshow(pairwise_color_dist)
    cbar = ax.figure.colorbar(im)
    cbar.ax.set_ylabel("color distance", rotation=-90, va="bottom")
    ax.set_xticks(np.arange(pairwise_color_dist.shape[0]))
    ax.set_yticks(np.arange(pairwise_color_dist.shape[1]))
    short_names = [
        f"{x[0][:2]}_{x[2][:2]}_{x[1]}_{x[3]}"
        for x in pairwise_color_dist.index]
    ax.set_xticklabels(
        short_names, rotation=45, ha="right", rotation_mode="anchor")
    ax.set_yticklabels(short_names)
    return fig


def distance_dendrogram(pairwise_color_df):
    from scipy.cluster import hierarchy
    from matplotlib import cm
    pairwise_color_dist = pairwise_color_df.copy().sort_index()
    fig = plt.figure(figsize=(8, 6))
    ax = plt.subplot(111)

    dendrogram_df = hierarchy.linkage(pairwise_color_dist, 'ward')
    short_names = [
        f"{x[0][:2]}_{x[2][:2]}_{x[1]}_{x[3]}"
        for x in pairwise_color_dist.index]
    hierarchy.dendrogram(
        dendrogram_df, labels=short_names, ax=ax, orientation='bottom',
        leaf_rotation=90, color_threshold=0)
    for edge in ["left", "right", "top", "bottom"]:
        ax.spines[edge].set_color('none')
    ax.set_xlabel("Clustering method: Ward's D")
    taxon_short_names = set([y[:5] for y in short_names])
    cmap = cm.get_cmap("Dark2", 256)
    color_dict = {
        x: cmap(i+1/(len(taxon_short_names)))
        for i, x in enumerate(taxon_short_names)}
    labels = ax.get_xmajorticklabels()
    for label in labels:
        color = color_dict[label.get_text()[:5]]
        label.set_color(color)
    return fig


def pca_snsplot(color_data, pcomp_a=1, pcomp_b=2):
    import numpy as np
    import matplotlib.pyplot as plt
    from matplotlib import cm
    import seaborn as sns
    import pandas as pd
    from sklearn.decomposition import PCA
    from sklearn.impute import SimpleImputer

    imp = SimpleImputer(missing_values=np.nan, strategy='mean')
    imp = imp.fit(color_data)
    imp_color_data = imp.transform(color_data)
    taxon_assignments = [
        f"{x[0]}_{x[2]}".replace("_", " ") for x in color_data.columns]

    pca = PCA(svd_solver='full').fit(imp_color_data.transpose())

    COLOR_MAP = 'spring'  # 'tab20'

    pc_a = f'PC{pcomp_a}'
    pc_b = f'PC{pcomp_b}'
    try:
        all_markers = np.array(
            range(len(taxon_assignments))) / len(list(set(taxon_assignments)))
    except ZeroDivisionError:
        all_markers = np.array([0])
    given_markers = {group: all_markers[i] for i, group in enumerate(
        list(set(taxon_assignments)))}
    colors = [given_markers[el] for el in taxon_assignments]

    col_dic = {asg: co for asg, co in zip(taxon_assignments, colors)}
    cmap = cm.get_cmap(COLOR_MAP, 256)

    cmap_dic = {asg: cmap(co) for asg, co in zip(taxon_assignments, colors)}
    color_order = [cl for i, cl in enumerate(colors) if cl not in colors[:i]]

    pca_color_data = pca.transform(imp_color_data.transpose())

    pca_plotting_data = pd.DataFrame(
        zip(pca_color_data[:, pcomp_a-1], pca_color_data[:, pcomp_b-1],
            taxon_assignments, color_data.columns.get_level_values(1)),
        columns=[pc_a, pc_b, 'taxon', 'area'],
        index=color_data.columns)

    fig = plt.figure(figsize=(8, 6))
    ax = plt.subplot(1, 1, 1)
    try:
        ax = sns.kdeplot(
            data=pca_plotting_data, x=pc_a, y=pc_b, hue='taxon', palette=cmap_dic,
            ax=ax, fill=True, levels=2, thresh=.1, alpha=.2)
    except:
        pass
    S = 40
    sns.scatterplot(data=pca_plotting_data, x=pc_a, y=pc_b,
                    hue="taxon", alpha=.9, s=S, style="area",
                    edgecolor='#333333', lw=.3, ax=ax)
    ax.spines['left'].set_position('zero')
    ax.spines['left'].set_color('gray')
    ax.spines['bottom'].set_position('zero')
    ax.spines['bottom'].set_color('gray')
    ax.spines['right'].set_color('none')
    ax.spines['top'].set_color('none')
    ax.legend(loc="upper left", bbox_to_anchor=(1, 1))
    return fig


# main

default()