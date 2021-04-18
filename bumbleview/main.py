#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar 12 08:12:30 2021

@author: Thomsn
"""
import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

os.chdir("/Users/Thomsn/Desktop/island_in_the_sun/python/2021_03_bumbleview/bumbleview")
from bumbleview import bumblecore as bc

wl_df = pd.read_csv("data/xmpl_data.csv", sep=",", header=None)
meta_df = pd.read_csv("data/xmpl_meta.csv", sep=",", header=None)

# build new object
flowers = bc.new_floral_spectra(wl_df, meta_df, colab=True)
flowers.erg = bc.load_data_colab(bc.get_file_name("lucilia"), bc.get_header("bombus"))
flowers.bombus_vision()

# fig1 = flowers.plot_triangle(genus="Rhododendron")
# fig2 = flowers.plot_hexagon(genus="Gentiana", area="ventmed")
fig2 = flowers.plot_tetrachromate(genus="Gentiana", area="ventmed", show_fig=True)

# flowers.bombus_vision()
flowers.plot_pca(genus="Rhododendron", area="ventr", pc_a=1, pc_b=2,
                 data_type="insect_vision", show_fig=True)

flowers.pairwise_color_dist
flowers.data

flowers.triangle_df.transpose().copy()

