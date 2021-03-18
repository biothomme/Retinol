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

os.chdir("/Users/Thomsn/Desktop/island_in_the_sun/python/2021_03_bumbleview/bumbleview/bumbleview")
import convwale
from convwale import Perceived_Signals

wl_df = pd.read_csv("data/xmpl_data.csv", header=None)
meta_df = pd.read_csv("data/xmpl_meta.csv", header=None)

# build new object
flowers = convwale.new_floral_spectra(wl_df, meta_df)

fig1 = flowers.plot_triangle(genus="Rhododendron")
fig2 = flowers.plot_hexagon(genus="Gentiana", area="ventmed")

# flowers.bombus_vision()
flowers.plot_pca(genus="Rhododendron", area="ventr", pc_a=1, pc_b=2,
                 data_type="insect_vision", show_fig=True)

flowers.pairwise_color_dist
flowers.data

flowers.triangle_df.transpose().copy()
