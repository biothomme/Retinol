#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 10 08:28:10 2021

@author: Thomsn
"""
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

def load_data(file_name:str, header=None):
    from importlib import resources
    try:
        with resources.open_text("data", file_name) as fid:
            df = pd.read_csv(fid, header=None, index_col=0)
            if header != None:
                df.columns = header
    except FileNotFoundError:
        print(f"There was no file {file_name}. Try again.")
        return None
    return df

def get_header(name:str):
    FILES_DICT = {
        "bombus": ["uv", "blue", "green"],
        "apis": ["uv", "blue", "green"],
        "d65": ["reflectance"],
        "green_leaf_std": ["reflectance"]
        }
    name = name.lower()
    return FILES_DICT[name] if name in FILES_DICT.keys() else None

def get_file_name(name:str):
    FILES_DICT = {
        "bombus": "bombus_sensoring.csv",
        "apis": "apis_sensoring.csv",
        "d65": "d65_standards.csv",
        "green_leaf_std": "background_green_leaf.csv"
        }
    name = name.lower()
    return FILES_DICT[name] if name in FILES_DICT.keys() else None

def input_flowers():
    from IPython.display import display
    import ipywidgets as widgets
    uploader = widgets.FileUpload(accept='.csv', multiple=False)
    display(uploader)
    return uploader

def parse_flowers(uploader):
    import codecs
    from io import StringIO
    try:
        data = list(uploader.value.values())[0]
    except IndexError:
        print("Please upload a csv file first.")
    csv_data = codecs.decode(data["content"], encoding="utf-8")
    # check if comma or semicolon separated
    newline_count = csv_data.count('\n')
    for seperator in [",", ";", "\t"]:
        seperator_counts = [line.count(seperator) for line in csv_data.split("\n")]
        if (seperator_counts[0] != 0 and len(set(seperator_counts[:-1])) == 1):
            df = pd.read_csv(StringIO(csv_data), sep=seperator, header=None)
            return df

    print("csv file was neither seperated by comma, semicolon nor tabs. Please check and try again.")
    return

def new_floral_spectra(wl_df:pd.DataFrame, meta_df:pd.DataFrame):
    floral_spectra = Floral_Spectra(wl_df,
                                    genus_names=meta_df.iloc[:,0],
                                    species_names=meta_df.iloc[:,1],
                                    area_names=meta_df.iloc[:,2],
                                    additional=meta_df.iloc[:,3])
    return floral_spectra

def get_dropdown_value(key_choice):
    return key_choice.options[key_choice.value][0]

class Floral_Spectra:
    '''This implements floral spectra to convert those to insect vision.
    '''
    def __init__(self, floral_spectra_data, genus_names=None,
                  species_names=None, area_names=None, additional=None):
        self.data = floral_spectra_data.iloc[:, 1:]
        self.data.index = floral_spectra_data.iloc[:, 0]

        df_mid =  pd.MultiIndex.from_arrays(
            [genus_names, area_names, species_names, additional],
            names=("genus", "area", "species", "specimen"))
        self.data.columns = df_mid
        self.normalized = False
        self.converted_to_iv = False
        self.hexagon_df = None
        self.triangle_df = None
        self.pairwise_color_dist = None
        self.make_directory()
        return

    def make_directory(self):
        import tempfile
        import os
        import datetime
        import re
        self.temp = tempfile.TemporaryDirectory(dir=os.getcwd())
        self.directory = re.sub(
            r"[-: \.]", "",
            f"{self.temp.name}/bumble_view_{datetime.datetime.now()}")
        os.makedirs(self.directory)
        return

    def normalize(self):
        if (not self.normalized):
            df = self.data - self.data.apply(min)
            self.data = df / df.apply(max)
            self.normalized = True
        return

    def select_key(self, key="genus", genus_choice=None):
        from IPython.display import display
        import ipywidgets as widgets
        KEY_DICT = {"genus": False, "area": True}
        level_index = 0
        if (KEY_DICT[key] and genus_choice):
            df = self.data[get_dropdown_value(genus_choice)]
        else:
            df = self.data
            level_index += 1 if KEY_DICT[key] else 0
        options = [(k, i) for i, k in enumerate(set(
            df.keys().get_level_values(level_index)))]
        options += [(None, len(options))]
        key_choice = widgets.Dropdown(
            options=options,
            value=0,
            description=f"{key.capitalize()}:")
        display(key_choice)
        return key_choice

    def bombus_vision(self):
        if not self.converted_to_iv:
            self.normalize()
            bombus_df = load_data(get_file_name("bombus"), get_header("bombus"))
            d65_df = load_data(get_file_name("d65"), get_header("d65"))
            green_leaf_std_df = load_data(get_file_name("green_leaf_std"),
                                          get_header("green_leaf_std"))
            df = self.data.loc[self.get_wavelength_index(), :]
            df.index = np.round(df.index)
    
            # build arrays of minimal wl range for computation
            minimal_wl = max(
                [min(d65_df.index), min(bombus_df.index),
                 min(green_leaf_std_df.index), min(df.index)])
            maximal_wl = min(
                [max(d65_df.index), max(bombus_df.index),
                 max(green_leaf_std_df.index), max(df.index)])
            minimal_range_index = bombus_df.loc[minimal_wl:maximal_wl, :].index
            bombus_array = np.asarray(bombus_df.loc[minimal_wl:maximal_wl, :])
            green_leaf_std_array = np.asarray(
                green_leaf_std_df.loc[minimal_wl:maximal_wl, :]).flatten()
            d65_array = np.asarray(d65_df.loc[minimal_wl:maximal_wl, :]).flatten()
            data_array = np.asarray(df.loc[minimal_wl:maximal_wl, :])

            # quantum catch values for background adaptation (with green leaf std)
            # non adapted values for use in triangle spectrum loci
            qc_non_adapted = pd.DataFrame(
                np.apply_along_axis(
                    lambda x: x * d65_array, 0, bombus_array),
                index=minimal_range_index,
                columns=bombus_df.columns)
            qc_general = pd.DataFrame(
                np.apply_along_axis(
                    lambda x: x * green_leaf_std_array, 0,
                    np.asarray(qc_non_adapted)),
                index=minimal_range_index,
                columns=bombus_df.columns)

            # quantum catch values in all different receptors for all samples
            qc_specific = {}
            for i, receptor in enumerate(bombus_df.columns):
                qc_recepetor = pd.DataFrame(
                    np.apply_along_axis(lambda x: x *
                                        bombus_array[:,i]*
                                        d65_array, 0, data_array),
                    index=minimal_range_index,
                    columns=df.columns)
                qc_specific[receptor] = qc_recepetor

            # integral of quantum catch general and 
            # receptor specific sensitivity factors R
            qc_general_integral = np.sum(qc_general, axis=0)
            sensitivity_factors = qc_general_integral.apply(np.reciprocal, 0)

            # integral of quantum catch general,
            # amount of light absorbed by receptor type, called P and
            # the estimated excitation E of the non-linear phototransduction.
            qc_specific_integral = {}
            absorbed_lights_p = {}
            excitations_e = {}
            for i, receptor in enumerate(bombus_df.columns):
                qc_recepetor_integral = np.sum(qc_specific[receptor], axis=0)
                absorbed_light_p = qc_recepetor_integral *\
                    sensitivity_factors[receptor]
                excitation_e = absorbed_light_p.apply(lambda x: x/(x+1))
                qc_specific_integral[receptor] = qc_recepetor_integral
                absorbed_lights_p[receptor] = absorbed_light_p
                excitations_e[receptor] = excitation_e
            self.hexagon_df = pd.DataFrame(excitations_e)

            # relative light absorptions (following chittka u, b, g) and
            relative_absorptions = {
                key: v.divide(pd.DataFrame(absorbed_lights_p).sum(axis=1))
                for key, v in absorbed_lights_p.items()}
            self.triangle_df = pd.DataFrame(relative_absorptions)

            # receptor potential sensitivity called spectrum locus - hexagon
            intermediate_1 = qc_general*sensitivity_factors
            intermediate_2 = intermediate_1.apply(
                lambda x: np.divide(3*x, intermediate_1.sum(axis=1)),
                axis=0)
            spectrum_locus = intermediate_2 / (intermediate_2+1)
            self.hexagon_sl = spectrum_locus

            # receptor potential sensitivity called spectrum locus - triangle
            qc_spectrum_loci_triangle = qc_non_adapted.apply(
                lambda x: np.divide(x, qc_general_integral), axis=1)
            qc_spectrum_loci_triangle_rel = qc_spectrum_loci_triangle.apply(
                lambda x: np.divide(x, qc_spectrum_loci_triangle.sum(axis=1)),
                axis=0)
            self.triangle_sl = qc_spectrum_loci_triangle_rel

            # compute pairwise color distance between all samples
            excitation_signals = Perceived_Signals(self.hexagon_df)
            excitation_signals.get_x()
            excitation_signals.get_y()
            color_distance = excitation_signals.data.apply(
                lambda x: np.sqrt(np.subtract(
                    excitation_signals.data["x"], x["x"])**2 + np.subtract(
                        excitation_signals.data["y"], x["y"])**2), axis=1)
            self.pairwise_color_dist = color_distance

            self.converted_to_iv = True
        return

    def plot_wl_spectra(self, genus, area, p_value_threshold=.05,
                        show_fig=False):
        from plotting import single_plot
        self.normalize()
        if (genus.lower() in [k.lower() for k in self.data.keys().levels[0]]):
            if area is None:
                fig = single_plot(
                    self.data[genus].swaplevel(0, 1, axis=1),
                    p_value_threshold)
            else:
                fig = single_plot(self.data[genus][area], p_value_threshold)
            checkmake_dir_existence(f"{self.directory}/wl_spectra")
            fig.savefig(f"{self.directory}/wl_spectra/wl_spectra_{genus}_{area if area != None else 'all'}.pdf")
            if show_fig:
                return fig
            plt.close(fig)
            return
        print("Did not find given genus. Please check and try again.")
        return

    def plot_hexagon(self, genus=None, area=None,
                     axis_label=True, spectrum_loci_annotations=True,
                     show_fig=False):
        from plotting import polygon_plot
        self.bombus_vision()
        plotting_hex_df = self.subset_plotting_frame(
            self.hexagon_df, genus=genus, area=area)
        plotting_hex_sl = self.subset_plotting_frame(
            self.hexagon_sl, genus=genus, area=area)
        fig = polygon_plot(
            plotting_hex_df, plotting_hex_sl, axis_label=axis_label,
            spectrum_loci_annotations=spectrum_loci_annotations)
        checkmake_dir_existence(f"{self.directory}/hexagon")
        fig.savefig(f"{self.directory}/hexagon/ins_vis_hex_{genus}_{area if area != None else 'all'}.pdf")
        if show_fig:
            return fig
        plt.close(fig)
        return

    def plot_triangle(self, genus=None, area=None,
                      axis_label=True, spectrum_loci_annotations=True,
                      show_fig=False):
        from plotting import polygon_plot
        self.bombus_vision()
        plotting_tri_df = self.subset_plotting_frame(
            self.triangle_df, genus=genus, area=area)
        plotting_tri_sl = self.subset_plotting_frame(
            self.triangle_sl, genus=genus, area=area)
        fig = polygon_plot(
            plotting_tri_df, plotting_tri_sl, plot_type="triangle",
            axis_label=axis_label,
            spectrum_loci_annotations=spectrum_loci_annotations)
        checkmake_dir_existence(f"{self.directory}/triangle")
        fig.savefig(
            f"{self.directory}/triangle/ins_vis_tri_{genus}_{area if area != None else 'all'}.pdf")
        if show_fig:
            return fig
        plt.close(fig)
        return

    def plot_pca(self, genus=None, area=None, pc_a=1, pc_b=2,
                 data_type="physical", show_fig=False):
        from plotting import pca_snsplot
        self.bombus_vision()
        if data_type == "physical":
            df = self.data
            axis = 1
        elif data_type == "insect_vision":
            df = self.triangle_df.transpose().copy()
            axis = 1
        else:
            print(f"""It is not possible to plot the data of type {data_type}.
                  Try to set 'data_type' either to 'physical' or 'insect_vision'.""")
            return
        if (genus.lower() in [k.lower() for k in self.data.keys().levels[0]]):
            if area is None:
                fig = pca_snsplot(self.subset_plotting_frame(
                    df, genus=genus, axis=axis), pcomp_a=pc_a, pcomp_b=pc_b)
            else:
                fig = pca_snsplot(self.subset_plotting_frame(
                    df, genus=genus, area=area, axis=axis), pcomp_a=pc_a,
                    pcomp_b=pc_b)
            checkmake_dir_existence(f"{self.directory}/pca_{data_type}")
            fig.savefig(f"{self.directory}/pca_{data_type}/pca_{data_type}_{genus}_{area if area != None else 'all'}.pdf")
            if show_fig:
                return fig
            plt.close(fig)
            return
        print("Did not find given genus. Please check and try again.")
        return

    def plot_distances(self, genus=None, area=None, plot_type="heatmap",
                       show_fig=False):
        from plotting import distance_dendrogram
        from plotting import distance_heatmap
        PLOT_TYPE_DICT = {
            "dendrogram": distance_dendrogram,
            "heatmap": distance_heatmap
            }
        self.bombus_vision()
        if plot_type in PLOT_TYPE_DICT.keys():
            if (genus.lower() in [
                    k.lower() for k in self.data.keys().levels[0]]):
                if area is None:
                    df = self.subset_plotting_frame(
                        self.pairwise_color_dist, genus=genus)
                    df = self.subset_plotting_frame(df.transpose(), genus=genus)
                else:
                    df = self.subset_plotting_frame(
                        self.pairwise_color_dist, genus=genus, area=area)
                    df = self.subset_plotting_frame(
                        df.transpose(), genus=genus, area=area)
                fig = PLOT_TYPE_DICT[plot_type](df)
                checkmake_dir_existence(
                    f"{self.directory}/color_dist_{plot_type}")
                fig.savefig(f"{self.directory}/color_dist_{plot_type}/cd_{plot_type}_{genus}_{area if area != None else 'all'}.pdf")
                if show_fig:
                    return fig
                plt.close(fig)
                return
            print("Did not find given genus. Please check and try again.")
        print(f"""There is no such plot type as {plot_type} available. Try to
               set 'data_type' either to 'heatmap' or 'dendrogram'.""")
        return

    def subset_plotting_frame(self, df, genus=None, area=None, axis=0):
        plotting_frame = df
        if axis == 0:
            genus_present = str(genus) in df.index.get_level_values(0)
        else:
            genus_present = str(genus) in df.columns.get_level_values(0)
        if genus_present:
            plotting_frame = df.xs(
                    genus, level="genus", axis=axis, drop_level=False)
            if axis == 0:
                area_present = str(
                    area) in plotting_frame.index.get_level_values(1)
            else:
                area_present = str(
                    area) in plotting_frame.columns.get_level_values(1)
            if area_present:
                plotting_frame = plotting_frame.xs(
                        area, level="area", axis=axis, drop_level=False)
        return plotting_frame

    def get_wavelength_index(self):
        wavelength_index = []
        for wavelength in range(300, round(max(self.data.index)),5):
            min_wavelength = min(abs(self.data.index-wavelength))
            if wavelength+min_wavelength in self.data.index:
                wavelength_index.append(wavelength+min_wavelength)
            elif wavelength-min_wavelength in self.data.index:
                wavelength_index.append(wavelength-min_wavelength)
            else:
                print("There was a missing wavelength. Please check.")
        return wavelength_index

    def save_data(self, data_file):
        FILE_DICT = {"wl_spectra": (self.data, "wl_spectra_normalized.csv"),
                     "hexagon": (
                         self.hexagon_df, "ins_vis_hex_excitations.csv"),
                     "triangle": (
                         self.triangle_df, "ins_vis_tri_rel_absorptions.csv"),
                     "pca_physical": (self.data, "wl_spectra_normalized.csv"),
                     "pca_insect_vision": (
                         self.hexagon_df, "ins_vis_hex_excitations.csv"),
                     "heatmap": (self.pairwise_color_dist,
                         "ins_vis_pairwise_color_dist.csv"),
                     "dendrogram": (self.pairwise_color_dist,
                         "ins_vis_pairwise_color_dist.csv")}
        if data_file in ["wl_spectra", "pca_physical", "pca_insect_vision"]:
            FILE_DICT[data_file][0].to_csv(
                f"{self.directory}/{data_file}/{FILE_DICT[data_file][1]}")
        elif data_file in ["heatmap", "dendrogram"]:
            FILE_DICT[data_file][0].to_csv(
                f"{self.directory}/color_dist_{data_file}/{FILE_DICT[data_file][1]}")
        elif data_file in FILE_DICT.keys():
            df_signal = Perceived_Signals(FILE_DICT[data_file][0])
            df_signal.get_x()
            df_signal.get_y()
            df_signal.data.to_csv(
                f"{self.directory}/{data_file}/{FILE_DICT[data_file][1]}")
        else:
            print("There is no such data to save: {data_file}")
        return

    def plot_all_inclusive(self, plot_type="wl_spectra"):
        for genus in set(self.data.columns.get_level_values(0)):
            areas = list(
                set(self.data[genus].columns.get_level_values(0))) + [None]
            for area in areas:
                if plot_type == "wl_spectra":
                    self.plot_wl_spectra(genus, area)
                elif plot_type == "hexagon":
                    self.plot_hexagon(genus=genus, area=area)
                elif plot_type == "triangle":
                    self.plot_triangle(genus=genus, area=area)
                elif plot_type == "pca_physical":
                    self.plot_pca(genus=genus, area=area)
                elif plot_type == "pca_insect_vision":
                    self.plot_pca(
                        genus=genus, area=area, data_type="insect_vision")
                elif plot_type == "heatmap":
                    self.plot_distances(genus=genus, area=area)
                elif plot_type == "dendrogram":
                    self.plot_distances(
                        genus=genus, area=area, plot_type="dendrogram")
                else:
                    print(f"Plotting type {plot_type} was not recognized.")
                    return
        self.save_data(plot_type)
        return

    def download_data(self):
        import shutil
        from IPython.display import FileLink
        from IPython.display import display
        output_zip = self.directory.split("/")[-1]
        shutil.make_archive(output_zip, "zip", self.temp.name)
        display(FileLink(f"{output_zip}.zip"))

    def close_temporary_dir(self):
        self.temp.cleanup()


def checkmake_dir_existence(directory):
    import os
    if not os.path.exists(directory):
        os.makedirs(directory)
    return


class Perceived_Signals:
    '''Visual signals, which are present in a receptor specific dataframe
    can be stored in this class, to simplify its transformation to x, y values.
    In addition this class contains the shapes of hexagon and triangle.'''
    TRIANGLE_HEIGHT = np.sqrt(3/4)
    TRIANGLE_COORDINATES = [[-TRIANGLE_HEIGHT, 0, TRIANGLE_HEIGHT, -TRIANGLE_HEIGHT],
                            [-.5, 1, -.5, -.5]]
    HEXAGON_COORDINATES = [[-TRIANGLE_HEIGHT,
                            -TRIANGLE_HEIGHT,
                            0,
                            TRIANGLE_HEIGHT,
                            TRIANGLE_HEIGHT,
                            0,
                            -TRIANGLE_HEIGHT],
                           [-.5, .5, 1, .5, -.5, -1, -.5]]

    def __init__(self, signals_df):
        self.data = signals_df.copy()
        self.x = False
        self.y = False
        self.taxa = np.array([])
        return

    def get_x(self):
        if not self.x:
            self.data.loc[:, "x"] = (
                self.data.iloc[:, 2]-self.data.iloc[:, 0]) * self.TRIANGLE_HEIGHT
            self.x = True
        return self.data["x"]

    def get_y(self):
        if not self.y:
            self.data.loc[:, "y"] = (
                self.data.iloc[:, 1]) - (
                    self.data.iloc[:, 2]+self.data.iloc[:, 0])/2
            self.y = True
        return self.data["y"]

    def get_taxa(self):
        if self.taxa.size == 0:
            self.taxa = np.asarray([f"{x[0]}_{x[2]}".replace("_", " ")
                                    for x in self.data.index])
        return self.taxa

def reset_directory():
    import os
    import shutil
    temporaries = [x for x in os.listdir() if (".zip" in x) | ("tmp" in x)]
    print(temporaries)
    [shutil.rmtree(tmp)  if os.path.isdir(tmp) else os.remove(tmp)
        for tmp in temporaries]
    return

reset_directory()