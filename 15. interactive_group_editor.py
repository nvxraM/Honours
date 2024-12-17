import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.backend_bases import MouseButton
import matplotlib
import customtkinter as ctk
import tkinter as tk
import sys
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from matplotlib.widgets import RectangleSelector
import copy
import os
import time
import webbrowser
from Bio import SeqIO
from geopy.geocoders import Nominatim
from geopy.distance import geodesic
import logging
import shutil
from tkinter import messagebox

# Set up logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

# Use TkAgg backend for GUI functionality
matplotlib.use('TkAgg')


class InteractiveGroupEditor:
    """
    A GUI application to interactively edit group assignments of MDS data points.
    If no CSV files exist, tries to create them from MDS files.
    After editing, copies all files (CSV, PNG, MDS) to:
    sequences/Interactive_group_editor/[Genus]/[Species]/.mds
    """

    def __init__(self, base_directory):
        self.base_directory = base_directory
        # Attempt to get CSV files; if none exist, create them from MDS files
        self.csv_files = self._get_or_create_csv_files()

        self.csv_file = self.csv_files[0][1][0] if self.csv_files else None
        if not self.csv_file:
            raise FileNotFoundError("No CSV files found or created in the specified directory.")

        # Load CSV
        try:
            self.df = pd.read_csv(self.csv_file)
            self._ensure_columns_exist()
        except (FileNotFoundError, pd.errors.EmptyDataError) as e:
            raise ValueError(f"Error reading the CSV file: {e}")

        self.df_original = copy.deepcopy(self.df)

        # Configure customtkinter appearance
        ctk.set_appearance_mode("System")
        ctk.set_default_color_theme("blue")

        # Create main application window
        self.root = ctk.CTk()
        self.root.title("Interactive Group Editor")
        self.root.state('zoomed')

        # Set up matplotlib figure
        self.fig, self.ax = plt.subplots(figsize=(12, 10))
        self.selected_points = []
        self.highlight_scatter = None

        # Define group colors
        self.group_colors = {
            1: '#e6194b', 2: '#3cb44b', 3: '#ffe119',
            4: '#4363d8', 5: '#f58231', 6: '#911eb4',
            7: '#42d4f4'
        }

        # Species file map for updating dropdown menus
        self.species_file_map = {}

        # Load data, set up interactions, add UI elements
        self._load_data()
        self._connect_events()
        self._add_command_buttons()

        # Handle window close event
        self.root.protocol("WM_DELETE_WINDOW", self._on_close)

        # Rectangle selector for selecting points by drag
        self.rectangle_selector = RectangleSelector(
            self.ax, self._on_select,
            interactive=True, useblit=True,
            button={MouseButton.LEFT}, minspanx=5, minspany=5, spancoords='pixels'
        )

        # Geolocator for location lookups
        self.geolocator = Nominatim(user_agent="interactive_group_editor")

        # Basic oceans list, add more if needed
        self.oceans = [
            ("Arctic Ocean", "Arctic Ocean", (82.5, -45)),
            ("Atlantic Ocean", "Atlantic Ocean", (0, -30)),
            ("Indian Ocean", "Indian Ocean", (-20, 80)),
            ("Southern Ocean", "Southern Ocean", (-65, 140)),
            ("Pacific Ocean", "Pacific Ocean", (0, -140))
        ]

    def _get_or_create_csv_files(self):
        csv_files = self._get_csv_files()
        if not csv_files:
            # No CSV found, try to create from MDS files
            self._create_csv_from_mds()
            csv_files = self._get_csv_files()  # Check again after creation
        return csv_files

    def _get_csv_files(self):
        """
        Find all '_clustered.csv' files under the base_directory.
        Group them by genus (folder name) to populate dropdown menus.
        """
        csv_files = {}
        for root, dirs, files in os.walk(self.base_directory):
            for file in files:
                if file.endswith("_clustered.csv"):
                    family = os.path.basename(root)
                    if family not in csv_files:
                        csv_files[family] = []
                    csv_files[family].append(os.path.join(root, file))
        return sorted(csv_files.items())

    def _create_csv_from_mds(self):
        """
        If no CSV files exist, try to create them from MDS files.
        MDS files should be named '[Species]_mds.mds' in:
        sequences/CDS_Genus/[Genus]/CDS_nucleotide_gapped/MDS/[Species]/[Species]_mds.mds
        """
        for root, dirs, files in os.walk(self.base_directory):
            for file in files:
                if file.endswith("_mds.mds"):
                    mds_path = os.path.join(root, file)
                    species_name = file.replace("_mds.mds", "")
                    csv_name = f"{species_name}_clustered.csv"
                    csv_path = os.path.join(root, csv_name)

                    if not os.path.exists(csv_path):
                        try:
                            df_mds = pd.read_csv(mds_path, delim_whitespace=True)
                            if not {'FID', 'IID', 'C1', 'C2'}.issubset(df_mds.columns):
                                logging.error(f"MDS file {mds_path} missing required columns (FID,IID,C1,C2).")
                                continue

                            df_mds['Group'] = 1
                            df_mds['lat_lon'] = 'N/A'
                            df_mds['location'] = 'N/A'
                            df_mds['Location Source'] = 'N/A'

                            df_mds = df_mds[['FID', 'IID', 'C1', 'C2', 'Group', 'lat_lon', 'location', 'Location Source']]

                            df_mds.to_csv(csv_path, index=False)
                            logging.info(f"Created CSV from MDS file: {csv_path}")
                        except Exception as e:
                            logging.error(f"Error creating CSV from MDS {mds_path}: {e}")

    def _ensure_columns_exist(self):
        required_columns = ['lat_lon', 'location', 'Location Source']
        for col in required_columns:
            if col not in self.df.columns:
                self.df[col] = 'N/A'

    def _load_data(self):
        """
        Load data from self.df into the scatter plot.
        """
        self.x = self.df['C1'].values
        self.y = self.df['C2'].values
        self.groups = self.df['Group'].values
        self.iids = self.df['IID'].values

        colors = [self.group_colors.get(g, 'gray') for g in self.groups]
        if hasattr(self, 'scatters'):
            self.scatters.set_offsets(list(zip(self.x, self.y)))
            self.scatters.set_facecolors(colors)
        else:
            self.scatters = self.ax.scatter(self.x, self.y, c=colors, alpha=0.8, s=50, picker=True)

        self.ax.set_xlabel('C1', fontsize=12, labelpad=10)
        self.ax.set_ylabel('C2', fontsize=12, labelpad=10)
        species_name = os.path.basename(self.csv_file).replace('_clustered.csv', '')
        self.ax.set_title(f'{species_name} - 2D Scatter Plot', fontsize=18, fontweight='bold', pad=20)
        self.ax.set_xlim(self.x.min() - 0.05, self.x.max() + 0.05)
        self.ax.set_ylim(self.y.min() - 0.05, self.y.max() + 0.05)
        self.ax.grid(True, linestyle='--', alpha=0.5)
        self._add_legend()

    def _connect_events(self):
        self.cid_click = self.fig.canvas.mpl_connect('button_press_event', self._on_click)
        self.cid_pick = self.fig.canvas.mpl_connect('pick_event', self._on_pick)

    def _add_command_buttons(self):
        button_frame = ctk.CTkScrollableFrame(self.root)
        button_frame.pack(side='right', fill='both', padx=20, pady=20, expand=True)

        ctk.CTkLabel(button_frame, text="Select Family and Species", font=('Helvetica', 18, 'bold'), anchor='center').pack(pady=(0, 15), fill='x')

        self.selected_family = ctk.StringVar(self.root)
        self.selected_species = ctk.StringVar(self.root)

        self.family_menu = ctk.CTkOptionMenu(button_frame, variable=self.selected_family,
                                             values=[family for family, _ in self.csv_files],
                                             command=self._update_species_menu)
        self.family_menu.pack(pady=5, fill='x', ipady=5, ipadx=5)

        self.species_menu = ctk.CTkOptionMenu(button_frame, variable=self.selected_species, values=[""],
                                              command=lambda val: self._change_csv_file(val))
        self.species_menu.pack(pady=5, fill='x', ipady=5, ipadx=5)

        if self.csv_files:
            self.selected_family.set(self.csv_files[0][0])
            self._update_species_menu(self.csv_files[0][0])

        ctk.CTkLabel(button_frame, text="Edit Functions", font=('Helvetica', 18, 'bold'), anchor='center').pack(pady=(20, 15), fill='x')

        ctk.CTkButton(button_frame, text="Deselect All Points", command=self._deselect_all_points, width=200, height=40, font=('Helvetica', 16)).pack(pady=5, fill='x')
        ctk.CTkButton(button_frame, text="Delete Selected Points", command=self._delete_selected_points, width=200, height=40, font=('Helvetica', 16)).pack(pady=5, fill='x')
        ctk.CTkButton(button_frame, text="Reassign Group", command=self._reassign_group, width=200, height=40, font=('Helvetica', 16)).pack(pady=5, fill='x')
        ctk.CTkButton(button_frame, text="Undo Changes", command=self._undo_changes, width=200, height=40, font=('Helvetica', 16)).pack(pady=5, fill='x')
        ctk.CTkButton(button_frame, text="Show Selected Point Values", command=self._show_selected_point_values, width=200, height=40, font=('Helvetica', 16)).pack(pady=5, fill='x')
        ctk.CTkButton(button_frame, text="Toggle All Point Values", command=self._toggle_all_point_values, width=200, height=40, font=('Helvetica', 16)).pack(pady=5, fill='x')
        ctk.CTkButton(button_frame, text="Search Selected Points on NCBI", command=self._search_selected_points_ncbi, width=200, height=40, font=('Helvetica', 16)).pack(pady=5, fill='x')

        ctk.CTkLabel(button_frame, text="Assign Groups", font=('Helvetica', 18, 'bold'), anchor='center').pack(pady=(20, 15), fill='x')

        for group_num in range(1, 8):
            ctk.CTkButton(button_frame, text=f"Group {group_num}",
                          command=lambda g=group_num: self._assign_to_group(g),
                          fg_color=self.group_colors.get(group_num, '#e0e0e0'),
                          width=200, height=40, font=('Helvetica', 16)).pack(pady=5, fill='x')

        ctk.CTkButton(button_frame, text="Save CSV and PNG", command=self._save_csv_and_png, width=200, height=40,
                      font=('Helvetica', 16, 'bold'), fg_color='#34a853').pack(pady=20, fill='x')

        canvas_plot = FigureCanvasTkAgg(self.fig, master=self.root)
        canvas_plot.get_tk_widget().pack(side='left', fill='both', expand=True, padx=20, pady=20)
        canvas_plot.draw()

    def _update_species_menu(self, selected_family):
        species_files = self.csv_files[[family for family, _ in self.csv_files].index(selected_family)][1]
        self.species_file_map = {os.path.basename(species): species for species in species_files}
        self.species_menu.configure(values=list(self.species_file_map.keys()))
        if species_files:
            first_species = os.path.basename(species_files[0])
            self.selected_species.set(first_species)
            self._change_csv_file(first_species)

    def _add_legend(self):
        unique_groups = sorted(self.df['Group'].unique())
        group_counts = self.df['Group'].value_counts()
        handles = []
        for group in unique_groups:
            count = group_counts.get(group, 0)
            label = f'Group {group} (n={count})'
            line = plt.Line2D([0], [0], marker='o', color='w',
                              markerfacecolor=self.group_colors.get(group, 'gray'),
                              markersize=10, label=label)
            handles.append(line)
        self.legend = self.ax.legend(handles=handles, title='Groups',
                                     loc='best', framealpha=0.6)

    def _on_click(self, event):
        if event.button == MouseButton.LEFT and event.dblclick:
            if event.inaxes != self.ax:
                return
            clicked_point = (event.xdata, event.ydata)
            distances = ((self.x - clicked_point[0])**2 + (self.y - clicked_point[1])**2)**0.5
            closest_index = distances.argmin()
            selected_group = self.groups[closest_index]
            self.selected_points = [i for i, g in enumerate(self.groups) if g == selected_group]
            self._update_selected_points()
            if hasattr(self, 'rectangle_selector'):
                self.rectangle_selector.set_active(False)
                self.rectangle_selector.set_active(True)

    def _on_pick(self, event):
        ind = event.ind
        if ind is not None:
            self.selected_points.extend(ind)
            self.selected_points = list(set(self.selected_points))
            self._update_selected_points()

    def _on_select(self, eclick, erelease):
        x_min, x_max = sorted([eclick.xdata, erelease.xdata])
        y_min, y_max = sorted([eclick.ydata, erelease.ydata])
        self.selected_points = [i for i in range(len(self.df))
                                if x_min <= self.df['C1'].iloc[i] <= x_max and y_min <= self.df['C2'].iloc[i] <= y_max]
        self._update_selected_points()

    def _deselect_all_points(self):
        self.selected_points.clear()
        self._update_selected_points()
        if hasattr(self, 'rectangle_selector'):
            self.rectangle_selector.set_active(False)
            self.rectangle_selector.set_active(True)

    def _delete_selected_points(self):
        if self.selected_points:
            self.df.drop(self.df.index[self.selected_points], inplace=True)
            self.df.reset_index(drop=True, inplace=True)
            self._load_data()
            self.fig.canvas.draw()

    def _reassign_group(self):
        if self.selected_points:
            new_group = self._prompt_for_group()
            if new_group is not None:
                for i in self.selected_points:
                    if i < len(self.df):
                        self.df.at[self.df.index[i], 'Group'] = new_group
                self._update_scatter()
                self.selected_points.clear()
                self._update_selected_points()
                self.legend.remove()
                self._add_legend()

    def _assign_to_group(self, group):
        if self.selected_points:
            for i in self.selected_points:
                if i < len(self.df):
                    self.df.at[self.df.index[i], 'Group'] = group
            self._update_scatter()
            self.selected_points.clear()
            self._update_selected_points()
            self.legend.remove()
            self._add_legend()

    def _undo_changes(self):
        self.df = copy.deepcopy(self.df_original)
        self._load_data()
        self.fig.canvas.draw()

    def _show_selected_point_values(self):
        if not self.selected_points:
            messagebox.showwarning("Warning", "Please select at least one point to view its values.")
            return

        loading_label = ctk.CTkLabel(self.root, text="Loading, please wait...", font=('Helvetica', 16, 'bold'))
        loading_label.pack(pady=10)
        self.root.update_idletasks()

        values = []
        for idx in self.selected_points:
            if idx < len(self.df):
                c1_value = self.df['C1'].iloc[idx]
                c2_value = self.df['C2'].iloc[idx]
                iid_value = self.df['IID'].iloc[idx]
                lat_lon_value = self.df['lat_lon'].iloc[idx]
                location_value = self.df['location'].iloc[idx]
                values.append(f"IID: {iid_value} | ({c1_value:.2f}, {c2_value:.2f}) | lat_lon={lat_lon_value} | location={location_value}")

        loading_label.destroy()
        value_window = tk.Toplevel(self.root)
        value_window.title("Selected Point Values")
        value_window.geometry('1200x800')
        text_box = tk.Text(value_window, wrap='word', font=('Helvetica', 12))
        text_box.insert('1.0', "\n".join(values))
        text_box.config(state='normal')
        text_box.see('1.0')
        text_box.pack(expand=True, fill='both')
        ok_button = ctk.CTkButton(value_window, text="OK", command=value_window.destroy)
        ok_button.pack(pady=5)

    def _get_lat_lon_and_location(self, iid):
        # Stub function: returns N/A. Implement logic as needed.
        return "N/A", "N/A"

    def _parse_lat_lon(self, lat_lon_str):
        if isinstance(lat_lon_str, float):
            return None
        try:
            parts = lat_lon_str.split()
            lat = float(parts[0])
            if parts[1].upper() == 'S':
                lat = -lat
            lon = float(parts[2])
            if parts[3].upper() == 'W':
                lon = -lon
            return (lat, lon)
        except (ValueError, IndexError):
            logging.error(f"Error parsing lat_lon string: {lat_lon_str}")
            return None

    def _get_location_name(self, lat_lon):
        lat_lon_tuple = self._parse_lat_lon(lat_lon)
        if lat_lon_tuple is None:
            return "N/A"
        try:
            time.sleep(1)
            location = self.geolocator.reverse(lat_lon_tuple, language='en')
            if location:
                return location.address
        except Exception as e:
            logging.error(f"Error getting location name for {lat_lon_tuple}: {e}")
            return "N/A"
        return "N/A"

    def _get_closest_ocean(self, lat_lon):
        lat_lon_tuple = self._parse_lat_lon(lat_lon)
        if lat_lon_tuple is None:
            return ("N/A", "N/A")
        closest_ocean = ("N/A", "N/A")
        min_distance = float('inf')
        for ocean_name, major_ocean, ocean_coords in self.oceans:
            distance = geodesic(lat_lon_tuple, ocean_coords).kilometers
            if distance < min_distance:
                min_distance = distance
                closest_ocean = (major_ocean, ocean_name)
        return closest_ocean

    def _get_combined_location(self, lat_lon):
        if lat_lon == "N/A":
            return "N/A"
        location_name = self._get_location_name(lat_lon)
        major_ocean, specific_ocean = self._get_closest_ocean(lat_lon)
        return f"{location_name} ({major_ocean}, {specific_ocean})"

    def _toggle_all_point_values(self):
        if hasattr(self, 'all_point_annotations') and self.all_point_annotations:
            for annotation in self.all_point_annotations:
                annotation.remove()
            self.all_point_annotations = []
        else:
            self.all_point_annotations = []
            for i in range(len(self.df)):
                lat_lon_value = self.df['lat_lon'].iloc[i]
                location_value = self.df['location'].iloc[i]
                annotation = self.ax.annotate(
                    f"IID: {self.df['IID'].iloc[i]} | ({self.df['C1'].iloc[i]:.2f}, {self.df['C2'].iloc[i]:.2f}) | lat_lon={lat_lon_value} | location={location_value}",
                    (self.df['C1'].iloc[i], self.df['C2'].iloc[i]),
                    textcoords="offset points", xytext=(5, 5), ha='center', fontsize=10,
                    color='#333333')
                self.all_point_annotations.append(annotation)
        self.fig.canvas.draw()

    def _prompt_for_group(self):
        new_group = ctk.CTkInputDialog(title="Input", text="Enter new group for selected points:")
        try:
            return int(new_group.get_input())
        except ValueError:
            return None

    def _update_scatter(self):
        colors = [self.group_colors.get(g, 'gray') for g in self.df['Group'].values]
        self.scatters.set_facecolors(colors)
        self.fig.canvas.draw()

    def _update_selected_points(self):
        if self.highlight_scatter:
            self.highlight_scatter.remove()
            self.highlight_scatter = None
        if self.selected_points:
            selected_x = [self.df['C1'].iloc[i] for i in self.selected_points if i < len(self.df)]
            selected_y = [self.df['C2'].iloc[i] for i in self.selected_points if i < len(self.df)]
            self.highlight_scatter = self.ax.scatter(selected_x, selected_y, facecolors='yellow', edgecolors='black',
                                                     s=100, linewidths=2)
        self.fig.canvas.draw()

    def _change_csv_file(self, selected_file):
        if not selected_file:
            return
        full_path = self.species_file_map.get(selected_file)
        if not full_path:
            messagebox.showerror("Error", f"File {selected_file} not found.")
            return
        self.csv_file = full_path
        try:
            self.df = pd.read_csv(self.csv_file)
            self._ensure_columns_exist()
        except (FileNotFoundError, pd.errors.EmptyDataError) as e:
            raise ValueError(f"Error reading the CSV file: {e}")
        self.df_original = copy.deepcopy(self.df)
        self._load_data()
        self.fig.canvas.draw()

    def _save_csv_and_png(self):
        """
        Save CSV and PNG, then copy CSV, PNG, and MDS to:
        sequences/Interactive_group_editor/[Genus]/[Species]/.mds
        """
        self._deselect_all_points()

        if hasattr(self, 'rectangle_selector'):
            self.rectangle_selector.set_active(False)

        time.sleep(1)

        if hasattr(self, 'rectangle_selector'):
            self.rectangle_selector.set_visible(False)
            self.fig.canvas.draw()

        def update_lat_lon_and_location(row):
            if pd.isna(row['lat_lon']) or row['lat_lon'] == "N/A":
                lat_lon, location = self._get_lat_lon_and_location(row['IID'])
                row['lat_lon'] = lat_lon
                if lat_lon != "N/A":
                    row['location'] = self._get_combined_location(lat_lon)
                elif location != "N/A":
                    row['location'] = location
            return row

        self.df = self.df.apply(update_lat_lon_and_location, axis=1)

        self.df['Location Source'] = self.df.apply(
            lambda r: r['Location Source'] if pd.notna(r['Location Source']) and r['Location Source'] not in ["", "N/A"]
            else ("gb file" if r['lat_lon'] != "N/A" else r['Location Source']),
            axis=1
        )

        csv_path = self.csv_file
        species_name = os.path.basename(csv_path).replace("_clustered.csv", "")

        # Create PNG
        fig_save, ax_save = plt.subplots(figsize=(12, 10))
        ax_save.scatter(self.df['C1'], self.df['C2'],
                        c=[self.group_colors.get(g, 'gray') for g in self.df['Group']],
                        alpha=0.8, s=50)
        ax_save.set_xlabel('C1', fontsize=12, labelpad=10)
        ax_save.set_ylabel('C2', fontsize=12, labelpad=10)
        ax_save.set_title(f'{species_name}', fontsize=25, fontweight='bold', pad=20)
        ax_save.set_xlim(self.df['C1'].min() - 0.05, self.df['C1'].max() + 0.05)
        ax_save.set_ylim(self.df['C2'].min() - 0.05, self.df['C2'].max() + 0.05)
        ax_save.grid(True, linestyle='--', alpha=0.5)

        unique_groups = sorted(self.df['Group'].unique())
        group_counts = self.df['Group'].value_counts()
        handles = []
        for g in unique_groups:
            count = group_counts.get(g, 0)
            label = f'Group {g} (n={count})'
            line = plt.Line2D([0], [0], marker='o', color='w',
                              markerfacecolor=self.group_colors.get(g, 'gray'),
                              markersize=10, label=label)
            handles.append(line)
        ax_save.legend(handles=handles, title='Groups', loc='best', framealpha=0.6)

        png_path = csv_path.replace("_clustered.csv", "_2d_scatterplot.png")
        fig_save.savefig(png_path, bbox_inches='tight', facecolor=fig_save.get_facecolor(), transparent=True)
        plt.close(fig_save)

        self.df.to_csv(csv_path, index=False)

        # Determine genus and species from file path
        relative_path = os.path.relpath(csv_path, self.base_directory)
        parts = relative_path.split(os.sep)
        # parts[0] = Genus
        genus = parts[0]

        # MDS file: [Species]_mds.mds
        mds_file_name = f"{species_name}_mds.mds"
        mds_file_path = os.path.join(os.path.dirname(csv_path), mds_file_name)

        # New directory: sequences/Interactive_group_editor/[Genus]/[Species]/.mds
        interactive_dir = os.path.join("sequences", "Interactive_group_editor", genus, species_name, ".mds")
        if not os.path.exists(interactive_dir):
            os.makedirs(interactive_dir)

        # Copy files
        dest_csv_path = os.path.join(interactive_dir, os.path.basename(csv_path))
        dest_png_path = os.path.join(interactive_dir, os.path.basename(png_path))
        shutil.copy2(csv_path, dest_csv_path)
        shutil.copy2(png_path, dest_png_path)
        if os.path.exists(mds_file_path):
            dest_mds_path = os.path.join(interactive_dir, os.path.basename(mds_file_path))
            shutil.copy2(mds_file_path, dest_mds_path)

        if hasattr(self, 'rectangle_selector'):
            self.rectangle_selector.set_active(True)
            self.rectangle_selector.set_visible(True)

        messagebox.showinfo("Save Successful", f"CSV, PNG, and MDS (if available) have been saved to:\n{interactive_dir}")

    def _open_ncbi_search(self, iids):
        modified_iids = [f"NC_{iid}" if len(str(iid)) == 6 else iid for iid in iids]
        search_query = ' OR '.join([f"{iid}[All Fields]" for iid in modified_iids])
        url = f"https://www.ncbi.nlm.nih.gov/pmc?term={search_query}"
        try:
            webbrowser.open(url, new=2)
        except Exception as e:
            messagebox.showerror("Error", f"Failed to open browser: {e}")

    def _search_selected_points_ncbi(self):
        if self.selected_points:
            iids = [self.df['IID'].iloc[idx] for idx in self.selected_points if idx < len(self.df)]
            self._open_ncbi_search(iids)
        else:
            messagebox.showwarning("Warning", "Please select at least one point to search on NCBI.")

    def _on_close(self):
        if messagebox.askyesno("Exit", "Do you want to save changes before exiting?"):
            self._save_csv_and_png()
        self.root.destroy()
        sys.exit()

    def show(self):
        self.root.mainloop()


# Example usage:
if __name__ == '__main__':
    # Adjust the base_directory as needed
    base_directory = 'sequences/CDS_Genus'
    editor = InteractiveGroupEditor(base_directory)
    editor.show()
