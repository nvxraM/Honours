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
    If no CSV files exist, tries to create them from MDS files and directly places them in:
    sequences/Interactive_group_editor/[Genus]/[Species].
    After editing, copies all files (CSV, PNG, MDS, FASTA) to:
    sequences/Interactive_group_editor/[Genus]/[Species].
    """

    def __init__(self, base_directory):
        self.base_directory = base_directory
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

        # Parse genus and species from current CSV file path
        self.genus, self.species_name = self._parse_genus_species_from_path(self.csv_file)

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

        self.species_file_map = {}

        # Extended oceans and seas list
        self.oceans = [
            ("Arctic Ocean", "Arctic Ocean", (82.5, -45)),
            ("Atlantic Ocean", "Atlantic Ocean", (0, -30)),
            ("Indian Ocean", "Indian Ocean", (-20, 80)),
            ("Southern Ocean", "Southern Ocean", (-65, 140)),
            ("Pacific Ocean", "Pacific Ocean", (0, -140)),

            # Atlantic Ocean Regions
            ("Mediterranean Sea", "Atlantic Ocean", (35, 18)),
            ("Caribbean Sea", "Atlantic Ocean", (15, -75)),
            ("Gulf of Mexico", "Atlantic Ocean", (25, -90)),
            ("Baltic Sea", "Atlantic Ocean", (58, 20)),
            ("North Sea", "Atlantic Ocean", (56, 3)),
            ("Norwegian Sea", "Atlantic Ocean", (68, 5)),
            ("Sargasso Sea", "Atlantic Ocean", (30, -60)),
            ("Labrador Sea", "Atlantic Ocean", (58, -55)),
            ("Gulf of Guinea", "Atlantic Ocean", (0, 5)),
            ("English Channel", "Atlantic Ocean", (50, -1)),
            ("Bay of Biscay", "Atlantic Ocean", (45, -5)),
            ("Irish Sea", "Atlantic Ocean", (54, -5)),
            ("Gulf of St. Lawrence", "Atlantic Ocean", (48, -62)),

            # Pacific Ocean Regions
            ("South China Sea", "Pacific Ocean", (15, 115)),
            ("East China Sea", "Pacific Ocean", (28, 125)),
            ("Sea of Japan", "Pacific Ocean", (40, 135)),
            ("Philippine Sea", "Pacific Ocean", (20, 130)),
            ("Coral Sea", "Pacific Ocean", (-18, 152)),
            ("Tasman Sea", "Pacific Ocean", (-40, 160)),
            ("Bering Sea", "Pacific Ocean", (58, -175)),
            ("Gulf of Alaska", "Pacific Ocean", (55, -145)),
            ("Sea of Okhotsk", "Pacific Ocean", (55, 150)),
            ("Yellow Sea", "Pacific Ocean", (35, 123)),
            ("Gulf of California", "Pacific Ocean", (25, -110)),
            ("Arafura Sea", "Pacific Ocean", (-10, 135)),
            ("Timor Sea", "Indian Ocean", (-10, 127)),
            ("Gulf of Thailand", "Pacific Ocean", (10, 100)),
            ("Java Sea", "Pacific Ocean", (-5, 110)),
            ("Andaman Sea", "Indian Ocean", (10, 95)),
            ("Banda Sea", "Pacific Ocean", (-5, 128)),
            ("Celebes Sea", "Pacific Ocean", (5, 125)),
            ("Bismarck Sea", "Pacific Ocean", (-3, 146)),
            ("Solomon Sea", "Pacific Ocean", (-8, 154)),
            ("Gulf of Carpentaria", "Pacific Ocean", (-14, 137)),
            ("Bering Strait", "Pacific Ocean", (66, -169)),
            ("Cook Strait", "Pacific Ocean", (-41, 174)),
            ("Strait of Malacca", "Indian Ocean", (4, 99)),
            ("Gulf of Tonkin", "Pacific Ocean", (20, 108)),

            # Indian Ocean Regions
            ("Red Sea", "Indian Ocean", (20, 38)),
            ("Arabian Sea", "Indian Ocean", (15, 65)),
            ("Bay of Bengal", "Indian Ocean", (15, 90)),
            ("Persian Gulf", "Indian Ocean", (26, 52)),
            ("Gulf of Aden", "Indian Ocean", (12, 48)),
            ("Mozambique Channel", "Indian Ocean", (-17, 42)),
            ("Laccadive Sea", "Indian Ocean", (10, 75)),
            ("Great Australian Bight", "Indian Ocean", (-35, 130)),

            # Arctic Ocean Regions
            ("Barents Sea", "Arctic Ocean", (75, 45)),
            ("Kara Sea", "Arctic Ocean", (75, 70)),
            ("Laptev Sea", "Arctic Ocean", (75, 130)),
            ("East Siberian Sea", "Arctic Ocean", (72, 165)),
            ("Chukchi Sea", "Arctic Ocean", (70, -165)),
            ("Beaufort Sea", "Arctic Ocean", (72, -140)),
            ("Greenland Sea", "Arctic Ocean", (75, -5)),
            ("Lincoln Sea", "Arctic Ocean", (83, -50)),
            ("Hudson Bay", "Arctic Ocean", (60, -85)),

            # Southern Ocean Regions
            ("Weddell Sea", "Southern Ocean", (-73, -45)),
            ("Ross Sea", "Southern Ocean", (-75, 175)),
            ("Amundsen Sea", "Southern Ocean", (-73, -115)),
            ("Bellingshausen Sea", "Southern Ocean", (-70, -90)),
            ("Scotia Sea", "Southern Ocean", (-55, -45)),

            # Additional Seas and Straits
            ("Black Sea", "Atlantic Ocean", (43, 35)),
            ("Caspian Sea", "Atlantic Ocean", (41, 50)),
            ("Azov Sea", "Atlantic Ocean", (46, 36)),
            ("Marmara Sea", "Atlantic Ocean", (40.8, 28)),
            ("Bosporus Strait", "Atlantic Ocean", (41, 29)),
            ("Dardanelles Strait", "Atlantic Ocean", (40, 26)),
            ("Strait of Gibraltar", "Atlantic Ocean", (36, -5)),
            ("Strait of Hormuz", "Indian Ocean", (26, 56)),
            ("Davis Strait", "Atlantic Ocean", (66, -58)),
            ("Skagerrak", "Atlantic Ocean", (57, 10)),
            ("Kattegat", "Atlantic Ocean", (56, 12)),
            ("Gulf of Bothnia", "Atlantic Ocean", (63, 20)),
            ("Gulf of Finland", "Atlantic Ocean", (60, 25)),
            ("Sulu Sea", "Pacific Ocean", (8, 120)),
            ("Gulf of Oman", "Indian Ocean", (24, 58)),
            ("North Channel", "Atlantic Ocean", (55, -6)),
            ("Norfolk Strait", "Pacific Ocean", (-39, 174)),
            ("Bass Strait", "Pacific Ocean", (-40, 146)),
            ("Korea Strait", "Pacific Ocean", (34, 129)),
            ("Denmark Strait", "Atlantic Ocean", (66, -30)),
            ("Fram Strait", "Arctic Ocean", (79, 0)),
            ("Yucatan Channel", "Atlantic Ocean", (21.5, -86)),
            ("Taiwan Strait", "Pacific Ocean", (24, 119))
        ]

        # Geolocator for location lookups
        self.geolocator = Nominatim(user_agent="interactive_group_editor")

        # Before plotting, try to update lat_lon and location from gb files
        self._update_all_locations_from_gb()

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

    def _parse_genus_species_from_path(self, csv_path):
        relative_path = os.path.relpath(csv_path, 'sequences/Interactive_group_editor')
        parts = relative_path.split(os.sep)
        genus = parts[0]
        species_name = os.path.basename(csv_path).replace("_clustered.csv", "")
        return genus, species_name

    def _get_or_create_csv_files(self):
        csv_files = self._get_csv_files()
        if not csv_files:
            self._create_csv_from_mds()
            csv_files = self._get_csv_files()
        return csv_files

    def _get_csv_files(self):
        base_dir = os.path.join("sequences", "Interactive_group_editor")
        csv_files = {}
        for root, dirs, files in os.walk(base_dir):
            for file in files:
                if file.endswith("_clustered.csv"):
                    family = os.path.basename(root)
                    if family not in csv_files:
                        csv_files[family] = []
                    csv_files[family].append(os.path.join(root, file))
        return sorted(csv_files.items())

    def _create_csv_from_mds(self):
        for root, dirs, files in os.walk(self.base_directory):
            for file in files:
                if file.endswith("_mds.mds"):
                    mds_path = os.path.join(root, file)
                    species_name = file.replace("_mds.mds", "")
                    rel_path = os.path.relpath(mds_path, self.base_directory)
                    parts = rel_path.split(os.sep)
                    genus = parts[0]

                    interactive_dir = os.path.join("sequences", "Interactive_group_editor", genus, species_name)
                    if not os.path.exists(interactive_dir):
                        os.makedirs(interactive_dir)

                    csv_path = os.path.join(interactive_dir, f"{species_name}_clustered.csv")

                    if not os.path.exists(csv_path):
                        try:
                            df_mds = pd.read_csv(mds_path, delim_whitespace=True)
                            if not {'FID', 'IID', 'C1', 'C2'}.issubset(df_mds.columns):
                                logging.error(f"MDS file {mds_path} missing required columns.")
                                continue

                            df_mds['Group'] = 1
                            df_mds['lat_lon'] = 'N/A'
                            df_mds['location'] = 'N/A'
                            df_mds['Location Source'] = 'N/A'

                            df_mds = df_mds[['FID', 'IID', 'C1', 'C2', 'Group', 'lat_lon', 'location', 'Location Source']]
                            df_mds.to_csv(csv_path, index=False)
                            logging.info(f"Created CSV from MDS file: {csv_path}")

                            fasta_path = os.path.join(os.path.dirname(mds_path), f"{species_name}_concatenated_gapped_sequences.fasta")
                            if os.path.exists(fasta_path):
                                dest_fasta_path = os.path.join(interactive_dir, f"{species_name}_concatenated_gapped_sequences.fasta")
                                shutil.copy2(fasta_path, dest_fasta_path)
                                logging.info(f"Copied FASTA file: {fasta_path} to {dest_fasta_path}")
                            else:
                                logging.warning(f"No FASTA file found for {species_name}.")
                        except Exception as e:
                            logging.error(f"Error creating CSV from MDS {mds_path}: {e}")

    def _ensure_columns_exist(self):
        required_columns = ['lat_lon', 'location', 'Location Source']
        for col in required_columns:
            if col not in self.df.columns:
                self.df[col] = 'N/A'

    def _load_data(self):
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

            interactive_dir = os.path.join("sequences", "Interactive_group_editor", self.genus, self.species_name)
            interactive_fasta_path = os.path.join(interactive_dir, f"{self.species_name}_concatenated_gapped_sequences.fasta")
            if os.path.exists(interactive_fasta_path):
                from Bio import SeqIO
                records = list(SeqIO.parse(interactive_fasta_path, "fasta"))
                df_iids = set(self.df['IID'].tolist())
                filtered_records = [r for r in records if r.id in df_iids]
                with open(interactive_fasta_path, "w") as handle:
                    SeqIO.write(filtered_records, handle, "fasta")

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

        # Update genus and species based on new file
        self.genus, self.species_name = self._parse_genus_species_from_path(self.csv_file)

        # Update locations from gb again for the new file
        self._update_all_locations_from_gb()

        self._load_data()
        self.fig.canvas.draw()

    def _save_csv_and_png(self):
        self._deselect_all_points()

        if hasattr(self, 'rectangle_selector'):
            self.rectangle_selector.set_active(False)

        time.sleep(1)

        if hasattr(self, 'rectangle_selector'):
            self.rectangle_selector.set_visible(False)
            self.fig.canvas.draw()

        # Re-apply lat_lon and location updates before saving
        self.df = self.df.apply(self._final_location_updates, axis=1)

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

        genus, species_name = self._parse_genus_species_from_path(csv_path)
        mds_file_name = f"{species_name}_mds.mds"
        fasta_file_name = f"{species_name}_concatenated_gapped_sequences.fasta"

        interactive_dir = os.path.join("sequences", "Interactive_group_editor", genus, species_name)
        if not os.path.exists(interactive_dir):
            os.makedirs(interactive_dir)

        def safe_copy(src, dst):
            if os.path.abspath(src) != os.path.abspath(dst):
                shutil.copy2(src, dst)

        # CSV
        dest_csv_path = os.path.join(interactive_dir, os.path.basename(csv_path))
        safe_copy(csv_path, dest_csv_path)

        # PNG
        dest_png_path = os.path.join(interactive_dir, os.path.basename(png_path))
        safe_copy(png_path, dest_png_path)

        # MDS
        mds_file_path = os.path.join(interactive_dir, mds_file_name)
        if os.path.exists(mds_file_path):
            dest_mds_path = os.path.join(interactive_dir, os.path.basename(mds_file_path))
            safe_copy(mds_file_path, dest_mds_path)

        # FASTA
        fasta_file_path = os.path.join(interactive_dir, fasta_file_name)
        if os.path.exists(fasta_file_path):
            dest_fasta_path = os.path.join(interactive_dir, fasta_file_name)
            safe_copy(fasta_file_path, dest_fasta_path)

        if hasattr(self, 'rectangle_selector'):
            self.rectangle_selector.set_active(True)
            self.rectangle_selector.set_visible(True)

        messagebox.showinfo("Save Successful", f"Files have been saved to:\n{interactive_dir}")

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

    def _parse_lat_lon(self, lat_lon_str):
        if isinstance(lat_lon_str, float):
            return None
        try:
            parts = lat_lon_str.strip().split()
            if len(parts) < 2:
                return None
            lat = float(parts[0])
            # Check next part for hemisphere
            if parts[1].upper().startswith('S'):
                lat = -lat
            elif parts[1].upper().startswith('N'):
                pass
            else:
                # no hemisphere info, assume first coordinate is latitude anyway
                # if no hemisphere given, just trust sign
                pass

            if len(parts) >= 4:
                lon = float(parts[2])
                if parts[3].upper().startswith('W'):
                    lon = -lon
            elif len(parts) == 2:
                # just lat, lon with no hemisphere (unlikely)
                lon = float(parts[1])
            else:
                return None
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

    def _get_lat_lon_and_location(self, iid):
        parts = iid.split("_", 1)
        if len(parts) > 1:
            acc_id = parts[0]
            full_species = parts[1]
        else:
            acc_id = iid.split("_")[0]
            full_species = self.species_name

        gb_dir = os.path.join("sequences", "gb", full_species)
        lat_lon = "N/A"
        location = "N/A"
        location_source = "N/A"

        possible_files = [
            os.path.join(gb_dir, acc_id + ".gb"),
            os.path.join(gb_dir, acc_id + ".gbk")
        ]

        gb_file_path = None
        for p in possible_files:
            if os.path.exists(p):
                gb_file_path = p
                break

        if not gb_file_path:
            logging.debug(f"No gb files found for IID {iid} in {gb_dir}.")
            return lat_lon, location, location_source

        logging.debug(f"Parsing GenBank file {gb_file_path} for IID {iid}")
        try:
            with open(gb_file_path, "r") as gb_file:
                for record in SeqIO.parse(gb_file, "genbank"):
                    found_geo_loc = None
                    found_pop_variant = None
                    found_country = None
                    found_altitude = None

                    if "lat_lon" in record.annotations:
                        lat_lon = record.annotations["lat_lon"]

                    for feature in record.features:
                        if "lat_lon" in feature.qualifiers:
                            lat_lon = feature.qualifiers["lat_lon"][0]
                        if "geo_loc_name" in feature.qualifiers:
                            found_geo_loc = feature.qualifiers["geo_loc_name"][0]
                        if "pop_variant" in feature.qualifiers:
                            found_pop_variant = feature.qualifiers["pop_variant"][0]
                        if "country" in feature.qualifiers:
                            found_country = feature.qualifiers["country"][0]
                        if "altitude" in feature.qualifiers:
                            found_altitude = feature.qualifiers["altitude"][0]

                    if found_geo_loc and found_geo_loc.strip():
                        location = found_geo_loc
                        location_source = "gb file"
                    elif found_pop_variant and found_pop_variant.strip():
                        location = found_pop_variant
                        location_source = "gb file"
                    elif found_country and found_country.strip():
                        location = found_country
                        location_source = "gb file"
                    elif found_altitude and found_altitude.strip():
                        location = found_altitude
                        location_source = "gb file"

                    # If lat_lon found but no location from gb file, try geolocator
                    if lat_lon != "N/A" and (location == "N/A" or location.strip() == ""):
                        resolved_name = self._get_location_name(lat_lon)
                        if resolved_name and resolved_name != "N/A":
                            location = resolved_name
                            location_source = "lat_lon lookup"
                        else:
                            # No geolocator location, try oceans
                            major_ocean, ocean_name = self._get_closest_ocean(lat_lon)
                            if major_ocean != "N/A":
                                location = f"{major_ocean}, {ocean_name}"
                                location_source = "lat_lon lookup"

                    break
        except Exception as e:
            logging.error(f"Error reading GenBank file {gb_file_path}: {e}")

        logging.debug(f"For IID {iid}, extracted lat_lon={lat_lon}, location={location}, source={location_source}")
        return lat_lon, location, location_source

    def _update_all_locations_from_gb(self):
        logging.debug("Updating all locations from gb files at initialization.")

        def update_lat_lon_and_location(row):
            if ((pd.isna(row['lat_lon']) or row['lat_lon'] == "N/A")
                and (pd.isna(row['location']) or row['location'] == "N/A")):
                lat_lon, loc, loc_source = self._get_lat_lon_and_location(row['IID'])
                row['lat_lon'] = lat_lon
                row['location'] = loc
                row['Location Source'] = loc_source
            else:
                # If we have lat_lon but no location, try reverse lookup and ocean fallback
                if row['lat_lon'] != "N/A" and (row['location'] == "N/A" or row['location'].strip() == ""):
                    resolved_name = self._get_location_name(row['lat_lon'])
                    if resolved_name and resolved_name != "N/A":
                        row['location'] = resolved_name
                        row['Location Source'] = "lat_lon lookup"
                    else:
                        # fallback to ocean
                        major_ocean, ocean_name = self._get_closest_ocean(row['lat_lon'])
                        if major_ocean != "N/A":
                            row['location'] = f"{major_ocean}, {ocean_name}"
                            row['Location Source'] = "lat_lon lookup"

                if row['location'] != "N/A" and (row['Location Source'] == "N/A" or row['Location Source'].strip() == ""):
                    row['Location Source'] = "gb file"
            return row

        self.df = self.df.apply(update_lat_lon_and_location, axis=1)

    def _final_location_updates(self, row):
        if row['lat_lon'] != "N/A" and (row['location'] == "N/A" or row['location'].strip() == ""):
            resolved_name = self._get_location_name(row['lat_lon'])
            if resolved_name and resolved_name != "N/A":
                row['location'] = resolved_name
                row['Location Source'] = "lat_lon lookup"
            else:
                # fallback to ocean
                major_ocean, ocean_name = self._get_closest_ocean(row['lat_lon'])
                if major_ocean != "N/A":
                    row['location'] = f"{major_ocean}, {ocean_name}"
                    row['Location Source'] = "lat_lon lookup"

        if row['location'] != "N/A" and (row['Location Source'] == "N/A" or row['Location Source'].strip() == ""):
            row['Location Source'] = "gb file"
        return row


# Example usage:
if __name__ == '__main__':
    base_directory = 'sequences/CDS_Genus'
    editor = InteractiveGroupEditor(base_directory)
    editor.show()
