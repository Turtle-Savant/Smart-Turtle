
# -*- coding: utf-8 -*-
"""
Created on Mon Jul  1 20:19:49 2024

@author: Coleman Nielsen, Ph.D Candidate 

Copyright Turtle-Savant 2024

License:  GNU General Public License (GPL) v3.0.
"""

import tkinter as tk
import numpy as np
from tkinter import filedialog
import pandas as pd
import re
import random
from collections import Counter
import os
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from tkinter import ttk
from decimal import Decimal, getcontext
import traceback
from tkinter import PhotoImage
from PIL import ImageTk, Image
getcontext().prec = 10  # Adjust the precision as needed

class IsotopeDistributionSimulator:
    def __init__(self, master):
        self.master = master
        master.title("Smart Turtle: Isotope Distribution and Mass Spectrometry Simulator (1.18)")
        
        # Set the icon
        self.icon = PhotoImage(file='smart_turtle_cropped_soft_edges.png')
        master.iconphoto(False, self.icon)
        
        self.colors = ['b', 'r', 'c', 'm', 'y']
        self.current_color = 'b'
        self.current_index = 0
        self.run_counter = 0
        self.max = 0
        self.mz = 0
    
        # Configure the grid layout
        master.grid_rowconfigure(0, weight=1)
        master.grid_columnconfigure(0, weight=1)
        master.grid_columnconfigure(1, weight=1)
        master.grid_rowconfigure(1, weight=0)
        
        # Create a frame to contain the text output and the plot
        self.output_frame = ttk.Frame(master)
        self.output_frame.grid(row=1, column=0, sticky="nsew")
    
        # Create a text widget to display the isotopic distribution
        self.text_output = tk.Text(self.output_frame, height=10, width=35)
        self.text_output.pack(fill=tk.BOTH, expand=True, padx=10, pady=10)
        
    
        # Create a notebook to organize the widgets
        self.notebook = ttk.Notebook(master)
        self.notebook.grid(row=1, column=1, sticky="nsew")
    
        # Create a frame for basic settings
        self.basic_frame = ttk.Frame(self.notebook)
        self.notebook.add(self.basic_frame, text='Basic Settings')
    
        # Create a frame for advanced settings
        self.advanced_frame = ttk.Frame(self.notebook)
        self.notebook.add(self.advanced_frame, text='Advanced Settings')
        
        # Create a frame for isotopic enrichment settings
        self.enrichment_frame = ttk.Frame(self.notebook)
        self.notebook.add(self.enrichment_frame, text='Isotope Enrichment')
    
        # Create labels and entry widgets for user inputs in the basic frame
        self.label_formula = tk.Label(self.basic_frame, text="Chemical Formula:")
        self.label_formula.pack()
        
        self.entry_formula = tk.Entry(self.basic_frame)
        self.entry_formula.pack()
        
        # Create a checkbox for showing Gaussian data in the basic frame
        self.gaussian_var = tk.BooleanVar(value=True)  # Default checked
        self.gaussian_checkbox = tk.Checkbutton(self.basic_frame, text="Show Profile Data", variable=self.gaussian_var)
        self.gaussian_checkbox.pack()
        
        # Create a checkbox for showing centroided data in the basic frame
        self.centroid_var = tk.BooleanVar(value=False)  # Default checked
        self.centroid_checkbox = tk.Checkbutton(self.basic_frame, text="Show Centroided Data", variable=self.centroid_var)
        self.centroid_checkbox.pack()
    
        
    
        # Show Centroided or not
        self.show_centroid = self.centroid_var.get()
    
        # Show Gaussian or not
        self.show_gaussian = self.gaussian_var.get()
    
       # Create a frame to hold the buttons side by side
        self.button_frame = ttk.Frame(master)
        self.button_frame.grid(row=2, column=0, columnspan=2)
        
        # Load your image
        image_path = "smart_turtle_cropped_soft_edges.png"  
        img = Image.open(image_path)
        img = img.resize((70, 70), Image.LANCZOS)
        self.img = ImageTk.PhotoImage(img)  
        
        # Create a label to display the image
        self.image_label = tk.Label(master, image=self.img)
        self.image_label.grid(row=2, column=0, sticky='sw', padx=5, pady=5)  # Place the image in the bottom left corner
        
        # Copyright
        cr_text = "© 2024 Coleman Nielsen. All rights reserved."
        self.cr_text_label = tk.Label(master, text=cr_text)
        self.cr_text_label.grid(row=2, column=1, sticky='se', padx=5, pady=5) 
        
        
        # Create a button to trigger the simulation in the basic frame
        self.simulate_button_basic = tk.Button(self.button_frame, text="Simulate", command=self.run_simulation_and_next_value, bg="#9CCC65", fg="black")
        self.simulate_button_basic.pack(side=tk.LEFT, padx=5, pady=5)
        
        # Create a button to reset the plots
        self.reset_button = tk.Button(self.button_frame, text="Reset Graph", command=self.reset_plots, bg="#FFAB91", fg="black")
        self.reset_button.pack(side=tk.LEFT, padx=5, pady=5)
        
        # Add a button to save the graph
        self.save_button = tk.Button(self.button_frame, text="Save Graph", command=self.save_graph, bg="#90CAF9", fg="black")
        self.save_button.pack(side=tk.LEFT, padx=5, pady=5)

    
        # Create labels and entry widgets for advanced settings with default values
        # Left column
        self.label_x = tk.Label(self.advanced_frame, text="Number of Molecules:")
        self.label_x.grid(row=0, column=0, sticky="e", padx=10, pady=5)
    
        self.entry_x = tk.Entry(self.advanced_frame)
        self.entry_x.grid(row=0, column=1, sticky="w", padx=10, pady=5)
        self.entry_x.insert(tk.END, "10000")  # Set the default value to 10000
    
        self.label_resolving_power = tk.Label(self.advanced_frame, text="Resolving Power:")
        self.label_resolving_power.grid(row=1, column=0, sticky="w", padx=10, pady=5)
    
        self.entry_resolving_power = tk.Entry(self.advanced_frame)
        self.entry_resolving_power.grid(row=1, column=1, sticky="w", padx=10, pady=5)
        self.entry_resolving_power.insert(tk.END, '10000')
        
        self.label_noise = tk.Label(self.advanced_frame, text="% Baseline:")
        self.label_noise.grid(row=2, column=0, sticky="w", padx=10, pady=5)
    
        self.entry_noise = tk.Entry(self.advanced_frame)
        self.entry_noise.grid(row=2, column=1, sticky="w", padx=10, pady=5)
        self.entry_noise.insert(tk.END, '0')
    
        self.label_adduct = tk.Label(self.advanced_frame, text="Select Adduct Type:")
        self.label_adduct.grid(row=3, column=0, sticky="w", padx=10, pady=5)
    
        adduct_options = [
    "+H⁺",
    "+Na⁺",
    "+K⁺",
    "+NH₄⁺",
    "+2H²⁺",
    "+2Na²⁺",
    "++2K²⁺",
    "+2NH₄²⁺",
    "-H⁻",
    "+Cl⁻",
    "+Br⁻",
    "+I⁻",
    "+F⁻",
    "+HCOO⁻",
    "+CH₃COO⁻",
    "-2H²⁻",
    "+2Cl²⁻",
    "+2Br²⁻",
    "+2I²⁻",
    "+2F²⁻",
    "+2HCOO²⁻",
    "+2CH₃COO²⁻"
]
        self.combobox_adduct = ttk.Combobox(self.advanced_frame, values=adduct_options)
        self.combobox_adduct.set(adduct_options[0])  # Set the default value
        self.combobox_adduct.grid(row=3, column=1, sticky="w", padx=10, pady=5)
    
        # Right column
        #self.compress_var = tk.BooleanVar(value=True)  # Default checked
        #self.compress_checkbox = tk.Checkbutton(self.advanced_frame, text="Compress Elemental Isotopomers", variable=self.compress_var)
        #self.compress_checkbox.grid(row=0, column=2, sticky="w", padx=10, pady=5)
        
        self.fix_graph_var = tk.BooleanVar(value=False)  # Default unchecked
        self.fix_graph_checkbox = tk.Checkbutton(self.advanced_frame, text="Fixate Graph", variable=self.fix_graph_var, command=self.toggle_mz_entry)
        self.fix_graph_checkbox.grid(row=1, column=2, sticky="w", padx=10, pady=5)
        
        self.label_mz_range = tk.Label(self.advanced_frame, text="m/z range:")
        self.label_mz_range.grid(row=2, column=2, sticky="w", padx=10, pady=5)
    
        self.entry_mz_range = tk.Entry(self.advanced_frame, state=tk.DISABLED)  # Greyed out initially
        self.entry_mz_range.grid(row=2, column=3, sticky="w", padx=0, pady=5)
        
        self.help_tools_var = tk.BooleanVar(value=False)  # Default unchecked
        self.help_tools_checkbox = tk.Checkbutton(self.advanced_frame, text="Activate Help Tools?", variable=self.help_tools_var)
        self.help_tools_checkbox.grid(row=3, column=2, sticky="w", padx=10, pady=5)
        
    
        # Create labels and entry widgets for isotopic enrichment settings in a 2x2 grid
        self.label_deuterium = tk.Label(self.enrichment_frame, text="% Deuterium")
        self.label_deuterium.grid(row=0, column=0, padx=10, pady=10, sticky="e")
        
        self.entry_deuterium = tk.Entry(self.enrichment_frame)
        self.entry_deuterium.grid(row=0, column=1, padx=10, pady=10, sticky="w")
        
        self.label_c13 = tk.Label(self.enrichment_frame, text="% Carbon-13")
        self.label_c13.grid(row=1, column=0, padx=10, pady=10, sticky="e")
        
        self.entry_c13 = tk.Entry(self.enrichment_frame)
        self.entry_c13.grid(row=1, column=1, padx=10, pady=10, sticky="w")
        
        self.label_n15 = tk.Label(self.enrichment_frame, text="% Nitrogen-15")
        self.label_n15.grid(row=0, column=2, padx=10, pady=10, sticky="e")
        
        self.entry_n15 = tk.Entry(self.enrichment_frame)
        self.entry_n15.grid(row=0, column=3, padx=10, pady=10, sticky="w")
        
        self.label_o18 = tk.Label(self.enrichment_frame, text="% Oxygen-18")
        self.label_o18.grid(row=1, column=2, padx=10, pady=10, sticky="e")
        
        self.entry_o18 = tk.Entry(self.enrichment_frame)
        self.entry_o18.grid(row=1, column=3, padx=10, pady=10, sticky="w")

        # Create a frame to contain the plot
        self.plot_frame = ttk.Frame(master)
        self.plot_frame.grid(row=0, column=0, columnspan=2, sticky="nsew")
    
        # Create a figure for the plot
        self.fig, self.ax = plt.subplots()
        self.ax.set_xlabel('m/z')
        self.ax.set_ylabel('Normalized Intensity')
        self.ax.set_ylim(0, 110)
    
        # Adjust the bottom margin to make the x-axis label more visible
        self.fig.subplots_adjust(bottom=0.2)
    
        # Create a canvas for the plot and pack it into the GUI
        self.canvas = FigureCanvasTkAgg(self.fig, master=self.plot_frame)
        self.canvas.get_tk_widget().pack(fill=tk.BOTH, expand=True, padx=10, pady=10)
    
        # Compress data or not
        #self.compress = self.compress_var.get()
    
        # Add hover tooltips to all labels and checkboxes
        self.add_tooltips()
        
    def toggle_mz_entry(self):
        if self.fix_graph_var.get():
            self.entry_mz_range.config(state=tk.NORMAL)
        else:
            self.entry_mz_range.config(state=tk.DISABLED)
    
    
    

    def show_tooltip(self, event, text):
        if self.tooltip is None:
            self.tooltip = tk.Toplevel(self.master)
            self.tooltip.wm_overrideredirect(True)
            self.tooltip.wm_geometry(f"+{event.x_root + 20}+{event.y_root + 10}")
            label = tk.Label(self.tooltip, text=text, background="lightblue", relief="flat", borderwidth=0.5)
            label.pack()

    def hide_tooltip(self, event):
        if self.tooltip:
            self.tooltip.destroy()
            self.tooltip = None

    def add_tooltips(self):
        self.tooltip = None
        
        def update_tooltips():
            if self.help_tools_var.get():
                # Bind tooltips
                self.label_formula.bind("<Enter>", lambda event: self.show_tooltip(event, "Enter the chemical formula for your molecule\n Example Formats: C6H12O6, NaCOOH"))
                self.label_formula.bind("<Leave>", self.hide_tooltip)
                
                self.fix_graph_checkbox.bind("<Enter>", lambda event: self.show_tooltip(event, "Check to fixate the graph to a specific m/z range"))
                self.fix_graph_checkbox.bind("<Leave>", self.hide_tooltip)
        
                self.label_mz_range.bind("<Enter>", lambda event: self.show_tooltip(event, "Enter the desired m/z range for the fixed graph"))
                self.label_mz_range.bind("<Leave>", self.hide_tooltip)
        
            
                self.centroid_checkbox.bind("<Enter>", lambda event: self.show_tooltip(event, "Select to show a centroided mass spectrum"))
                self.centroid_checkbox.bind("<Leave>", self.hide_tooltip)
        
                self.gaussian_checkbox.bind("<Enter>", lambda event: self.show_tooltip(event, "Select to show a profile (guassian) mass spectrum"))
                self.gaussian_checkbox.bind("<Leave>", self.hide_tooltip)
        
                self.label_x.bind("<Enter>", lambda event: self.show_tooltip(event, "Enter the number of molecules you would like to simulate.\n More molecules results in greater accuracy,\n with increased processing time"))
                self.label_x.bind("<Leave>", self.hide_tooltip)
        
                self.label_resolving_power.bind("<Enter>", lambda event: self.show_tooltip(event, "Resolving power (R) represents the ratio of the mass/charge (m/z) value of an ion,\n versus the smallest differentiable mass/charge value difference (Δm/z)\n for that mass spectrometer, such that R=(m/z)/(Δm/z).\nFor example if our ion of interest has an m/z of 100 and we are able to distinguish\n the difference between an ion of 99.9 and 100.0 (Δm/z= 0.1),\nthen our mass resolution is 100/0.1 = 1000.\n A mass resolution of around 10 thousand is quite common for mass spectrometers,\n while the best can achieve mass resolutions of over 1 million"))
                self.label_resolving_power.bind("<Leave>", self.hide_tooltip)
        
                self.label_noise.bind("<Enter>", lambda event: self.show_tooltip(event, "Enter the baseline noise percentage"))
                self.label_noise.bind("<Leave>", self.hide_tooltip)
        
                self.label_adduct.bind("<Enter>", lambda event: self.show_tooltip(event, "In order for a molecule to be measured by a mass spectrometer, it must first be ionized.\nThis process most commonly happens during a process called electrospray ionization.\nAn electric field is used to isolate small, charged particles at the end of a stream of solvent.\nThis clustering of charged particles causes increased interactions of these species\nwith analytes of interest, causing those analytes to become charged as well.\nThe interaction forms an ionized 'adduct', which has a unique isotopic distribution."))
                self.label_adduct.bind("<Leave>", self.hide_tooltip)
        
                #self.compress_checkbox.bind("<Enter>", lambda event: self.show_tooltip(event, "Isotopomers are molecules which have the same number of neutrons,\n yet differ in the location of those neutrons.\n In the case where a molecule contains an extra hydrogen neutron,\nand an isomeric molecule has an extra carbon neutron,\n a difference in nuclear binding energy will cause the two\n to have slightly different masses.\nMass spectrometers with high mass resolution can actually\n resolve the mass difference between these isotopomers.\nFor the sake of simplicity, this setting collapses such isotopomers\n into a single data point."))
                #self.compress_checkbox.bind("<Leave>", self.hide_tooltip)
        
                self.label_deuterium.bind("<Enter>", lambda event: self.show_tooltip(event, "Isotopic labeling is a common practice\nwhere researchers intentionally increase the percentage of heavy isotopes\nbeyond their natural occurance.\nThe result has significant impacts on the isotopic distributions\nof molecules in the affected system.\nEnter the percentage of artifical deuterium enrichment"))
                self.label_deuterium.bind("<Leave>", self.hide_tooltip)
        
                self.label_c13.bind("<Enter>", lambda event: self.show_tooltip(event, "Isotopic labeling is a common practice\nwhere researchers intentionally increase the percentage of heavy isotopes\nbeyond their natural occurance.\nThe result has significant impacts on the isotopic distributions\nof molecules in the affected system.\nEnter the percentage of artificial carbon-13 enrichment"))
                self.label_c13.bind("<Leave>", self.hide_tooltip)
        
                self.label_n15.bind("<Enter>", lambda event: self.show_tooltip(event, "Isotopic labeling is a common practice\nwhere researchers intentionally increase the percentage of heavy isotopes\nbeyond their natural occurance.\nThe result has significant impacts on the isotopic distributions\nof molecules in the affected system.\nEnter the percentage of artificial nitrogen-15 enrichment"))
                self.label_n15.bind("<Leave>", self.hide_tooltip)
        
                self.label_o18.bind("<Enter>", lambda event: self.show_tooltip(event, "Isotopic labeling is a common practice\nwhere researchers intentionally increase the percentage of heavy isotopes\nbeyond their natural occurance.\nThe result has significant impacts on the isotopic distributions\nof molecules in the affected system.\nEnter the percentage of artificial oxygen-18 enrichment"))
                self.label_o18.bind("<Leave>", self.hide_tooltip)
                
            else:
                # Unbind tooltips
                self.label_formula.unbind("<Enter>")
                self.label_formula.unbind("<Leave>")
                
                self.fix_graph_checkbox.unbind("<Enter>")
                self.fix_graph_checkbox.unbind("<Leave>")
        
                self.label_mz_range.unbind("<Enter>")
                self.label_mz_range.unbind("<Leave>")
        
        
                self.centroid_checkbox.unbind("<Enter>")
                self.centroid_checkbox.unbind("<Leave>")
        
                self.gaussian_checkbox.unbind("<Enter>")
                self.gaussian_checkbox.unbind("<Leave>")
        
                self.label_x.unbind("<Enter>")
                self.label_x.unbind("<Leave>")
        
                self.label_resolving_power.unbind("<Enter>")
                self.label_resolving_power.unbind("<Leave>")
        
                self.label_noise.unbind("<Enter>")
                self.label_noise.unbind("<Leave>")
        
                self.label_adduct.unbind("<Enter>")
                self.label_adduct.unbind("<Leave>")
        
                #self.compress_checkbox.unbind("<Enter>")
                #self.compress_checkbox.unbind("<Leave>")
        
                self.label_deuterium.unbind("<Enter>")
                self.label_deuterium.unbind("<Leave>")
        
                self.label_c13.unbind("<Enter>")
                self.label_c13.unbind("<Leave>")
        
                self.label_n15.unbind("<Enter>")
                self.label_n15.unbind("<Leave>")
        
                self.label_o18.unbind("<Enter>")
                self.label_o18.unbind("<Leave>")
        
        update_tooltips()  # Initial setup
    
        # Monitor changes to self.help_tools_var and update tooltips accordingly
        self.help_tools_var.trace_add('write', lambda *args: update_tooltips())
        
    def save_graph(self):
        file_path = filedialog.asksaveasfilename(defaultextension=".png", filetypes=[("PNG files", "*.png"), ("All files", "*.*")])
        if file_path:
            self.fig.savefig(file_path,format='png', dpi=300)
            
    def reset_plots(self):
        # Clear the plot
        self.ax.clear()
        self.ax.set_xlabel('m/z')
        self.ax.set_ylabel('Normalized Intensity')
        self.canvas.draw_idle()
        self.ax.set_ylim(0, 110)
        
        # Reset the colors
        self.current_color = 'b'
        self.current_index = 0
        self.run_counter = 0
        
    
    def update_boxes(self):
       # Update the value of self.compress when the checkbox is clicked
       #self.compress = self.compress_var.get()
       self.show_centroid = self.centroid_var.get()
       self.show_gaussian = self.gaussian_var.get()
       #print(self.show_gaussian)
       #print(self.compress)


    def next_value(self):
        # Display the next value from the list
        self.current_index = (self.current_index + 1) % len(self.colors)
        self.current_color = self.colors[self.current_index]
        self.run_counter = self.run_counter + 1
        
        
    def modify_cf_for_adducts(self): # Not in use yet, plug into version 16
        adduct = self.combobox_adduct.get()
        if adduct.startswith('+'):
            adduct_type = 'positive'
        elif adduct.startswith('-'):
            adduct_type = 'negative'
        adducts = {
    "+H⁺": {'H': 1},
    "+Na⁺": {'Na': 1},
    "+K⁺": {'K': 1},
    "+NH₄⁺": {'N': 1, 'H': 4},
    "+2H²⁺": {'H': 2},
    "+2Na²⁺": {'Na': 2},
    "++2K²⁺ ": {'K': 2},
    "+2NH₄²⁺": {'N': 2, 'H': 8},
    "-H⁻": {'H': -1},
    "+Cl⁻": {'Cl': 1}, #something is weird here
    "+Br⁻": {'Br': 1},
    "+I⁻": {'I': 1},
    "+F⁻": {'F': 1},
    "+HCOO⁻": {'H': 1, 'C': 1, 'O':2},  # Formate
    "+CH₃COO⁻": {'H': 3, 'C': 2, 'O':2},  # Acetate
    "-2H²⁻": {'H': -2},
    "+2Cl²⁻": {'Cl': 2},
    "+2Br²⁻": {'Br': 2},
    "+2I²⁻": {'I': 2},
    "+2F²⁻": {'F': 2},
    "+2HCOO²⁻": {'H': 2, 'C': 2, 'O':4},
    "+2CH₃COO²⁻": {'H': 6, 'C': 4, 'O':4}
}
        
        
        if adduct_type == 'positive':
            adduct_dict = {k: v for k, v in adducts.items() if k.startswith('+')}
        elif adduct_type == 'negative':
            adduct_dict = {k: v for k, v in adducts.items() if k.startswith('-')}
        else:
            raise ValueError("Invalid adduct type. Must be 'positive' or 'negative'.")
        
        if adduct in adduct_dict:
            for element, count in adduct_dict[adduct].items():
                if element in self.cf_dictionary:
                    self.cf_dictionary[element] += count
                    # Check if any element count goes below zero
                    if self.cf_dictionary[element] < 0:
                        raise ValueError(f"Element count for {element} cannot be negative.")
                else:
                    # If the element does not exist in the dictionary, add it
                    self.cf_dictionary[element] = count
                    # Check if any element count goes below zero
                    if self.cf_dictionary[element] < 0:
                        raise f"Element count for {element} cannot be negative."
        else:
            print(f"{adduct} not recognized for type {adduct_type}")
            
        if adduct in  {
    "+H⁺": {'H': 1},
    "+Na⁺": {'Na': 1},
    "+K⁺": {'K': 1},
    "+NH₄⁺": {'N': 1, 'H': 4},
    "-H⁻": {'H': -1},
    "+Cl⁻": {'Cl': 1},
    "+Br⁻": {'Br': 1},
    "+I⁻": {'I': 1},
    "+F⁻": {'F': 1},
    "+HCOO⁻": {'H': 1, 'C': 1, 'O':2},  # Formate
    "+CH₃COO⁻": {'H': 3, 'C': 2, 'O':2}}.keys():
            self.charge = 1
        else: 
            self.charge = 2
        #print(self.charge)
        #print(self.cf_dictionary)

    def create_chemical_formula_dictionary(self):
            
        # Create an empty dictionary to store DataFrames
        element_dfs = {}

        # Group the DataFrame by the 'Name' column
        grouped = self.abundances_and_exact_masses.groupby('Code')

        # Iterate through groups and create a dictionary entry for each group
        for code, group_df in grouped:
            element_dfs[code] = group_df

        # Create the atomic pools from which to sample
        self.atomic_pools = {}
        for key in element_dfs:
            pool = []
            for row in element_dfs[key].iterrows():
                mass = row[1]['Mass']
                abundance = int((row[1]['Abund.']) / 100 * 100000)
                pool.extend([mass] * abundance)
            self.atomic_pools[key] = pool

        # Turn chemical formula into a dictionary of atoms and associated amounts
        chemical_formula_list = self.split_string(self.chemical_formula)
        result_list = [re.split(r'(\d+)', s, 1)[:2] for s in chemical_formula_list]
        
        #Handle multiple instances of a particular element
        unique_first_values = []
        for item in result_list:
            if item[0] not in unique_first_values:
                unique_first_values.append(item[0])

        chemical_formula_dictionary = {}
        for i in unique_first_values:
            occurrences = 0
            for j in result_list:
                if i == j[0]:
                    occurrences += float(j[1])
            chemical_formula_dictionary[i] = occurrences
        self.cf_dictionary = chemical_formula_dictionary
                

    def run_simulation(self):
        # Read the CSV file into a DataFrame
        script_directory = os.path.dirname(os.path.abspath(__file__))
        file_path = ('SCIENTIFIC INSTRUMENT SERVICES_df.csv')
        self.abundances_and_exact_masses = pd.read_csv(file_path)
        self.abundances_and_exact_masses['Code'] = self.abundances_and_exact_masses['Symbol'].apply(
            lambda symbol: self.extract_left_of_parenthesis(symbol))

        
        self.update_boxes()
        # Get user inputs from entry widgets
        self.chemical_formula = self.entry_formula.get()
        if not self.chemical_formula:
            raise ValueError("Chemical formula cannot be empty")
            
        
        
            
        
            
        x = int(self.entry_x.get())
        
        # Read the CSV file into a DataFrame
        script_directory = os.path.dirname(os.path.abspath(__file__))
        file_path = ('SCIENTIFIC INSTRUMENT SERVICES_df.csv')
        abundances_and_exact_masses = pd.read_csv(file_path)
        #print(abundances_and_exact_masses)
        abundances_and_exact_masses['Code'] = abundances_and_exact_masses['Symbol'].apply(
            lambda symbol: self.extract_left_of_parenthesis(symbol))
        
        def get_valid_float(entry_widget):
            try:
                value = float(entry_widget.get())
                return value, True
            except ValueError:
                return None, False

        abundances_and_exact_masses = self.abundances_and_exact_masses.copy()
        # Deuterium enrichment
        deuterium_enrichment, valid_deuterium = get_valid_float(self.entry_deuterium)
        if valid_deuterium:
    
            abundances_and_exact_masses.loc[abundances_and_exact_masses['Symbol'] == 'H(2)', 'Abund.'] = deuterium_enrichment
            abundances_and_exact_masses.loc[abundances_and_exact_masses['Symbol'] == 'H(1)', 'Abund.'] = 100 - deuterium_enrichment
            #print(abundances_and_exact_masses.loc[abundances_and_exact_masses['Symbol'] == 'H(2)', 'Abund.'].values)
        
        # Carbon-13 enrichment
        carbon13_enrichment, valid_c13 = get_valid_float(self.entry_c13)
        if valid_c13:
            abundances_and_exact_masses.loc[abundances_and_exact_masses['Symbol'] == 'C(13)', 'Abund.'] = carbon13_enrichment
            abundances_and_exact_masses.loc[abundances_and_exact_masses['Symbol'] == 'C(12)', 'Abund.'] = 100 - carbon13_enrichment
        
        # Nitrogen-15 enrichment
        nitrogen15_enrichment, valid_n15 = get_valid_float(self.entry_n15)
        if valid_n15:
            abundances_and_exact_masses.loc[abundances_and_exact_masses['Symbol'] == 'N(15)', 'Abund.'] = nitrogen15_enrichment
            abundances_and_exact_masses.loc[abundances_and_exact_masses['Symbol'] == 'N(14)', 'Abund.'] = 100 - nitrogen15_enrichment
        
        # Oxygen-18 enrichment
        oxygen18_enrichment, valid_o18 = get_valid_float(self.entry_o18)
        if valid_o18:
            abundances_and_exact_masses.loc[abundances_and_exact_masses['Symbol'] == 'O(18)', 'Abund.'] = oxygen18_enrichment
            abundances_and_exact_masses.loc[abundances_and_exact_masses['Symbol'] == 'O(16)', 'Abund.'] = 100 - oxygen18_enrichment
        self.abundances_and_exact_masses = abundances_and_exact_masses
        
        self.create_chemical_formula_dictionary()
        self.modify_cf_for_adducts()
            

        # Original simulation code...
        '''Sample each molecule x times in order to make a predicted isotopic distribution'''
        '''The higher the number of reps, the more precise your isotopic distribution will be'''
        individual_molecules = []
        #print(self.cf_dictionary)
        for molecule in range(0, x):
            mass = Decimal('0.0')
            #print(mass)
            for atom in self.cf_dictionary:
                #print(atom)
                for i in range(int(self.cf_dictionary[atom])):
                    random_atom_mass = Decimal(random.choice(self.atomic_pools[atom]))
                    #print(random_atom_mass)
                    mass += random_atom_mass
            #print(mass)
            individual_molecules.append(round(float(mass), 6))
            #print(individual_molecules)


        # Use Counter to count occurrences of each value
        value_counts_dictionary = Counter(sorted(individual_molecules))
        #print(value_counts_dictionary)
        
        # Print unique values and their counts
        isotopic_distribution = []
        isotopic_mz = []
        for value, count in value_counts_dictionary.items():
            isotopic_distribution.append(count)
            isotopic_mz.append(value/float(self.charge))
        #print(isotopic_distribution)
        #print(isotopic_mz)
        self.mz = isotopic_mz[0]
        #print(isotopic_mz[0])
        
        if self.run_counter == 0:
                self.normalization_value = max(isotopic_distribution)    
        
        normalized_distribution = [x/self.normalization_value *100 for x in isotopic_distribution]
        #print(f'isotope_distribution_len: {len()}')
        
        def centroided_float_values(float_values):
            
            def gaussian(x, area, mean, resolving_power):
                # Calculate FWHM based on resolving power
                fwhm = mean / resolving_power
                sigma = fwhm / (2 * np.sqrt(2 * np.log(2)))
                
                # Calculate the amplitude from the given area
                amplitude = area / (sigma * np.sqrt(2 * np.pi))
                
                return amplitude * np.exp(-(x - mean)**2 / (2 * sigma**2))
            
            def add_noise(input_list):
                # Get the noise range from the Entry widget
                noise_level = float(self.entry_noise.get())
            
                # Generate the modified list using a list comprehension
                modified_numbers = [num + random.uniform(-noise_level, noise_level) for num in input_list]
            
                return modified_numbers
            
            def normalize_list_by_max(input_list, total=100):
                if self.max < max(input_list):  # Step 1: Find the maximum value in the list
                    self.max = max(input_list)
                if self.max == 0:
                    return [0] * len(input_list)  # Handle case where all values are zero
                
                normalized_values = [value / self.max * total for value in input_list]  # Step 2: Normalize each value
                
                return normalized_values  # Step 3: Return normalized list
            
            # Initialize an array to store the cumulative spectrum
            common_x = np.arange(min(isotopic_mz) - 3 * max(isotopic_mz) / float(self.entry_resolving_power.get()), max(isotopic_mz) + 3 * max(isotopic_mz) / float(self.entry_resolving_power.get()) + 0.00001, 0.00001)
            total_spectrum = np.zeros_like(common_x)

            # Apply the method
            for i in range(len(isotopic_mz)):
                #print(f'float_val_len: {len(float_values)}')
                #print(f'isotope_mz_len: {len(isotopic_mz)}')
                y = gaussian(common_x, float_values[i], isotopic_mz[i], float(self.entry_resolving_power.get()))
                total_spectrum += y
            #print(max(total_spectrum))
            total_spectrum = normalize_list_by_max(total_spectrum)
            total_spectrum = add_noise(total_spectrum)
            #print(max(total_spectrum))
            if self.show_gaussian:
                plt.plot(common_x, total_spectrum, color = self.current_color)
            if self.show_centroid:
                plt.vlines(isotopic_mz, 0, float_values, linestyles='dashed', color = self.current_color, alpha = 0.5)
                print(isotopic_mz)
                
             #Set limits to avoid autoscaling too far out
            #print(self.fix_graph_var.get())
            try:
                if self.fix_graph_var.get() and float(self.entry_mz_range.get()) > 0.1:
                    self.ax.set_xlim(min(isotopic_mz) - 0.5, min(isotopic_mz) + float(self.entry_mz_range.get()))
            except:
                print('That did not work... Oh well.')
                
            
            #print(isotopic_mz)
            #print(float_values)
            
        centroided_float_values(normalized_distribution)
        
        def generate_neutromers(masses):
            # Create a dictionary to store Neutromers with the same x-value
            neutromers_dict = {}
        
            # Sort the masses in ascending order
            sorted_masses = sorted(enumerate(masses), key=lambda x: x[1])
            #print(sorted_masses)
        
            # Initialize variables to track the current x-value and letter index
            current_x_value = 0
            letter_index = 0
        
            # Iterate through the sorted masses
            #print(sorted_masses)
            for idx, mass in sorted_masses:
                x_value = round(float(str(mass - sorted_masses[0][1])[:3])*self.charge)
                #print(x_value)
        
                # Reset letter index when x-value changes
                if x_value != current_x_value:
                    current_x_value = x_value
                    letter_index = 0
        
                # Check if the current letter index is already used for the current x-value
                while f'M{x_value}{chr(ord("a") + letter_index)}' in neutromers_dict.get(x_value, []):
                    letter_index += 1
        
                x_label = chr(ord('a') + letter_index)
                letter_index += 1
        
                neutromer_label = f'M{x_value}{x_label}' if x_value >= 0 else f'M{x_label}'
        
                if x_value not in neutromers_dict:
                    neutromers_dict[x_value] = [neutromer_label]
                else:
                    neutromers_dict[x_value].append(neutromer_label)
        
            # Remove the lowercase letter if there is only one Neutromer for a particular x-value
            neutromers_list = [label if len(labels) > 1 else f'M{x}' for x, labels in neutromers_dict.items() for label in labels]
            #print(neutromers_list)
            return neutromers_list
            
        
        neutromers_list = generate_neutromers(isotopic_mz)
        #print(neutromers_list)
        
        # Display the normalized distribution in the text widget
        self.text_output.delete(1.0, tk.END)  # Clear previous content
        self.text_output.insert(tk.END, "Normalized Isotopic Distribution:\n")
        #print(normalized_distribution)
        #print(neutromers_list)

        for value, count in zip(range(len(normalized_distribution)), normalized_distribution):
            self.text_output.insert(tk.END, f"{neutromers_list[value]}, Intensity: {round(count, 3)}, m/z: {round(isotopic_mz[value],4)}\n")

        # Plot the isotopic distribution on the existing figure
        # Plot the isotopic distribution on the existing figure without clearing
        centroided_float_values(normalized_distribution)
        self.canvas.draw_idle()  # Update the canvas immediately
        
        

    
    def run_simulation_and_next_value(self):
        try:
            self.run_simulation()
            self.next_value()
        except Exception as e:
            tb = traceback.extract_tb(e.__traceback__)
            line_number = tb[-1].lineno
            self.text_output.delete(1.0, tk.END)
            self.text_output.insert(tk.END, f"Analysis failed. Error on line {line_number}: {e}")

    def extract_left_of_parenthesis(self, input_string):
         index = input_string.find('(')
         if index == -1:
             return input_string
         result = input_string[:index].strip()
         return result
 
 
    def separate_string(self, input_string):
         letters = ""
         numbers = ""
         for char in input_string:
             if char.isalpha():
                 letters += char
             elif char.isdigit():
                 numbers += char
         return (letters, numbers)
     
    def append_ones(self, string_list):
         result_list = []
 
         for value in string_list:
             if any(char.isdigit() for char in value):
                 result_list.append(value)
             else:
                 result_list.append(value + "1")
 
         return result_list
 
    def split_string(self,input_string):
         result_list = []
         current_word = ""
 
         for char in input_string:
             if char.isupper():
                 if current_word:
                     result_list.append(current_word)
                 current_word = char
             else:
                 current_word += char
 
         if current_word:
             result_list.append(current_word)
 
         return self.append_ones(result_list)
# Main application loop

# Main application loop
if __name__ == "__main__":
    root = tk.Tk()
    app = IsotopeDistributionSimulator(root)
    root.mainloop()
