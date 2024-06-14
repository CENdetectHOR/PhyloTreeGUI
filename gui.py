# -*- coding: utf-8 -*-
# @author Paolo Pagliuca <paolo.pagliuca@istc.cnr.it>

import math, copy, matplotlib
import matplotlib.pyplot as plt
from colour import Color
from matplotlib import cm
import tkinter as tk
import tkinter.font as font
from tkinter import ttk
import numpy as np
import json
from Bio import Phylo
from Bio.Phylo.PAML import codeml
from Bio.Phylo import PhyloXML as PX
from matplotlib.figure import Figure 
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
import sys
import os
import re
import matplotlib.patches as patches
import matplotlib.colors as mcolors
import matplotlib.gridspec as gridspec
import matplotlib.ticker as ticker
import matplotlib.collections as mpcollections
import copy
import string
import pandas as pd
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from matplotlib.widgets import Slider

def atoi(text):
    return int(text) if text.isdigit() else text

def natural_keys(text):
    return [ atoi(c) for c in re.split(r'(\d+)', text) ]

########################################################################################
## Pysage GUI
########################################################################################

##########################################################################
# factory to dynamically create the gui
class GUIFactory:
    factories = {}
    def add_factory(idx, gui_factory):
        GUIFactory.factories[idx] = gui_factory
    add_factory = staticmethod(add_factory)

    def create_gui(master, folder=None):
        return PysageGUI.Factory().create(master, folder=folder)
    create_gui = staticmethod(create_gui)


##########################################################################
# GUI main class
class PysageGUI(object):

    class Factory:
        def create(self, master, folder=None): return PysageGUI(master, folder=folder)

    ##########################################################################
    # standart class init
    def __init__(self, master, folder=None):

        self.master = master
        # Default folder is current directory
        folder = os.getcwd()
        if folder is not None:
            # Use passed folder
            self.folder = folder
        
        # Phylogenetic tree
        self.tree = None
        self.tree_backup = None # Backup for monomers' tree
        self.actual_root = None
        # HOR tree
        self.hor_tree = None
        self.hor_tree_backup = None # Backup for HORs' tree
        self.hor_root = None
        self.hor_dict = None
        self.hor_lengths = None
        # Family name counter
        self.fcnt = None
        # Chromosome sequence (start and end)
        self.chr_seq = None
        self.seq_name = None
        self.seq_len = 1.0
        # List of monomers to be plotted
        self.monomers = None
        self.monomer_colors = None
        # Location of the HORS to be plotted inside the whole sequence
        self.hors = None
        self.hor_locations = None
        self.hor_len = 0
        # Canvas where plots are shown
        self.canvas = None
        self.tree_canvas = None
        self.other_canvas = None
        # Data for visualization (collapse of some nodes and other info)
        self.clade_coords = {}
        self.clade_ids = {}
        self.clicked = None
        self.clicked_colors = None
        self.num_clicked = 0
        self.patches = None
        self.shown = False
        # Nodes to be collapsed and relative information
        self.clades_to_collapse = {}
        self.collapsed_clades = {}
        self.collapsed_indices = None
        self.collapsed_colors = None
        self.collapsed_patches = None
        # Zoom
        self.zoomed = False
        self.zoomed_nodes = []
        self.zoomed_coords = []
        # Object used to plot the highlight and scroll
        self.treeForSubPlot = None
        # Tree axis
        self.ax_tree = None
        self.ax_tree_ins = None
        self.axes_ins = None
        # Data to scroll
        self.scroll_y = 0.0
        # Colors to highlight branches
        self.colors = ['red', 'green', 'blue', 'orange', 'yellow', 'purple', 'grey', 'brown', 'cyan', 'magenta', 'pink', 'gold', 'salmon', 'lime', 'teal', 'silver', 'fuchsia', 'aqua', 'maroon', 'navy', 'olive', 'gray']#mcolors.TABLEAU_COLORS
        self.hor_colors = ['cyan', 'magenta', 'orange', 'purple', 'pink', 'yellow', 'brown', 'blue', 'green', 'red', 'lime', 'navy', 'gold', 'salmon']
        # Directory containing files (json, xml and others)
        self.filedir = os.getcwd() # The tool assumes that the file are located in current directory!!!
        self.filename = None
        # Start the GUI
        self.initialize()
            
    ##########################################################################
    # Method that generates a pop-up with an error message
    def popupMsg(self, msg):
        popup = tk.Toplevel(self.w)
        popup.wm_title("Error")
        popup.tkraise(self.w) # This just tells the message to be on top of the root window.
        tk.Label(popup, text=msg).pack(side="top", fill="x", pady=10)
        tk.Button(popup, text="Ok", command = popup.destroy).pack()
        
    ##########################################################################
    # Method computing a dictionary for each HOR containing monomers and locations
    def calcHorsMonomersList(self, old_names, new_names):
        all_clades = self.hor_tree.find_clades()
        # Build dict associating to each HOR the corresponding list of monomer(s)
        self.hor_dict = {}
        self.hor_lengths = {}
        cnt = 0
        for clade in all_clades:
            if cnt == 0:
                # Store the sequence start and end + the length
                if clade.sequences:
                    # Extract sequence start and end
                    loc_str = clade.sequences[0].location
                    substr = loc_str.split('[')[0]
                    substr2 = substr.split(':')[1]
                    seq_start = substr2.split('-')[0]
                    seq_end = substr2.split('-')[1]
                    # Store chromosome name
                    self.seq_name = substr.split(':')[0]
                    # Store sequence start and end
                    self.chr_seq = [seq_start, seq_end]
                    # Compute sequence length
                    self.seq_len = int(self.chr_seq[1]) - int(self.chr_seq[0])
                    cnt += 1
            if clade.properties:
                for prop in clade.properties:
                    monomers_seq = prop.value
                    try:
                        # Monomers in a HOR are separated by commas
                        monomers = monomers_seq.split(",")
                    except:
                        # Single monomer
                        monomers = monomers_seq
                    # Before storing data, we change name of monomers (and HORs)
                    new_monos = []
                    for mono in monomers:
                        idx = old_names.index(mono)
                        new_monos.append(new_names[idx])
                    # Change HOR name
                    hor_name = str(len(new_monos))
                    for mono in new_monos:
                        hor_name += mono
                    clade.name = hor_name
                    # Change clade properties
                    new_value = ""
                    nmonos = len(new_monos)
                    cnt = 0
                    for mono in new_monos:
                        new_value += mono
                        cnt += 1
                        if cnt < nmonos:
                            new_value += ","
                    prop.value = new_value
                    # Store monomers and locations for this HOR
                    self.hor_dict[clade.name] = [new_monos, [seq.location for seq in clade.sequences if clade.sequences is not None]]
                    self.hor_lengths[clade.name] = len(new_monos)
                    if len(new_monos) > self.hor_len:
                        self.hor_len = len(new_monos)
     
    ##########################################################################
    # Method that loads json and xml files with basename <filename>
    def loadFile(self, filename=None):
        self.fcnt = 1
        if filename is not None:
            res = Phylo.parse(filename, "phyloxml")
            trees = []
            for i, elem in enumerate(res):
                trees.append(elem)
            assert len(trees) == 2, "Invalid number of data contained in phyloXML file!!!"
            # Monomers
            self.tree = trees[0]
            self.tree_backup = copy.deepcopy(self.tree)
            self.tree.rooted = True
            # Look for the actual root of the tree
            for clade in self.tree.find_clades():
                try:
                    self.actual_root = clade.common_ancestor(self.tree.root)
                    break
                except:
                    pass
            assert self.actual_root is not None, "Monomers tree: root not found!!!"
            # Change name of clades
            old_names = []
            new_names = []
            for clade in self.tree.find_clades():
                 if clade.name and "chr" not in clade.name:
                     old_names.append(clade.name)
                     cname = "F" + str(self.fcnt)
                     new_names.append(cname)
                     # Change clade name
                     clade.name = cname
                     self.fcnt += 1

            assert len(old_names) == len(new_names), "Weird error in list append"
            # HORs
            self.hor_tree = trees[1]
            self.hot_tree_backup = copy.deepcopy(self.hor_tree)
            self.hor_tree.rooted = True
            # Look for the actual root of the tree
            for clade in self.hor_tree.find_clades():
                try:
                    self.hor_root = clade.common_ancestor(self.hor_tree.root)
                    break
                except:
                    pass
            assert self.hor_root is not None, "HORs tree: root not found!!!"
            # Get list of monomers for each hor (useful to print data)
            self.calcHorsMonomersList(old_names, new_names)
            # Save CSV file containing associations
            d = {"old_name": old_names, "new_name": new_names}
            df = pd.DataFrame(data=d)
            df.to_csv(os.path.join(self.folder, self.seq_name + '_name_association.csv'), sep=',', index=False) 
            # Save new XML file
            new_trees = PX.Phyloxml(phylogenies=[self.tree, self.hor_tree], attributes={'xsd': 'http://www.w3.org/2001/XMLSchema'})
            fname_no_ext = os.path.splitext(self.filename)[0]
            fname_no_ext += "_new.xml"
            outfile = os.path.join(self.folder, fname_no_ext)
            Phylo.write(new_trees, outfile, format='phyloxml')
        else:
            self.popupMsg("You must select the file before loading!!!")
            return
            
    ##########################################################################
    # Method that checks if a monomer belongs also to other HORs
    def checkMonomerInOtherHORs(self, mono, hor):
        isIn = False
        omono = None
        # Loop over clicked nodes
        for elem in self.clicked:
            if elem != hor:
                mono_and_locs = self.hor_dict[elem]
                # Loop over monomers of this HOR
                monomers = mono_and_locs[0]
                for cmono in monomers:
                    clades = self.tree.find_clades(cmono)
                    for clade in clades:
                        subclades = clade.find_clades()
                        for sclade in subclades:
                            if sclade.name and sclade.name == mono:
                                # The monomer has been found in another HOR
                                isIn = True
                                omono = cmono
                                break
            if isIn:
                break
        return isIn, omono
        
    ##########################################################################
    # Method that returns the HORs
    def extractHORs(self):
        # Loop over the HORs
        self.monomers = []
        self.hors = []
        self.locations = []
        self.monomer_colors = []
        elem = 0
        # Dict of checked monomers
        checked = {}
        for hor in self.clicked:
            # Extract monomers and locations
            mono_and_locs = self.hor_dict[hor]
            monomers = mono_and_locs[0]
            mono_locs = mono_and_locs[1]
            mono_colors = []
            new_monomers = []
            mod = False
            # Color clades associated to monomers
            for mono in monomers:
                # Check if current monomer is in more than one HOR
                found = True
                check_vals = None
                try:
                    # Check the existence of the monomer in the dict
                    check_vals = checked[mono]
                except:
                    found = False
                if found:
                    # Monomer already analyzed, get info
                    isInOtherHor = check_vals[0]
                    otherMono = check_vals[1]
                else:
                    # Monomer not yet analyzed, check if it belongs to another HOR
                    isInOtherHor, otherMono = self.checkMonomerInOtherHORs(mono, hor)
                    # Save the information
                    checked[mono] = [isInOtherHor, otherMono]
                if isInOtherHor:
                    # If the monomer belongs to another HOR, we append the parent monomer
                    new_monomers.append(otherMono)
                    mod = True
                else:
                    # Append current monomer
                    new_monomers.append(mono)
                clades = self.tree.find_clades(mono)
                for clade in clades:
                    if clade.color.to_hex() == "#000000":
                        clade.color = PX.BranchColor.from_name(self.colors[elem % len(self.colors)])
                        elem += 1
                    mono_colors.append(clade.color)
                    mono_clades = clade.find_clades()
                    for mono_clade in mono_clades:
                        # We use the color list
                        if mono_clade.color.to_hex() == "#000000":
                            mono_clade.color = clade.color
            # Extract rel_start, rel_end and strand for all the locations to be plotted
            locations = []
            for loc in mono_locs:
                substr = loc.split('[')[1]
                substr2 = substr.split(']')[0]
                substr3 = substr.split(']')[1]
                rel_start = substr2.split(':')[0]
                rel_end = substr2.split(':')[1]
                strand = substr3[1]
                locations.append([int(rel_start), int(rel_end), strand])
            # Sort locations
            locations.sort()
            # Insert information in the corresponding lists
            new_hor = hor
            if mod:
                # Modify HOR name for visualization
                hor_len = self.hor_lengths[hor]
                new_hor = str(hor_len)
                for mono in new_monomers:
                    new_hor += mono
                # Change the name of the HOR in clicked
                hor_id = self.clicked.index(hor)
                self.clicked[hor_id] = new_hor
                # Add new entry in the HOR dict
                self.hor_dict[new_hor] = [new_monomers, mono_locs]
            self.hors.append(new_hor)
            self.monomers.append(new_monomers)
            self.locations.append(locations)
            self.monomer_colors.append(mono_colors)
        
    ##########################################################################
    # Method that allows to click on the plot and do something
    def on_click(self, event):
        # Find the HOR within a range of 0.1 (the radius of the circle)
        found = False
        cid = -1
        ccoord = None
        clade_keys = self.clade_coords.keys()
        # Initialize lists of clicked nodes and colors
        if self.clicked is None:
            self.clicked = []
        if self.clicked_colors is None:
            self.clicked_colors = []
        # Loop over patches
        found = False
        for patch in self.patches:
            center = patch.center
            cx, cy = tuple(center)
            ox, oy = event.xdata, event.ydata
            d = math.sqrt(math.pow((cx - ox), 2.0) + math.pow((cy - oy), 2.0))
            click_all = False
            unclick_all = False
            if d <= 0.1:
                ccolor = None
                ckey = None
                for key in clade_keys:
                    coords = self.clade_coords[key]
                    if isinstance(coords, tuple):
                        if tuple(coords) == (cx, cy):
                            if key not in self.clicked:
                                ccolor = self.hor_colors[self.num_clicked % len(self.hor_colors)]
                                patch.set_color(ccolor)
                                self.clicked.append(key)
                                self.clicked_colors.append(ccolor)
                                self.num_clicked += 1
                            else:
                                # HOR already selected, unselect
                                # Before changing color, we must remove the right entry in the <clicked_colors> list
                                pcolor = patch.get_facecolor()
                                for strcolor in self.clicked_colors:
                                    clcolor = mcolors.to_rgba(strcolor)
                                    if clcolor == pcolor:
                                        self.clicked_colors.remove(strcolor)
                                        break
                                patch.set_color('black')
                                self.clicked.remove(key)
                    elif isinstance(coords, list):
                        for coord in coords:
                            if tuple(coord) == (cx, cy):
                                if key not in self.clicked:
                                    ckey = key
                                    ccolor = self.hor_colors[self.num_clicked % len(self.hor_colors)]
                                    patch.set_color(ccolor)
                                    self.clicked.append(key)
                                    self.num_clicked += 1
                                    click_all = True
                                else:
                                    # HOR already selected, unselect
                                    ckey = key
                                    # Before changing color, we must remove the right entry in the <clicked_colors> list
                                    pcolor = patch.get_facecolor()
                                    for strcolor in self.clicked_colors:
                                        clcolor = mcolors.to_rgba(strcolor)
                                        if clcolor == pcolor:
                                            self.clicked_colors.remove(strcolor)
                                            break
                                    patch.set_color('black')
                                    self.clicked.remove(key)
                                    unclick_all = True
                # Only for HORs appearing in different branches in the HOR tree we check if they have to be all clicked/unclicked
                if click_all:
                    all_coords = self.clade_coords[ckey]
                    for coord in all_coords:
                        ccx, ccy = tuple(coord)
                        if ccx != cx or ccy != cy:
                            # Found new coordinates for the HOR in the HOR tree
                            # Look for the patch
                            for mp in self.patches:
                                px, py = tuple(mp.center)
                                if ccx == px and ccy == py:
                                    mp.set_color(ccolor)
                if unclick_all:
                    all_coords = self.clade_coords[ckey]
                    for coord in all_coords:
                        ccx, ccy = tuple(coord)
                        if ccx != cx or ccy != cy:
                            # Found new coordinates for the HOR in the HOR tree
                            # Look for the patch
                            for mp in self.patches:
                                px, py = tuple(mp.center)
                                if ccx == px and ccy == py:
                                    mp.set_color('black')
                found = True
                break
                
        if found:
            self.canvas.draw()
            
    ##########################################################################
    def expandSubTree(self, event):
        found = False
        coords = self.collapsed_clades.keys()
        elem = None
        mindist = 9999.0
        for coord in coords:
            cx, cy = tuple(coord)
            ox, oy = event.xdata, event.ydata
            d = math.sqrt(math.pow((cx - ox), 2.0) + math.pow((cy - oy), 2.0))
            if d < mindist:
                mindist = d
                elem = coord
        if mindist <= 50.0: # It is constant now, maybe it could be changed...
            # Found patch
            found = True
        if found:
            # Append the element to the list of nodes to zoom
            self.zoomed_nodes.append(self.collapsed_clades[elem])
            self.zoomed_coords.append(elem)
            
    ##########################################################################    
    # Zoom in
    def zoomIn(self):
        if len(self.zoomed_nodes) == 0:
            self.popupMsg("You must click a node before zooming!!!")
            return
        # Clear figure
        self.ax_tree.clear()
        
        self.ax_tree.get_xaxis().set_visible(False)
        self.ax_tree.get_yaxis().set_visible(False)
        
        node_colors = []
        start_x = 0.0
        start_y = 0.0
        width = 1.0
        height = 1.0 / len(self.zoomed_nodes)
        self.ax_tree_ins = []
        for elem in self.zoomed_nodes:
            node_colors.append(self.collapsed_colors[elem.root])
            
        # Create two sub-plots, one for the highlight and one for the scrollbar
        self.ax_tree_ins = []
        ax_tree_ins = self.ax_tree.inset_axes([0.0, 0.0, 0.95, 1.0])
        self.ax_tree_ins.append(ax_tree_ins)
        ax_slider_ins = self.ax_tree.inset_axes([0.95, 0.0, 0.05, 1.0])
        self.ax_tree_ins.append(ax_slider_ins)
            
        # Copy of original tree
        self.treeForSubPlot = copy.deepcopy(self.tree)
        all_clades = self.treeForSubPlot.find_clades()
        found = False
        # Look for the clades to be colored, all others will be painted in black
        for clade in all_clades:
            toColor = False
            found = False
            for node_color in node_colors:
                if clade.color.to_hex() == node_color.to_hex():
                    found = True
                    break
            if not found:
                if clade.color.to_hex() != "#000000":
                    toColor = True
            if toColor:
                clade.color = PX.BranchColor.from_name('black')
            
        # Create list of nodes to expand
        nodes_to_highlight = []
        nodes_to_collapse = []
        coords_to_collapse = []
        keys = list(self.collapsed_clades.keys())
        for elem in keys:
            cclade = self.collapsed_clades[elem]
            if cclade not in self.zoomed_nodes:
                nodes_to_collapse.append(cclade)
                coords_to_collapse.append(elem)
            else:
                nodes_to_highlight.append(cclade)

        # Compute the min and max values for y
        ymax = self.treeForSubPlot.count_terminals() + 0.8 # Derived from original draw() method
        ymin = 0.2 # Derived from original draw() method
        # Compute differece between max and min
        ydiff = ymax - ymin
        # Compute the amount of plot to be shown
        yrange = ydiff / 100.0
        
        # Choose the Slider color
        slider_color = 'White'
        # Create the slider
        slider_position = Slider(ax_slider_ins, label='', valmin=0.0, valmax=1.0, valinit=1.0, orientation="vertical")

        # Draw part of the tree from ymin to ymax
        self.drawTreePart(self.treeForSubPlot, nodes_to_highlight=nodes_to_highlight, nodes_to_collapse=nodes_to_collapse, coords_to_collapse=coords_to_collapse, y_min=ymin, y_max=ymin + yrange, axes=ax_tree_ins)
        
        # update() function to change the graph when the slider is in use
        def update(val):
            pos = slider_position.val
            
            # Round slider position to second decimal digit
            posy = round(pos, 2)
            
            # Compute the delta y
            dy = 1.0 - posy
            # Compute new plot y ranges (minimum and maximum)
            bottom_y = 0.2 + ydiff * dy
            top_y = bottom_y + yrange
            if top_y > ymax:
                top_y = ymax
                bottom_y = ymax - yrange
            # Draw updated part of the tree from ymin to ymax
            self.drawTreePart(self.treeForSubPlot, nodes_to_highlight=nodes_to_highlight, nodes_to_collapse=nodes_to_collapse, coords_to_collapse=coords_to_collapse, y_min=bottom_y, y_max=top_y, axes=ax_tree_ins)
        
        # update function called using on_changed() function
        slider_position.on_changed(update)
        # make slider text not visible
        slider_position.valtext.set_visible(False)
            
        self.tree_canvas.draw()
            
        # Set zoom flag to true
        self.zoomed = True

    ##########################################################################
    # Copy of Phylo.draw() used to store coordinates of clades
    def drawTree(self, 
        tree,
        label_func=str,
        do_show=True,
        show_confidence=True,
        # For power users
        axes=None,
        branch_labels=None,
        label_colors=None,
        *args,
        **kwargs,
    ):
    
        # Dict: clade name (key) and x,y (value)
        self.clade_coords = {}
    
        # Arrays that store lines for the plot of clades
        horizontal_linecollections = []
        vertical_linecollections = []

        # Options for displaying branch labels / confidence
        def conf2str(conf):
            if int(conf) == conf:
                return str(int(conf))
            return str(conf)

        if not branch_labels:
            if show_confidence:

                def format_branch_label(clade):
                    try:
                        confidences = clade.confidencesPhylo.draw
                        # phyloXML supports multiple confidences
                    except AttributeError:
                        pass
                    else:
                        return "/".join(conf2str(cnf.value) for cnf in confidences)
                    if clade.confidence is not None:
                        return conf2str(clade.confidence)
                    return None

            else:

                def format_branch_label(clade):
                    return None

        elif isinstance(branch_labels, dict):

            def format_branch_label(clade):
                return branch_labels.get(clade)

        else:
            if not callable(branch_labels):
                raise TypeError(
                    "branch_labels must be either a dict or a callable (function)"
                )
            format_branch_label = branch_labels

        # options for displaying label colors.
        if label_colors:
            if callable(label_colors):

                def get_label_color(label):
                    return label_colors(label)

            else:
                # label_colors is presumed to be a dictPhylo.draw(self.tree, axes=ax_tree)
                def get_label_color(label):
                    return label_colors.get(label, "black")

        else:

            def get_label_color(label):
                # if label_colors is not specified, use black
                return "black"

        # Layout

        def get_x_positions(tree):
            depths = tree.depths()
            # If there are no branch lengths, assume unit branch lengths
            if not max(depths.values()):
                depths = tree.depths(unit_branch_lengths=True)
            return depths

        def get_y_positions(tree):
            maxheight = tree.count_terminals()
            # Rows are defined by the tips
            heights = {tip: maxheight - i for i, tip in enumerate(reversed(tree.get_terminals()))}

            # Internal nodes: place at midpoint of children
            def calc_row(clade):
                for subclade in clade:
                    if subclade not in heights:
                        calc_row(subclade)
                # Closure over heights
                heights[clade] = (heights[clade.clades[0]] + heights[clade.clades[-1]]) / 2.0

            if tree.root.clades:
                calc_row(tree.root)
            return heights

        x_posns = get_x_positions(tree)
        y_posns = get_y_positions(tree)
        # The function draw_clade closes over the axes object
        if axes is None:
            fig = plt.figure()
            axes = fig.add_subplot(1, 1, 1)
        elif not isinstance(axes, plt.matplotlib.axes.Axes):
            raise ValueError(f"Invalid argument for axes: {axes}")

        def draw_clade_lines(
            use_linecollection=False,
            orientation="horizontal",
            y_here=0,
            x_start=0,
            x_here=0,
            y_bot=0,
            y_top=0,
            color="black",
            lw=".1",
        ):

            if not use_linecollection and orientation == "horizontal":
                axes.hlines(y_here, x_start, x_here, color=color, lw=lw)
            elif use_linecollection and orientation == "horizontal":
                horizontal_linecollections.append(
                    mpcollections.LineCollection(
                        [[(x_start, y_here), (x_here, y_here)]], color=color, lw=lw
                    )
                )
            elif not use_linecollection and orientation == "vertical":
                axes.vlines(x_here, y_bot, y_top, color=color)
            elif use_linecollection and orientation == "vertical":
                vertical_linecollections.append(
                    mpcollections.LineCollection(
                        [[(x_here, y_bot), (x_here, y_top)]], color=color, lw=lw
                    )
                )

        def draw_clade(clade, x_start, color, lw):
            """Recursively draw a tree, down from the given clade."""
            x_here = x_posns[clade]
            y_here = y_posns[clade]
            if clade.name:
                if clade.name not in self.clade_coords:
                    self.clade_coords[clade.name] = (x_here, y_here)
                else:
                    val = self.clade_coords[clade.name]
                    if isinstance(val, tuple):
                        cx, cy = tuple(val)
                        self.clade_coords[clade.name] = [(cx, cy), (x_here, y_here)]
                    elif isinstance(val, list):
                        val.append((x_here, y_here))
                        self.clade_coords[clade.name] = val
                
            # phyloXML-only graphics annotations
            if hasattr(clade, "color") and clade.color is not None:
                color = clade.color.to_hex()
            if hasattr(clade, "width") and clade.width is not None:
                lw = clade.width * plt.rcParams["lines.linewidth"]
            # Draw a horizontal line from start to here
            draw_clade_lines(
                use_linecollection=True,
                orientation="horizontal",
                y_here=y_here,
                x_start=x_start,
                x_here=x_here,
                color=color,
                lw=lw,
            )
            if clade.clades:
                # Draw a vertical line connecting all children
                y_top = y_posns[clade.clades[0]]
                y_bot = y_posns[clade.clades[-1]]
                # Only apply widths to horizontal lines, like Archaeopteryx
                draw_clade_lines(
                    use_linecollection=True,
                    orientation="vertical",
                    x_here=x_here,
                    y_bot=y_bot,
                    y_top=y_top,
                    color=color,
                    lw=lw,
                )
                # Draw descendents
                for child in clade:
                    draw_clade(child, x_here, color, lw)

        draw_clade(tree.root, 0, "k", plt.rcParams["lines.linewidth"])

        # If line collections were used to create clade lines, here they are added
        # to the pyplot plot.
        for i in horizontal_linecollections:
            axes.add_collection(i)
        for i in vertical_linecollections:
            axes.add_collection(i)

        # Aesthetics

        try:
            name = tree.name
        except AttributeError:
            pass
        else:
            if name:
                axes.set_title(name)
        axes.set_xlabel("branch length")
        axes.set_ylabel("taxa")
        # Add margins around the tree to prevent overlapping the axes
        xmax = max(x_posns.values())
        axes.set_xlim(-0.05 * xmax, 1.25 * xmax)
        # Also invert the y-axis (origin at the top)
        # Add a small vertical margin, but avoid including 0 and N+1 on the y axis
        axes.set_ylim(max(y_posns.values()) + 0.8, 0.2)

        # Parse and process key word arguments as pyplot options
        for key, value in kwargs.items():
            try:
                # Check that the pyplot option input is iterable, as required
                list(value)
            except TypeError:
                raise ValueError(
                    'Keyword argument "%s=%s" is not in the format '
                    "pyplot_option_name=(tuple), pyplot_option_name=(tuple, dict),"
                    " or pyplot_option_name=(dict) " % (key, value)
                ) from None
            if isinstance(value, dict):
                getattr(plt, str(key))(**dict(value))
            elif not (isinstance(value[0], tuple)):
                getattr(plt, str(key))(*value)
            elif isinstance(value[0], tuple):
                getattr(plt, str(key))(*value[0], **dict(value[1]))

        if do_show:
            plt.show()
        
    ##########################################################################
    # Copy of Phylo.draw() used to make some clades collapsing based on the selected HORs and the number of sub-clades
    def drawCollapsedTree(self, 
        tree,
        label_func=str,
        do_show=True,
        show_confidence=True,
        # For power users
        axes=None,
        branch_labels=None,
        label_colors=None,
        *args,
        **kwargs,
    ):
    
        # Collapsed_clades is a dict: key: coordinate of the patch; value: the phylogenetic tree with the collased clade as root
        self.clades_to_collapse = {}
        self.collapsed_clades = {}
        self.collapsed_indices = {}
        self.collapsed_colors = {}
        self.collapsed_patches = {}
    
        # Arrays that store lines for the plot of clades
        horizontal_linecollections = []
        vertical_linecollections = []

        # Options for displaying branch labels / confidence
        def conf2str(conf):
            if int(conf) == conf:
                return str(int(conf))
            return str(conf)

        if not branch_labels:
            if show_confidence:

                def format_branch_label(clade):
                    try:
                        confidences = clade.confidences
                        # phyloXML supports multiple confidences
                    except AttributeError:
                        pass
                    else:
                        return "/".join(conf2str(cnf.value) for cnf in confidences)
                    if clade.confidence is not None:
                        return conf2str(clade.confidence)
                    return None

            else:

                def format_branch_label(clade):
                    return None

        elif isinstance(branch_labels, dict):

            def format_branch_label(clade):
                return branch_labels.get(clade)

        else:
            if not callable(branch_labels):
                raise TypeError(
                    "branch_labels must be either a dict or a callable (function)"
                )
            format_branch_label = branch_labels

        # options for displaying label colors.
        if label_colors:
            if callable(label_colors):

                def get_label_color(label):
                    return label_colors(label)

            else:
                # label_colors is presumed to be a dictPhylo.draw(self.tree, axes=ax_tree)
                def get_label_color(label):
                    return label_colors.get(label, "black")

        else:

            def get_label_color(label):
                # if label_colors is not specified, use black
                return "black"
                
        # Check clades to collapse
        def checkCollapsedClades(clade):
            collapsed = False
            check = False
            found = False
            if clade.name:
                # Check whether clade in at least one of the HORs
                for monomer in self.monomers:
                    if clade.name in monomer:
                        found = True
                        break
                if not found:
                    # Clade is not contained in any of the HORs -> check whether at least one subclade is contained in the HORs
                    subclades = clade.find_clades()
                    for sclade in subclades:
                        for monomer in self.monomers:
                            if sclade.name in monomer:
                                found = True
                                break
                    if not found:
                        # If neither the clade is in the HOR list, nor any of its subclades is, we must check the size of this portion of the tree
                        check = True
            else:
                # Check whether at least one subclade is contained in the HORs
                subclades = clade.find_clades()
                for sclade in subclades:
                    for monomer in self.monomers:
                        if sclade.name in monomer:
                            found = True
                            break
                if not found:
                    # If none of the subclades is in the list, we must check the size of this portion of the tree
                    check = True
            if check:
                # Check sub-clades
                nsubclades = sum([1 for elem in clade.find_clades()])
                if nsubclades >= 10: # It is constant, maybe it could be changed...
                    collapsed = True
            self.clades_to_collapse[clade] = [collapsed, None]
            if collapsed:
                # Collapse all sub-clades
                for elem in clade.find_clades():
                    if elem not in self.clades_to_collapse:
                        self.clades_to_collapse[elem] = [True, clade]

        # Dict containing a flag for each clade indicating whether or not it will be collapsed
        for clade in tree.find_clades():
            if clade not in self.clades_to_collapse:
                checkCollapsedClades(clade)

        # Layout
        
        def get_xy_positions(tree):
            x = {}
            y = {}
            # Compute the list of visible clades and the list of visible collapsed clades
            visible_collapsed = []
            visible = []
            for elem in self.clades_to_collapse:
                flag, parent = self.clades_to_collapse[elem]
                if parent is None:
                    visible.append(elem)
                    if flag:
                        visible_collapsed.append(elem)
            
            depths = tree.depths()
            # If there are no branch lengths, assume unit branch lengths
            if not max(depths.values()):
                depths = tree.depths(unit_branch_lengths=True)
            for elem in depths:
                if elem in visible:
                    x[elem] = depths[elem]
            
            # Compute the list of visible leaves in the tree
            visible_leaves = []
            leaves = tree.get_terminals()
            for leaf in leaves:
                if leaf in visible:
                    visible_leaves.append(leaf)
            
            # Height of the collapsed tree depends on the number of visible leaves and visible collapsed clades
            maxheight = len(visible_leaves) + len(visible_collapsed)
            visible_nodes = []
            all_clades = tree.find_clades()
            for elem in all_clades:
                if elem in visible_leaves or elem in visible_collapsed:
                    visible_nodes.append(elem)
                    
            # Rows are defined by the tips
            heights = {tip: maxheight - i for i, tip in enumerate(reversed(visible_nodes))}#tree.get_terminals()))}
            
            # Compute visible sub-clades for each visible clade
            visible_subclades = {}
            for elem in visible:
                subclades = []
                for sclade in elem.clades:
                    if sclade in visible:
                        subclades.append(sclade)
                visible_subclades[elem] = subclades

            # Internal nodes: place at midpoint of children
            def calc_row(clade):
                try:
                    subclades = visible_subclades[clade]
                except:
                    clade_leaves = clade.get_terminals()
                    subclades = []
                    for leaf in clade_leaves:
                        if leaf in visible_nodes:
                            subclades.append(leaf)
                for subclade in subclades:
                    if subclade not in heights:
                        calc_row(subclade)
                # Closure over heights
                if len(subclades) > 1:
                    heights[clade] = (heights[subclades[0]] + heights[subclades[-1]]) / 2.0
                else:
                    heights[clade] = heights[subclades[0]]
                
            calc_row(tree.root)
            
            y = heights

            return x, y, visible

        x_posns, y_posns, visible_clades = get_xy_positions(tree)
        # The function draw_clade closes over the axes object
        if axes is None:
            fig = plt.figure()
            axes = fig.add_subplot(1, 1, 1)
        elif not isinstance(axes, plt.matplotlib.axes.Axes):
            raise ValueError(f"Invalid argument for axes: {axes}")

        def draw_clade_lines(
            use_linecollection=False,
            orientation="horizontal",
            y_here=0,
            x_start=0,
            x_here=0,
            y_bot=0,
            y_top=0,
            color="black",
            lw=".1",
        ):

            if not use_linecollection and orientation == "horizontal":
                axes.hlines(y_here, x_start, x_here, color=color, lw=lw)
            elif use_linecollection and orientation == "horizontal":
                horizontal_linecollections.append(
                    mpcollections.LineCollection(
                        [[(x_start, y_here), (x_here, y_here)]], color=color, lw=lw
                    )
                )
            elif not use_linecollection and orientation == "vertical":
                axes.vlines(x_here, y_bot, y_top, color=color)
            elif use_linecollection and orientation == "vertical":
                vertical_linecollections.append(
                    mpcollections.LineCollection(
                        [[(x_here, y_bot), (x_here, y_top)]], color=color, lw=lw
                    )
                )
                        
        def draw_clade(clade, x_start, color, lw):
            """Recursively draw a tree, down from the given clade."""
            x_here = x_posns[clade]
            y_here = y_posns[clade]
            draw = True
            flag, _ = self.clades_to_collapse[clade]
            if flag:
                # Add a circle (clickable)
                circle = plt.Circle((x_here, y_here), 0.25, color=clade.color.to_hex())
                axes.add_patch(circle)
                # Draw a horizontal line from start to here
                draw_clade_lines(
                    use_linecollection=True,
                    orientation="horizontal",
                    y_here=y_here,
                    x_start=x_start,
                    x_here=x_here,
                    color=clade.color.to_hex(),
                    lw=lw,
                )
                # Add the "new" tree to the dict
                self.collapsed_clades[(x_here, y_here)] = PX.Phylogeny(root=clade, name=clade.name)
                idx = None
                cnt = 0
                for elem in tree.find_clades():
                    if elem == clade:
                        idx = cnt
                        break
                    cnt += 1
                self.collapsed_indices[clade] = idx
                self.collapsed_colors[clade] = clade.color
                draw = False
            if draw:
                # phyloXML-only graphics annotations
                if hasattr(clade, "color") and clade.color is not None:
                    color = clade.color.to_hex()
                if hasattr(clade, "width") and clade.width is not None:
                    lw = clade.width * plt.rcParams["lines.linewidth"]
                # Draw a horizontal line from start to here
                draw_clade_lines(
                    use_linecollection=True,
                    orientation="horizontal",
                    y_here=y_here,
                    x_start=x_start,
                    x_here=x_here,
                    color=color,
                    lw=lw,
                )
                if clade.clades:
                    subclades = []
                    for subclade in clade.clades:
                        if subclade in visible_clades:
                            subclades.append(subclade)
                    # Draw a vertical line connecting all children
                    y_top = y_posns[subclades[0]]
                    y_bot = y_posns[subclades[-1]]
                    # Only apply widths to horizontal lines, like Archaeopteryx
                    draw_clade_lines(
                        use_linecollection=True,
                        orientation="vertical",
                        x_here=x_here,
                        y_bot=y_bot,
                        y_top=y_top,
                        color=color,
                        lw=lw,
                    )
                    # Draw descendents
                    for child in subclades:
                        draw_clade(child, x_here, color, lw)

        draw_clade(tree.root, 0, "k", plt.rcParams["lines.linewidth"])

        # If line collections were used to create clade lines, here they are added
        # to the pyplot plot.
        for i in horizontal_linecollections:
            axes.add_collection(i)
        for i in vertical_linecollections:
            axes.add_collection(i)

        # Aesthetics

        try:
            name = tree.name
        except AttributeError:
            pass
        else:
            if name:
                axes.set_title(name)
        axes.set_xlabel("branch length")
        axes.set_ylabel("taxa")
        # Add margins around the tree to prevent overlapping the axes
        xmax = max(x_posns.values())
        axes.set_xlim(-0.05 * xmax, 1.25 * xmax)
        # Also invert the y-axis (origin at the top)
        # Add a small vertical margin, but avoid including 0 and N+1 on the y axis
        axes.set_ylim(max(y_posns.values()) + 0.8, 0.2)

        # Parse and process key word arguments as pyplot options
        for key, value in kwargs.items():
            try:
                # Check that the pyplot option input is iterable, as required
                list(value)
            except TypeError:
                raise ValueError(
                    'Keyword argument "%s=%s" is not in the format '
                    "pyplot_option_name=(tuple), pyplot_option_name=(tuple, dict),"
                    " or pyplot_option_name=(dict) " % (key, value)
                ) from None
            if isinstance(value, dict):
                getattr(plt, str(key))(**dict(value))
            elif not (isinstance(value[0], tuple)):
                getattr(plt, str(key))(*value)
            elif isinstance(value[0], tuple):
                getattr(plt, str(key))(*value[0], **dict(value[1]))

        if do_show:
            plt.show()
            
        self.collapsed_patches = axes.patches
            
    ##########################################################################
    # Copy of Phylo.draw() used to show part of the tree (from y_max to y_min) and to expand some nodes, while others are still collapsed
    def drawTreePart(self, 
        tree,
        nodes_to_highlight = None,
        nodes_to_collapse = None,
        coords_to_collapse = None,
        y_min = None,
        y_max = None,
        label_func=str,
        do_show=True,
        show_confidence=True,
        # For power users
        axes=None,
        branch_labels=None,
        label_colors=None,
        *args,
        **kwargs,
    ):
    
        # Collapsed_clades is a dict: key: coordinate of the patch; value: the phylogenetic tree with the collased clade as root
        clades_to_collapse = {}
    
        # Arrays that store lines for the plot of clades
        horizontal_linecollections = []
        vertical_linecollections = []

        # Options for displaying branch labels / confidence
        def conf2str(conf):
            if int(conf) == conf:
                return str(int(conf))
            return str(conf)

        if not branch_labels:
            if show_confidence:

                def format_branch_label(clade):
                    try:
                        confidences = clade.confidences
                        # phyloXML supports multiple confidences
                    except AttributeError:
                        pass
                    else:
                        return "/".join(conf2str(cnf.value) for cnf in confidences)
                    if clade.confidence is not None:
                        return conf2str(clade.confidence)
                    return None

            else:

                def format_branch_label(clade):
                    return None

        elif isinstance(branch_labels, dict):

            def format_branch_label(clade):
                return branch_labels.get(clade)

        else:
            if not callable(branch_labels):
                raise TypeError(
                    "branch_labels must be either a dict or a callable (function)"
                )
            format_branch_label = branch_labels

        # options for displaying label colors.
        if label_colors:
            if callable(label_colors):

                def get_label_color(label):
                    return label_colors(label)

            else:
                # label_colors is presumed to be a dictPhylo.draw(self.tree, axes=ax_tree)
                def get_label_color(label):
                    return label_colors.get(label, "black")

        else:

            def get_label_color(label):
                # if label_colors is not specified, use black
                return "black"
                
        def getCladesToCollapse():
            cladesToCollapse = []
            cnt = 0
            for clade in tree.find_clades():
                for elem in nodes_to_collapse:
                    idx = self.collapsed_indices[elem.root]
                    if cnt == idx:
                        cladesToCollapse.append(clade)
                        break
                cnt += 1
            return cladesToCollapse
            
        def getCladesToHighlight():
            cladesToHighlight = []
            cnt = 0
            for clade in tree.find_clades():
                for elem in nodes_to_highlight:
                    idx = self.collapsed_indices[elem.root]
                    if cnt == idx:
                        cladesToHighlight.append(clade)
                        break
                cnt += 1
            return cladesToHighlight
            
        cladesToCollapse = getCladesToCollapse()
        cladesToHighlight = getCladesToHighlight()
                  
        # Check clades to collapse
        def checkCollapsedClades(clade):
            collapsed = False
            check = False
            found = False
            if clade in cladesToCollapse:
                collapsed = True
            clades_to_collapse[clade] = [collapsed, None]
            if collapsed:
                # Collapse all sub-clades
                for elem in clade.find_clades():
                    if elem not in clades_to_collapse:
                        clades_to_collapse[elem] = [True, clade]

        # Dict containing a flag for each clade indicating whether or not it will be collapsed
        cnt = 0
        for clade in tree.find_clades():
            if clade not in clades_to_collapse:
                checkCollapsedClades(clade)
            cnt += 1

        # Layout
        
        def get_xy_positions(tree):
            x = {}
            y = {}
            # Compute the list of visible clades and the list of visible collapsed clades
            visible_collapsed = []
            visible = []
            for elem in clades_to_collapse:
                flag, parent = clades_to_collapse[elem]
                if parent is None:
                    visible.append(elem)
                    if flag:
                        visible_collapsed.append(elem)
            
            depths = tree.depths()
            # If there are no branch lengths, assume unit branch lengths
            if not max(depths.values()):
                depths = tree.depths(unit_branch_lengths=True)
            for elem in depths:
                if elem in visible:
                    x[elem] = depths[elem]
            
            # Compute the list of visible leaves in the tree
            visible_leaves = []
            leaves = tree.get_terminals()
            for leaf in leaves:
                if leaf in visible:
                    visible_leaves.append(leaf)
            
            # Height of the collapsed tree depends on the number of visible leaves and visible collapsed clades
            maxheight = len(visible_leaves) + len(visible_collapsed)
            visible_nodes = []
            all_clades = tree.find_clades()
            for elem in all_clades:
                if elem in visible_leaves or elem in visible_collapsed:
                    visible_nodes.append(elem)
                    
            # Rows are defined by the tips
            heights = {tip: maxheight - i for i, tip in enumerate(reversed(visible_nodes))}#tree.get_terminals()))}
            
            # Compute visible sub-clades for each visible clade
            visible_subclades = {}
            for elem in visible:
                subclades = []
                for sclade in elem.clades:
                    if sclade in visible:
                        subclades.append(sclade)
                visible_subclades[elem] = subclades

            # Internal nodes: place at midpoint of children
            def calc_row(clade):
                try:
                    subclades = visible_subclades[clade]
                except:
                    clade_leaves = clade.get_terminals()
                    subclades = []
                    for leaf in clade_leaves:
                        if leaf in visible_nodes:
                            subclades.append(leaf)
                for subclade in subclades:
                    if subclade not in heights:
                        calc_row(subclade)
                # Closure over heights
                if len(subclades) > 1:
                    heights[clade] = (heights[subclades[0]] + heights[subclades[-1]]) / 2.0
                else:
                    heights[clade] = heights[subclades[0]]
                
            calc_row(tree.root)
            
            y = heights

            return x, y, visible, maxheight
            
        x_posns, y_posns, visible_clades, max_height = get_xy_positions(tree)
        
        # Now compute depth and height of all sub-trees to highlight
        subtreesToHighlight = {}
        for clade in cladesToHighlight:
            xmax = -1e10
            subtree = PX.Phylogeny(root=clade, name=clade.name)
            xpos, _, _, _ = get_xy_positions(subtree)
            for xelem in xpos:
                x = xpos[xelem]
                if x > xmax:
                    xmax = x
            ymin = 1e10
            ymax = -1e10
            subclades = clade.find_clades()
            for subclade in subclades:
                ypos = y_posns[subclade]
                if ypos < ymin:
                    ymin = ypos
                if ypos > ymax:
                    ymax = ypos
            subtreesToHighlight[clade] = (xmax, ymin, ymax)

        # The function draw_clade closes over the axes object
        if axes is None:
            fig = plt.figure()
            axes = fig.add_subplot(1, 1, 1)
        elif not isinstance(axes, plt.matplotlib.axes.Axes):
            raise ValueError(f"Invalid argument for axes: {axes}")

        def draw_clade_lines(
            use_linecollection=False,
            orientation="horizontal",
            y_here=0,
            x_start=0,
            x_here=0,
            y_bot=0,
            y_top=0,
            color="black",
            lw=".1",
        ):

            if not use_linecollection and orientation == "horizontal":
                axes.hlines(y_here, x_start, x_here, color=color, lw=lw)
            elif use_linecollection and orientation == "horizontal":
                horizontal_linecollections.append(
                    mpcollections.LineCollection(
                        [[(x_start, y_here), (x_here, y_here)]], color=color, lw=lw
                    )
                )
            elif not use_linecollection and orientation == "vertical":
                axes.vlines(x_here, y_bot, y_top, color=color)
            elif use_linecollection and orientation == "vertical":
                vertical_linecollections.append(
                    mpcollections.LineCollection(
                        [[(x_here, y_bot), (x_here, y_top)]], color=color, lw=lw
                    )
                )

        def draw_clade(clade, x_start, color, lw):
            """Recursively draw a tree, down from the given clade."""
            x_here = x_posns[clade]
            y_here = y_posns[clade]
            draw = True
            flag, _ = clades_to_collapse[clade]
            if flag:
                # Add a circle (clickable)
                circle = plt.Circle((x_here, y_here), 0.25, color=clade.color.to_hex())
                axes.add_patch(circle)
                # Draw a horwidthizontal line from start to here
                draw_clade_lines(
                    use_linecollection=True,
                    orientation="horizontal",
                    y_here=y_here,
                    x_start=x_start,
                    x_here=x_here,
                    color=clade.color.to_hex(),
                    lw=lw,
                )
                draw = False
            if draw:
                # phyloXML-only graphics annotations
                if hasattr(clade, "color") and clade.color is not None:
                    color = clade.color.to_hex()
                if hasattr(clade, "width") and clade.width is not None:
                    lw = clade.width * plt.rcParams["lines.linewidth"]
                # Draw a horizontal line from start to here
                draw_clade_lines(
                    use_linecollection=True,
                    orientation="horizontal",
                    y_here=y_here,
                    x_start=x_start,
                    x_here=x_here,
                    color=color,
                    lw=lw,
                )
                if clade.clades:
                    subclades = []
                    for subclade in clade.clades:
                        if subclade in visible_clades:
                            subclades.append(subclade)
                    # Draw a vertical line connecting all children
                    y_top = y_posns[subclades[0]]
                    y_bot = y_posns[subclades[-1]]
                    # Only apply widths to horizontal lines, like Archaeopteryx
                    draw_clade_lines(
                        use_linecollection=True,
                        orientation="vertical",
                        x_here=x_here,
                        y_bot=y_bot,
                        y_top=y_top,
                        color=color,
                        lw=lw,
                    )
                    # Draw descendents
                    for child in subclades:
                        draw_clade(child, x_here, color, lw)
                        
                if clade in subtreesToHighlight:
                    sizes = subtreesToHighlight[clade]
                    width, height_min, height_max = tuple(sizes)
                    rect = patches.Rectangle((x_start, height_min), width, height_max - height_min, color='#D3D3D3') # +1 helps with colored branches???
                    axes.add_patch(rect)

        draw_clade(tree.root, 0, "k", plt.rcParams["lines.linewidth"])

        # If line collections were used to create clade lines, here they are added
        # to the pyplot plot.
        for i in horizontal_linecollections:
            axes.add_collection(i)
        for i in vertical_linecollections:
            axes.add_collection(i)

        # Aesthetics

        try:
            name = tree.name
        except AttributeError:
            pass
        else:
            if name:
                axes.set_title(name)
        axes.set_xlabel("branch length")
        axes.set_ylabel("taxa")
        # Add margins around the tree to prevent overlapping the axes
        xmax = max(x_posns.values())
        axes.set_xlim(-0.05 * xmax, 1.25 * xmax)
        # Also invert the y-axis (origin at the top)
        # Add a small vertical margin, but avoid including 0 and N+1 on the y axis
        # Check difference between max height of the tree and passed y_max
        y_range = y_min - y_max
        if y_max > max_height or abs(max_height - y_max) < 100.0:
            y_max = max_height + 0.8
        # Check whether y_min is higher than max_height
        if y_min >= max_height:
            y_min = max_height - 150.0 # Currently constant, this can be changed...
        axes.set_ylim(y_max, y_min)

        # Parse and process key word arguments as pyplot options
        for key, value in kwargs.items():
            try:
                # Check that the pyplot option input is iterable, as required
                list(value)
            except TypeError:
                raise ValueError(
                    'Keyword argument "%s=%s" is not in the format '
                    "pyplot_option_name=(tuple), pyplot_option_name=(tuple, dict),"
                    " or pyplot_option_name=(dict) " % (key, value)
                ) from None
            if isinstance(value, dict):
                getattr(plt, str(key))(**dict(value))
            elif not (isinstance(value[0], tuple)):
                getattr(plt, str(key))(*value)
            elif isinstance(value[0], tuple):
                getattr(plt, str(key))(*value[0], **dict(value[1]))

        if do_show:
            plt.show()

    ##########################################################################
    # GUI plotTree: plot tree (detailed)
    def plotTree(self):
        if self.hor_tree is None:
            self.popupMsg("You must load the phyloXML file before plotting!!!")
            return
        if self.canvas is not None:
            self.canvas.get_tk_widget().destroy()
            
        # Draw HOR tree
        fig = Figure(figsize = (6, 6), dpi = 100, constrained_layout=True)
        fig.canvas.mpl_connect('button_press_event', self.on_click)
        matplotlib.rc('font', size=6)
        ax = fig.add_subplot(1, 1, 1)
        
        self.drawTree(self.hor_tree, axes=ax)
        self.clade_ids = {}
        self.clicked = []
        self.clicked_colors = []
        self.num_clicked = 0 # Index to iterate over color palette (it can be only increased or decreased)
        # Now we add clickable circles where HORs are placed
        cid = 0
        clade_keys = self.clade_coords.keys()
        for key in clade_keys:
            coords = self.clade_coords[key]
            if isinstance(coords, tuple):
                # Single location
                x, y = tuple(coords)
                circle = plt.Circle((x, y), 0.1, color='black', zorder=5)
                ax.add_patch(circle)
                self.clade_ids[key] = cid
                cid += 1
            elif isinstance(coords, list):
                cids = []
                for coord in coords:
                    x, y = tuple(coord)
                    circle = plt.Circle((x, y), 0.1, color='black', zorder=5)
                    ax.add_patch(circle)
                    cids.append(cid)
                    cid += 1
                self.clade_ids[key] = cids
                
        self.patches = ax.patches
        self.canvas = FigureCanvasTkAgg(fig, master=self.v)
        self.canvas.draw()
        
        # placing the canvas on the Tkinter window 
        self.canvas.get_tk_widget().pack() 
        self.master.update()
        self.master.update_idletasks()
        
    ##########################################################################
    # Check whether HORs fully overlap and manage corresponding locations
    def checkFullOverlap(self):
        # Check consistency of the data
        assert len(self.hors) == len(self.monomers) == len(self.locations), "Inconsistent data!"
        # Loop over HORs to see whether some of them overlap
        i = 0
        while i < len(self.hors) - 1:
            # Get current monomer
            cmono = self.monomers[i]
            clocs = self.locations[i]
            # List of locations to remove
            clocs_to_remove = []
            j = i + 1
            while j < len(self.hors):
                # Get other monomer
                omono = self.monomers[j]
                olocs = self.locations[j]
                # List of locations to remove
                olocs_to_remove = []
                for oloc in olocs:
                    oloc_start = oloc[0]
                    oloc_end = oloc[1]
                    for cloc in clocs:
                        cloc_start = cloc[0]
                        cloc_end = cloc[1]
                        # Check whether locations overlap
                        if cloc_start <= oloc_start and cloc_end >= oloc_end:
                            # cloc contains oloc -> remove oloc
                            if oloc not in olocs_to_remove:
                                olocs_to_remove.append(oloc)
                        else:
                            # Check whether oloc contains cloc
                            if oloc_start <= cloc_start and oloc_end >= cloc_end:
                                # cloc contains oloc -> remove oloc
                                if cloc not in clocs_to_remove:
                                    clocs_to_remove.append(cloc)
                # Remove locations overlapping
                for oloc in olocs_to_remove:
                    self.locations[j].remove(oloc)
                j += 1
            # Remove locations overlapping
            for cloc in clocs_to_remove:
                self.locations[i].remove(cloc)
            i += 1
            
    ##########################################################################
    # Check whether HORs partially overlap and manage corresponding locations
    def checkPartialOverlap(self):
        # Check consistency of the data
        assert len(self.hors) == len(self.monomers) == len(self.locations), "Inconsistent data!"
        # Compute coverage for each HOR
        self.coverage = {}
        for i, hor in enumerate(self.hors):
            locs = self.locations[i]
            coverage = 0
            for loc in locs:
                coverage += (loc[1] - loc[0])
            self.coverage[hor] = coverage
        # Loop over HORs to see whether some of them overlap
        for i in range(len(self.hors)):
            chor = self.hors[i]
            # Get current locations
            clocs = self.locations[i]
            for j in range(len(self.hors)):
                if i != j:
                    ohor = self.hors[j]
                    # Get other locations
                    olocs = self.locations[j]
                    for oloc in olocs:
                        oloc_start = oloc[0]
                        oloc_end = oloc[1]
                        for cloc in clocs:
                            cloc_start = cloc[0]
                            cloc_end = cloc[1]
                            # Check whether locations partially overlap
                            if cloc_end > oloc_start and cloc_end < oloc_end:
                                if self.coverage[chor] >= self.coverage[ohor]:
                                    oloc[0] = cloc[1]
                                else:
                                    cloc[1] = oloc[0]
                j += 1
            i += 1
            
    ##########################################################################    
    # Show data
    def showData(self):
        if self.tree is None:
            self.popupMsg("You must load the phyloXML file before plotting!!!")
            return
        if self.shown:
            self.popupMsg("You have already shown data!!!")
            return
        if self.tree_canvas is not None:
            self.tree_canvas.get_tk_widget().destroy()
        if self.other_canvas is not None:
            self.other_canvas.get_tk_widget().destroy()
        
        # Color the tree with black
        for clade in self.tree.find_clades():
             clade.color = PX.BranchColor.from_name('black')
             
        # Check whether at least one of the HORs has been clicked
        if len(self.clicked) < 1:
            self.popupMsg("You can see monomers' tree only after choosing at least one of the HORs!!!")
            return
            
        # Extract HORs
        self.extractHORs()
        # Check whether HORs fully overlap
        self.checkFullOverlap()
        # Check whether HORs partially overlap
        self.checkPartialOverlap()
        # Set zoom flag to false
        self.zoomed = False
        
        # Generate BED file
        monomer_string = []
        abs_start = int(self.chr_seq[0])
        abs_end = int(self.chr_seq[1])
        # self.seq_name self.monomers (numero di HORS) self.locations (numero di HORS)
        bed_data = []
        for monos, locs in zip(self.monomers, self.locations):
            # Re-build monomers' string (i.e., HOR)
            mono_str = ""
            # Add monomers in the HOR
            for mono in monos:
                mono_str += str(mono)
                mono_str += ","
            # Remove last comma
            mono_str = mono_str[:-1]
            monomer_string.append(mono_str)
            # Append start, end, HOR and strand
            for loc in locs:
                bed_data.append([loc[0], loc[1], mono_str, loc[2]])
        # Retrieve all other data (those not related to HORs) from monomers' tree
        all_clades = self.tree.find_clades()
        for clade in all_clades:
            if clade.name:
                found = False
                for mono in self.monomers:
                    if clade.name in mono:
                        found = True
                if not found:
                    # Not in HORs, check if it contains sequences
                    if clade.sequences:                    
                        seqs = clade.sequences
                        # Loop over sequences
                        for seq in seqs:
                            # Extract relative start, end and strand
                            loc = seq.location
                            substr = loc.split('[')[1]
                            substr2 = substr.split(']')[0]
                            substr3 = substr.split(']')[1]
                            rel_start = substr2.split(':')[0]
                            rel_start = int(rel_start)
                            rel_end = substr2.split(':')[1]
                            rel_end = int(rel_end)
                            strand = substr3[1]
                            ndata = len(bed_data)
                            found = False
                            i = 0
                            while i < ndata and not found:
                                cdata = bed_data[i]
                                cstart = cdata[0]
                                cend = cdata[1]
                                if rel_start >= cstart and rel_end <= cend:
                                    found = True
                                i += 1
                            if not found: 
                                bed_data.append([rel_start, rel_end, "mono", strand])
        # Sort data based on locations
        bed_data.sort()
        # Build actual data to be stored (i.e., we collapse info when needed)
        bdata = []
        ndata = len(bed_data)
        # First entry is treated separatedly
        cdata = bed_data[0]
        if cdata[0] != 0:
            # First part is a monomer organization
            # Check if the first data to be stored is a monomer
            if cdata[2] != "mono":
                # HOR, we fill mono until the HOR
                bdata.append([0, cdata[0], "mono", "+"])
                bdata.append([cdata[0], cdata[1], cdata[2], cdata[3]])
            else:
                # Monomer organization, fill until the end
                bdata.append([0, cdata[1], cdata[2], cdata[3]])
        else:
            bdata.append([cdata[0], cdata[1], cdata[2], cdata[3]])
        # Store previous data
        prev_start = cdata[0]
        prev_end = cdata[1]
        prev_mono = cdata[2]
        prev_strand = cdata[3]
        # Loop over remaining data
        i = 1
        while i < ndata:
            # Get current entry
            cdata = bed_data[i]
            curr_start = cdata[0]
            curr_end = cdata[1]
            curr_mono = cdata[2]
            curr_strand = cdata[3]
            # Check if there is an overlap
            if curr_start >= prev_start and curr_end <= prev_end:
                # Overlap -> Ignore current row
                pass
            else:
                if curr_start == prev_end:
                    if curr_mono == prev_mono:
                        if curr_strand == prev_strand:
                            # We simply need to modify end of previous entry
                            pdata = bdata[-1]
                            pdata[1] = curr_end
                            prev_end = curr_end
                        else:
                            # Add new line
                            bdata.append([curr_start, curr_end, curr_mono, curr_strand])
                            prev_start = curr_start
                            prev_end = curr_end
                            prev_strand = curr_strand
                    else:
                        # Add new line
                        bdata.append([curr_start, curr_end, curr_mono, curr_strand])
                        prev_start = curr_start
                        prev_end = curr_end
                        prev_mono = curr_mono
                        prev_strand = curr_strand
                else:
                    # There is a gap, fill with a mono
                    if curr_mono != "mono":
                        if prev_mono != "mono":
                            bdata.append([prev_end, curr_start, "mono", "+"])
                        else:
                            pdata = bdata[-1]
                            pdata[1] = curr_start
                        bdata.append([curr_start, curr_end, curr_mono, curr_strand])
                        prev_start = curr_start
                        prev_end = curr_end
                        prev_mono = curr_mono
                        prev_strand = curr_strand
                    else:
                        pdata = bdata[-1]
                        pdata[1] = curr_end
                        prev_end = curr_end
            i += 1
               
        # Return the list of selected HORs
        horfile = "selected_hor_list.txt"
        fp = open(os.path.join(self.folder, horfile), "w")
        for hor in self.hors:
            fp.write("%s\n" % hor)
        fp.close()       
        
        mono_to_save = []
        # Extract monomers belonging to families (i.e., leaves of the tree)
        for monos in self.monomers:
            for mono in monos:
                if not mono in mono_to_save:
                    for clade in self.tree.find_clades():
                        if clade.name == mono:
                            tmp_tree = PX.Phylogeny(root=clade, name=clade.name)
                            leaves = tmp_tree.get_terminals()
                            monofile = self.seq_name + "_" + mono + ".txt"
                            fp = open(os.path.join(self.folder, monofile), "w")
                            for leaf in leaves:
                                fp.write("%s\n" % leaf.name)
                            fp.close()
                            mono_to_save.append(mono)
               
        # Build output BED filename (there is one file for each selection of the HORs)
        filename = os.path.splitext(self.filename)[0]
        outfile = filename + "_HORs.bed"
        """
        for hor in self.hors:
            outfile += hor + "-"
        outfile = outfile[:-1]
        outfile += ".bed"
        """
        fp = open(os.path.join(self.folder, outfile), "w")
        # Write data
        rows = len(bdata)
        # First row
        cdata = bdata[0]
        cols = len(cdata)
        assert cols == 4, "Inconsistent data length {%d}!".format(cols)
        # Header
        fp.write("track name=\"ItemRGBDemo\" description=\"Item RGB demonstration\" itemRgb=\"On\"\n")
        if cdata[0] != 0:
            fp.write("%s\t%d\t%d\tmono\t0\t+\t%d\t%d\t0,0,0\n" % (self.seq_name, abs_start, abs_start + cdata[0], abs_start, abs_start + cdata[0]))
        # Get index of HOR in list in order to retrieve the corresponding color
        idx = -1
        cnt = 0
        found = False
        for mono_str in monomer_string:
            if cdata[2] == mono_str:
                found = True
                idx = cnt
                break
            cnt += 1
        if found:
            # Extract the color
            ccolor = self.hor_colors[idx % len(self.hor_colors)]
            # Convert the color string into RGB (values bounded in the range [0,1])
            red, green, blue = mcolors.to_rgb(ccolor)
            # Adjust red, green and blue in the range [0,255]
            red = int(red * 255)
            green = int(green * 255)
            blue = int(blue * 255)
            fp.write("%s\t%d\t%d\t%s\t0\t%s\t%d\t%d\t%d,%d,%d\n" % (self.seq_name, abs_start + cdata[0], abs_start + cdata[1], cdata[2], cdata[3], abs_start + cdata[0], abs_start + cdata[1], red, green, blue))
        else:
            fp.write("%s\t%d\t%d\t%s\t0\t%s\t%d\t%d\t0,0,0\n" % (self.seq_name, abs_start + cdata[0], abs_start + cdata[1], cdata[2], cdata[3], abs_start + cdata[0], abs_start + cdata[1]))
        # Other rows
        row = 1
        while row < rows:
            # Get current data
            cdata = bdata[row]
            cols = len(cdata)
            assert cols == 4, "Inconsistent data length {%d}!".format(cols)
            # Get index of HOR in list in order to retrieve the corresponding color
            idx = -1
            cnt = 0
            found = False
            for mono_str in monomer_string:
                if cdata[2] == mono_str:
                    found = True
                    idx = cnt
                    break
                cnt += 1
            if found:
                # Extract the color
                ccolor = self.hor_colors[idx % len(self.hor_colors)]
                # Convert the color string into RGB (values bounded in the range [0,1])
                red, green, blue = mcolors.to_rgb(ccolor)
                # Adjust red, green and blue in the range [0,255]
                red = int(red * 255)
                green = int(green * 255)
                blue = int(blue * 255)
                fp.write("%s\t%d\t%d\t%s\t0\t%s\t%d\t%d\t%d,%d,%d\n" % (self.seq_name, abs_start + cdata[0], abs_start + cdata[1], cdata[2], cdata[3], abs_start + cdata[0], abs_start + cdata[1], red, green, blue))
            else:
                fp.write("%s\t%d\t%d\t%s\t0\t%s\t%d\t%d\t0,0,0\n" % (self.seq_name, abs_start + cdata[0], abs_start + cdata[1], cdata[2], cdata[3], abs_start + cdata[0], abs_start + cdata[1]))
            row += 1
        # Additional information to complete the sequence
        if abs_start + cdata[1] != abs_end:
            fp.write("%s\t%d\t%d\tmono\t0\t+\t%d\t%d\t0,0,0\n" % (self.seq_name, abs_start + cdata[1], abs_end, abs_start + cdata[1], abs_end))
        fp.close()

        # Plot data
        nmonos = len(self.monomers)
        assert nmonos >= 1, "Invalid number of HORs!!!"
        # Determine the width of the bar depending on the HOR with the highest number of monomers (only in case of multiple HORs to be plotted)
        mono_size = len(self.monomers[0])
        if nmonos > 1:
            i = 1
            while i < nmonos:
                if len(self.monomers[i]) > mono_size:
                    mono_size = len(self.monomers[i])
                i += 1
        
        # the figure that will contain the plot 
        self.fig = Figure(figsize = (10, 10), dpi = 100, constrained_layout=True)
        matplotlib.rc('font', size=6)
        self.fig.canvas.mpl_connect('button_press_event', self.expandSubTree)

        # adding the subplot 
        self.ax_tree = self.fig.add_subplot(1, 1, 1)
        # containing the Matplotlib figure 
        matplotlib.rcParams["lines.linewidth"] = 0.5
        # Create a copy of the tree
        treeToPlot = copy.deepcopy(self.tree)
        self.drawCollapsedTree(treeToPlot, axes=self.ax_tree)
        self.ax_tree.get_yaxis().set_visible(False)
        
        # the figure that will contain the plot 
        other_fig = Figure(figsize = (4, 4), dpi = 100, constrained_layout=True)
        H = 1.0 / (nmonos + 1)
        gs = other_fig.add_gridspec(nmonos + 1, 1, height_ratios=[H for _ in range(nmonos + 1)], hspace=0.1)
            
        # Monomers in HORs
        for j in range(nmonos):
            hor_size = len(self.monomers[j])
            # Workaround to make the HORs (with monomers) be displayed in an acceptable fashion
            gs_hor = gridspec.GridSpecFromSubplotSpec(2, 2, subplot_spec=gs[j, 0], width_ratios=[hor_size, self.hor_len - hor_size], height_ratios=[1, 1], hspace=0.1)
            cmono = self.monomers[j]
            cmono_colors = self.monomer_colors[j]
            ax_hor = other_fig.add_subplot(gs_hor[0,0])
            ax_hor.set_xlim([0, len(cmono)])
            ax_hor.set_ylim([0, 1])
            ax_hor.set_xticks(np.arange(0, len(cmono) + 1, 1))
            ax_hor.get_yaxis().set_visible(False)
            ax_hor.get_xaxis().set_visible(False)
            for i in range(len(cmono)):
                rect = patches.Rectangle((i, 0), 1, 1, facecolor=cmono_colors[i].to_hex(), edgecolor='black')
                ax_hor.add_patch(rect)
                ax_hor.text(i, 1.25, str(cmono[i]))
            hor_rect = patches.Rectangle((-4.5, 0.25), 3, 0.5, color=self.clicked_colors[j], clip_on=False)
            ax_hor.add_patch(hor_rect)
            N = len(self.hors[j])
            ax_hor.text(-(2 + 0.1 * N / 2), 1, self.hors[j], fontsize='xx-small')
           
        # Locations of HORs in sequence
        nlocs = len(self.locations)
        assert nmonos == nlocs, "Inconsistent data sizes: {nmonos} vs {nlocs}".format(nmonos, nlocs)
        
        # Workaround to make the chromosome sequence be displayed in an acceptable fashion
        gs_seq = gridspec.GridSpecFromSubplotSpec(2, 1, subplot_spec=gs[nmonos, 0], width_ratios=[self.seq_len], height_ratios=[1, 1], hspace=0.1)
        ax_seq = other_fig.add_subplot(gs_seq[0,:])
        ax_seq.set_xlim([0, self.seq_len])
        ax_seq.set_ylim([0, 1])
        ax_seq.set_yticks([])
        ax_seq.set_xticks([0, self.seq_len])
        ax_seq.set_xticklabels(self.chr_seq)
        ax_seq.add_patch(patches.Rectangle((0, 0), self.seq_len, 1, facecolor='grey'))
        ax_seq.text(-500000.5, 0.5, self.seq_name) # Maybe the x-coordinate of this text can be changed...
        plts = []
        for j in range(nlocs):
            # Compute the ratio between actual length and plot
            # Color locations of the HORS with black
            trans = ax_seq.get_xaxis_transform()
            cplt = None
            for loc in self.locations[j]:
                cplt = ax_seq.add_patch(patches.Rectangle((loc[0], 0), (loc[1] - loc[0]), 1, facecolor=self.clicked_colors[j]))
                # Check whether the strand is inverted (i.e., '-'). If negative, add an arrow over the bar corresponding to the location
                if "-" in loc[2]:
                    ax_seq.annotate("", (loc[1], 1.25), (loc[0], 1.25), xycoords=trans, arrowprops=dict(arrowstyle='<|-'))#width=1))
            if cplt is not None:
                plts.append(cplt)

        # Create the canvas
        self.tree_canvas = FigureCanvasTkAgg(self.fig, master=self.w)
        self.tree_canvas.draw()
        self.other_canvas = FigureCanvasTkAgg(other_fig, master=self.z)
        self.other_canvas.draw()
        # placing the canvas on the Tkinter window 
        self.tree_canvas.get_tk_widget().pack() 
        self.other_canvas.get_tk_widget().pack() 
        self.master.update()
        self.master.update_idletasks()
        # Set flag to true
        self.shown = True
        
    ##########################################################################    
    # Zoom out
    def zoomOut(self):
        # It works only if a clickable patch has been pressed
        if self.zoomed:
            # Clear the window
            self.ax_tree.clear()
            if self.axes_ins is not None:
                self.axes_ins.remove()
            if self.ax_tree_ins is not None:
                for elem in self.ax_tree_ins:
                    elem.remove()
                self.ax_tree_ins = []
            treeToPlot = copy.deepcopy(self.tree)
            # Plot the collapsed tree
            self.drawCollapsedTree(treeToPlot, axes=self.ax_tree)
            self.ax_tree.get_yaxis().set_visible(False)
            self.tree_canvas.draw()
            # Unset list
            self.zoomed_nodes = []
            self.zoomed_coords = []
            # Set zoom flag to false
            self.zoomed = False
        else:
            self.popupMsg("No active zoom found!!!")
        return
        
    ##########################################################################    
    # Reset
    def reset(self):
        # Set color of nodes in the HOR tree to black
        for patch in self.patches:
            patch.set_color('black')
        # Reset number of clicked items
        self.clicked = []
        self.num_clicked = 0
        # Reset flag
        self.shown = False
        self.canvas.draw()
        # Delete visualization of other figures (i.e., HORs, sequence and monomers' tree)
        if self.tree_canvas is not None:
            self.tree_canvas.get_tk_widget().destroy()
        if self.other_canvas is not None:
            self.other_canvas.get_tk_widget().destroy()
        
    ##########################################################################    
    # Select the file to be loaded
    def chooseFile(self, event):
        self.filename = self.combo['file'].get()
        # We add the suffix previously removed for visualization purposes
        self.filename += ".tree.xml"

    ##########################################################################
    # GUI intialize function: setup the tk environment
    def initialize(self):
        #print(self.master.winfo_screenwidth(), self.master.winfo_screenheight())
        self.toolbar = tk.Frame(self.master, relief='raised', bd=2)
        self.toolbar.pack(side='top', fill='x')
        self.font=("Arial", 25)
        self.fontA=("Arial", 14)
        self.fontB=("Arial", 12)
        # Buttons (and associated commands)
        self.load_file = tk.Button(self.toolbar, text="LoadFile", command=lambda: self.loadFile(filename=self.filename))
        self.plot_tree = tk.Button(self.toolbar, text="PlotTree", command=lambda: self.plotTree())
        self.show_data = tk.Button(self.toolbar, text="ShowData", command=lambda: self.showData())
        self.zoom_in = tk.Button(self.toolbar, text="ZoomIn", command=lambda: self.zoomIn())
        self.zoom_out = tk.Button(self.toolbar, text="ZoomOut", command=lambda: self.zoomOut())
        self.reset_win = tk.Button(self.toolbar, text="Reset", command=lambda: self.reset())
        self.load_file.pack(side='left')
        self.plot_tree.pack(side='left')
        self.show_data.pack(side='left')
        self.zoom_in.pack(side='left')
        self.zoom_out.pack(side='left')
        self.reset_win.pack(side='left')
        # Create frame where the monomers' tree will be displayed
        self.w = tk.Frame(self.master, background="dimgray")
        self.w.pack(side='left',fill='both',expand='True')
        # Create frame containing the visualization of the HORs' tree
        self.v = tk.Frame(self.master, background="white")
        self.v.pack(side='top',fill='both',expand='True')
        # Create frame hosting the plot of selected HORs and their locations in the chromosome sequence
        self.z = tk.Frame(self.master, background="white")
        self.z.pack(side='bottom',fill='both',expand='True')
        # Create combobox for files
        self.combo_var = {}
        self.combo = {}
        # Combobox creation 
        xml_file_values = [os.path.splitext(os.path.splitext(filename)[0])[0] for filename in os.listdir(self.filedir) if filename.endswith(".xml")]
        # Make visible only files without the suffix "monomers" or "hors"
        file_values = []
        for elem in xml_file_values:
            if not "monomers" in elem and not "hors" in elem:
                file_values.append(elem)
        file_values.sort(key=natural_keys)
        # Extract a set from the list (unique values)
        file_values = set(file_values)
        # Sort files alphabetically
        file_values = sorted(file_values)
        self.combo_var['file'] = tk.StringVar() 
        self.combo['file'] = ttk.Combobox(self.toolbar, values = tuple(file_values), width = 15, textvariable = self.combo_var['file'])
        # Adding combobox drop down list 
        self.combo['file'].pack(side='right') 
        self.combo['file'].current()
        self.combo['file'].bind("<<ComboboxSelected>>", self.chooseFile)
        # label 
        file_label = tk.Label(self.toolbar, text = "Load file(s) :",  font = ("Times New Roman", 10))
        file_label.pack(side='right')
        
