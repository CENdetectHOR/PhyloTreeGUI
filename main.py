# -*- coding: utf-8 -*-
# @author Paolo Pagliuca <paolo.pagliuca@istc.cnr.it>

import os, sys, getopt, importlib
import argparse
from gui import GUIFactory
########################################################################################
## main functions
########################################################################################

def start(path=os.getcwd(), unit_length=0, folder=None, fullscreen=None):
    
    global tk
    import tkinter as tk
    root = tk.Tk()
    root.title('PhyloTreeGUI')
    root.lift()
    if fullscreen:
        width, height = root.winfo_screenwidth(), root.winfo_screenheight()
        root.geometry('%dx%d+0+0' % (width, height))
    
    # Window can be closed by pressing "ESC" keyboard button
    def close_escape(event=None):
        root.destroy()
    
    root.bind("<Escape>", close_escape)
    GUIFactory.create_gui(root, path=path, unit_length=unit_length, folder=folder)
    root.mainloop()

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='PhyloTreeGUIv1.0')
    parser.add_argument('--path', type=str, default=os.getcwd())
    parser.add_argument('--unit_length', type=int, default=0)
    parser.add_argument('--folder', type=str, default=None)
    parser.add_argument('--fullscreen', action=argparse.BooleanOptionalAction)
    arg = parser.parse_args()
    assert arg.unit_length == 0 or arg.unit_length == 1, "Invalid argument for parameter <unit_length>"
    start(path=arg.path, unit_length=arg.unit_length, folder=arg.folder, fullscreen=arg.fullscreen)

#sys.argv is the list of commandline arguments passed to the Python program
