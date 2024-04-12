# -*- coding: utf-8 -*-
# @author Paolo Pagliuca <paolo.pagliuca@istc.cnr.it>

import os, sys, getopt, importlib
import argparse
from gui import GUIFactory
########################################################################################
## main functions
########################################################################################

def start(argv, fullscreen=None):
    
    global tk
    import tkinter as tk
    root = tk.Tk()
    root.title('PhyloTreeGUIv1.0')
    root.lift()
    if fullscreen:
        width, height = root.winfo_screenwidth(), root.winfo_screenheight()
        root.geometry('%dx%d+0+0' % (width, height))
    
    # Window can be closed by pressing "ESC" keyboard button
    def close_escape(event=None):
        root.destroy()
    
    root.bind("<Escape>", close_escape)
    GUIFactory.create_gui(root)
    root.mainloop()

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='PhyloTreeGUIv1.0')
    parser.add_argument('--fullscreen', action=argparse.BooleanOptionalAction)
    arg = parser.parse_args()
    start(sys.argv[1:], fullscreen=arg.fullscreen)

#sys.argv is the list of commandline arguments passed to the Python program
