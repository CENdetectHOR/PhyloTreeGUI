# -*- coding: utf-8 -*-
# @author Paolo Pagliuca <paolo.pagliuca@istc.cnr.it>

import os, sys, getopt, importlib
from gui import GUIFactory
########################################################################################
## main functions
########################################################################################

def start(argv):
    
    global tk
    import tkinter as tk
    root = tk.Tk()
    root.title('PhyloTreeGUIv1.0')
    root.lift()
    GUIFactory.create_gui(root)
    root.mainloop()

if __name__ == "__main__":
    start(sys.argv[1:])

#sys.argv is the list of commandline arguments passed to the Python program
