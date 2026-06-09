# -*- coding: utf-8 -*-
# @author Paolo Pagliuca <paolo.pagliuca@istc.cnr.it>

import os
import traceback
import tkinter as tk
from tkinter import messagebox

import click

from phylotreegui.gui import GUIFactory


def _report_callback_exception(exc_type, exc_value, exc_tb):
    """Surface uncaught Tk callback errors instead of letting Tk swallow them
    into a blank dialog. The full traceback goes to the terminal; a concise,
    never-empty message is shown in a dialog."""
    traceback.print_exception(exc_type, exc_value, exc_tb)
    message = "".join(traceback.format_exception_only(exc_type, exc_value)).strip()
    if not message:
        message = f"{exc_type.__name__} (no message)"
    messagebox.showerror("PhyloTreeGUI error", message)


def start(path=None, unit_length=0, folder=None, fullscreen=False):
    if path is None:
        path = os.getcwd()
    root = tk.Tk()
    root.report_callback_exception = _report_callback_exception
    root.title("PhyloTreeGUI")
    root.lift()
    if fullscreen:
        width, height = root.winfo_screenwidth(), root.winfo_screenheight()
        root.geometry("%dx%d+0+0" % (width, height))

    # Window can be closed by pressing "ESC" keyboard button
    def close_escape(event=None):
        root.destroy()

    root.bind("<Escape>", close_escape)
    GUIFactory.create_gui(root, path=path, unit_length=unit_length, folder=folder)
    root.mainloop()


@click.command()
@click.version_option(message="PhyloTreeGUI %(version)s")
@click.option(
    "--path",
    default=None,
    type=click.Path(exists=True, file_okay=False, dir_okay=True),
    help="Directory containing the phyloXML files (defaults to the current directory).",
)
@click.option(
    "--unit-length",
    "unit_length",
    default=0,
    type=click.IntRange(0, 1),
    help="HOR branch lengths: 0 = floating-point lengths, 1 = unit lengths.",
)
@click.option(
    "--folder",
    default=None,
    type=click.Path(file_okay=False, dir_okay=True),
    help="Output directory for BED/TXT/PNG/CSV files (created if missing).",
)
@click.option(
    "--fullscreen/--no-fullscreen",
    default=False,
    help="Launch the GUI at full screen size.",
)
def main(path, unit_length, folder, fullscreen):
    """Launch the PhyloTreeGUI desktop application."""
    start(path=path, unit_length=unit_length, folder=folder, fullscreen=fullscreen)


if __name__ == "__main__":
    main()
