import os
import argparse

filedir = os.path.dirname(__file__)

def get_results_dir():
    resdir = os.path.join(filedir, "results")
    if os.path.isfile(resdir):
        raise OSError(f"Default results directory {resdir} is already a file")
    elif not os.path.isdir(resdir):
        os.mkdir(resdir)
    return resdir

def get_image_dir():
    plotdir = os.path.join(filedir, "images")
    if os.path.isfile(plotdir):
        raise OSError(f"Default image directory {plotdir} is already a file")
    elif not os.path.isdir(plotdir):
        os.mkdir(plotdir)
    return plotdir

def get_argparser():
    argparser = argparse.ArgumentParser()
    argparser.add_argument("--results-dir", default=get_results_dir())
    return argparser

def get_plot_argparser():
    argparser = get_argparser()
    argparser.add_argument("--image-dir", default=get_image_dir())
    argparser.add_argument("--show", action="store_true")
    argparser.add_argument("--no-save", action="store_true")
    argparser.add_argument("--opaque", action="store_true")
    argparser.add_argument("--title", default=None)
    return argparser
