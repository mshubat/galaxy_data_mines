# Galaxy Agent - User Interface
from .data_controller import DataController
from .comparison_tree import ComparisonTree
import sys
import click
import os
import logging
import datetime
import pandas as pd

this_dir, this_filename = os.path.split(__file__)
test_dir = os.path.dirname(__file__)


@click.group()
@click.option('--log',
              type=click.Choice(['debug', 'info', 'warning', 'error', 'critical']),
              help="Use this option to display messages at log level of choice.")
@click.option('--glossary',
              is_flag=True,
              help="Shows the glossary of terms used by gdmines tool.")
@click.option('--showtree',
              is_flag=True,
              help="Shows the simbad tree for relationship comparisons.")
@click.option('--showplot',
              is_flag=True,
              help="Shows graphical plot of comparison results.")
@click.option('--saveplot',
              is_flag=True,
              help="Shows graphical plot of comparison results.")
@click.option('--savetable',
              type=click.Choice(['csv', 'fits']),
              help="Use this option to save the comparison table to the output\
              folder in the format of your choice.")
@click.option('--savestats',
              is_flag=True,
              help="Saves statistics of comparison results.")
@click.pass_context
def main(ctx, log, glossary,
         showtree, showplot,
         saveplot, savetable, savestats):
    '''
    Compares object classifications between NED and SIMBAD.
    '''

    # Environment setup and checks.
    ctx.ensure_object(dict)  # Click convention.

    # Attempt to create directory if a relevant option is provided.
    dir_exists = False
    if log or savetable or saveplot or savestats:
        try:
            os.makedirs("gdm_output", exist_ok=True)  # Succeeds even if directory exists.
        except OSError:
            logging.error("Creation of the directory /gdm_output failed")
        else:
            dir_exists = True
            logging.info("Successfully created 'gdm_output' directory.")

    # Deal with option(s) which don't need any additional information.
    if log:

        currentDT = datetime.datetime.now()

        filename = (currentDT.strftime("%Y-%m-%d|%Hhr-%Mm-%Ss")) + "-gdm.log"
        numeric_level = getattr(logging, log.upper(), None)

        if dir_exists:
            logging.basicConfig(filename='gdm_output/'+filename,
                                level=numeric_level)
        else:
            logging.basicConfig(filename=filename,
                                level=numeric_level)

        root = logging.getLogger()
        root.setLevel(numeric_level)

        handler = logging.StreamHandler(sys.stdout)
        handler.setLevel(logging.DEBUG)
        formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
        handler.setFormatter(formatter)
        root.addHandler(handler)

    if showtree:
        treefile = "/data/savedTree.txt"
        safelyopenfile(treefile)

    # Create important objects
    ct = ComparisonTree(run_mode=True)
    dc = DataController()

    # Store reused objects in context.
    ctx.obj['ct'] = ct
    ctx.obj['dc'] = dc

    # Store option values in context.
    ctx.obj['glossary'] = glossary
    ctx.obj['showplot'] = showplot
    ctx.obj['savetable'] = savetable
    ctx.obj['saveplot'] = saveplot
    ctx.obj['savestats'] = savestats

    # Store additional variables in context.
    ctx.obj['dir_exists'] = dir_exists


@main.command()
@click.argument('name', type=str)
@click.option('-match-tol', type=float,
              help='Optional 2D match tolerance (in arc seconds). Default value is 1.0 arcsec.')
@click.option('-obj-radius', type=float,
              help='Optional radius to search around object (in arc minutes). Default value is 1.0 arcmin.')
@click.pass_context
def byname(ctx, name, match_tol, obj_radius):
    '''
    Downloads objects, via NED and SIMBAD, from region around object name

    Arguments:
        NAME - the name of the object to be searched around
    '''
    dc = ctx.obj['dc']
    ct = ctx.obj['ct']

    logging.info("Query requested: {}".format(name))

    # Option values are assumed to be valid by default.
    match_tol_valid = True
    obj_radius_valid = True

    # If supplied, confirm match_tol and obj_radius are valid values.
    if match_tol:
        match_tol_valid = withinbounds(match_tol, lower=0, upper=60)

    if obj_radius:
        obj_radius_valid = withinbounds(obj_radius, lower=0, upper=60)

    # Proceed based on validity of any options passed.
    if match_tol_valid == False:
        logging.error("Invalid match_tol entered.")
        logging.error("mmatch_tol must be between 0 and 60 arcsecs")
    elif obj_radius_valid == False:
        logging.error("Invalid obj_radius entered.")
        logging.error("obj_radius must be between 0 and 60 arcmins")
    else:
        if match_tol and obj_radius:
            logging.info("Confirm: match-tol & obj-radius passed.")
            dc.query_region_by_name(name, match_tol=match_tol, obj_radius=obj_radius)
        elif match_tol:
            logging.info("Confirm: match-tol passed.")
            dc.query_region_by_name(name, match_tol=match_tol)
        elif obj_radius:
            logging.info("Confirm: obj-radius passed.")
            dc.query_region_by_name(name, obj_radius=obj_radius)
        else:
            logging.info("Default settings used for byname query.")
            dc.query_region_by_name(name)

        if dc.combined_table is not None:
            # Pass table to comparison tree to compare each object
            ct.compare_objects(dc.combined_table)
            common_option_handler(ctx)


@main.command()
@click.pass_context
def bycoord(ctx):
    '''
    Downloads objects, via NED and SIMBAD, from region described by coordinates
    '''
    pass


@main.command()
@click.argument('ned-name', type=str)
@click.argument('simbad-name', type=str)
@click.pass_context
def local(ctx, ned_name, simbad_name):
    '''
    Uses local data provided to compare objects in NED and SIMBAD
    '''
    combined_table = os.path.join(this_dir, "data", "M83_NScomb_mar1.fits")

    dc = ctx.obj['dc']
    ct = ctx.obj['ct']
    # ctx test

    dc.load_data(combined_table)

    # Pass table to comparison tree to compare each object
    ct.compare_objects(dc.combined_table)

    if dc.combined_table is not None:
        dc.combined_table.show_in_browser(jsviewer=True)


# ----------- Helper Functions ----------- #

def withinbounds(val, lower, upper):
    '''
    Function created to ensure float is within bounds.
    '''
    if val > lower and val <= upper:
        return True
    else:
        return False


def safelyopenfile(filepath):
    '''open file on any platform'''

    if (sys.platform == "darwin"):  # MacOS
        os.system("open " + test_dir + filepath)
    elif (sys.platform == "cygwin" or sys.platform == "win32"):  # Windows
        os.system("start " + test_dir + filepath)
    elif (sys.platform == "linux"):
        os.system("xdg-open " + test_dir + filepath)


def common_option_handler(ctx):
    '''
    Handle options/operations relavant to all query types.
    '''
    dc = ctx.obj['dc']

    dc.combined_table.show_in_browser(jsviewer=True)
    # dc.stats.derive_table_stats(dc.combined_table)

    # If show plot option present.
    if ctx.obj['showplot']:
        DataController.plot_match_table(dc.combined_table)

    # If glossary option present.
    if ctx.obj['glossary']:
        safelyopenfile('/data/glossary.txt')


if __name__ == "__main__":
    main(obj={})
