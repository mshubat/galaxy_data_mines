# Galaxy Agent - User Interface
from .data_controller import DataController
from .comparison_tree import ComparisonTree
import sys
import click
import os
import logging
import datetime
this_dir, this_filename = os.path.split(__file__)
test_dir = os.path.dirname(__file__)

import pandas as pd


@click.group()
@click.option('--log',
              type=click.Choice(['debug', 'info', 'warning', 'error', 'critical']),
              help="Use this option to display messages at log level of choice.")
@click.option('--savetable',
              is_flag=True,
              help="Use this option to save comparison table to output folder.")
@click.option('--showplot',
              is_flag=True,
              help="Shows graphical plot of comparison results.")
@click.option('--showtree',
              is_flag=True,
              help="Shows the simbad tree for relationship comparisons.")
@click.pass_context
def main(ctx, log, savetable, showplot, showtree):
    '''
    Compares object classifications between NED and SIMBAD.
    '''

    # Environment setup and checks.
    ctx.ensure_object(dict)  # Click convention.

    # Attempt to create directory if relevant option provided
    dir_exists = False
    if log or savetable:
        try:
            os.makedirs("gdm_output", exist_ok=True)  # Succeeds even if directory exists.
        except OSError:
            logging.error("Creation of the directory /gdm_output failed")
        else:
            dir_exists = True
            logging.info("Successfully created 'gdm_output' directory.")

    # Deal with relevant options
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

    if savetable:
        if dir_exists:
            dc.saveTable(fileName="alteredTable.fits", file_format="fits")
        else:
            dc.saveTable(fileName="gdm_output/alteredTable.fits", file_format="fits")

    if showtree:
        treefile = "data/savedTree.txt"
        if (sys.platform == "darwin"):  # MacOS
            os.system("open " + test_dir + "/data/savedTree.txt")
        elif (sys.platform == "cygwin" or sys.platform == "win32"):  # Windows
            os.system("start "+treefile)
        elif (sys.platform == "linux"):
            os.system("xdg-open "+treefile)

    # Compares object and handles/retrieves data respectively.
    ct = ComparisonTree(run_mode=True)
    dc = DataController()

    # Pass ct and dc to context object for click commands.
    ctx.obj['ct'] = ct
    ctx.obj['dc'] = dc
    ctx.obj['sp'] = showplot
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
        NAME - the name of the object to    be searched around
    '''
    dc = ctx.obj['dc']
    ct = ctx.obj['ct']

    logging.info("Querying region by name: {}".format(name))

    # Confirm match_tol and obj_radius are valid values.
    if match_tol:
        if match_tol > 0 and match_tol <= 60:
            match_tol_status = "valid"
        else:
            match_tol_status = "invalid"
    else:
        match_tol_status = "valid"

    if obj_radius:
        if obj_radius > 0 and obj_radius <= 60:
            obj_radius_status = "valid"
        else:
            obj_radius_status = "invalid"
    else:
        obj_radius_status = "valid"

    if match_tol_status == "invalid":
        logging.error("Invalid match_tol entered.")
        logging.error("mmatch_tol must be between 0 and 60 arcsecs")
    elif obj_radius_status == "invalid":
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

        # Pass table to comparison tree to compare each object
        ct.compare_objects(dc.combined_table)

        if dc.combined_table is not None:
            dc.combined_table.show_in_browser(jsviewer=True)
            # dc.get_table_stats().to_csv(r"gdm_output/result_stats.csv")

            # If show plot option is present.
            if ctx.obj['sp'] == True:
                DataController.plot_match_table(dc.combined_table)


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

    dc.load_data(combined_table)

    # Pass table to comparison tree to compare each object
    ct.compare_objects(dc.combined_table)

    if dc.combined_table is not None:
        dc.combined_table.show_in_browser(jsviewer=True)


if __name__ == "__main__":
    main(obj={})
