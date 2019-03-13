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

    # Attempt to create directory
    dir_exists = False
    try:
        os.makedirs("gdm_output", exist_ok=True)  # Succeeds even if directory exists.
    except OSError:
        print("Creation of the directory /gdm_output failed")
    else:
        dir_exists = True
        print("Successfully created 'gdm_output' directory.")

    # Deal with relevant options
    if log:

        currentDT = datetime.datetime.now()

        filename = (currentDT.strftime("%Y-%m-%d|%Hhr-%Mm-%Ss")) + "-gdm.log"
        numeric_level = getattr(logging, log.upper(), None)

        if dir_exists:
            print("log is set and dir_exists")
            logging.basicConfig(filename='gdm_output/'+filename,
                                level=numeric_level)
        else:
            logging.basicConfig(filename=filename,
                                level=numeric_level)

    if savetable:
        if dir_exists:
            dc.saveTable(fileName="alteredTable.fits", file_format="fits")
        else:
            dc.saveTable(fileName="gdm_output/alteredTable.fits", file_format="fits")

    if showtree:
        treefile = "data/savedTree.txt"
        if (sys.platform == "darwin"):  # MacOS
            print("test")
            os.system("open " + test_dir + "/data/savedTree.txt")
        elif (sys.platform == "cygwin" or sys.platform == "win32"):  # Windows
            os.system("start "+treefile)
        elif (sys.platform == "linux"):
            os.system("xdg-open "+treefile)

    # Compares object and handles/retrieves data respectively.
    ct = ComparisonTree(run_mode=True)
    dc = DataController()

    # Pass ct and dc to context object for click commands
    ctx.obj['ct'] = ct
    ctx.obj['dc'] = dc
    ctx.obj['sp'] = showplot


@main.command()
@click.argument('name', type=str)
@click.option('-match-tol', type=float,
              help='Optional 2D match tolerance (in arc seconds). Default value is 1.0 arcsec')
@click.pass_context
def byname(ctx, name, match_tol):
    '''
    Downloads objects, via NED and SIMBAD, from region around object name

    Arguments:
        NAME - the name of the object to be searched around
    '''
    dc = ctx.obj['dc']
    ct = ctx.obj['ct']

    logging.info("Querying region by name: {}".format(name))

    if match_tol:
        logging.info("Confirm: match-toll passed")
        dc.query_region_by_name(name, match_tol=match_tol)
    else:
        logging.info("No match-tol passed")
        dc.query_region_by_name(name)

    # Pass table to comparison tree to compare each object
    ct.compare_objects(dc.combined_table)

    if dc.combined_table is not None:
        dc.combined_table.show_in_browser(jsviewer=True)
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
