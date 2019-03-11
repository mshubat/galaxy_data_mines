# Galaxy Agent - User Interface
from .data_controller import DataController
from .comparison_tree import ComparisonTree
import sys
import click
import os
import logging
this_dir, this_filename = os.path.split(__file__)
test_dir = os.path.dirname(__file__)


@click.group()
@click.option('-log',
              type=click.Choice(['debug', 'info', 'warning', 'error', 'critical']),
              help="Use this option to display messages at log level of choice.")
@click.option('-savetable',
              is_flag=True,
              help="Use this option to save comparison table to output folder.")
@click.pass_context
def main(ctx, log, savetable):
    '''
    Compares object classifications between NED and SIMBAD.
    '''
    # Ensure that ctx.obj exists and is a dict (in case `cli()` is called
    # by means other than the `if __name__ == "__main__"` block below
    ctx.ensure_object(dict)

    # Create comparison tree for object comparisons
    ct = ComparisonTree(run_mode=True)

    # Create data controller to handle and retrieve data
    dc = DataController()

    dir_exists = False
    # If log or savetable options set, attempt to create directory
    if log or savetable:
        try:
            os.makedirs("output", exist_ok=True)  # succeeds even if directory exists.
        except OSError:
            logging.warning("Creation of the directory /output failed")
            dc.saveTable(fileName="alteredTable.fits", file_format="fits")
        else:
            dir_exists = True
            log.info("Successfully created the directory /output")
            dc.saveTable(fileName="output/alteredTable.fits", file_format="fits")

    if log:
        if dir_exists:
            logging.basicConfig(filename='data/output.log', level=logging.DEBUG)
        else:
            logging.basicConfig(filename='output.log', level=logging.DEBUG)

    if savetable:
        if dir_exists:
            dc.saveTable(fileName="alteredTable.fits", file_format="fits")
        else:
            dc.saveTable(fileName="data/alteredTable.fits", file_format="fits")

    # Pass ct and dc to context object for click commands
    ctx.obj['ct'] = ct
    ctx.obj['dc'] = dc


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
