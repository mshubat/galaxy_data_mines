# Galaxy Agent - User Interface
from .data_controller import DataController
from .comparison_tree import ComparisonTree
import sys
import click
import os
this_dir, this_filename = os.path.split(__file__)
test_dir = os.path.dirname(__file__)


@click.group()
@click.option('-log', is_flag=True, help="Save console output to log file.")
@click.option('-savetable', is_flag=True, help="Use this option to save comparison table to output folder")
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

    if log:
        original = sys.stdout
        try:
            os.makedirs("output", exist_ok=True)  # succeeds even if directory exists.
        except OSError:
            print("Creation of the directory /data failed")
            sys.stdout = open('log.txt', 'w+')
        else:
            print("Successfully created the directory /data")
            sys.stdout = open('data/log.txt', 'w+')

    if savetable:
        dc.saveTable(fileName="alteredTable.fits", file_format="fits")

    ctx.obj['ct'] = ct
    ctx.obj['dc'] = dc


@main.command()
@click.argument('name', type=str)
@click.pass_context
def byname(ctx, name):
    '''
    Downloads objects, via NED and SIMBAD, from region around object name

    Arguments:
        NAME - the name of the object to be searched around
    '''
    dc = ctx.obj['dc']

    print("Querying region by name: {}".format(name))
    dc.query_region_by_name(name)

    # Pass table to comparison tree to compare each object
    # ct.compare_objects(dc.combined_table)

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
