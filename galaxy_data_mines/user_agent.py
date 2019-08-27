# Galaxy Agent - User Interface
from .data_controller import DataController
from .comparison_tree import ComparisonTree
import sys
import click
import os
import logging
import datetime
import pandas as pd

package_dir = os.path.dirname(__file__)
working_dir = os.getcwd()


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
@click.option('--showtable',
              is_flag=True,
              help="Shows the table of object comparisons.")
@click.option('--savelog',
              type=click.Choice(['debug', 'info', 'warning', 'error', 'critical']),
              help="Use this option to save logs at the level of choice.")
@click.option('--saveplot',
              is_flag=True,
              help="Saves graphical plot of comparison results.")
@click.option('--savetable',
              type=click.Choice(['csv', 'fits']),
              help="Use this option to save the comparison table to the output\
              folder in the format of your choice.")
@click.option('--savestats',
              is_flag=True,
              help="Saves statistics of comparison results.")
@click.option('--shortstats',
              is_flag=True,
              help="Show statistics of comparison results in short format.")
@click.pass_context
def main(ctx, log, glossary,
         showtree, showplot, showtable,
         savelog, saveplot, savetable, savestats, shortstats):
    '''
    Compares object classifications between NED and SIMBAD.
    '''

    ctx.ensure_object(dict)  # Click convention.

    # "Save" options create a directory.
    dir_exists = False
    if savelog or savetable or saveplot or savestats:
        currentDT = datetime.datetime.now()
        filename = (currentDT.strftime("%Y-%m-%d|%Hhr-%Mm-%Ss")) + "-gdm"
        ctx.obj['filename'] = filename

        try:
            # Succeeds even if directory exists.
            os.makedirs(working_dir+"/gdm_output", exist_ok=True)
        except OSError:
            logging.error("Creation of the directory 'gdm_output' failed")
        else:
            dir_exists = True
            logging.info("Successfully created 'gdm_output' directory.")

    if log or savelog:
        logFormatter = logging.Formatter(
            "%(asctime)s [%(threadName)-12.12s] [%(levelname)-5.5s]  %(message)s")
        rootLogger = logging.getLogger()
        rootLogger.handlers = []
        rootLogger.setLevel(logging.DEBUG)

    # OPTION: log
    if log:
        numeric_level = getattr(logging, log.upper(), None)

        consoleHandler = logging.StreamHandler()
        consoleHandler.setFormatter(logFormatter)
        consoleHandler.setLevel(numeric_level)
        rootLogger.addHandler(consoleHandler)

    # OPTION: Save Log
    if savelog:
        numeric_level = getattr(logging, savelog.upper(), None)

        if dir_exists:
            logPath = working_dir+"/gdm_output"
        else:
            logPath = working_dir

        fileHandler = logging.FileHandler("{0}/{1}.log".format(logPath, filename))
        fileHandler.setFormatter(logFormatter)
        fileHandler.setLevel(numeric_level)
        rootLogger.addHandler(fileHandler)

    if showtree:
        treefile = "/data/savedTree.txt"
        safelyopenfile(treefile)

    # OPTION: shortstats
    if shortstats:
      stats_template = "statsSchema_short.txt"
    else:
      stats_template = "statsSchema_long.txt"

    # Create important objects
    ct = ComparisonTree(run_mode=True)
    dc = DataController(st=stats_template)

    # Store reused objects in context.
    ctx.obj['ct'] = ct
    ctx.obj['dc'] = dc
    ctx.obj['stats_template'] = stats_template

    # Store option values in context.
    ctx.obj['glossary'] = glossary
    ctx.obj['showplot'] = showplot
    ctx.obj['showtable'] = showtable
    ctx.obj['savetable'] = savetable
    ctx.obj['saveplot'] = saveplot
    ctx.obj['savestats'] = savestats
    ctx.obj['shortstats'] = shortstats

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
    ctx.obj['name'] = name

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
            dc.query_region(name, match_tol=match_tol, obj_radius=obj_radius)
        elif match_tol:
            logging.info("Confirm: match-tol passed.")
            dc.query_region(name, match_tol=match_tol)
        elif obj_radius:
            logging.info("Confirm: obj-radius passed.")
            dc.query_region(name, obj_radius=obj_radius)
        else:
            logging.info("Default settings used for byname query.")
            dc.query_region(name)

        if dc.combined_table is not None:
            # Pass table to comparison tree to compare each object
            ct.compare_objects(dc.combined_table)
            common_option_handler(ctx, dc)


@main.command()
@click.argument('coord', type=str)
@click.option('-match-tol', type=float,
              help='Optional 2D match tolerance (in arc seconds). Default value is 1.0 arcsec.')
@click.option('-obj-radius', type=float,
              help='Optional radius to search around object (in arc minutes). Default value is 1.0 arcmin.')
@click.pass_context
def bycoord(ctx, coord, match_tol, obj_radius):
    '''
    Downloads objects, via NED and SIMBAD, from region described by coordinates

    Arguments:
        COORD - the location of coordinates to be searched around
    '''
    ctx.obj['name'] =  coord

    dc = ctx.obj['dc']
    ct = ctx.obj['ct']

    logging.info("Query requested: {}".format(coord))

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
        logging.error("mmatch_tol must be between 0 and 60 arcsec")
    elif obj_radius_valid == False:
        logging.error("Invalid obj_radius entered.")
        logging.error("obj_radius must be between 0 and 60 arcmin")
    else:
        if match_tol and obj_radius:
            logging.info("Confirm: match-tol & obj-radius passed.")
            dc.query_region(coord, match_tol=match_tol, obj_radius=obj_radius, bycoord=True)
        elif match_tol:
            logging.info("Confirm: match-tol passed.")
            dc.query_region(coord, match_tol=match_tol, bycoord=True)
        elif obj_radius:
            logging.info("Confirm: obj-radius passed.")
            dc.query_region(coord, obj_radius=obj_radius, bycoord=True)
        else:
            logging.info("Default settings used for bycoord query.")
            dc.query_region(coord, bycoord=True)

        if dc.combined_table is not None:
            # Pass table to comparison tree to compare each object
            ct.compare_objects(dc.combined_table)
            common_option_handler(ctx, dc)


'''
Not yet implemented. Local data less of a priority. Perhaps future versions.


@main.command()
@click.argument('ned-name', type=str)
@click.argument('simbad-name', type=str)
@click.pass_context
def local(ctx, ned_name, simbad_name):
    # Uses local data provided to compare objects in NED and SIMBAD
    combined_table = os.path.join(package_dir, "data", "M83_NScomb_mar1.fits")

    dc = ctx.obj['dc']
    ct = ctx.obj['ct']

    dc.load_data(combined_table)

    # Pass table to comparison tree to compare each object
    ct.compare_objects(dc.combined_table)

    if dc.combined_table is not None:
        dc.combined_table.show_in_browser(jsviewer=True)
'''

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
        os.system("open " + package_dir + filepath)
    elif (sys.platform == "cygwin" or sys.platform == "win32"):  # Windows
        os.system("start " + package_dir + filepath)
    elif (sys.platform == "linux"):
        os.system("xdg-open " + package_dir + filepath)


def common_option_handler(ctx, dc):
    '''
    Handle options/operations relevant to all query types.
    '''
    # Stats are always shown.
    dc.stats.derive_table_stats(dc.combined_table)
    stats_output = dc.stats.generateStatTemplate()

    objname = ctx.obj['name']

    # Save file variables.
    if ctx.obj["savetable"] or ctx.obj['saveplot'] or ctx.obj['savestats']:
        filename = ctx.obj['filename']

    # If table option(s) present.
    if ctx.obj["showtable"]:
        dc.combined_table.show_in_browser(jsviewer=True)
    if ctx.obj["savetable"]:
        dc.saveTable(fileName=working_dir+'/gdm_output'+'/'+objname+'-'+filename,
                     file_format=ctx.obj['savetable'])

    # If plot option(s) present.
    if ctx.obj['showplot'] or ctx.obj['saveplot']:
        plot = DataController.plot_match_table(dc.combined_table, name=objname)
        if ctx.obj['saveplot']:
            plot.savefig(working_dir+'/gdm_output'+'/'+objname+'-'+filename+'.png')
        if ctx.obj['showplot']:
            plot.show()

    if ctx.obj['savestats']:
        with open(working_dir+'/gdm_output'+'/'+objname+"-stats-"+filename+".txt", "w") as text_file:
            text_file.write("{}".format(stats_output))

    # If glossary option present.
    if ctx.obj['glossary']:
        safelyopenfile('/data/glossary.txt')

# ----------- main entry point | script run directly ----------- #


if __name__ == "__main__":
    main(obj={})
