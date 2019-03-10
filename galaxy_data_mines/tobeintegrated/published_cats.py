from astropy.table import Table, Column, hstack, vstack
from astropy.coordinates import SkyCoord, match_coordinates_sky
from astropy import units as u
import numpy as np
import os

# published_cats:
#  tools for manipulating published catalogs of m83 objects, for comparison
#  with Chandar et al catalog
#     - reformat and combine NED and SIMBAD data tables
#     - add in tables of data not included in NED or SIMBAD
#     - match to Chandar catalog to produce a list of best-matches with object types
#

# usage:
# published_cats.ns_combine('ned-20160629.fits','simbad-20160629.fits','M83_NScomb.fits','M83_NSall.fits')
# published_cats.add_tables('M83_NSall.fits',['williams15_rsg.fits','kim12_wr.fits'],'M83_final.fits')
# published_cats.catalog_match('M83_final.fits', 'hlsp_wfc3ers_hst_wfc3_m83_cat_all_v2-corrected.txt','M83_ers_pubcat.txt')

ned_rename = [('Name_N', 'Name'), ('RA(deg)', 'RA'), ('DEC(deg)', 'Dec'), ('Type_N', 'Type')]
sim_rename = [('Name_S', 'Name'), ('RA_d', 'RA'), ('DEC_d', 'Dec'), ('Type_S', 'Type')]

# ---------------- ns_combine AND helper functions ---------------- #


def ns_combine(ned_name, simbad_name, ns_combine, final_tab, match_tol=1.0):  # match_tol in arcsec

    ned_in = Table.read(ned_name)
    simbad_in = Table.read(simbad_name)

    # prelim processing
    ned_proc = reformat_cat(ned_in, old_name='Object Name', new_name='Name_N', old_type='Type', new_type='Type_N',
                            keepcols=['Object Name', 'RA(deg)', 'DEC(deg)', 'Type'])
    sim_proc = reformat_cat(simbad_in, old_name='MAIN_ID', new_name='Name_S',
                            old_type='OTYPE', new_type='Type_S')

    # construct coordinates needed for matching # ***** MATT - created SkyCoord's w correct unit columns
    ned_coo = SkyCoord(ra=ned_proc['RA(deg)'], dec=ned_proc['DEC(deg)'])
    sim_coo = SkyCoord(ra=sim_proc['RA_d'], dec=sim_proc['DEC_d'])

    # do the matching # ***** MATT - Returns indices of matched col's for ned+sim tables
    matched_ned, matched_sim, ned_only, sim_only = symmetric_match_sky_coords(
        ned_coo, sim_coo, match_tol*u.arcsec)

    print("Matched NED column:")
    print(ned_proc[matched_ned])
    print("Matched SIMBAD column:")
    print(sim_proc[matched_sim])
    print("Unmatched NED:")
    print(ned_proc[ned_only])
    print("Unmatched SIMBAD:")
    print(sim_proc[sim_only])

    # generate the matched table
    matchtab = hstack([ned_proc[matched_ned], sim_proc[matched_sim]], join_type='outer')
    # mark the really good matches
    matchtab2 = process_match(matchtab)
    matchtab2.write(ns_combine, format='fits')
    # rename some columns
    nedcat = process_unmatch(ned_proc[ned_only], src='N', rename_cols=ned_rename)
    simcat = process_unmatch(sim_proc[sim_only], src='S', rename_cols=sim_rename)
    keeplist = ['Name_N', 'RA(deg)', 'DEC(deg)', 'Type_N']
    matchtab3 = process_unmatch(Table(matchtab2[keeplist]), src='NS', rename_cols=ned_rename)

    # add on the unmatched objects
    finaltab = vstack([matchtab3, nedcat, simcat], join_type='outer')

    # save the result
    finaltab.write(final_tab, format='fits')

    return

# ---------------- ns_combine's helpers ---------------- #


ns_replace_names = [(" ", ""), ("MESSIER083:", ""), ("NGC5236:", ""), ("M83-", ""), ("NGC5236", "")]
ns_replace_types = [('*Cl', 'Cl*'), ('PofG', 'Galaxy'), ('X', 'XrayS'), ('Radio', 'RadioS')]
ns_remove_ids = ['NAMENGC5236Group', 'M83', 'MESSIER083', 'NGC5236GROUP']
ra_dec_cols = ['RA(deg)', 'DEC(deg)', 'RA_d', 'DEC_d']


def reformat_cat(in_tab, old_name, new_name, old_type, new_type, replace_names=ns_replace_names, replace_types=ns_replace_types, remove_id=ns_remove_ids, keepcols=None):
    ''' reformat NED or SIMBAD catalog to make more intercompatible'''
    # just keep selected columns
    if keepcols != None:
        in_tab = in_tab[keepcols]

    # change units for RA/Dec
    for col in ra_dec_cols:
        if col in in_tab.colnames:
            in_tab[col].unit = u.degree

    # change ID for name & type columns
    in_tab.rename_column(old_name, new_name)
    in_tab.rename_column(old_type, new_type)

    '''
    # reformat some object names
    for pair in replace_names:
        in_tab[new_name] = np.char.replace(in_tab[new_name], pair[0], pair[1])
        '''
    # reformat some object types # ***** MATT S. Can remove this code, will be handled by dictionary
    #in_tab[new_type] = np.char.replace(in_tab[new_type], "?", "")
    #in_tab[new_type] = np.char.replace(in_tab[new_type], " ", "")
    # for pair in replace_types:
    #    in_tab[new_type][in_tab[new_type] == pair[0]] = pair[1]

    # delete rows whose names are in remove_id
    # there's a non-loopy way to do this but I can't remember it
    remove_idx = []
    for i in range(0, len(in_tab)):
        if in_tab[i][new_name] in remove_id:
            remove_idx.append(i)
    in_tab.remove_rows(remove_idx)

    # all done
    return(in_tab)


def symmetric_match_sky_coords(coord1, coord2, tolerance):
    '''produce the symmetric match of coord1 to coord2
       output:
       index1_matched: index into coord1 for matched objects
       index2_matched: index into coord2 for matches of objects in index1_matched
       index1_unmatch: indices for unmatched objects in coord1
       index2_unmatch: indices for unmatched objects in coord2
    '''
    closest_2to1, sep2d_2to1, sep3d = match_coordinates_sky(
        coord1, coord2)  # location in coord2 for closest match to each coord1. len = len(coord1)
    # location in coord1 for closest match to each coord2. len = len(coord2)
    closest_1to2, sep2d_1to2, sep3d = match_coordinates_sky(coord2, coord1)

    index1_matched = []
    index2_matched = []
    index1_unmatched = []
    index2_unmatched = []

    for i in range(0, len(coord1)):  # doubtless there is a more Pythonic way to do this..
        # not sure this condition covers all of the possible corner cases. But good enough.
        if sep2d_1to2[i] < tolerance and i == closest_2to1[closest_1to2[i]]:
            index1_matched.append(i)
            index2_matched.append(closest_2to1[i])
        else:
            index1_unmatched.append(i)

    for j in range(0, len(coord2)):
        if j not in index2_matched:
            index2_unmatched.append(j)

    return(index1_matched, index2_matched, index1_unmatched, index2_unmatched)


def process_match(matched_tab_in):
    '''find secure matches btw NED and SIMBAD'''
    goodmatch = np.logical_or(matched_tab_in['Name_S'] == matched_tab_in['Name_N'],
                              matched_tab_in['Type_S'] == matched_tab_in['Type_N'])
    matched_tab_in.add_column(Column(goodmatch, name='Secure'))
    return(matched_tab_in)


def process_unmatch(tab_in, src, rename_cols):
    '''find secure matches btw NED and SIMBAD'''
    for pair in rename_cols:
        tab_in.rename_column(pair[0], pair[1])
    tab_in.add_column(Column(name='Source', length=len(tab_in), dtype='S2'))
    tab_in['Source'] = src
    return(tab_in)

# ------------------------- Add tables ------------------------- #


def add_tables(basetab_file, tab_file_list, outfile, jt='outer'):
    basetab = Table.read(basetab_file)
    tablist = [basetab]
    for filename in tab_file_list:
        tab = Table.read(filename)
        tablist.append(tab)
    stack = vstack(tablist, join_type=jt)
    stack.write(outfile)
    return

# ------------------------ match up objects ------------------------ #


def catalog_match(pubcat_file, erscat_file, match_out_file, match_tol=1.0):
    pubcat = Table.read(pubcat_file, format='fits')
    erscat = Table.read(erscat_file, format='ascii.commented_header')

    # construct coordinates needed for matching
    pub_coo = SkyCoord(ra=pubcat['RA'], dec=pubcat['Dec'])
    ers_coo = SkyCoord(ra=erscat['ra']*u.degree, dec=erscat['dec']*u.degree)

    # do the matching
#    closest_2to1, sep2d_2to1, sep3d = match_coordinates_sky(coord1, coord2) # location in coord2 for closest match to each coord1. len = len(coord1)
    # location in coord2 for closest match to each coord1. len = len(coord1)
    closest, sep2d, sep3d = match_coordinates_sky(pub_coo, ers_coo)
    matched = sep2d < match_tol*u.arcsec
#    matched_ers, matched_pub, ers_only, pub_only = symmetric_match_sky_coords(ers_coo, pub_coo, match_tol*u.arcsec)

    # generate the matched table
    keeplist = ['id_', 'ra', 'dec']
    tmpcat = Table(erscat[keeplist])
    matchtab = hstack([tmpcat[closest][matched], pubcat[matched]], join_type='outer')

    # write the matched catalog to a file
    matchtab.write(match_out_file, format='ascii.commented_header')

    return


# BELOW HERE IS OLD STUFF, not used
# --------------------------------------------------------
within_galaxy_types_ned = ['*Cl', 'HII', 'PofG', 'Neb', 'SN', 'SNR', 'V*', 'WR*']
background_types_ned = ['G', 'GClstr', 'GGroup', 'QSO']
obs_types_ned = ['IrS', 'RadioS', 'UvES', 'UvS', 'VisS', 'XrayS']

within_galaxy_types_simbad = ['**', 'Assoc*', 'Candidate_SN*', 'Candidate_WR*', 'Cepheid', 'Cl*', 'HII', 'ISM',
                              'LMXB', 'MolCld', 'Nova', 'PN', 'PartofG', 'SN', 'SNR', 'SNR?', 'Star', 'ULX', 'V*', 'WR*', 'semi-regV*']
background_types_simbad = ['BLLac', 'ClG', 'EmG', 'Galaxy',
                           'GinCl', 'GroupG', 'Possible_G', 'StarburstG']
obs_types_simbad = ['EmObj', 'IR', 'Radio', 'Radio(sub-mm)', 'UV', 'X']


def name_match(simbad_name, ned_name):
    matched = np.zeros(len(simbad_name), dtype='bool')
    matched[np.char.replace(simbad_name, " ", "") == np.char.replace(ned_name, " ", "")] = True
    # TODO: account for cases where one name has leading zero in an ID (eg [WLC2001]03) and other doesn't

    return(matched)

# not using this


def process_tab(tab_in, tab_out, type_col, select_list=within_galaxy_types_ned, rfmt_fn=None):
    '''select specific object types from a table'''
    tab = Table.read(tab_in)
    if type_col != 'Type':
        tab.rename_column(type_col, 'Type')
    tab['Type'] = np.char.strip(tab['Type'])  # remove whitespace -- helps with matching
    tg = tab.group_by('Type')
    mask = np.in1d(tg.groups.keys['Type'], select_list)  # create mask for only wanted types
    wanted_types = tg.groups[mask]

    # do some reformatting
    if rfmt_fn != None:
        wanted_types = rfmt_fn(wanted_types)

    if tab_out != None:  # write to file
        if os.path.exists(tab_out):
            os.unlink(tab_out)
        wanted_types.write(tab_out, format='fits')
    return(wanted_types)
