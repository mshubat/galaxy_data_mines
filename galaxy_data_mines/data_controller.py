# Data Storage and Organization Class
#
# This class serves as an import, organization, and processing tool for NED
# and SIMBAD data.
#
# Possible data sources:
# - local - on machine datasets
# - remote - via astroquery quieries to NED and SIMBAD
#

from .match_stats import MatchStats
from astropy.table import Table, Column, hstack, vstack
from astroquery.simbad import Simbad
from astroquery.ned import Ned
from astropy import units as u
from astropy.coordinates import SkyCoord, match_coordinates_sky, get_icrs_coordinates
from astropy.coordinates.name_resolve import sesame_database, NameResolveError
from astropy.io import fits
from requests.exceptions import Timeout

import matplotlib.pyplot as plt
import logging
import warnings
# from astropy.table import Table, Column, hstack, vstack


class DataController:

    ned_to_simbad_cond_dict = {
        "*": "*",
        "**": "**",
        "*Ass": "As*",
        "*Cl": "Cl*",
        "AbLS": "ALS",
        "Blue*": ("BS*", "s*b"),
        "C*": "C*",
        "EmLS": ("Em*", "EmG"),  # Emission line source could be galaxy or star...
        "EmObj": "EmO",
        "exG*": "*",
        "Flare*": "Fl*",
        "G": "G",
        "GammaS": "gam",
        "GClstr": "ClG",
        "GGroup": "GrG",
        "GPair": "PaG",
        "GTrpl": "GrG",  # Triple not exactly the same as a group..
        "G_Lens": "LeG",
        "HII": "HII",
        "IrS": "IR",
        "MCld": "MoC",
        "Neb": "Cld",
        "Nova": "No*",
        "Other": ("?", "err"),  # unknown or non-existent
        "PN": "PN",
        "PofG": "PoG",
        "Psr": "Psr",
        "QGroup": "",  # not too common
        "QSO": "QSO",
        "Q_Lens": "LeQ",
        "RadioS": "Rad",
        "Red*": ("RG*", "s*r"),
        "RfN": "RNe",
        "SN": "SN*",
        "SNR": "SNR",
        "UvES": "UV",
        "UvS": "UV",
        "V*": "V*",
        "VisS": "",  # no direct match found
        "WD*": "WD*",
        "WR*": "WR*",
        "XrayS": "X",
        # Generally, SIMBAD does not distinguish between "Galactic"
        # variants of objects
        "!*": "*",
        "!**": "**",
        "!*Ass": "As*",
        "!*Cl": "Cl*",
        "!Blue*": ("BS*", "s*b"),
        "!C*": "C*",
        "!EmObj": "EmO",
        "!Flar*": "Fl*",
        "!HII": "HII",
        "!MCld": "MoC",
        "!Neb": "GNe",
        "!Nova": "No*",
        "!PN": "PN",
        "!Psr": "Psr",
        "!RfN": "RNe",
        "!Red*": ("RG*", "s*r"),
        "!SN": "SN*",
        "!SNR": "SNR",
        "!V*": "V*",
        "!WD*": "WD*",
        "!WR*": "WR*"}

    simbad_std_to_cond = {
        "Unknown": "?",
        "Transient": "ev",
        "Radio": "Rad",
        "Radio(m)": "mR",
        "Radio(cm)": "cm",
        "Radio(mm)": "mm",
        "Radio(sub-mm)": "smm",
        "HI": "HI",
        "radioBurst": "rB",
        "Maser": "Mas",
        "IR": "IR",
        "IR>30um": "FIR",
        "IR<10um": "NIR",
        "Red": "red",
        "RedExtreme": "ERO",
        "Blue": "blu",
        "UV": "UV",
        "X": "X",
        "ULX?": "UX?",
        "ULX": "ULX",
        "gamma": "gam",
        "gammaBurst": "gB",
        "Inexistent": "err",
        "Gravitation": "grv",
        "LensingEv": "Lev",
        "Candidate_LensSystem": "LS?",
        "Candidate_Lens": "Le?",
        "Possible_lensImage": "LI?",
        "GravLens": "gLe",
        "GravLensSystem": "gLS",
        "Candidates": "..?",
        "Possible_G": "G?",
        "Possible_SClG": "SC?",
        "Possible_ClG": "C?G",
        "Possible_GrG": "Gr?",
        "Possible_As*": "As?",
        "Candidate_**": "**?",
        "Candidate_EB*": "EB?",
        "Candidate_Symb*": "Sy?",
        "Candidate_CV*": "CV?",
        "Candidate_Nova": "No?",
        "Candidate_XB*": "XB?",
        "Candidate_LMXB": "LX?",
        "Candidate_HMXB": "HX?",
        "Candidate_Pec*": "Pec?",
        "Candidate_YSO": "Y*?",
        "Candidate_pMS*": "pr?",
        "Candidate_TTau*": "TT?",
        "Candidate_C*": "C*?",
        "Candidate_S*": "S*?",
        "Candidate_OH": "OH?",
        "Candidate_CH": "CH?",
        "Candidate_WR*": "WR?",
        "Candidate_Be*": "Be?",
        "Candidate_Ae*": "Ae?",
        "Candidate_HB*": "HB?",
        "Candidate_RRLyr": "RR?",
        "Candidate_Cepheid": "Ce?",
        "Candidate_RGB*": "RB?",
        "Candidate_SG*": "sg?",
        "Candidate_RSG*": "s?r",
        "Candidate_YSG*": "s?y",
        "Candidate_BSG*": "s?b",
        "Candidate_AGB*": "AB?",
        "Candidate_LP*": "LP?",
        "Candidate_Mi*": "Mi?",
        "Candiate_sr*": "sv?",
        "Candidate_post-AGB*": "pA?",
        "Candidate_BSS": "BS?",
        "Candidate_Hsd": "HS?",
        "Candidate_WD*": "WD?",
        "Candidate_NS": "N*?",
        "Candidate_BH": "BH?",
        "Candidate_SN*": "SN?",
        "Candidate_low-mass*": "LM?",
        "Candidate_brownD*": "BD?",
        "multiple_object": "mul",
        "Region": "reg",
        "Void": "vid",
        "SuperClG": "SCG",
        "ClG": "ClG",
        "GroupG": "GrG",
        "Compact_Gr_G": "CGG",
        "PairG": "PaG",
        "IG": "IG",
        "Cl*?": "C?*",
        "GlCl?": "Gl?",
        "Cl*": "Cl*",
        "GlCl": "GlC",
        "OpCl": "OpC",
        "Assoc*": "As*",
        "Stream*": "St*",
        "MouvGroup": "MGr",
        "**": "**",
        "EB*": "EB*",
        "EB*Algol": "Al*",
        "EB*betLyr": "bL*",
        "EB*WUMa": "WU*",
        "EB*Planet": "EP*",
        "SB*": "SB*",
        "EllipVar": "El*",
        "Symbiotic*": "Sy*",
        "CataclyV*": "CV*",
        "DQHer": "DQ*",
        "AMHer": "AM*",
        "Nova-like": "NL*",
        "Nova": "No*",
        "DwarfNova": "DN*",
        "XB": "XB*",
        "LMXB": "LXB",
        "HMXB": "HXB",
        "ISM": "ISM",
        "PartofCloud": "PoC",
        "PN?": "PN?",
        "ComGlob": "CGb",
        "Bubble": "bub",
        "EmObj": "EmO",
        "Cloud": "Cld",
        "GalNeb": "GNe",
        "BrNeb": "BNe",
        "DkNeb": "DNe",
        "RfNeb": "RNe",
        "MolCld": "MoC",
        "Globule": "glb",
        "denseCore": "cor",
        "SFregion": "SFR",
        "HVCld": "HVC",
        "HII": "HII",
        "PN": "PN",
        "HIshell": "sh",
        "SNR?": "SR?",
        "SNR": "SNR",
        "Circumstellar": "cir",
        "outflow?": "of?",
        "Outflow": "out",
        "HH": "HH",
        "Star": "*",
        "*inCl": "*iC",
        "*inNeb": "*iN",
        "*inAssoc": "*iA",
        "*in**": "*i*",
        "V*?": "V*?",
        "Pec*": "Pe*",
        "HB*": "HB*",
        "YSO": "Y*O",
        "Ae*": "Ae*",
        "Em*": "Em*",
        "Be*": "Be*",
        "BlueStraggler": "BS*",
        "RGB*": "RG*",
        "AGB*": "AB*",
        "C*": "C*",
        "S*": "S*",
        "SG*": "sg*",
        "RedSG*": "s*r",
        "YellowSG*": "s*y",
        "BlueSG*": "s*b",
        "HotSubdwarf": "HS*",
        "post-AGB*": "pA*",
        "WD*": "WD*",
        "pulsWD*": "ZZ*",
        "low-mass*": "LM*",
        "brownD*": "BD*",
        "Neutron*": "N*",
        "OH/IR": "OH*",
        "CH": "CH*",
        "pMS*": "pr*",
        "TTau*": "TT*",
        "WR*": "WR*",
        "PM*": "PM*",
        "HV*": "HV*",
        "V*": "V*",
        "Irregular_V*": "Ir*",
        "Orion_V*": "Or*",
        "Rapid_Irreg_V*": "RI*",
        "Eruptive*": "Er*",
        "Flare*": "Fl*",
        "FUOr": "FU*",
        "Erupt*RCrB": "RC*",
        "RCrB_Candidate": "RC?",
        "RotV*": "Ro*",
        "RotV*alf2CVn": "a2*",
        "Pulsar": "Psr",
        "BYDra": "BY*",
        "RSCVn": "RS*",
        "PulsV*": "Pu*",
        "RRLyr": "RR*",
        "Cepheid": "Ce*",
        "PulsV*delSct": "dS*",
        "PulsV*RVTau": "RV*",
        "PulsV*WVir": "WV*",
        "PulsV*bCep": "bC*",
        "deltaCep": "cC*",
        "gammaDor": "gD*",
        "pulsV*SX": "SX*",
        "LPV*": "LP*",
        "Mira": "Mi*",
        "semi-regV*": "sr*",
        "SN": "SN*",
        "Sub-stellar": "su*",
        "Planet?": "Pl?",
        "Planet": "Pl",
        "Galaxy": "G",
        "PartofG": "PoG",
        "GinCl": "GiC",
        "BClG": "BiC",
        "GinGroup": "GiG",
        "GinPair": "GiP",
        "High_z_G": "HzG",
        "AbsLineSystem": "ALS",
        "Ly-alpha_ALS": "LyA",
        "DLy-alpha_ALS": "DLA",
        "metal_ALS": "mAL",
        "Ly-limit_ALS": "LLS",
        "Broad_ALS": "BAL",
        "RadioG": "rG",
        "HII_G": "H2G",
        "LSB_G": "LSB",
        "AGN_Candidate": "AG?",
        "QSO_Candidate": "Q?",
        "Blazar_Candidate": "Bz?",
        "BLLac_Candidate": "BL?",
        "EmG": "EmG",
        "StarburstG": "SBG",
        "BlueCompG": "bCG",
        "LensedImage": "LeI",
        "LensedG": "LeG",
        "LensedQ": "LeQ",
        "AGN": "AGN",
        "LINER": "LIN",
        "Seyfert": "SyG",
        "Seyfert_1": "Sy1",
        "Seyfert_2": "Sy2",
        "Blazar": "Bla",
        "BLLac": "BLL",
        "OVV": "OVV",
        "QSO": "QSO",
    }

    candidate_dict = {
        "ULX": "UX?",
        "gLS": "LS?",
        "gLe": "Le?",
        "LeI": "LI?",
        "G": "G?",
        "SCG": "SC?",
        "ClG": "C?G",
        "GrG": "Gr?",
        "**": "**?",
        "EB*": "EB?",
        "Sy*": "Sy?",
        "CV*": "CV?",
        "No*": "No?",
        "XB*": "XB?",
        "LXB": "LX?",
        "HXB": "HX?",
        "Pe*": "Pec?",
        "Y*O": "Y*?",
        "pr*": "pr?",
        "TT*": "TT?",
        "C*": "C*?",
        "S*": "S*?",
        "OH*": "OH?",
        "CH*": "CH?",
        "WR*": "WR?",
        "Be*": "Be?",
        "Ae*": "Ae?",
        "HB*": "HB?",
        "RR*": "RR?",
        "Ce*": "Ce?",
        "RG*": "RB?",
        "sg*": "sg?",
        "s*r": "s?r",
        "s*y": "s?y",
        "s*b": "s?b",
        "AB*": "AB?",
        "LP*": "LP?",
        "Mi*": "Mi?",
        "sr*": "sv?",
        "pA*": "pA?",
        "BS*": "BS?",
        "WD*": "WD?",
        "N*": "N*?",
        "SN*": "SN?",
        "LM*": "LM?",
        "BD*": "BD?",
        "Cl*": "C?*",
        "GlC": "Gl?",
        "PN": "PN?",
        "SNR": "SR?",
        "out": "of?",
        "V*": "V*?",
        "RC*": "RC?",
        "Pl": "Pl?",
        "AGN": "AG?",
        "QSO": "Q?",
        "Bla": "Bz?",
        "BLL": "BL?",
    }

    simbad_table_cols = ["MAIN_ID", "OTYPE", "RA_d", "DEC_d"]
    simbad_table_dtypes = ["object","object","float64","float64"]

    def __init__(self, st):
        self.combined_table = None
        self.stats = MatchStats(template_file=st)

    @staticmethod
    def ned_to_simbad(ned_entry):
        if type(ned_entry) is bytes:
            ned_entry = ned_entry.decode("utf-8")

        if ned_entry == "":
            logging.info("Empty string was passed as NED object type.")
            return ""

        return DataController.ned_to_simbad_cond_dict[ned_entry]

    @staticmethod
    def simbad_long_to_small(simbad_std):
        if type(simbad_std) is bytes:
            simbad_std = simbad_std.decode("utf-8")

        if simbad_std == "":
            logging.info("Empty string was passed as SIMBAD object type.\
            SIMBAD must not have given this object a classification.")
            return ""

        return DataController.simbad_std_to_cond[simbad_std]

    @staticmethod
    def candidate_match(non_candidate):
        if type(non_candidate) is bytes:
            non_candidate = non_candidate.decode("utf-8")

        if non_candidate == "":
            logging.info("Empty string was passed during candidate lookup.")
            return ""

        return DataController.candidate_dict[non_candidate]

    def query_region(self, objectname, match_tol=1.0, obj_radius=1.0, bycoord=False):
        '''
        Fetch remote data from NED and SIMBAD matching coordinates and build table.
        '''
        # Create custom query objects.
        customSimbad = Simbad()
        customNed = Ned()

        # Log SIMBAD votable (changeable) fields.
        logging.debug("SIMBAD votable fields")
        logging.debug(customSimbad.get_votable_fields())

        customSimbad.remove_votable_fields('coordinates')
        # customSimbad.add_votable_fields("otype(3)", "ra(d)", "dec(d)")
        customSimbad.add_votable_fields("otype", "ra(d)", "dec(d)")

        # Download object data from both SIMBAD and NED.
        logging.info("Querying SIMBAD and NED for region {}".format(objectname))

        if bycoord:
            try:
                objectcoords = SkyCoord(objectname)
            except (ValueError, u.UnitsError):
                 logging.info("Invalid coordinates.")
                 return
        else:
        # Resolve the object name into sky coordinate using NED
        # ensures that NED and SIMBAD searches are using the same position
            sesame_database.set('ned')
            try:
                objectcoords = get_icrs_coordinates(objectname)
            except NameResolveError:
                 logging.info("Name resolution failed.")
                 return

        logging.info("Name resolved to coordinates {}".format(objectcoords))

        # SIMBAD
        logging.info("SIMBAD is currently being queried...")
        try:
            with warnings.catch_warnings(): # suppress warnings generated by SIMBAD query
                warnings.simplefilter("ignore")
                simbad_table = customSimbad.query_region(objectcoords, radius=obj_radius*u.arcmin)
        # workaround. If SIMBAD query finds nothing, returns None but we want a zero-length table
            if type(simbad_table) is not Table:
                logging.debug("No SIMBAD objects")
                simbad_table = Table(data=None,
                                   names=DataController.simbad_table_cols,
                                   dtype=DataController.simbad_table_dtypes,
                                   masked=True)
            logging.info("SUCCESS: SIMBAD Data retrieved.")
        except Timeout:
            logging.debug("SIMBAD timeout error")
            return

       # NED
        logging.info("NED is currently being queried...")
        try:
            ned_table = Ned.query_region(objectcoords, radius=obj_radius*u.arcmin)
            logging.info("SUCCESS: NED Data retrieved.")
        except Timeout:
            logging.debug("NED timeout error")
            return

# Save some query stats.
        self.stats.query_name = objectname
        self.stats.sim_count = len(simbad_table)
        self.stats.ned_count = len(ned_table)

        # process tables
        ned_table = self.reformat_table(ned_table,
                                        keepcolsifpresent=['Object Name',
                                                           # cover NED changing names of cols
                                                           'RA(deg)',
                                                           'RA',
                                                           'DEC(deg)',
                                                           'DEC',
                                                           'Type'],
                                        old_name='Object Name', new_name='Name_N',
                                        old_type='Type', new_type='Type_N')

        logging.info("Reformatting tables.")
        simbad_table = self.reformat_table(simbad_table,
                                           keepcolsifpresent=["MAIN_ID", "RA_d", "DEC_d", "OTYPE"],
                                           old_name='MAIN_ID', new_name='Name_S',
                                           old_type='OTYPE', new_type='Type_S')

        logging.info("Building sky coordinates.")
        # Build SkyCoord from appropriate ned and simbad col's with matching units
        ned_coo = SkyCoord(ra=ned_table['RA(deg)'], dec=ned_table['DEC(deg)'])
        sim_coo = SkyCoord(ra=simbad_table['RA_d'], dec=simbad_table['DEC_d'])

        logging.info("Finding object matches.")
        # Find object matches
        if len(ned_coo)> 0 and len(sim_coo) > 0:
            matched_ned, matched_sim, ned_only, sim_only = self.symmetric_match_sky_coords_v2(
                ned_coo, sim_coo, match_tol*u.arcsec)
        else:
             matched_ned = []
             matched_sim = []
             ned_only = []
             sim_only = []
        logging.debug("")
        logging.debug("Matched NED rows:")
        logging.debug(ned_table[matched_ned])
        logging.debug("Matched SIMBAD rows:")
        logging.debug(simbad_table[matched_sim])
        logging.debug("")

        self.stats.overlap_count = len(matched_ned)

        # Explore results
        logging.debug("Matched NED:")
        logging.debug(matched_ned)
        logging.debug("Matched SIMBAD")
        logging.debug(matched_sim)
        logging.debug("NED ONLY")
        logging.debug(ned_only)
        logging.debug("SIMBAD ONLY")
        logging.debug(sim_only)

        # Generate the matched table and save the result.
        logging.info("Building combined table.")
        matched_table = hstack(
            [ned_table[matched_ned], simbad_table[matched_sim]],
            join_type='outer',
            metadata_conflicts='silent')  # Hide the metadata warning.
        self.combined_table = matched_table

    def query_region_by_coord(self, coord_type, RA, DEC):
        pass

    def reformat_table(self, table, keepcolsifpresent, old_name, new_name, old_type, new_type):
        '''Reformat NED or SIMBAD catalog to make more intercompatible'''

        ra_dec_cols = ['RA(deg)', 'DEC(deg)', 'RA', 'DEC', 'RA_d', 'DEC_d']

        # Just keep selected columns.
        keepcols = []
        if keepcolsifpresent != None:
            for col in keepcolsifpresent:
                if col in table.colnames:
                    keepcols.append(col)
            table = table[keepcols]

        # Change units for RA/DEC.
        for col in ra_dec_cols:
            if col in table.colnames:
                table[col].unit = u.degree

        # Change ID for name & type columns.
        table.rename_column(old_name, new_name)
        table.rename_column(old_type, new_type)

        # If ned name changes, revert to consistent name.
        if 'DEC' in table.colnames:
            table.rename_column('DEC', 'DEC(deg)')
        if 'RA' in table.colnames:
            table.rename_column('RA', 'RA(deg)')

        return(table)

    def symmetric_match_sky_coords(self, coord1, coord2, tolerance):
        '''produce the symmetric match of coord1 to coord2
           output:
           index1_matched: index into coord1 for matched objects
           index2_matched: index into coord2 for matches of objects in index1_matched
           index1_unmatch: indices for unmatched objects in coord1
           index2_unmatch: indices for unmatched objects in coord2
        '''
        closest_2to1, sep2d_2to1, sep3d = match_coordinates_sky(
            coord1, coord2)  # indices for "coord2" for closest match to each coord1. len = len(coord1)
        # location in coord1 for closest match to each coord2. len = len(coord2)
        closest_1to2, sep2d_1to2, sep3d = match_coordinates_sky(coord2, coord1)

        index1_matched = []
        index2_matched = []
        index1_unmatched = []
        index2_unmatched = []

        logging.debug("DEBUG STATEMENTS:")
        logging.debug("tolerance = {}".format(tolerance))
        logging.debug("len(sep2d_2to1) = {}".format(len(sep2d_2to1)))
        logging.debug("len(sep2d_1to2) = {}".format(len(sep2d_1to2)))
        logging.debug("len(closest_2to1) = {}".format(len(closest_2to1)))
        logging.debug("len(closest_1to2) = {}".format(len(closest_1to2)))
        logging.debug("len(coord1) = {}".format(len(coord1)))
        logging.debug("len(coord2) = {}".format(len(coord2)))

        for i in range(0, len(coord1)):  # doubtless there is a more Pythonic way to do this..
            # not sure this condition covers all of the possible corner cases. But good enough.
            logging.debug("-------------------")
            logging.debug("iteration i = {}".format(i))
            logging.debug("-------------------")
            # logging.debug("sep2d_1to2[i] = {}".format(sep2d_1to2[i]))
            # logging.debug("closest_1to2[i] = {}".format(closest_1to2[i]))

            if sep2d_1to2[i] < tolerance and i == closest_2to1[closest_1to2[i]]:
                index1_matched.append(i)
                index2_matched.append(closest_2to1[i])
            else:
                index1_unmatched.append(i)

        for j in range(0, len(coord2)):
            if j not in index2_matched:
                index2_unmatched.append(j)

        return(index1_matched, index2_matched, index1_unmatched, index2_unmatched)

    def symmetric_match_sky_coords_v2(self, coord1, coord2, tolerance):
        '''produce the symmetric match of coord1 to coord2
           output:
           index1_matched: index into coord1 for matched objects
           index2_matched: index into coord2 for matches of objects in index1_matched
           index1_unmatch: indices for unmatched objects in coord1
           index2_unmatch: indices for unmatched objects in coord2
        '''
        closest_2to1, sep2d_2to1, sep3d = match_coordinates_sky(
            coord1, coord2)  # indices for "coord2" for closest match to each coord1. len = len(coord1)
        # location in coord1 for closest match to each coord2. len = len(coord2)
        closest_1to2, sep2d_1to2, sep3d = match_coordinates_sky(coord2, coord1)

        index1_matched = []
        index2_matched = []
        index1_unmatched = []
        index2_unmatched = []

        logging.debug("DEBUG STATEMENTS:")
        logging.debug("tolerance = {}".format(tolerance))
        logging.debug("len(sep2d_2to1) = {}".format(len(sep2d_2to1)))
        logging.debug("len(sep2d_1to2) = {}".format(len(sep2d_1to2)))
        logging.debug("len(closest_2to1) = {}".format(len(closest_2to1)))
        logging.debug("len(closest_1to2) = {}".format(len(closest_1to2)))
        logging.debug("len(coord1) = {}".format(len(coord1)))
        logging.debug("len(coord2) = {}".format(len(coord2)))

        # ----------- Matt's attempt ----------- #
        # -------------------------------------- #

        if len(coord1) < len(coord2):
            shortest_len = len(coord1)
            longest_len = len(coord2)
            shortest = "coord1"
            longest = "coord2"
        else:
            shortest_len = len(coord2)
            longest_len = len(coord1)
            shortest = "coord2"
            longest = "coord1"

        i = 0
        while i < shortest_len:
            if shortest == "coord1":
                if sep2d_2to1[i] < tolerance and i == closest_1to2[closest_2to1[i]]:
                    index1_matched.append(i)
                    index2_matched.append(closest_2to1[i])
                else:
                    index1_unmatched.append(i)

            elif shortest == "coord2":
                if sep2d_1to2[i] < tolerance and i == closest_2to1[closest_1to2[i]]:
                    index2_matched.append(i)
                    index1_matched.append(closest_1to2[i])
            i += 1

        for j in range(longest_len):
            if longest == "coord1":
                if j not in index1_matched:
                    index1_unmatched.append(j)

            elif longest == "coord2":
                if j not in index2_matched:
                    index2_unmatched.append(j)

        return (index1_matched, index2_matched, index1_unmatched, index2_unmatched)

        # -------------------------------------- #
        # ----------- Matt's attempt ----------- #

    def load_data(self, first_file, sec_file=None):
        '''
        Load local data into table.

        Parameters:

        filename1 - required - if only argument passed, table will be loaded
        without processing or combining with another table (assumes table has
        already been processed).

        filename2 - optional - if given 2 tables will be processed and joined
        into common table. First file MUST be NED Second MUST be SIMBAD.

        '''
        if (sec_file == None):
            self.combined_table = Table.read(first_file)
            self.combined_table.remove_column("Secure")  # temporary for test data

        else:
            ned_in = Table.read(first_file)
            simbad_in = Table.read(sec_file)

    @staticmethod
    def plot_match_table(combtab, name, col_option=3):
        '''
        The more blue the closer the match.
        '''

        xmask = combtab['Exact Match'] == True
        cmask = combtab['Candidate Match'] == True
        otmask = combtab['ofType Match'] == True
        scatmask = combtab['Shared Category Match'] == True
        gmask = combtab['Generalization Match'] == True
        nomatchmask = combtab['Non Match'] == True

        xmatches = combtab[xmask]
        cmatches = combtab[cmask]
        otmatches = combtab[otmask]
        scatmatches = combtab[scatmask]
        gmatches = combtab[gmask]
        nonmatches = combtab[nomatchmask]

        if col_option == 1:
            cols = ['violet', 'blueviolet', 'blue', 'yellow', 'orange', 'red']
        elif col_option == 2:
            cols = ['blueviolet', 'royalblue', 'cyan', 'mediumseagreen', 'orange', 'red']
        elif col_option == 3:
            cols = ['blueviolet', 'dodgerblue', 'mediumseagreen', 'gold', 'orange', 'red']

        labels = ['Exact Match', 'Candidate Match', 'ofType Match',
                  'Shared Category Match', 'Generalization Match', 'Non Matches']
        matchtypes = [xmatches, cmatches, otmatches, scatmatches, gmatches, nonmatches]

        plt.figure(figsize=(11, 8))

        for i, m in enumerate(matchtypes):
            c = cols[i]
            l = labels[i]
            plt.scatter(m['RA(deg)'], m['DEC(deg)'], color=c, label=l, s=25)

        plt.legend(bbox_to_anchor=(0., 1.0, 1., .102), loc=3, ncol=6, mode='expand')
        plt.xlabel("RA (degrees)")
        plt.ylabel("DEC (degrees)")
        plt.title("{} Object Overlap coloured by Match Type".format(
            name.upper()), loc='center', pad=30.0)
        plt.tight_layout()  # make room for plot labels
        # plt.show()
        return plt

    def saveTable(self, *, fileName, file_format):
        if file_format == "csv":
            self.combined_table.write(fileName+".csv",
                                      format='ascii.csv',
                                      fast_writer=False)
        elif file_format == "fits":
            logging.warning("Sorry, saving as fits is not currently implemented.")
