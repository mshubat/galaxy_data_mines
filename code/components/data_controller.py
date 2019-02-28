# Data Storage and organization class
#
# This class serves as an import, organization, and processing tool for NED
# and SIMBAD data.
#
# Possible data sources:
# - local - on machine datasets
# - remote - via astroquery quieries to NED and SIMBAD
#

from astropy.table import Table, Column
# from astropy.table import Table, Column, hstack, vstack


class DataController:

    ned_to_simbad_dict = {
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
        "QGroup": "",  # # not too common
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
        "VisS": "",  # no match found
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
        "long_name": "tag",
        "Radio": "Rad",
        "Radio_m": "mR",
        "Radio_cm": "cm",
        "Radio_mm": "mm",
        "Radio_sub-mm": "smm",
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

    def __init__(self):
        self.combined_table = None

    @staticmethod
    def ned_to_simbad(ned_entry):
        return DataController.ned_to_simbad_dict[ned_entry]

    @staticmethod
    def simbad_long_to_small(simbad_std):
        return DataController.simbad_std_to_cond[simbad_std]

    @staticmethod
    def candidate_match(non_candidate):
        return DataController.candidate_dict[non_candidate]

    def download_data(self):
        '''
        Fetch remote data from NED and SIMBAD matching coordinates and build table.
        '''
        pass

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

        # rest of published_cats relevant functionality

        pass

    def build_joint_table(self):
        '''

        '''
        pass

    def match_by_objects_location(self):
        pass

    def saveTable(self, *, fileName, file_format):
        if file_format == "csv":
            self.combined_table.write(fileName,
                                      format='ascii.csv',
                                      fast_writer=False)
        elif file_format == "fits":
            self.combined_table.write(fileName)
