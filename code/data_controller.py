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

    def __init__(self):
        self.combined_table = None

    @staticmethod
    def ned_to_simbad(ned_entry):
        return DataController.ned_to_simbad_dict[ned_entry]

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
