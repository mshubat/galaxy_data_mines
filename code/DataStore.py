class DataStore:

    # NED(left) to SIMBAD(right) mapping.
    # - If no direct map or one-to-many possibility,
    #   then a tuple with all possibilities is used.
    # - Candidates are currently ignored

    ned_to_simbad_dict = {
        "*": "*",
        "**": "**",
        "*Ass": "As*",
        "*Cl": "Cl*",
        "AbLS": "ALS",
        "Blue*": ("BS*", "s*b"),  # ??
        "C*": "C*",
        "EmLS": "Em*",  # ??
        "EmObj": "EmO",
        "exG*": "",  # ??? Extragalactic star: not a member of a known galaxy.
        "Flare*": "Fl*",
        "G": "G",
        "GammaS": "gam",
        "GClstr": "ClG",
        "GGroup": "GrG",
        "GPair": "PaG",
        "GTrpl": "",  # ???
        "G_Lens": "LeG",
        "HII": "HII",
        "IrS": "IR",
        "MCld": "MoC",
        "Neb": "Cld",  # ??
        "Nova": "No*",
        "Other": "",  # ??? Perhaps multiple matches
        "PN": "PN",
        "PofG": "PoG",
        "Psr": "Psr",
        "QGroup": "",  # ???
        "QSO": "QSO",
        "Q_Lens": "LeQ",
        "RadioS": "Rad",
        "Red*": ("RG*", "s*r"),  # ??
        "RfN": "RNe",
        "SN": "SN*",
        "SNR": "SNR",
        "UvES": "UV",  # ???
        "UvS": "UV",
        "V*": "V*",
        "VisS": "",  # ???
        "WD*": "WD*",
        "WR*": "WR*",
        "XrayS": "X",
        # Generally, SIMBAD does not distinguish between "Galactic"
        # variants of objects
        "!*": "*",
        "!**": "**",
        "!*Ass": "As*",
        "!*Cl": "Cl*",
        "!Blue*": ("BS*", "s*b"),  # ??
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
        "!Red*": ("RG*", "s*r"),  # ??
        "!SN": "SN*",
        "!SNR": "SNR",
        "!V*": "V*",
        "!WD*": "WD*",
        "!WR*": "WR*"
    }


print(DataStore.ned_to_simbad_dict["Blue*"])
print(DataStore.ned_to_simbad_dict["*"])

print("*" in DataStore.ned_to_simbad_dict)
print(DataStore.ned_to_simbad_dict)

for key, entry in DataStore.ned_to_simbad_dict.items():
    print("{}, {}".format(key, entry))
