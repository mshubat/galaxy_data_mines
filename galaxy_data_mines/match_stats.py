# Statistics Class
#
# Class served to generate, save, and print statistics
# of object comparison # match findings.
#
#

import os
from string import Template
import pandas as pd

this_dir = os.path.dirname(__file__)


class MatchStats:
    def __init__(self):
        # Open statistics template file.
        stats_template = open(this_dir + "/data/statsSchema.txt", 'r')

        # Set attributes
        self.template = Template(stats_template.read())
        self.processed_template = None

        # ---- Stats variables ----- #
        self.tot_obj_count = None     # Derived
        self.ned_count = None
        self.sim_count = None
        self.ned_match_count = None
        self.sim_match_count = None
        self.ned_match_perc = None    # Derived
        self.sim_match_perc = None    # Derived
        # Object Comparison Stats
        self.comp_count = None        # Derived
        # High Level
        self.strong_count = None      # Derived
        self.weak_count = None        # Derived
        self.non_count = None         # Derived
        self.strong_perc = None       # Derived
        self.weak_perc = None         # Derived
        self.non_perc = None          # Derived
        # Detailed
        self.exact_count = None       # Derived
        self.cand_count = None        # Derived
        self.oftype_count = None      # Derived
        self.samecat_count = None     # Derived
        self.general_count = None     # Derived
        self.exact_perc = None        # Derived
        self.cand_perc = None         # Derived
        self.oftype_perc = None       # Derived
        self.samecat_perc = None      # Derived
        self.general_perc = None      # Derived
        self.table_describe = None    # Derived

    def template_test():
        s = "The cat stopped on the $word1. It then bathed in the $word2."
        s_template = Template(s)
        s = s_template.substitute(word1='grass', word2='sun')
        print(s)

    def generateStatTemplate(self):
        self.processed_template = self.template.substitute(
            # General Statistics
            tot_obj_count=self.tot_obj_count,
            ned_count=self.ned_count,
            sim_count=self.sim_count,
            ned_match_perc=self.ned_match_perc,
            sim_match_perc=self.sim_match_perc,
            # Object Comparison Stats
            comp_count=self.comp_count,
            # High Level
            strong_count=self.strong_count,
            weak_count=self.weak_count,
            non_count=self.non_count,
            strong_perc=self.strong_perc,
            weak_perc=self.weak_perc,
            non_perc=self.non_perc,
            # Detailed
            exact_count=self.exact_count,
            cand_count=self.cand_count,
            oftype_count=self.oftype_count,
            samecat_count=self.samecat_count,
            general_count=self.general_count,
            exact_perc=self.exact_perc,
            cand_perc=self.cand_perc,
            oftype_perc=self.oftype_perc,
            samecat_perc=self.samecat_perc,
            general_perc=self.general_perc,
            table_describe=self.table_describe
        )
        print(self.processed_template)

    def derive_table_stats(self, table):
        '''
        Use the Pandas package to get some stats about the generated
        match table and derive others.
        '''
        # Calculate known stats.
        self.tot_obj_count = self.sim_count + self.ned_count

        df = table.to_pandas()
        matchcols = df.columns[10:]

        print(df["Exact Match"].describe(include="all"))
        exact_col_counts = df["Exact Match"].value_counts(
            sort=False)  # Ensures False is always index 0

        print(exact_col_counts)
        print("False count: ".format(exact_col_counts[0]))
        print("True count: ".format(exact_col_counts[1]))

        self.table_describe = df.describe(include="all")

    def saveStats(self):
        stats_output = open(this_dir + "data/stats_output.txt", 'r')
        # dc.get_table_stats().to_csv(r"gdm_output/result_stats.csv")

        pass


if __name__ == "__main__":
    cs = MatchStats()
    cs.generateStatTemplate()
    # MatchStats.template_test()
