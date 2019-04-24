# Statistics Class
#
# Class served to generate, save, and print statistics
# of object comparison # match findings.
#
#

import os
from string import Template
import pandas as pd
import logging

this_dir = os.path.dirname(__file__)


class MatchStats:
    def __init__(self):
        # Open statistics template file.
        stats_template = open(this_dir + "/data/statsSchema.txt", 'r')

        # Set attributes
        self.template = Template(stats_template.read())
        self.processed_template = None

        # ---- Stats variables ----- #
        self.query_name = None        # Set in data_controller
        self.tot_obj_count = None     # Derived
        self.ned_count = None         # Set in data_controller
        self.sim_count = None         # Set in data_controller
        self.overlap_count = None   # Set in data_controller
        self.ned_match_perc = None    # Derived
        self.sim_match_perc = None    # Derived
        # Object Comparison Stats
        self.comp_count = None        # Derived
        # High Level
        self.strong_count = None      # Derived
        self.weak_count = None        # Derived
        self.combined_count = None    # Derived
        self.non_count = None         # Fetched
        self.strong_perc = None       # Derived
        self.weak_perc = None         # Derived
        self.combined_perc = None     # Derived
        self.non_perc = None          # Derived
        # Detailed
        self.exact_count = None       # Fetched
        self.cand_count = None        # Fetched
        self.oftype_count = None      # Fetched
        self.samecat_count = None     # Fetched
        self.general_count = None     # Fetched
        self.exact_perc = None        # Derived
        self.cand_perc = None         # Derived
        self.oftype_perc = None       # Derived
        self.samecat_perc = None      # Derived
        self.general_perc = None      # Derived

    def template_test():
        s = "The cat stopped on the $word1. It then bathed in the $word2."
        s_template = Template(s)
        s = s_template.substitute(word1='grass', word2='sun')
        print(s)

    def generateStatTemplate(self):
        self.processed_template = self.template.substitute(
            query_name=self.query_name,
            # General Statistics
            tot_obj_count=self.tot_obj_count,
            overlap_perc="{:5.2f}%".format((self.overlap_count/self.tot_obj_count)*100),
            ned_count=self.ned_count,
            sim_count=self.sim_count,
            overlap_count=self.overlap_count,
            ned_match_perc="{:5.2f}%".format(self.ned_match_perc),
            sim_match_perc="{:5.2f}%".format(self.sim_match_perc),
            # Object Comparison Stats
            comp_count=self.comp_count,
            # High Level
            strong_count=self.strong_count,
            weak_count=self.weak_count,
            combined_count=self.combined_count,
            non_count=self.non_count,
            strong_perc="{:5.2f}%".format(self.strong_perc),
            weak_perc="{:5.2f}%".format(self.weak_perc),
            combined_perc="{:5.2f}%".format(self.combined_perc),
            non_perc="{:5.2f}%".format(self.non_perc),
            # Detailed
            exact_count=self.exact_count,
            cand_count=self.cand_count,
            oftype_count=self.oftype_count,
            samecat_count=self.samecat_count,
            general_count=self.general_count,
            exact_perc="{:5.2f}%".format(self.exact_perc),
            cand_perc="{:5.2f}%".format(self.cand_perc),
            oftype_perc="{:5.2f}%".format(self.oftype_perc),
            samecat_perc="{:5.2f}%".format(self.samecat_perc),
            general_perc="{:5.2f}%".format(self.general_perc)
        )
        print(self.processed_template)
        return self.processed_template

    def derive_table_stats(self, table):
        '''
        Use the Pandas package to get some stats about the generated
        match table and derive others.
        '''
        # High Level Overview - Derived Values:
        self.tot_obj_count = self.sim_count + self.ned_count
        self.ned_match_perc = (self.overlap_count/self.ned_count)*100
        self.sim_match_perc = (self.overlap_count/self.sim_count)*100
        self.comp_count = self.overlap_count

        # Object Comparison Results - Fetched Values:
        df = table.to_pandas()

        # Count number of each match (sort=False, so index 1 is counts True)
        exact_dict = df["Exact Match"].value_counts(sort=False)
        cand_dict = df["Candidate Match"].value_counts(sort=False)
        oftype_dict = df["ofType Match"].value_counts(sort=False)
        samecat_dict = df["Shared Category Match"].value_counts(sort=False)
        gen_dict = df["Generalization Match"].value_counts().to_dict()
        non_dict = df["Non Match"].value_counts(sort=False)

        self.exact_count = 0 if True not in exact_dict else exact_dict[True]
        self.cand_count = 0 if True not in cand_dict else cand_dict[True]
        self.oftype_count = 0 if True not in oftype_dict else oftype_dict[True]
        self.samecat_count = 0 if True not in samecat_dict else samecat_dict[True]
        self.general_count = 0 if True not in gen_dict else gen_dict[True]
        self.non_count = 0 if True not in non_dict else non_dict[True]

        # Object Comparison Results - Derived Values:
        self.exact_perc = 0 if self.comp_count == 0 else (self.exact_count/self.comp_count)*100
        self.cand_perc = 0 if self.comp_count == 0 else (self.cand_count/self.comp_count)*100
        self.oftype_perc = 0 if self.comp_count == 0 else (self.oftype_count/self.comp_count)*100
        self.samecat_perc = 0 if self.comp_count == 0 else (self.samecat_count/self.comp_count)*100
        self.general_perc = 0 if self.comp_count == 0 else (self.general_count/self.comp_count)*100
        self.non_perc = 0 if self.comp_count == 0 else (self.non_count/self.comp_count)*100

        self.strong_count = self.exact_count + self.cand_count + self.oftype_count
        self.weak_count = self.samecat_count + self.general_count
        self.combined_count = self.strong_count+self.weak_count

        self.strong_perc = self.exact_perc + self.cand_perc + self.oftype_perc
        self.weak_perc = self.samecat_perc + self.general_perc
        self.combined_perc = self.strong_perc + self.weak_perc

    def saveStats(self):
        stats_output = open(this_dir + "data/stats_output.txt", 'r')
        # dc.get_table_stats().to_csv(r"gdm_output/result_stats.csv")

        pass


if __name__ == "__main__":
    cs = MatchStats()
    cs.generateStatTemplate()
    # MatchStats.template_test()
