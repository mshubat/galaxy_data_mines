# Galaxy Agent - User Interface
from .data_controller import DataController
from .comparison_tree import ComparisonTree
import sys
#sys.path.insert(0, '../Components/')
#original = sys.stdout
#sys.stdout = open('output_log.txt', 'w')


def main():
    combined_table = "data/M83_NScomb_mar1.fits"
    # ned_name = input("Please enter the NED filename: ")

    ct = ComparisonTree(run_mode=True)
    data = DataController()
    data.load_data(combined_table)

    ct.compare_objects(data.combined_table)
    # load the data from the file and
    #data.saveTable(fileName="alteredTable.fits", file_format="fits")
    data.combined_table.show_in_browser(jsviewer=True)
