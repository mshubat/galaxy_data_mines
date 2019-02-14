# Galaxy Agent - User Interface
from data_controller import DataController
from comparison_tree import ComparisonTree


combined_table = "M83_NScomb_myTest.fits"
# ned_name = input("Please enter the NED filename: ")

ct = ComparisonTree(run_mode=True)
data = DataController()
data.load_data(combined_table)

ct.compare_objects(data.combined_table)
# load the data from the file and
data.saveTable(fileName="alteredTable.csv")
data.combined_table.show_in_browser(jsviewer=True)
