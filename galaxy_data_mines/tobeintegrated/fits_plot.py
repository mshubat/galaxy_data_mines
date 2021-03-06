from astropy.io import fits
from astropy.table import Table
import numpy as np
import matplotlib.pyplot as plt

hdul = fits.open("data/output/alteredTable.fits")

data = hdul[1].data
# print(data[122])

for d in data:
    c = "black"
    if d["Exact Match"] == True:
        c = "blue"
    if d["Candidate Match"] == True:
        c = "green"
    plt.scatter(d["RA(deg)"], d["DEC(deg)"], color=c)

# plt.plot(object_data["RA(deg)"], object_data["DEC(deg)"], 'o', color="red")
plt.show()
