import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits

from astropy.utils.data import download_file


table_data = image_file = download_file(
    'http://data.astropy.org/tutorials/FITS-images/HorseHead.fits', cache=True)


hdu_list = fits.open(image_file)
hdu_list.info()

image_data = hdu_list[0].data
print(type(image_data))
print(image_data.shape)

# hdu_list.close()

plt.imshow(image_data, cmap="gray")
plt.colorbar()

plt.show()
