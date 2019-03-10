'''
from astropy.table import Table
from astroquery.ned import Ned


result_table = Ned.query_region("m83")
result_table.pprint(show_unit=True)

cols = result_table.columns
print(type(cols))

t = Table(result_table)
for col in cols:
    if not col == "RA(deg)":
        del(t[col])
result_table.show_in_browser(jsviewer=True)


for col in cols:
    if not col == "RA(deg)":
        del(result_table[col])


result_table.show_in_browser(jsviewer=True)
'''

from astroquery.simbad import Simbad

# Simbad.list_votable_fields()
# Simbad.get_field_description("otype")
customSimbad = Simbad()
# customSimbad.get_votable_fields()
customSimbad.add_votable_fields("otype")
result_table = customSimbad.query_region("m3")

result_table.pprint(show_unit=True)


result_table.show_in_browser(jsviewer=True)

'''
from astroquery.simbad import Simbad
from astropy import coordinates
import astropy.units as u

# works only for ICRS coordinates:
c = coordinates.SkyCoord("05h35m17.3s -05d23m28s", frame='icrs')
r = 5 * u.arcminute
result_table = Simbad.query_region(c, radius=r)
result_table.pprint(show_unit=True, max_width=80, max_lines=5)

result_table.show_in_browser()
'''
