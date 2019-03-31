## Galaxy Data Mines

A command line tool to compare object classifications given by [NED](https://ned.ipac.caltech.edu) and [SIMBAD](http://simbad.u-strasbg.fr/simbad/) in common regions. Objects are found in user provided searches: "cone search" or via region names ("M83").

## Motivation
When astronomers work between NED *(NASA extragalactic Database)* and SIMBAD *(Set of Identifications, Measurements, and Bibliography for Astronomical Data)* it becomes noticeable that both systems use different classification schemes for objects. Although many of these objects are identical, they may be given different names or symbols between the two systems.

As can be imagined, this situation produces some confusion, and causes a researcher (or curious mind) to go through more work in order to sort out whether or not the two systems do or do not agree about a particular object, or set of objects, in a given region of the sky.

The goal of this project is 2 fold:
1. Create a command line tool to automate some of this inter-database object comparison process.
2. Attempt to to quantify and infer overall similarity between the two systems: done through selective queries to the [Messier object set](https://en.wikipedia.org/wiki/Messier_object).

## Code style
[![PEP8](https://img.shields.io/badge/code%20style-pep8-orange.svg)](https://www.python.org/dev/peps/pep-0008/)  

Though there are portions of the project which do not yet conform perfectly to this standard. It has and is being used as the guiding style for the project. Progress will continue to be made to match the the PEP8 standard.

## Screenshots
A few screenshots of the gdmines tool.

**Help Screen** - Provides overview of commands and options  
![help-screen](/docs/screenshots/help-screen.png)

**Match Stats** - Provides a summary of the match likeness of given query  
*example: m83 galaxy*  
![m83-stats](/docs/screenshots/m83-stats.png)

**Plot Results 1** - Graphical display of the objects in the sky, coloured by match status  
*example: m83 galaxy, radius = 1 arcmin*  
![m83-plot-smallr](/docs/screenshots/m83-plot-smallr.png)

**Plot Results 2** - Graphical display of the objects in the sky, coloured by match status  
*example: m83 galaxy, radius = 9 arcmin*  
![m83-plot-bigr](/docs/screenshots/m83-plot-bigr.png)

## Tech/framework used
<b>Built with</b>
- [Python 3.7+](https://www.python.org/)
- [Click](https://palletsprojects.com/p/click/)
- [Astropy](http://www.astropy.org/index.html)
- [Astroquery](https://astroquery.readthedocs.io/en/latest/#)
- [Treelib](https://treelib.readthedocs.io/en/latest/)
- [SciPy](https://www.scipy.org)
   - [Matplotlib](https://matplotlib.org/)
   - [pandas](https://pandas.pydata.org)
   - [NumPy](https://www.numpy.org/)

## Features
* Query a region by name, M31 for example and get comparison results between NED and SIMBAD.
* Adjust parameters such as match-tolerance and obj-radius to fine tune your query.
* Do a cone search of a particular location for those who know exactly what they're looking for.
* View match statistics to get a better idea of the match breakdown.
* Show the match table for detailed match results for each object pair and their relationship.
* Display a 2D plot of overlapping objects in the sky, coloured by their computed match type.
* View the SIMBAD tree structure used to compare object classifications.
* Get a glossary of terms to define what each match type really means.

## Installation
Will be coming to [PyPI](https://pypi.org) soon...

```
pip install galaxy_data_mines
```

## Tests
There is a small test file tree_tester.py which tests the comparison_tree.py operations used to determine the relationships between NED and SIMBAD objects. This can be run by in the following way:
```
python3 tree_tester.py
```

## How to use?
Brief instructions coming soon...

For now try this:
```bash
gdmines --help
```

## Contribute

If you feel intrigued by the project and would like to improve upon it, please feel free to fork it and submit a pull request. Thank you.

## Credits
Thank you to my thesis advisor [Pauline Barmby](https://github.com/PBarmby), who developed the original idea for this project, contributed some foundational elements to the codebase, and provided guidance. The project would not exist without her contributions.

## License
[MIT](LICENSE) © [mshubat](https://github.com/mshubat)
