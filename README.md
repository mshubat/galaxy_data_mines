# Galaxy Data Mines
This project is focused on comparing object classifications between two popular astronomical datasets: [NED](https://ned.ipac.caltech.edu) and [SIMBAD](http://simbad.u-strasbg.fr/simbad/).

<img src="https://i.pinimg.com/originals/56/85/cd/5685cdbaacb0b043347c34113b206a92.jpg" alt="space"/>

This tool compares classifications given to common astronomical objects in NED and SIMBAD. Objects are found in user provided searches: from regions via coordinates (RA/DEC/rad) or via object names ("M83").

Once the regions are provided the software will query the particular areas for objects and compute matches. Once matches have been computed the class of each object is compared. Since SIMBAD uses a [hierarchical classification scheme](http://simbad.u-strasbg.fr/simbad/sim-display?data=otypes) and [NED a linear one](https://ned.ipac.caltech.edu/?q=help/srcnom/list-objecttypes&popup=1), the comparisons are done via a mapping from NED's scheme to SIMBAD's. (This produces few issues as the SIMBAD system is (mostly) a superset of NED.

Once mapping is completed the analysis can be done to gather the results in a few different ways. The user can then peruse the results to get a better sense of the object classification consistency between the two systems.

The workflow looks like:
1. Data
   * provided locally *or* via online search
2. Compute 
   * match objects in each dataset by coordinates *or* region
   * compare object classes between NED & SIMBAD
3. Display
   * output is provided in several ways:
      * as tabular data 
      * as plots
      * as a "match summary" giving the basic statistics of the matches
   
More project details will be added here soon. For now, you can read the initial project proposal below to get a sense of the impetus behind this tool.

## Project Plan 

[Original Project Proposal](./mshubat_cs4490_thesis_proposal.pdf)
