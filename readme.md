MAPEXT
======
*TJRennie* \
*Init:* 25<sup>th</sup> March 2022
Running default supplied MAPEXT scripts
---------------------------------------
1. Ensure all dependant moduleßs are installed
2. Setup MAPEXT using ```python setup.py install```.
3. Update run_mapext.py to match your run specifications
4. Run ```python run_mapext.py``` to run the program

DevLog
------
**0.0.5 Bug fixes for sample run files**

**0.0.4 Added beam correction algorithms and extended region backend**
 - Beam correction ialgorithms added to produce a model for and them account for a non-gaussian beam with significant sidelobes in the main pipeline
 - Implemented scientific colormaps as default option on some functions to ensure maps and plots are universally readable by all audiences

**0.0.3 Added region support**
- Allows sources (defined as gaussian/point like) and regions (extended emission parameterised by a boundary defined at a frequency)
- Added aperture photometry options for custom defined shapes

**0.0.2 Extended input/output options**
 - Added tab10 and tab20 to use later in plotting
 - Added spectrum function to plot SEDs

**0.0.1 Initial commit**
 - Main refactor of MAPEXT to run in a more object-oriented way

Acknowledgements
----------------
The Scientific colour maps implemented in this repository
[(Crameri 2018)](http://doi.org/10.5281/zenodo.1243862) are used to prevent visual distortion of the data and exclusion of users with colour­vision deficiencies [(Crameri et al., 2020)](https://doi.org/10.1038/s41467-020-19160-7).
