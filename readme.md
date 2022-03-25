MAPEXT
======
*TJRennie, Version 0.0.3*

Running default supplied MAPEXT scripts
---------------------------------------
1. Ensure all dependant modules are installed
1. Setup MAPEXT using ```python setup.py install```. If you are reinstalling MAPEXT use ```python setup.py clean --all install``` to update
1. Update run_mapext.py to match your run specifications
1. Run ```python run_mapext.py``` to run the program

DevLog
------
**0.0.1 Initial commit**
 - Main refactor of MAPEXT to run in a more object-oriented way
**0.0.2 Extended input/output options**
 - Added tab10 and tab20 to use later in plotting
 - Added spectrum function to plot SEDs
**0.0.3 Added region support**
 - Allows sources (defined as gaussian/point like) and regions (extended emission parameterised by a boundary defined at a frequency)
 - Added aperture photometry options for custom defined shapes
