#===============================================================================
# File output settings
#===============================================================================
import numpy as np

__all__ = []

# * Output file firectory
path = None

# * Author information
author = 'Thomas J. Rennie'
description = 'MAPEXT test run'

# * Main output file naming scheme
#    'auto'     for automatic (default)
#    '<name>'   for custom name given as <name>
default_naming_scheme = 'sequential'

# * Define naming schemes
naming_schemes = {
    'sequential' : {
        'main' : ['mapext-',['DateTime','%Y%m%dT%H:%M:%S']],
        'objs' : ['mapext-',['']],
    }
}

# * Class for output file
# Leave as None (updated by program during operation)
outfile = None

# * Active boolean flag
def active():
    print((path is not None) and isinstance(path, str),
        outfile is not None)
    isactive = np.all([
        (path is not None) and isinstance(path, str), # Check directory is given
        outfile is not None, # Check outfile has been initiated
    ])
    return isactive