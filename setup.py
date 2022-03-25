from distutils.core import setup

# Capture the current git commit to use as version
exec(open("mapext/version.py").read())

config = {
          'name':'mapext',
          'version':__version__,
          'packages':['mapext','mapext.analysis','mapext.core','mapext.io','mapext.catalogues'],
          'package_data':{'':["*.fit","*.fits"]},
          #'include_package_data':True,
          }

setup(**config)
print('Successfully installed version:'+__version__)
