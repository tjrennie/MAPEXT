from setuptools import setup, find_packages
import codecs
import os
import toml

from mapext._version import __version__, __date__

print(f'\nInstalling MAPEXT {__version__}\nLast updated on {__date__}\n')

with open("pyproject.toml", "r") as f:
    pyproject = toml.load(f)

here = os.path.abspath(os.path.dirname(__file__))

with codecs.open(os.path.join(here, pyproject["project"]["readme"]), encoding="utf-8") as fh:
    long_description = "\n" + fh.read()

# Setting up
setup(
    **pyproject,
    **{'packages':find_packages()},
)