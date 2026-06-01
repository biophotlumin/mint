from pathlib import Path
from setuptools.config.setupcfg import read_configuration

__version__ = read_configuration(
    Path(__file__).parent.parent.joinpath('setup.cfg'))['metadata']['version']
__author__ = 'Baptiste Grimaud'
__credits__ = 'LuMIn, École Normale Supérieure Paris-Saclay'