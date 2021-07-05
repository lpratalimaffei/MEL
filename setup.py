#from distutils.core import setup
from setuptools import setup

setup(name='MEL',
  version= '1.0.0',
  description='Master Equation based Lumping code for integrating single PESs into global kinetic schemes',
  author='Luna Pratali Maffei',
  author_email= 'luna.pratali@polimi.it',
  packages =['MEL'],
  entry_points={'console_scripts':['MEL = MEL.run:main']
  },
  keywords = ['Master Equation', 'Lumping', 'PES', 'Rate Constants'])

'''
$ python setup.py install
'''