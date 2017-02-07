#from distutils.core import setup
from setuptools import setup, find_packages

setup(
      name='rational_design',
      version='1.0',
      description='Rational peptide design for DMR lab',
      author='Anna Cunningham',
      author_email='annadcunningham@gmail.com',
      packages=find_packages(),
      scripts=['bin/design']
     )
