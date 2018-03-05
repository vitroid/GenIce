#!/usr/bin/env python3

from setuptools import setup
import os
import codecs
import re

#Copied from wheel package
here = os.path.abspath(os.path.dirname(__file__))
#README = codecs.open(os.path.join(here, 'README.txt'), encoding='utf8').read()
#CHANGES = codecs.open(os.path.join(here, 'CHANGES.txt'), encoding='utf8').read()

with codecs.open(os.path.join(os.path.dirname(__file__), 'genice', '__init__.py'),
                 encoding='utf8') as version_file:
    metadata = dict(re.findall(r"""__([a-z]+)__ = "([^"]+)""", version_file.read()))

setup(name='GenIce',
      version=metadata['version'],
      description='A Swiss army knife to generate hydrogen-disordered ice structures.',
      #long_description=README + '\n\n' +  CHANGES,
      classifiers=[
        "Development Status :: 4 - Beta",
        "Intended Audience :: Developers",
        "Programming Language :: Python",
        "Programming Language :: Python :: 3.5",
        ],
      author='Masakazu Matsumoto',
      author_email='vitroid@gmail.com',
      url='https://github.com/vitroid/GenIce/',
      keywords=['genice',],
      license='MIT',
      packages=['genice',
                'genice.molecules',
                'genice.lattices',
                'genice.formats',
                ],
      # install_requires=['numpy', 'networkx', 'countrings>=0.1.4', 'svgwrite', 'vpython'],
      install_requires=['numpy', 'networkx', 'countrings>=0.1.4', 'svgwrite', 'cif2ice>=0.2'],
      entry_points = {
              'console_scripts': [
                  'genice = genice.__main__:main'
              ]
          }
#      entry_points="""
#      [console_scripts]
#      genice = genice.__main__:main
#      """
      )

