#!/usr/bin/env python3

from setuptools import setup, find_packages
import os
import codecs
import re

# Copied from wheel package
here = os.path.abspath(os.path.dirname(__file__))
#README = codecs.open(os.path.join(here, 'README.txt'), encoding='utf8').read()
#CHANGES = codecs.open(os.path.join(here, 'CHANGES.txt'), encoding='utf8').read()

with codecs.open(os.path.join(os.path.dirname(__file__), 'genice2', '__init__.py'),
                 encoding='utf8') as version_file:
    metadata = dict(
        re.findall(
            r"""__([a-z]+)__ = "([^"]+)""",
            version_file.read()))

long_desc = "".join(open("README.md", encoding="utf-8").readlines())

with open("requirements.txt") as f:
    requires = [x.strip() for x in f.readlines()]

setup(name='GenIce2',
      python_requires='>=3.6',
      version=metadata['version'],
      description='A Swiss army knife to generate hydrogen-disordered ice structures.',
      long_description=long_desc,
      long_description_content_type="text/markdown",
      classifiers=[
          "Development Status :: 4 - Beta",
          "Intended Audience :: Developers",
          "Programming Language :: Python",
          "Programming Language :: Python :: 3.6",
      ],
      author='Masakazu Matsumoto',
      author_email='vitroid@gmail.com',
      url='https://github.com/vitroid/GenIce/',
      keywords=['genice2', ],
      license='MIT',
      packages=find_packages(),
      install_requires=requires,
      entry_points={
              'console_scripts': [
                  'genice2 = genice2.cli.genice:main',
                  'analice2 = genice2.cli.analice:main'
              ]
      }
      )
