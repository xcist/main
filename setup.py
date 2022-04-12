# Copyright 2020, General Electric Company. All rights reserved. See https://github.com/xcist/code/blob/master/LICENSE

# To install XCIST-CatSim, open Python console, run: pip install [folder name]
# e.g., you can navigate to this folder, run: pip install .

from setuptools import setup

setup(name='catsim',
      version='0.1.3',
      description='Simulation toolkit for X-ray based cancer imaging',
      url='https://github.com/xcist',
      author='Mingye Wu, Paul FitzGerald, Brion Sarachan, Bruno De Man',
      author_email='Mingye.Wu@ge.com',
      license='BSD 3-Clause License',
      install_requires=['numpy', 'scipy', 'matplotlib'],
      packages=['catsim'],
      zip_safe=False,
      package_data={'catsim':[r'lib/*.*', r'cfg/*.cfg', 
        r'data/bowtie/*.dat', r'data/material/*', r'data/material/edlp/*/*.dat', 
        r'data/phantom/*.*', r'data/scatter/*.dat', r'data/spectrum/*.dat']},
      include_package_data=True)
