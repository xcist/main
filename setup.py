# Copyright 2020, General Electric Company. All rights reserved. See https://github.com/xcist/code/blob/master/LICENSE

# To install XCIST-CatSim, open Python console, run: pip install [folder name]
# e.g., you can navigate to this folder, run: pip install .

from setuptools import setup

setup(name='catsim-xcist',
      version='1.0.0',
      description='Simulation toolkit for X-ray based cancer imaging',
      url='https://github.com/xcist',
      author='Mingye Wu, Paul FitzGerald, Brion Sarachan, Bruno De Man',
      author_email='Mingye.Wu@ge.com',
      license='BSD 3-Clause License',
      install_requires=['numpy', 'scipy', 'matplotlib'],
      packages=['gecatsim', 'gecatsim.pyfiles', 'gecatsim.reconstruction'],
      zip_safe=False,
      package_data={
          'gecatsim':[r'lib/*.*', r'cfg/*.cfg', 
                    r'pyfiles/doserecon/*.*',
                    r'dose_data/*.*',
                    r'bowtie/*.txt', r'material/*', r'material/edlp/*/*.dat',
                    r'phantom/*.*', r'phantom/CatSimLogo_1024/*.*', r'phantom/poly_bin/poly*',
                    r'scatter/*.dat', r'spectrum/*.dat', r'reconstruction/pyfiles/*.*',
                    r'reconstruction/lib/*.*']},
      include_package_data=True)
