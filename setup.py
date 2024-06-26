# Copyright 2020, General Electric Company. All rights reserved. See https://github.com/xcist/main/blob/master/license/LICENSE

# To install XCIST-CatSim, open Python console, run: pip install [folder name]
# e.g., you can navigate to this folder, run: pip install .

from setuptools import setup

setup(name='gecatsim',
      version='1.3.2',
      description='Simulation toolkit for X-ray based cancer imaging',
      long_description=open('README.md').read(),
      long_description_content_type='text/markdown',
      url='https://github.com/xcist/documentation/wiki',
      author='Mingye Wu, Paul FitzGerald, James Jobin, Jiayong Zhang, Nandan Reddy, Anupama Debnath, Brion Sarachan, Bruno De Man',
      author_email='Mingye.Wu@gehealthcare.com',
      license='BSD 3-Clause License',
      install_requires=['numpy', 'scipy', 'matplotlib', 'tqdm'],
      packages=['gecatsim', 'gecatsim.pyfiles', 'gecatsim.dose', 'gecatsim.reconstruction'],
      zip_safe=False,
      package_data={
          'gecatsim':[r'bowtie/*.*', 
                    r'cfg/*.*', 
                    r'dose/*.*', r'dose/data/*.*', r'dose/examples/*.*', r'dose/lib/*.*', r'dose/pyfiles/*.*', r'dose/src/*.*', 
                    r'examples/*.*', r'examples/cfg/*.*', r'examples/vct_examples/*.*',
                    r'focal_spot/*.*',
                    r'lib/*.*',
                    r'material/*', r'material/edlp/*/*.dat',
                    r'phantom/*.*', r'phantom/CatSimLogo_1024/*.*', r'phantom/poly_bin/poly*',
                    r'pyfiles/*.*',
                    r'reconstruction/*.*', r'reconstruction/lib/*.*', r'reconstruction/pyfiles/*.*', r'reconstruction/src/*.*',
                    r'response_matrix/*.*',
                    r'scatter/*.*',
                    r'spectrum/*.*']},
      include_package_data=True)
