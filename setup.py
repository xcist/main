# Copyright 2020, General Electric Company. All rights reserved. See https://github.com/xcist/main/blob/master/license/LICENSE

# To install XCIST-CatSim, open Python console, run: pip install [folder name]
# e.g., you can navigate to this folder, run: pip install .

from setuptools import setup, Extension
from setuptools.command.build_py import build_py as _build_py
from setuptools.command.build_ext import build_ext
import subprocess
import os
import sys


class BuildAssetsCommand(build_ext):
    """Run a bash script before building Python modules."""
    def run(self):
        # Only run if the script exists (useful when sdist excludes it)
        if sys.platform == 'linux':
            script = os.path.join(os.path.dirname(__file__), "gecatsim", "clib_build", "BuildLinux64")
            print(script)
            if os.path.exists(script):
                print("Building the catsim library with: ", script)
                if sys.platform.startswith("win"):
                    # Prefer cross-platform steps written in Python,
                    # but if you must, call via bash provided by Git-Bash/MSYS2.
                    bash = os.environ.get("BASH", "bash")
                    subprocess.check_call([bash, script])
                else:
                    subprocess.check_call(["bash", script], cwd = os.path.join(os.path.dirname(__file__), "gecatsim", "clib_build"))
            else:
                print("Error! script not found!")
            # Continue with normal build
            super().run()

setup(name='gecatsim',
      version='1.6.8',
      ext_modules = [Extension("dummy", sources=[])],
      cmdclass = {"build_ext": BuildAssetsCommand},
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
                    r'clib_build/*', r'clib_build/*/*', r'clib_build/*/*/*',
                    r'material/*', r'material/edlp/*/*.dat',
                    r'phantom/*.*', r'phantom/CatSimLogo_1024/*.*', r'phantom/poly_bin/poly*',
                    r'pyfiles/FlatPanel/*.*',
                    r'pyfiles/*.*',
                    r'reconstruction/*.*', r'reconstruction/lib/*.*', r'reconstruction/pyfiles/*.*', r'reconstruction/src/*.*',
                    r'response_matrix/*.*',
                    r'scatter/*.*',
                    r'spectrum/*.*']},
      include_package_data=True)
