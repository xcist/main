Paul FitzGerald, GE-GRC
May 4, 2022

*** IMPORTANT! *** These examples should be run by installing all the required files in your own folders, outside the XCIST repository.

You will also need to download the phantom files required for each experiment from the relevent repository. See the documentation at the top of each experiment_*.py file for the name of the phantom file that you will need.

Use the environment variable XCIST_UserPath to point to the top-level folder where the experiment and phantom files reside.

Be sure that Python can find all the python files in your copy of the examples. I do this using the environment variable PYTHONPATH and by running each experiment_*.py file while the current directory is the folder where that file resides.

EXAMPLE:
I created a folder C:\my_XCIST_files and copied the entire tree "examples" under that.
I then created C:\my_XCIST_files\my_phantoms and put the downloaded phantoms there.
I then assigned the environment variable XCIST_UserPath to C:\my_XCIST_files.
I also added C:\my_XCIST_files\examples\evaluation\pyfiles to the environment variable PYTHONPATH.