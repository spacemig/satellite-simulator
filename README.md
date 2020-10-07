# Satellite Simulator in Python (Position and Attitude) 
A simple satellite attitude simulator in python

Author: Miguel Nunes

For Documenting the code use:

http://www.python.org/dev/peps/pep-0008

example comment for functions
```
def square_and_rooter(x):
    """Returns the square root of self times self."""
    ...
```
## General Notes

This code using the following environment 
Mac OS 10.8.3 (Mountain Lion)
Enthought Canopy

dependencies
- numpy **
- matplotlib **
- mayavi (install via Package Manager)
- vtk (install via Package Manager)
- geomag (packaged in thirparty folder)
- sgp4 (packaged in thirparty folder)

** automatically installed on Enthought or Anaconda

## Instructions to install geomag

1. download from url
https://pypi.python.org/pypi/geomag/

2. unzip, open terminal and 'cd' to the unzipped folder

3. install: 
$python setup.py install 

4. copy 'WMM.COF' file
you must copy the file 'WMM.COF' inside the site-packages folder just created by the install. The installation step copies the file to another place that the Enthought canopy cannot recognize.

$ cp WMM.COF /Users/<your_user>/Library/Enthought/Canopy_64bit/User/lib/python2.7/site-packages/geomag/ 


## To install sgp4

1. download from url
https://pypi.python.org/pypi/sgp4/

2. unzip, open terminal and 'cd' to the unzipped folder

3. install: 
$python setup.py install 

Future development
- integrate pyEphem for the SGP4 integration
http://stackoverflow.com/questions/12093766/rectangular-coordinates-for-earth-satellites-in-pyephem
