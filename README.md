# Python Distlink - MOID computation
This is a Python wrapper for the C++ library `distlink` published by R.V. Baluev and D.V. Mikryukov used to 
compute the **minimum orbit intersection distance (MOID)** of two orbits described by classical Kepler elements.
Their original C++ library is unmodified and lies in the `distlink-2.0` subdirectory. 
Please refer to [the authors' README file](distlink-2.0/README.txt) for more information. 
The original code is published [here](https://sourceforge.net/projects/distlink/).

## Dependencies
The C++ library needs to be compiled for your OS and your architecture. Therfore, you need to be
able to compile Python extensions. That generally requires
* g++ compiler (**not** clang if you are using macOS)
* installed Python C headers and development libraries (e.g. Ubuntu `sudo apt-get install python3-dev`)

## Setup
* Clone the repository into the root of your project.
* Open a terminal inside `distlink`
* Run `make build` to compile the C++ part of this package
* Run `python test.py` and ensure that no error messages are raised.
* Run `make clean` to remove temporary build files

## Usage
```python
from distlink import COrbitData, MOID_fast

# Create two orbits
o1 = COrbitData(7149.23810, 0.01951, 98.34212, 163.46357, 40.32555)
o2 = COrbitData(7254.82582, 0.02770, 98.10792, 523.30395-360, 10.76767)

# Compute MOID between two orbits.
dist = MOID_fast(o1, o2, 2e-15, 1e-15)
```
`dist` is a `SMOIDResult` object with the following main attributes:
* `distance` and `distance_error`: Computed MOID between orbits and respective error
* `u1` and `u1_error`: Eccentric anomaly of orbit 1 with respective error. Refer to Kholshevnikov and Vassiliev (1999)
  for the exact equation.
* `u2` and `u2_error`: Eccentric anomaly of orbit 1 with respective error. Refer to Kholshevnikov and Vassiliev (1999)
  for the exact equation.

For more details refer to [distlink.h](distlink-2.0/distlink.h) or [the authors' README](distlink-2.0/README.txt) file.
See discussion in Baluev and Mikryukov (2019) on how to set the error limits of the `MOID_fast()` function.

## Learn more
* Python wrapping of C++ code with SWIG: http://www.swig.org/Doc4.0/Python.html
* C++ template functions with SWIG: http://www.swig.org/Doc4.0/SWIGDocumentation.html#SWIGPlus_nn30

## References
* Baluev, Roman V., and Denis V. Mikryukov. "Fast error-controlling MOID computation for confocal elliptic orbits." 
  Astronomy and Computing 27 (2019): 11-22.
* Kholshevnikov, Konstantin V., and Nikolay N. Vassiliev. "On the distance function between two Keplerian elliptic 
  orbits." Celestial Mechanics and Dynamical Astronomy 75.2 (1999): 75-83.