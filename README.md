# SatelliteCollisionAvoidance
This project looks for pairs of satellites that come within 10km of each other. For each pair, it writes the time, names, catalog numbers, positions, velocities, distances, and extrapolated sistances of the satellites to a file. It will start at a designated time and jumps forward every 10 seconds and adds it to the file. 

## Installation
To be able to run this code you will need the following Python libraries:
Matplotlib
```
$ python -m pip install -U pip
$ python -m pip install -U matplotlib
```
Scipy
```
$ python -m pip install --user numpy scipy matplotlib ipython jupyter pandas sympy nose
```
Skyfield
```
$ pip install skyfield
```
Numpy
```
$ pip install numpy
```
Astropy
```
$ pip install astropy
```
You will also need a full catalog TLE of satellites downloaded from [Space Track](https://www.space-track.org/) 

## Usage
To run the program, follow the Installation steps above, then run
```
$ python3 main.py
```
## Materials
If you would like to learn more about the plots made in this project, you can look at either the Satellite Collision Avoidance paper or Slideshow in the Presentation directory. 

## License
This project is liscenses under the GPL-3.0 License- see the LICENSE.md file for details.
