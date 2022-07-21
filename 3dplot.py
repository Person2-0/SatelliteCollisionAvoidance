import matplotlib.pyplot as plt
import numpy as np
import skyfield
from skyfield.api import load, wgs84
from skyfield.api import EarthSatellite


satellite_url = 'file:///home/bella/project/3le.txt'
satellites = load.tle_file(satellite_url)
for i in range(len(satellites)):
    print(satellites[i])



ts = load.timescale()
L1 = '1 00005U 58002B   22201.80587978  .00000357  00000-0  43440-3 0  9992'
L2 = '2 00005  34.2522  36.6233 1847469 187.6731 169.1333 10.84981346288165'
satellite = EarthSatellite(L1, L2, 'Vangaurd 1', ts)
print(satellite)

t = ts.now()

geocentric = satellite.at(t)
print(geocentric.position.km)


Roadster = EarthSatellite(L1, L2)

Rpos    = Roadster.at(t).position.km
Rposecl = Roadster.at(t).ecliptic_position().km

print(Rpos.shape)

halfpi, pi, twopi = [f*np.pi for f in [0.5, 1, 2]]
degs, rads = 180/pi, pi/180

re = 6378.

theta = np.linspace(0, twopi, 201)
cth, sth, zth = [f(theta) for f in [np.cos, np.sin, np.zeros_like]]
lon0 = re*np.vstack((cth, zth, sth))
lons = []
for phi in rads*np.arange(0, 180, 15):
    cph, sph = [f(phi) for f in [np.cos, np.sin]]
    lon = np.vstack((lon0[0]*cph - lon0[1]*sph,
                     lon0[1]*cph + lon0[0]*sph,
                     lon0[2]) )
    lons.append(lon)

lat0 = re*np.vstack((cth, sth, zth))
lats = []
for phi in rads*np.arange(-75, 90, 15):
    cph, sph = [f(phi) for f in [np.cos, np.sin]]
    lat = re*np.vstack((cth*cph, sth*cph, zth+sph))
    lats.append(lat)

if True:
    fig = plt.figure(figsize=[10, 8])  # [12, 10]

    ax  = fig.add_subplot(1, 1, 1, projection='3d')


    x, y, z = Rpos
    ax.plot(x, y, z)
    for x, y, z in lons:
        ax.plot(x, y, z, '-k')
    for x, y, z in lats:
        ax.plot(x, y, z, '-k')

    plt.show()

r_Roadster = np.sqrt((Rpos**2).sum(axis=0))
alt_roadster = r_Roadster - re
