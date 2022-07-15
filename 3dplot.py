import matplotlib.pyplot as plt
import numpy as np
import skyfield
from skyfield.api import load, wgs84
from skyfield.api import EarthSatellite



satellite_url = 'file:///home/bella/project/3le.txt'
satellites = load.tle_file(satellite_url)
print('Loaded', len(satellites), 'satellites')


ts = load.timescale()
line1 = '1 5U 58002B   22196.28486465  .00000351  00000-0  43722-3 0  9993'
line2 = '2 5  34.2492  53.6414 1847265 162.8155 204.4542 10.84977341287669'
satellite = EarthSatellite(line1, line2, 'Vangaurd 1', ts)
print(satellite)

t = ts.now()

geocentric = satellite.at(t)
print(geocentric.position.km)



plt.rcParams["figure.figsize"] = [7.00, 3.50]
plt.rcParams["figure.autolayout"] = True
fig = plt.figure()
ax = fig.add_subplot(projection='3d')
r = (0.5, -0.5)
u, v = np.mgrid[0:2 * np.pi:30j, 0:np.pi:20j]
x = np.cos(u) * np.sin(v)
y = np.sin(u) * np.sin(v)
z = np.cos(v)
ax.plot_surface(x, y, z, color='b')
plt.show()
