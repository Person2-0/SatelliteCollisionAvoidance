import matplotlib.pyplot as plt
import numpy as np
import skyfield
from skyfield.api import load, wgs84
from skyfield.api import EarthSatellite
from scipy.spatial import KDTree
from numpy.linalg import norm


def satellitepos_vels(satellites, time):
    pos_vels = []
    for satellite in satellites:
        geocentric = satellite.at(time)
        pos_vels.append((geocentric.position.km, geocentric.velocity.km_per_s))
    pos_vels = np.array(pos_vels)
    return pos_vels


def plotearth(ax):
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
    for x, y, z in lons:
        ax.plot(x, y, z, '-k')
    for x, y, z in lats:
        ax.plot(x, y, z, '-k')

    plt.show()


def findpairs(satellites, pos_vels, maxdistance=100):
    tree = KDTree(pos_vels[:, 0, :])
    pairs = tree.query_pairs(r=maxdistance)
    pairs = np.array(list(pairs))
    # make_array[pairnumber, sat0_or_sat1, pos_or_vel, x_y_z]
    pairpos_vels = pos_vels[pairs]
    appdist = approach_distance(pairpos_vels)
    pass


def approach_distance(pairpos_vels):
    # make_delta_array[pairnumber, pos_vels, x_y_z]
    deltas = (pairpos_vels[:, 1, :, :]
             -pairpos_vels[:, 0, :, :])
    delta_norms = np.sqrt(np.sum(deltas**2, axis=-1))
    # plt.hist(delta_norms[:, 0])
    # plt.show()
    distances = []
    for row in pairpos_vels:
        distance = mindist(row[0,0], row[1, 0], row[0,1], row[1,1])
        distances.append(distance)
    distances = np.array(distances)
    pass


def plotsats(ax, pos_vels):
    ax.plot(pos_vels[:, 0, 0], pos_vels[:, 0, 1], pos_vels[:, 0, 2], ",k")
    ax.set_box_aspect((1,1,1))
    r = 1e4
    ax.set(xlim = [-r, r], ylim = [-r, r], zlim = [-r, r])
    plt.show()
    pass


def mindist(pos0, pos1, vel0, vel1):
    r = pos1 - pos0
    v = vel1 - vel0
    costheta = costheta = r.dot(v)/(norm(r)*norm(v))
    sintheta = np.sqrt(1 - costheta ** 2)
    distance = norm(r) * sintheta
    return distance


satellite_url = 'file:///home/bella/project/3le.txt'
satellites = load.tle_file(satellite_url)

ts = load.timescale()

time = ts.now()
pos_vels = satellitepos_vels(satellites, time)
print(pos_vels.shape)
pairs = findpairs(satellites, pos_vels)
if True:
    fig = plt.figure(figsize=[10, 8])  # [12, 10]

    ax  = fig.add_subplot(1, 1, 1, projection='3d')
    #plotearth(ax)
    plotsats(ax, pos_vels)

    plt.show()
