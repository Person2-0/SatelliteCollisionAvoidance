import matplotlib.pyplot as plt
import numpy
import numpy as np
import skyfield
from skyfield.api import load, wgs84
from skyfield.api import EarthSatellite
from scipy.spatial import KDTree
from numpy.linalg import norm
import datetime
from astropy.io import fits
from pathlib import Path

seconds_per_day = 24 * 60 * 60
outdir_default = Path("/home/bella/Documents/Satellites")


def satellitepos_vels(satellites, time):
    pos_vels = []
    for satellite in satellites:
        geocentric = satellite.at(time)
        pos_vels.append((geocentric.position.km, geocentric.velocity.km_per_s))
    pos_vels = np.array(pos_vels)
    return pos_vels


def plotearth(ax):
    halfpi, pi, twopi = [f * np.pi for f in [0.5, 1, 2]]
    degs, rads = 180 / pi, pi / 180

    re = 6378.

    theta = np.linspace(0, twopi, 201)
    cth, sth, zth = [f(theta) for f in [np.cos, np.sin, np.zeros_like]]
    lon0 = re * np.vstack((cth, zth, sth))
    lons = []
    for phi in rads * np.arange(0, 180, 15):
        cph, sph = [f(phi) for f in [np.cos, np.sin]]
        lon = np.vstack((lon0[0] * cph - lon0[1] * sph,
                         lon0[1] * cph + lon0[0] * sph,
                         lon0[2]))
        lons.append(lon)

    # lat0 = re * np.vstack((cth, sth, zth))
    lats = []
    for phi in rads * np.arange(-75, 90, 15):
        cph, sph = [f(phi) for f in [np.cos, np.sin]]
        lat = re * np.vstack((cth * cph, sth * cph, zth + sph))
        lats.append(lat)
    for x, y, z in lons:
        ax.plot(x, y, z, '-k')
    for x, y, z in lats:
        ax.plot(x, y, z, '-k')

    plt.show()


def findpairs(pos_vels, search_radius=100, maxdistance=10, timestep=10):
    tree = KDTree(pos_vels[:, 0, :])
    pairs = tree.query_pairs(r=search_radius)
    pairs = np.array(list(pairs))
    # make_array[pairnumber, sat0_or_sat1, pos_or_vel, x_y_z]
    pairpos_vels = pos_vels[pairs]
    # plt.hist(pairpos_vels[:, 1]
    # plt.show
    distances, times = mindist_and_time(pairpos_vels)
    w = np.logical_and(~numpy.isnan(distances),
                       np.logical_and(times < timestep / 2,
                                      np.logical_and(times > -timestep / 2,
                                                     distances < maxdistance)))

    return pairs[w], distances[w], times[w]


# def approach_distance_and_time(pairpos_vels):
# make_delta_array[pairnumber, pos_vels, x_y_z]
# deltas = (pairpos_vels[:, 1, :, :]
#          - pairpos_vels[:, 0, :, :])
# delta_norms = np.sqrt(np.sum(deltas ** 2, axis=-1))
# plt.hist(delta_norms[:, 0])
# plt.show()
# distances = []
# for row in pairpos_vels:
#     distance = mindist(row[0, 0], row[1, 0], row[0, 1], row[1, 1])
#     distances.append(distance)
# slow_distances = np.array(distances)
# distances, times = mindist_and_time(pairpos_vels)
# return distances, times


def mindist_and_time(pairpos_vels):
    deltas = pairpos_vels[:, 1] - pairpos_vels[:, 0]
    r = deltas[:, 0, :]
    v = deltas[:, 1, :]
    rnorm = np.sqrt(np.sum(r ** 2, axis=-1))
    vnorm = np.sqrt(np.sum(v ** 2, axis=-1))
    rdotv = np.sum(r * v, axis=-1)
    costheta = rdotv / (rnorm * vnorm)
    sintheta = np.sqrt(1 - costheta ** 2)
    distance = rnorm * sintheta
    # the distance it has to travel to the minimum distance point
    travel = -rnorm * costheta
    time = travel / vnorm
    return distance, time


def mindist(pos0, pos1, vel0, vel1):
    r = pos1 - pos0
    v = vel1 - vel0
    costheta = r.dot(v) / (norm(r) * norm(v))
    sintheta = np.sqrt(1 - costheta ** 2)
    distance = norm(r) * sintheta
    return distance


def plotsats(ax, positions, ):
    ax.plot(positions[:, 0], positions[:, 1], positions[:, 2], ",k")
    ax.set_box_aspect((1, 1, 1))
    r = 1e4
    ax.set(xlim=[-r, r], ylim=[-r, r], zlim=[-r, r])
    plt.show()


def getpairs(satellite_url=None, tstart=None, tend=None, outdir=None):
    if satellite_url is None:
        satellite_url = 'file:///home/bella/project/3le.txt'
    if tstart is None:
        ts = load.timescale()
        tstart = ts.now()
    if tend is None:
        # tend = tstart.ts.tt_jd(tstart.tt+1)
        tend = tstart + 1
    if outdir is None:
        outdir = outdir_default
    outdir = Path(outdir)
    outdir.mkdir(exist_ok=True)
    outfile = outdir.joinpath(f"conjunctions {tstart.utc_datetime():%Y%m%d_%H%M}.fits")
    satellites = load.tle_file(satellite_url)
    results = []
    time = tstart
    timestep_s = 10  # seconds
    timestep = timestep_s / seconds_per_day  # days
    while time.tt < tend.tt:
        print(f"{time.utc_datetime()}")
        result = pairs_for_time(satellites, time, timestep=timestep_s)
        results.append(result)
        time = time + timestep
    results = np.concatenate(results)
    if outfile.exists():
        outfile.replace(outfile.with_suffix(".bak"))
    fits.writeto(outfile, results)
    read_back = fits.getdata(outfile)
    return results


def pairs_for_time(satellites, time, search_radius=100, maxdistance=10, timestep=10):
    pos_vels = satellitepos_vels(satellites, time)
    if np.any(np.isnan(pos_vels)):
        print(pos_vels.shape, np.count_nonzero(np.isnan(pos_vels[:, 0, 0])), "invalidpositions")
        pos_vels = np.nan_to_num(pos_vels)
    pairs, distances, delta_times = findpairs(pos_vels, search_radius=search_radius,
                                              maxdistance=maxdistance, timestep=timestep)
    times = time + delta_times / seconds_per_day
    pair_type = np.dtype([("time_jd", float), ("sat_names", "U30", 2), ("catnums", int, 2),
                          ("positions", float, (2, 3)), ("velocities", float, (2, 3)),
                          ("distance", float), ("extrapolated_distance", float)])
    result = np.zeros(len(distances), dtype=pair_type)
    result["extrapolated_distance"] = distances
    for i in range(len(distances)):
        result[i]["time_jd"] = times[i].tt
        for j, satnum in enumerate(pairs[i]):
            sat = satellites[satnum]
            geocentric = sat.at(times[i])
            result[i]["positions"][j] = geocentric.position.km
            result[i]["velocities"][j] = geocentric.velocity.km_per_s
            result[i]["sat_names"][j] = sat.name
            result[i]["catnums"][j] = sat.model.satnum

    delta_pos = np.diff(result["positions"], axis=1)[:, 0, :]
    new_distance = np.sqrt(np.sum(delta_pos ** 2, axis=-1))
    result["distance"] = new_distance
    return result


def plot_results(outdir=None):
    if outdir is None:
        outdir = outdir_default
    outdir = Path(outdir)
    fig = plt.figure(figsize=[10, 8])  # [12, 10]

    ax = fig.add_subplot(1, 1, 1, projection='3d')
    ax.set_xlabel('km')
    ax.set_ylabel('km')
    ax.set_zlabel('km')
    all_data = []
    for infile in outdir.glob("conjunction*.fits"):
        data = fits.getdata(infile)
        all_data.append(data)
    all_data = np.concatenate(all_data)
    # FIXME why is the distance different from extrapolated distance different for some satellites
    # quick patch
    all_data = all_data[all_data["distance"] < 10]
    # plotearth(ax)
    plotsats(ax, all_data["positions"].reshape(-1, 3))
    plt.show()

    fig.savefig(outdir.joinpath("conjunction_locations.png"))
    dright, dup = delta_in_vertical_coordinates(all_data)
    figvert = plt.figure(figsize=[10, 8])
    axvert = figvert.add_subplot(1, 1, 1)
    axvert.plot(dright, dup, ",")
    axvert.set(aspect=1)
    plt.show()
    pass


def cross(a, b):
    x = a[..., 1] * b[..., 2] - a[..., 2] * b[..., 1]
    y = a[..., 2] * b[..., 0] - a[..., 0] * b[..., 2]
    z = a[..., 0] * b[..., 1] - a[..., 1] * b[..., 0]
    result = np.vstack((x, y, z)).T
    return result


def delta_in_vertical_coordinates(data):
    positions = data["positions"]
    deltas = positions[:, 1] - positions[:, 0]
    velocities = data["velocities"]
    delta_velocities = velocities[:, 1] - velocities[:, 0]
    up = positions[:, 1]
    up = up / np.linalg.norm(up, axis=-1).reshape((-1, 1))
    right = cross(delta_velocities, up)
    right = right / np.linalg.norm(right, axis=-1).reshape((-1, 1))
    dright = np.sum(deltas * right, axis=-1)
    dup = np.sum(deltas * up, axis=-1)
    return dright, dup


def main():
    ts = load.timescale()
    outdir = outdir_default
    tfirst = ts.tt(2022, 8, 2, 0, 0, 0)
    tdelta = 0.001
    while True:
        result = getpairs(tstart=tfirst, tend=tfirst + tdelta, outdir=outdir)
        tfirst = tfirst + tdelta
        plot_results(outdir=outdir)
    # tstart = ts.tt(2022, 7, 25, 0, 0, 0)
    # tstart = tstart + 200/seconds_per_day
    # tend = tstart + 0.1
    # result = getpairs(tstart=tstart, tend=tend)


if __name__ == "__main__":
    # plt.ion()
    # main()
    plot_results()
