import numpy as np
from numpy.linalg import norm
import matplotlib.pyplot as plt

r = np.array([3, 4, 0])
v = np.array([-1, -4, 0])
print(r, v, r*v, sum(r*v))
print(r.dot(v))


costheta = r.dot(v)/(norm(r)*norm(v))
sintheta = np.sqrt(1 - costheta ** 2)
theta_deg = np.rad2deg(np.arccos(costheta))
theta_deg2 = np.rad2deg(np.arcsin(sintheta))
print(f"Theta is {theta_deg:0.2f} which should be equivalent to {theta_deg2:0.2f}")

distance = norm(r) * sintheta



def plotvec(v, start=[0,0]):
    plt.plot([start[0],start[0]+v[0]], [start[1],start[1]+v[1]])

plotvec(v)
plotvec(r)
plotvec(v, start=r)
plt.gca().set(aspect=1)


def mindist(pos0, pos1, vel0, vel1):
    r = pos1 - pos0
    v = vel1 - vel0
    costheta = costheta = r.dot(v)/(norm(r)*norm(v))
    sintheta = np.sqrt(1 - costheta ** 2)
    distance = norm(r) * sintheta
    return distance

print(mindist([0,0,0], r, [0,0,0], v))

pos0 = np.array([8372, 859, 4389])
pos1 = pos0 + r
v0 = [3.5, 7.1, 0.1]
v1 = v + v0

print(mindist(pos0, pos1, v0, v1))
