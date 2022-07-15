import matplotlib.pyplot as plt
import numpy as np
import skyfield
from skyfield.api import load_file
satellite=load_file('-/home/bella/project/3le.txt')



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
