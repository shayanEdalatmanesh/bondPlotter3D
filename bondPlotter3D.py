import numpy as np
import matplotlib.pyplot as plt
from matplotlib.cm import ScalarMappable
from matplotlib.colors import ListedColormap
from mpl_toolkits.mplot3d import Axes3D

# Print the bond lengths (True / False)
printBonds = False

# Read the XYZ file
filename = 'molecule.xyz'
with open(filename, 'r') as f:
    lines = f.readlines()

# Extract the atomic symbols and coordinates
n_atoms = int(lines[0])
symbols = []
coords = np.zeros((n_atoms, 3))
for i in range(2, n_atoms+2):
    symbol, x, y, z = lines[i].split()
    symbols.append(symbol)
    coords[i-2] = [float(x), float(y), float(z)]

# Define a dictionary that maps atomic symbols to colors
colors = {'H': 'grey', 'C': 'black', 'O': 'red', 'N': 'blue', 'F': 'green', 'Au':'gold'}

# Compute the pairwise distances between atoms
dists = np.zeros((n_atoms, n_atoms))
for i in range(n_atoms):
    for j in range(i+1, n_atoms):
        dists[i,j] = dists[j,i] = np.linalg.norm(coords[i]-coords[j])

# Define a range for the plotted bonds
min_dist = 1.3
max_dist = 1.55

# Define a color map that maps bond distances to colors within the range
#cmap = ListedColormap(['blue', 'green', 'yellow', 'orange', 'red'])
cmap = 'jet'
norm = plt.Normalize(min_dist, max_dist)
sm = ScalarMappable(cmap=cmap, norm=norm)
sm.set_array([])

# Define a range for the plotted bond thicknesses
min_linewidth = 1.0
max_linewidth = 10.0

# Define font size for labels
fontsize = 12

# Plot the atoms and bonds in 3D
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
for i in range(n_atoms):
    ax.scatter(coords[i, 0], coords[i, 1], coords[i, 2], marker='o', color=colors[symbols[i]])
for i in range(n_atoms):
    for j in range(i+1, n_atoms):
        if dists[i,j] < min_dist or dists[i,j] > max_dist:
            continue
        color = sm.to_rgba(dists[i,j])
        linewidth = (dists[i,j] - min_dist) / (max_dist - min_dist) * (max_linewidth - min_linewidth) + min_linewidth
        ax.plot([coords[i,0], coords[j,0]],
                [coords[i,1], coords[j,1]],
                [coords[i,2], coords[j,2]], color=color, linewidth=linewidth)
        if printBonds:
            x, y, z = (coords[i]+coords[j])/2
            ax.text(x, y, z, f'{dists[i,j]:.2f}', ha='center', va='center', fontsize=fontsize)
ax.set_xlabel('X', fontsize=fontsize)
ax.set_ylabel('Y', fontsize=fontsize)
ax.set_zlabel('Z', fontsize=fontsize)
cbar = fig.colorbar(sm)
cbar.set_label('Bond distance (Å)')
plt.show()
