# identify_cells

Algorithm that identifies clusters based on a set of criteria and calculates their height and width.

## Criterias
*  Minimum value of a grid cell in order to be part of a cluster
*  Minimum maximum value a cluster must have in order to be considered a cluster
*  Minimum height of grid cell with maximum value of the cluster for it to be considered
*  Maximum height of grid cell with maximum value of the cluster for it to be considered
*  Minimum height of the top of the cluster

**Remark:** These criterias are quite arbitrary and can be andjusted.

## Prerequisites to use it:
The algorithm is implemented in C++11. Mainly because native Python is very slow and the algorithm is a bit too complex for using pre-existing Python high-performance libraries such as NumPy. However, thanks to Pybind11 (https://pybind11.readthedocs.io) the algorithm can be used in Python and works with NumPy arrays.

What you need:
*  GCC
*  Python with NumPy
*  Pybind11 (`conda install -c conda-forge pybind11`)

## How to compile:
```
c++ -O3 -Wall -shared -std=c++11 -fPIC `python3 -m pybind11 --includes` identify.cpp -o identify`python3-config --extension-suffix`
```
Then you can copy the `.so` file into your directory with the Python code.

## How to use in Python:
For example, get clusters for convective updrafts:
```
import identify
from netCDF4 import Dataset

# constants and fields
dx = 2200
min_vel = 0.5               # if lower, it's not part of the cell
min_max_vel = 10.           # if lower, it's not real deep convection
min_height = 1250.          # no PBL stuff
max_height = 12500.         # no weird stuff in absorber
min_topheight = 2000.       # only deep convection
nc = Dataset('...')         # open the dataset - add path here
w = nc.variables['W']       # vertical wind velocity
lat = nc.variables['lat']   # latitude
lon = nc.variables['lon']   # longitude
hhl = nc.variables['HHL']   # height at half levels
hfl = 0.5 * (hhl[1:, :, :] + hhl[:-1, :, :])    # height at full levels
hsurf = nc.variables['HSURF']                   # surface height
height = hhl - hsurf        # height above surface

# Get clusters:
# Grid cells belonging to the same cluster will have the same number
# that serves as an identifier for the cluster.
clusters = identify.get_clusters(w, height, min_vel, min_max_vel,
                                 min_height, max_height, min_topheight)

# Count number of clusters and number of cells per cluster
# (grid cells with number 0 are not part of a cluster)
c_val, c_count = np.unique(clusters, return_counts=True)
c_count = c_count[c_val != 0]
c_val = c_val[c_val != 0]

# Get width and height of clusters
dict_dims = identify.get_cell_dimensions(hfl, lat, lon, clusters, c_val, dx)

# Store it in a dictionary
dict_cells = {'Height' : [],
              'Width' : [],
              'Ratio H/W' : [],
              'MaxVel' : [],
              'HeightMaxVel' : [],
              'ID' : []}

for key in dict_dims:                               
    # TODO: This if condition (width > 0 should not be necessary - check!)
    if dict_dims[key][3] > 0:
        dict_cells['Height'].append(dict_dims[key][2])
        dict_cells['Width'].append(dict_dims[key][3])
        dict_cells['Ratio H/W'].append(dict_dims[key][2] / dict_dims[key][3])
        dict_cells['MaxVel'].append(w[clusters == key].max())
        dict_cells['ID'].append(key)
```
