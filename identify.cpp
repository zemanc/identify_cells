/*******************************************************************************
 *
 * Identify clusters
 *
 * A cluster is a set of neighboring cells where all the cells have a value
 * which is above a certain treshold. This can for example be used to identify
 * convective cells with a vertical wind speed above a given velocity.
 *
 * The functions can be used with Python and NumPy and provide a significant
 * speedup compared to a Python implementation.
 * 
 * Compile it like this:
 * c++ -O3 -Wall -shared -std=c++11 -fPIC `python3 -m pybind11
 * --includes` identify.cpp -o identify`python3-config --extension-suffix`
 *
 * In Python:
 * import identify
 * clusters = identify.get_clusters(w, height, min_vel, min_max_vel, min_height,
 *                                  max_height, min_topheight)
 *
 * Christian Zeman, 2020
 *
 ******************************************************************************/

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>
#include <tuple>
#include <vector>
#include <map>
#include <cmath>
#include <iomanip>

namespace py = pybind11;

// get indices for a unique number
std::tuple<size_t, size_t, size_t> get_indices(size_t num, size_t ny, size_t nx)
{
    size_t k, j, i;
    k = size_t(num / (ny*nx));
    j = size_t((num - k*ny*nx) / nx);
    i = size_t(num - k*ny*nx - j*nx);
    return std::make_tuple(k, j, i);
}


// return identifier of subset
size_t find(const size_t *parents, size_t num, size_t ny, size_t nx)
{
    size_t k, j, i;
    std::tie(k, j, i) = get_indices(num, ny, nx);
    if (j >= ny || i >= nx)
        throw std::runtime_error("indices are wrong");
    size_t idx = k*ny*nx + j*nx + i;
    if (parents[idx] == num)
        return num;
    else {
        std::tie(k, j, i) = get_indices(parents[idx], ny, nx);
        idx = k*ny*nx + j*nx + i;
        return find(parents, parents[idx], ny, nx);
    }
}


// unite two sets
void unite(size_t *parents, size_t *ranks, size_t num_a, size_t num_b,
           size_t ny, size_t nx)
{
    size_t root_a = find(parents, num_a, ny, nx);
    size_t root_b = find(parents, num_b, ny, nx);
    size_t rank_a = ranks[root_a];
    size_t rank_b = ranks[root_b];
    if (root_a == root_b)
        return;

    if (rank_a > rank_b) {
        parents[root_b] = root_a;
    }
    else if (rank_a < rank_b) {
        parents[root_a] = root_b;
    }
    else if (rank_a == rank_b) {
        parents[root_b] = root_a;
        ranks[root_a]++;
    }

}


// calculate the distance in meter between two points on a sphere
double haversine(double lat1, double lon1, double lat2, double lon2)
{
    // distance between latitudes and longitudes
    double d_lat = (lat2 - lat1) * M_PI / 180.0;
    double d_lon = (lon2 - lon1) * M_PI / 180.0;

    // convert to radians
    lat1 = (lat1) * M_PI / 180.0;
    lat2 = (lat2) * M_PI / 180.0;

    // apply formula
    double a = std::pow(std::sin(d_lat / 2.0), 2) +
               std::pow(std::sin(d_lon / 2.0), 2) *
               std::cos(lat1) * std::cos(lat2);
    double rad = 6371000;
    double c = 2 * std::asin(std::sqrt(a));
    return rad * c;
}


// get the clusters
py::array_t<size_t> get_clusters(const py::array_t<double> values,
                                 const py::array_t<double> height_abv_sf,
                                 double min_val, double min_max_val,
                                 double min_height, double max_height,
                                 double min_top_height)
{
    // read array into buffer info and get pointer and shape
    py::buffer_info buf_values = values.request();
    py::buffer_info buf_height = height_abv_sf.request();
    double *ptr_values = (double *) buf_values.ptr;
    double *ptr_height = (double *) buf_height.ptr;
    size_t nz = buf_values.shape[0];
    size_t ny = buf_values.shape[1];
    size_t nx = buf_values.shape[2];

    // create arrays
    py::array_t<size_t> ids = py::array_t<size_t>(buf_values.size);
    py::array_t<size_t> points = py::array_t<size_t>(buf_values.size);
    py::array_t<size_t> ranks = py::array_t<size_t>(buf_values.size);
    py::array_t<size_t> parents = py::array_t<size_t>(buf_values.size);
    py::array_t<size_t> clusters = py::array_t<size_t>(buf_values.size);

    // create buffers and pointers
    py::buffer_info buf_ids = ids.request();
    py::buffer_info buf_points = points.request();
    py::buffer_info buf_ranks = ranks.request();
    py::buffer_info buf_parents = parents.request();
    py::buffer_info buf_clusters = clusters.request();
    size_t *ptr_ids = (size_t *) buf_ids.ptr;
    size_t *ptr_points = (size_t *) buf_points.ptr;
    size_t *ptr_ranks = (size_t *) buf_ranks.ptr;
    size_t *ptr_parents = (size_t *) buf_parents.ptr;
    size_t *ptr_clusters = (size_t *) buf_clusters.ptr;
   
    // Direct neighbors: In total there are 6 direct neighbors, but we only
    // have to check the ones we've already passed with the sweep. The other
    // neighbors will later do the check with the current grid cell.
    std::vector<std::tuple<int, int, int>> neighbors;
    neighbors.push_back(std::make_tuple(-1, 0, 0));
    neighbors.push_back(std::make_tuple(0, -1, 0));
    neighbors.push_back(std::make_tuple(0, 0, -1));

    // fill arrays
    size_t size = buf_ids.size;
    for (size_t id = 0; id < size; id++) {
        ptr_ids[id] = id;
        ptr_points[id] = 0;
        ptr_ranks[id] = 0;
        ptr_parents[id] = id;
        ptr_clusters[id] = id;
    }

    for (size_t k = 1; k < nz-1; k++)
        for (size_t j = 1; j < ny-1; j++)
            for (size_t i = 1; i < nx-1; i++) {
                size_t idx = k*ny*nx + j*nx + i;

                // check threshold criteria
                if (ptr_values[idx] >= min_val) {
                    ptr_points[idx] = 1;
                    ptr_parents[idx] = idx;
                    
                    // check the 6 direct neighbors
                    for (auto const& nb: neighbors) {
                        size_t pk, pj, pi, nk, nj, ni, idx_nb;
                        std::tie(pk, pj, pi) = nb; 
                        nk = k + pk;
                        nj = j + pj;
                        ni = i + pi;
                        idx_nb = nk*ny*nx + nj*nx + ni;

                        if (ptr_points[idx_nb] == 1)
                            unite(ptr_parents, ptr_ranks, idx_nb, idx, ny, nx);
                    }
                }
            }

    // flatten the thing -> set root as parent for every node
    // and while we're at it, also get the max. value and maximum height
    // of the respective cluster
    std::map<size_t, std::pair<double, double>> map_maxval;
    std::map<size_t, double> map_topheight;
    for (size_t k = 0; k < nz; k++)
        for (size_t j = 0; j < ny; j++)
            for (size_t i = 0; i < nx; i++)
            {
                size_t idx = k*ny*nx + j*nx + i;
                if (ptr_points[idx] == 1) {
                    size_t cluster_id = find(ptr_parents, idx, ny, nx);
                    ptr_clusters[idx] = cluster_id;
                    auto maxval = map_maxval.find(cluster_id);
                    if (maxval == map_maxval.end()) {
                        // insert new max value
                        map_maxval.insert(std::make_pair(cluster_id,
                                          std::make_pair(ptr_values[idx],
                                                         ptr_height[idx])));
                    } else {
                        // update max value
                        if (ptr_values[idx] > maxval->second.first)
                            maxval->second = std::make_pair(ptr_values[idx],
                                                            ptr_height[idx]);
                    }
                    auto topheight = map_topheight.find(cluster_id);
                    if (topheight == map_topheight.end()) {
                        // insert new top height
                        map_topheight.insert(
                                std::make_pair(cluster_id, ptr_height[idx]));
                    } else {
                        // update top height
                        if (ptr_height[idx] > topheight->second)
                            topheight->second = ptr_height[idx];
                    }
                } else {
                    ptr_clusters[idx] = 0;
                }
            }

    // only keep the clusters where the maximum value is above a certain
    // threshold, the height of the maximum value is within a certain range
    // and the top height of the cell is above a certain value
    for (size_t k = 1; k < nz-1; k++)
        for (size_t j = 1; j < ny-1; j++)
            for (size_t i = 1; i < nx-1; i++)
            {
                size_t idx = k*ny*nx + j*nx + i;
                if (ptr_points[idx] == 1) {
                    size_t cluster_id = ptr_clusters[idx];
                    auto maxval = map_maxval.find(cluster_id);
                    if (maxval == map_maxval.end())
                        throw std::runtime_error("cluster has no max. value");
                    auto topheight = map_topheight.find(cluster_id);
                    if (topheight == map_topheight.end())
                        throw std::runtime_error("cluster has no top height");
                    // throw the cell out if the value is below the threshold
                    if (maxval->second.first > min_max_val) {
                        ptr_points[idx] = 0;
                        ptr_clusters[idx] = 0;
                    }
                }
            }

    // reshape clusters to have the same shape as the values and return it
    clusters.resize({nz, ny, nx});
    return clusters;
}


// get the height and width of the clusters
std::map<size_t, std::tuple<double, double, double, double>>
get_cell_dimensions(const py::array_t<double> heights,
                    const py::array_t<double> lats,
                    const py::array_t<double> lons,
                    const py::array_t<size_t> clusters,
                    const py::array_t<size_t> ids,
                    double dx)
{
    // read array into buffer info
    py::buffer_info buf_heights = heights.request();
    py::buffer_info buf_lats = lats.request();
    py::buffer_info buf_lons = lons.request();
    py::buffer_info buf_clusters = clusters.request();
    py::buffer_info buf_ids = ids.request();

    // pointers
    double *ptr_heights = (double *) buf_heights.ptr;
    double *ptr_lats = (double *) buf_lats.ptr;
    double *ptr_lons = (double *) buf_lons.ptr;
    size_t *ptr_clusters = (size_t *) buf_clusters.ptr;
    size_t *ptr_ids = (size_t *) buf_ids.ptr;

    // dimensions
    size_t nz = buf_heights.shape[0];
    size_t ny = buf_heights.shape[1];
    size_t nx = buf_heights.shape[2];
    size_t nid = buf_ids.size;

    // create a map for the ids and vertical bounds
    double lb_default = 25000;
    double ub_default = 0;
    double height_default = 0;
    double width_default = 0;
    std::map<size_t, std::tuple<double, double, double, double>> map_dims;
    for (size_t i = 0; i < nid; i++)
        map_dims.insert(std::make_pair(ptr_ids[i],
                        std::make_tuple(lb_default, ub_default,
                                        height_default, width_default)));
    // loop over all values
    for (size_t k = 1; k < nz-1; k++) {
        // coordinates for calculating the width
        std::map<size_t, std::vector<std::pair<double, double>>> map_coords;

        for (size_t j = 1; j < ny-1; j++) {
            for (size_t i = 1; i < nx-1; i++) {
                // check if grid point is in a cluster
                size_t idx = k*ny*nx + j*nx + i;
                size_t idx_2d = j*nx + i;
                auto search = map_dims.find(ptr_clusters[idx]);
                if (search != map_dims.end()) {
                    // the grid point is in a cluster
                    double lb, ub, height, width;
                    size_t id = search->first;
                    std::tie(lb, ub, height, width) = search->second;
                    
                    // update lower and upper bound and height
                    size_t idx_lb = (k+1)*ny*nx + j*nx + i;
                    lb = std::min(lb, ptr_heights[idx_lb]);
                    ub = std::max(ub, ptr_heights[idx]);
                    height = ub - lb;
                    std::get<0>(search->second) = lb;
                    std::get<1>(search->second) = ub;
                    std::get<2>(search->second) = height;

                    // insert coordinates
                    map_coords[id].push_back(std::make_pair(ptr_lats[idx_2d],
                                                            ptr_lons[idx_2d]));
                }
            }
        }

        // go through all coordinates on this level to get the maximum width
        for (const auto& coords : map_coords) {
            size_t id = coords.first;
            auto& vec_coords = coords.second;
            if (vec_coords.size() == 0)
                throw std::runtime_error("coordinates vector is empty");
            else {
                double lb, ub, height, width;
                std::tie(lb, ub, height, width) = map_dims[id];
                if (vec_coords.size() == 1) {
                    // only one cell
                    width = std::max(dx, width);
                } else {
                    // calculate maximum distance between any two points
                    // of this cell on this level
                    double lat1, lat2, lon1, lon2;
                    for (auto it1 = vec_coords.begin();
                         it1 != vec_coords.end(); it1++)
                    {
                        lat1 = it1->first;
                        lon1 = it1->second;
                        for (auto it2 = std::next(it1, 1);
                             it2 != vec_coords.end(); it2++)
                        {
                            lat2 = it2->first;
                            lon2 = it2->second;
                            width = std::max(haversine(lat1, lon1, lat2, lon2)
                                             + dx, width);
                        }
                    }
                }
                // update width
                std::get<3>(map_dims[id]) = width;
            }
        }
    }
   
    // return map
    return map_dims;
}


PYBIND11_MODULE(identify, m) {
    m.doc() = "identify cells";
    m.def("get_clusters", &get_clusters, "get clusters");
    m.def("get_cell_dimensions", &get_cell_dimensions, "get height and width");
}
