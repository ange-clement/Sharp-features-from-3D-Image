#ifndef SURFACE_GRAPH_H
#define SURFACE_GRAPH_H

#include "Surfel.h"
#include "../helper/image_helper.h"

template <typename Vector>
struct Surface_graph
{
    typedef Surfel<Vector> Surfel;

    std::map<std::size_t, std::list<Surfel>> surface; // map of voxel to surfels (maximum 3 surfels per voxel)
    std::size_t xdim;
    std::size_t ydim;
    std::size_t zdim;
    std::size_t timestamp;
    std::size_t number_of_surfels = 0;

    Surface_graph() {}

    Surface_graph(const std::size_t &xdim, const std::size_t &ydim, const std::size_t &zdim)
        : xdim(xdim), ydim(ydim), zdim(zdim), timestamp(0)
    {
    }

    Surfel *add_surfel(const std::size_t &voxel_index, const typename Surfel::Orientation &orientation, const bool &inverse_direction = false)
    {
        return add_surfel(Surfel(voxel_index, orientation, inverse_direction));
    }

    Surfel *add_surfel(Surfel surfel)
    {
        std::list<Surfel> &voxel_surfaces = surface[surfel.voxel_index];
        voxel_surfaces.push_back(surfel);
        Surfel *surfel_ptr = &(voxel_surfaces.back());
        surfel_ptr->timestamp = timestamp;
        timestamp++;
        for (Surfel *neighbor : surfel.neighbors)
        {
            surfel_ptr->neighbors.insert(neighbor);
            neighbor->neighbors.insert(surfel_ptr);
        }
        number_of_surfels++;
        return surfel_ptr;
    }

    void remove_surfel(std::list<Surfel> &voxel_surfels, typename std::list<Surfel>::iterator &voxel_surfel_it)
    {
        // delete surface :
        //  we need to delete it from its neighbors's neighbor list
        Surfel &surfel = *voxel_surfel_it;
        for (Surfel *surfel_neigh : surfel.neighbors)
        {
            surfel_neigh->neighbors.erase(&surfel);
        }
        //  then delete the surfel
        voxel_surfel_it = voxel_surfels.erase(voxel_surfel_it);
        number_of_surfels--;
    }

    Surfel *get_first_surfel()
    {
        return &(*surface.begin()->second.begin());
    }
    const Surfel *get_first_surfel() const
    {
        return &(*surface.begin()->second.begin());
    }

    Surfel *get_surfel(const std::size_t &voxel_index, const typename Surfel::Orientation &orientation)
    {
        typename std::map<std::size_t, std::list<Surfel>>::iterator find_it = surface.find(voxel_index);
        if (find_it == surface.end())
            return nullptr;
        std::list<Surfel> &surfels = find_it->second;
        for (Surfel &s : surfels)
        {
            if (s.orientation == orientation)
                return &s;
        }
        return nullptr;
    }

    const Surfel *get_surfel(const std::size_t &voxel_index, const typename Surfel::Orientation &orientation) const
    {
        typename std::map<std::size_t, std::list<Surfel>>::const_iterator find_it = surface.find(voxel_index);
        if (find_it == surface.end())
            return nullptr;
        const std::list<Surfel> &surfels = find_it->second;
        for (const Surfel &s : surfels)
        {
            if (s.orientation == orientation)
                return &s;
        }
        return nullptr;
    }

    std::array<std::array<int, 3>, 4> get_surfel_points(const Surfel &surfel) const
    {
        return surfel.get_pointels(xdim, ydim);
    }

    template <typename Output_iterator>
    void get_surfels_adjacent_to_point(const std::size_t &x, const std::size_t &y, const std::size_t &z, Output_iterator out_iterator)
    {
        using Orientation = typename Surfel::Orientation;
        // maximum of 12 adjacent surfels contained on 7 voxels (surfel = voxel + orientation)
        std::array<std::pair<std::size_t, std::vector<Orientation>>, 7> point_surfels = {{{helper::coords_to_voxel_index(x - 1, y - 1, z - 1, xdim, ydim), {Orientation::X, Orientation::Y, Orientation::Z}},
                                                                                          {helper::coords_to_voxel_index(x, y - 1, z - 1, xdim, ydim), {Orientation::Y, Orientation::Z}},
                                                                                          {helper::coords_to_voxel_index(x - 1, y, z - 1, xdim, ydim), {Orientation::X, Orientation::Z}},
                                                                                          {helper::coords_to_voxel_index(x, y, z - 1, xdim, ydim), {Orientation::Z}},
                                                                                          {helper::coords_to_voxel_index(x - 1, y - 1, z, xdim, ydim), {Orientation::X, Orientation::Y}},
                                                                                          {helper::coords_to_voxel_index(x, y - 1, z, xdim, ydim), {Orientation::Y}},
                                                                                          {helper::coords_to_voxel_index(x - 1, y, z, xdim, ydim), {Orientation::X}}}};
        typename std::map<std::size_t, std::list<Surfel>>::const_iterator surface_it;
        for (const auto &[voxel_id, orientations] : point_surfels)
        {
            for (const Orientation &orientation : orientations)
            {
                auto surfel = get_surfel(voxel_id, orientation);
                if (surfel != nullptr)
                    *out_iterator++ = surfel;
            }
        }
    }

    template <typename Output_iterator>
    void get_surfels_adjacent_to_point(const std::size_t &x, const std::size_t &y, const std::size_t &z, Output_iterator out_iterator) const
    {
        return const_cast<Surface_graph *>(this)->get_surfels_adjacent_to_point(x, y, z, out_iterator);
    }

    template <typename Output_iterator>
    void get_surfels_adjacent_to_linel(
        const std::size_t &xs, const std::size_t &ys, const std::size_t &zs,
        const std::size_t &xe, const std::size_t &ye, const std::size_t &ze,
        Output_iterator out_iterator)
    {
        using Orientation = typename Surfel::Orientation;
        const std::uint8_t linel_axis = (xs != xe) ? 0 : ((ys != ye) ? 1 : 2);
        const std::uint8_t linel_axis_next = (linel_axis + 1) % 3;
        const std::uint8_t linel_axis_next_next = (linel_axis + 2) % 3;
        const long xmin = std::min(xs, xe);
        const long ymin = std::min(ys, ye);
        const long zmin = std::min(zs, ze);
        std::array<long, 3> coord_voxel_1 = {xmin - 1, ymin - 1, zmin - 1};
        coord_voxel_1[linel_axis]++;
        std::array<long, 3> coord_voxel_2 = {xmin, ymin, zmin};
        coord_voxel_2[linel_axis_next]--;
        std::array<long, 3> coord_voxel_3 = {xmin, ymin, zmin};
        coord_voxel_3[linel_axis_next_next]--;
        // maximum of 4 adjacent surfels contained on 3 voxels (surfel = voxel + orientation)
        std::array<std::pair<std::size_t, std::vector<Orientation>>, 3> point_surfels = {{
            {static_cast<std::size_t>(helper::coords_to_voxel_index(coord_voxel_1[0], coord_voxel_1[1], coord_voxel_1[2], static_cast<long>(xdim), static_cast<long>(ydim))),
             {static_cast<Orientation>(linel_axis_next), static_cast<Orientation>(linel_axis_next_next)}},
            {static_cast<std::size_t>(helper::coords_to_voxel_index(coord_voxel_2[0], coord_voxel_2[1], coord_voxel_2[2], static_cast<long>(xdim), static_cast<long>(ydim))),
             {static_cast<Orientation>(linel_axis_next)}},
            {static_cast<std::size_t>(helper::coords_to_voxel_index(coord_voxel_3[0], coord_voxel_3[1], coord_voxel_3[2], static_cast<long>(xdim), static_cast<long>(ydim))),
             {static_cast<Orientation>(linel_axis_next_next)}},
        }};
        typename std::map<std::size_t, std::list<Surfel>>::const_iterator surface_it;
        for (const auto &[voxel_id, orientations] : point_surfels)
        {
            for (const Orientation &orientation : orientations)
            {
                auto surfel = get_surfel(voxel_id, orientation);
                if (surfel != nullptr)
                    *out_iterator++ = surfel;
            }
        }
    }
    template <typename Output_iterator>
    void get_surfels_adjacent_to_linel(
        const std::size_t &xs, const std::size_t &ys, const std::size_t &zs,
        const std::size_t &xe, const std::size_t &ye, const std::size_t &ze,
        Output_iterator out_iterator) const
    {
        return const_cast<Surface_graph *>(this)->get_surfels_adjacent_to_linel(xs, ys, zs, xe, ye, ze, out_iterator);
    }

    template <typename Point>
    Surfel *get_surfel_from_midpoint(const Point &p)
    {
        using Orientation = typename Surfel::Orientation;
        std::size_t i2 = static_cast<std::size_t>(p[0] * 2);
        std::size_t j2 = static_cast<std::size_t>(p[1] * 2);
        std::size_t k2 = static_cast<std::size_t>(p[2] * 2);
        std::size_t voxel_index = helper::coords_to_voxel_index<std::size_t>((i2 - 1) / 2, (j2 - 1) / 2, (k2 - 1) / 2, xdim, ydim); // -1 ?
        Orientation orientation;
        if ((i2 & 1) == 0)
        {
            orientation = Orientation::X;
        }
        else if ((j2 & 1) == 0)
        {
            orientation = Orientation::Y;
        }
        else if ((k2 & 1) == 0)
        {
            orientation = Orientation::Z;
        }
        return get_surfel(voxel_index, orientation);
    }

    template <typename Point>
    const Surfel *get_surfel_from_midpoint(const Point &p) const
    {
        return const_cast<Surface_graph *>(this)->get_surfel_from_midpoint(p);
    }

    template <typename Point, typename OutputIterator>
    void get_surfels_around(const Point &p, const std::size_t &range, OutputIterator surfels_iterator)
    {
        const int x = static_cast<int>(p[0]);
        const int y = static_cast<int>(p[1]);
        const int z = static_cast<int>(p[2]);
        const std::size_t start_search_x = static_cast<std::size_t>(std::max<int>(x - range + 1, 0));
        const std::size_t start_search_y = static_cast<std::size_t>(std::max<int>(y - range + 1, 0));
        const std::size_t start_search_z = static_cast<std::size_t>(std::max<int>(z - range + 1, 0));
        const std::size_t end_search_x = std::min<std::size_t>(x + range + 1, xdim);
        const std::size_t end_search_y = std::min<std::size_t>(y + range + 1, ydim);
        const std::size_t end_search_z = std::min<std::size_t>(z + range + 1, zdim);
        for (std::size_t sz = start_search_z; sz < end_search_z; sz++)
            for (std::size_t sy = start_search_y; sy < end_search_y; sy++)
                for (std::size_t sx = start_search_x; sx < end_search_x; sx++)
                {
                    get_surfels_adjacent_to_point(sx, sy, sz, surfels_iterator);
                }
    }

    // works only if the point is close enough to surfels (under search_range voxel away)
    template <typename K>
    std::pair<Surfel *, typename K::Point_3> project_on_surfels(const typename K::Point_3 &p, const std::size_t &search_range = 10)
    {
        using Point = typename K::Point_3;
        const typename K::Compute_squared_distance_3 sq_distance_functor = K().compute_squared_distance_3_object();
        std::set<Surfel *> close_surfels;
        Surfel *closest_surfel = nullptr;
        Point projection_point;
        std::size_t curent_search_range = 1;
        while (closest_surfel == nullptr && curent_search_range <= search_range)
        {
            get_surfels_around(p, curent_search_range, std::inserter(close_surfels, close_surfels.end()));
            double min_sq_distance = std::numeric_limits<double>::max();
            for (Surfel *surfel : close_surfels)
            {
                const Point &surfel_projection_point = surfel->template project<K>(xdim, ydim, p);
                double sq_distance = sq_distance_functor(p, surfel_projection_point);
                if (sq_distance < min_sq_distance)
                {
                    min_sq_distance = sq_distance;
                    closest_surfel = surfel;
                    projection_point = surfel_projection_point;
                }
            }
            curent_search_range++;
        }
        return {closest_surfel, projection_point};
    }

    template <typename K>
    std::pair<const Surfel *, typename K::Point_3> project_on_surfels(const typename K::Point_3 &p, const std::size_t &search_range = 10) const
    {
        return const_cast<Surface_graph *>(this)->project_on_surfels<K>(p, search_range);
    }

    // works only if the point is close enough to surfels (under search_range voxel away)
    // compute closest to surfels center
    template <typename Point>
    Surfel *get_closest_surfel(const Point &p, const std::size_t &search_range = 10)
    {
        std::set<Surfel *> close_surfels;
        Surfel *closest_surfel = nullptr;
        std::size_t curent_search_range = 1;
        while (closest_surfel == nullptr && curent_search_range <= search_range)
        {
            get_surfels_around(p, curent_search_range, std::inserter(close_surfels, close_surfels.end()));
            double min_sq_distance = std::numeric_limits<double>::max();
            for (Surfel *surfel : close_surfels)
            {
                const Point &surfel_point = surfel->template get_center_point<Point>(xdim, ydim);
                double sq_distance = CGAL::squared_distance(p, surfel_point);
                if (sq_distance < min_sq_distance)
                {
                    min_sq_distance = sq_distance;
                    closest_surfel = surfel;
                }
            }
            curent_search_range++;
        }
        return closest_surfel;
    }

    template <typename Point>
    const Surfel *get_closest_surfel(const Point &p, const std::size_t &search_range = 10) const
    {
        return const_cast<Surface_graph *>(this)->get_closest_surfel<Point>(p, search_range);
    }
};

#endif // SURFACE_GRAPH_H