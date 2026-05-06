#ifndef SURFEL_H
#define SURFEL_H

#include "../helper/image_helper.h"

template <typename Vector>
struct Surfel
{
    enum Orientation
    {
        X = 0,
        Y = 1,
        Z = 2
    };

    Surfel(const std::size_t &voxel_index, const Orientation &orientation, const bool &inverse_direction = false)
        : voxel_index(voxel_index), orientation(orientation), inverse_direction(inverse_direction)
    {
    }

    void get_orientation_vectors(std::array<int, 3> &start, std::array<int, 3> &u_direction, std::array<int, 3> &v_direction) const
    {
        if (orientation == Orientation::X) // X ->YZ
        {
            start = {1, 0, 0};
            u_direction = {0, 1, 0};
            v_direction = {0, 0, 1};
        }
        if (orientation == Orientation::Y) // Y -> ZX
        {
            start = {0, 1, 0};
            u_direction = {0, 0, 1};
            v_direction = {1, 0, 0};
        }
        if (orientation == Orientation::Z) // Z -> XY
        {
            start = {0, 0, 1};
            u_direction = {1, 0, 0};
            v_direction = {0, 1, 0};
        }
    }

    std::array<int, 3> get_lowest_point(const std::size_t &xdim, const std::size_t &ydim) const
    {
        std::array<int, 3> s, u, v;
        get_orientation_vectors(s, u, v);
        const auto &[i, j, k] = helper::voxel_index_to_coord(voxel_index, xdim, ydim);
        return {static_cast<int>(i) + s[0], static_cast<int>(j) + s[1], static_cast<int>(k) + s[2]};
    }

    template <typename Point_3>
    Point_3 get_center_point(const std::size_t &xdim, const std::size_t &ydim) const
    {
        std::array<int, 3> s, u, v;
        get_orientation_vectors(s, u, v);
        const auto &[i, j, k] = helper::voxel_index_to_coord(voxel_index, xdim, ydim);
        return Point_3(i + s[0] + u[0] * 0.5 + v[0] * 0.5, j + s[1] + u[1] * 0.5 + v[1] * 0.5, k + s[2] + u[2] * 0.5 + v[2] * 0.5);
    }

    std::array<int, 3> get_normal() const
    {
        std::array<int, 3> normal;
        int v = inverse_direction ? -1 : 1;
        if (orientation == Orientation::X)
        {
            normal = {v, 0, 0};
        }
        if (orientation == Orientation::Y)
        {
            normal = {0, v, 0};
        }
        if (orientation == Orientation::Z)
        {
            normal = {0, 0, v};
        }
        return normal;
    }

    // get the four pointels with order
    //
    //          3 --u--> 2
    //          ^        ^
    //          v        v
    //          |        |
    // start -> 0 --u--> 1
    std::array<std::array<int, 3>, 4> get_pointels(const std::size_t &xdim, const std::size_t &ydim) const
    {
        std::array<int, 3> s, u, v;
        get_orientation_vectors(s, u, v);
        const auto &[i, j, k] = helper::voxel_index_to_coord(voxel_index, xdim, ydim);
        std::array<int, 3> p = {static_cast<int>(i) + s[0], static_cast<int>(j) + s[1], static_cast<int>(k) + s[2]};
        std::array<std::array<int, 3>, 4> surfel_points{
            p,
            {p[0] + u[0], p[1] + u[1], p[2] + u[2]},
            {p[0] + u[0] + v[0], p[1] + u[1] + v[1], p[2] + u[2] + v[2]},
            {p[0] + v[0], p[1] + v[1], p[2] + v[2]}};
        if (inverse_direction)
            std::swap(surfel_points[1], surfel_points[3]);
        return surfel_points;
    }

    template <typename K>
    typename K::Point_3 project(const std::size_t &xdim, const std::size_t &ydim, const typename K::Point_3 &point) const
    {
        using Point_3 = typename K::Point_3;
        using Vector_3 = typename K::Vector_3;
        using Plane_3 = typename K::Plane_3;
        // get lowest point on surfel
        std::array<int, 3> s, u, v;
        get_orientation_vectors(s, u, v);
        const auto &[i, j, k] = helper::voxel_index_to_coord(voxel_index, xdim, ydim);
        const Point_3 surfel_origin(i + s[0], j + s[1], k + s[2]);
        // Project on plane
        const std::array<int, 3> &normal_int = get_normal();
        const Vector_3 normal(normal_int[0], normal_int[1], normal_int[2]);
        const Plane_3 plane(surfel_origin, normal);
        const Point_3 &plane_projection = plane.projection(point);
        // Clamp projection in surfel (works because surfel is axis-aligned and have edge-length 1)
        const Vector_3 origin_to_projection(surfel_origin, plane_projection);
        const Vector_3 clamped_origin_to_projection(
            std::max(std::min(origin_to_projection[0], 1.0), 0.0),
            std::max(std::min(origin_to_projection[1], 1.0), 0.0),
            std::max(std::min(origin_to_projection[2], 1.0), 0.0));
        return surfel_origin + clamped_origin_to_projection;
    }

    std::size_t timestamp = static_cast<std::size_t>(-1);
    std::size_t voxel_index;
    Orientation orientation;
    bool inverse_direction = false;

    std::set<Surfel *> neighbors;

    Vector normal;
    double sharpness;
};

#endif // SURFEL_H