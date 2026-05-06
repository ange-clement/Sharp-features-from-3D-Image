#ifndef SURFEL_POINT_H
#define SURFEL_POINT_H

#include "../helper/image_helper.h"
#include "../helper/geometry_helper.h"
#include "../surface/Surface_graph.h"

template <typename K, typename Vector>
struct Surfel_point
{
    using Point_2 = typename K::Point_2;
    using Vector_2 = typename K::Vector_2;
    using Ray_2 = typename K::Ray_2;
    using Point_3 = typename K::Point_3;
    using Vector_3 = typename K::Vector_3;

    using Surface_graph = Surface_graph<Vector>;
    using Surfel = typename Surface_graph::Surfel;

    const Surfel *surfel;
    Point_2 point; // coordinates withing the surfel (from 0 to 1 on x and y)

    static Point_3 to_point_3(
        const std::array<int, 3> &s, const std::array<int, 3> &u, const std::array<int, 3> &v,
        const std::size_t &voxel_index,
        const Point_2 &point, const std::size_t &xdim, const std::size_t &ydim)
    {
        const auto &[i, j, k] = helper::voxel_index_to_coord(voxel_index, xdim, ydim);
        return Point_3(i + s[0] + point[0] * u[0] + point[1] * v[0],
                       j + s[1] + point[0] * u[1] + point[1] * v[1],
                       k + s[2] + point[0] * u[2] + point[1] * v[2]);
    }

    static Point_3 to_point_3(const Surfel *surfel, const Point_2 &point, const std::size_t &xdim, const std::size_t &ydim)
    {
        std::array<int, 3> s, u, v;
        surfel->get_orientation_vectors(s, u, v);
        return to_point_3(s, u, v, surfel->voxel_index, point, xdim, ydim);
    }

    // point_w is within surfel
    static Point_2 get_surfel_point(const Surfel *surfel, const Point_3 &point_w, const std::size_t &xdim, const std::size_t &ydim)
    {
        std::array<int, 3> s, u, v;
        surfel->get_orientation_vectors(s, u, v);
        const auto &[i, j, k] = helper::voxel_index_to_coord(surfel->voxel_index, xdim, ydim);
        //   surfel_vector_w vector from surfel lowest point to intersection point
        const Vector_3 surfel_vector_w(point_w[0] - i - s[0], point_w[1] - j - s[1], point_w[2] - k - s[2]);
        const std::uint8_t u_axis = u[0] == 1 ? 0 : u[1] == 1 ? 1
                                                              : 2;
        const std::uint8_t v_axis = v[0] == 1 ? 0 : v[1] == 1 ? 1
                                                              : 2;
        return Point_2(surfel_vector_w[u_axis], surfel_vector_w[v_axis]);
    }

    Surfel_point(const Surfel *surfel, const Point_2 &point) : surfel(surfel), point(point) {}
    Surfel_point(const Surfel_point &other) : surfel(other.surfel), point(other.point) {}
    Surfel_point(const Surfel *surfel, const Point_3 &point_w, const Surface_graph &surface_graph) // point_w is within surfel
        : surfel(surfel)
    {
        point = Surfel_point::get_surfel_point(surfel, point_w, surface_graph.xdim, surface_graph.ydim);
    }

    Point_3 get_point_3(const std::size_t &xdim, const std::size_t &ydim) const
    {
        return Surfel_point::to_point_3(surfel, point, xdim, ydim);
    }

    // move the point on the surface by delta. Handles the change of underlying surfel.
    void move(const Vector_2 &delta, const Surface_graph &surface_graph, const Surfel *previous_surfel = nullptr)
    {
        constexpr double epsilon = 1.e-3;
        const Point_2 new_point = point + delta;
        // do not update surfel if the point has not changed surfel
        // do not update surfel if delta is too small (within +or- epsilon). This avoids infinite loops due to float precisions.
        const bool is_new_point_in_surfel = new_point[0] >= 0.0 && new_point[0] <= 1.0 && new_point[1] >= 0.0 && new_point[1] <= 1.0;
        const bool is_delta_too_small = delta[0] >= 0.0 - epsilon && delta[0] <= 0.0 + epsilon && delta[1] >= 0.0 - epsilon && delta[1] <= 0.0 + epsilon;
        if (is_new_point_in_surfel || is_delta_too_small)
        {
            point = new_point;
            return;
        }
        // point+delta is outside the surfel : change surfel and update position

        const std::size_t xdim = surface_graph.xdim;
        const std::size_t ydim = surface_graph.ydim;
        std::array<int, 3> s, u, v;
        surfel->get_orientation_vectors(s, u, v);
        const std::size_t &surfel_voxel_index = surfel->voxel_index;
        const auto &[i, j, k] = helper::voxel_index_to_coord(surfel_voxel_index, xdim, ydim);

        //  get point at intersection of the direction and the surfel box
        const Ray_2 ray(new_point, -delta);
        //   The value of ray_bbox_intersection is not initialized if the ray hit perfecly the box corner. In such cases, the starting point is already at the intersection
        std::optional<Point_2> surfel_intersection_opt = helper::get_ray_bbox_intersection_2d<K>(ray, Point_2(0, 0), Point_2(1, 1));
        if (!surfel_intersection_opt)
        {
            const Point_2 goal_project_surfel_vector(
                std::max(std::min(new_point[0], 1.0), 0.0),
                std::max(std::min(new_point[1], 1.0), 0.0));
            surfel_intersection_opt = goal_project_surfel_vector;
        }
        const Point_2 surfel_intersection = surfel_intersection_opt.value_or(new_point);
        const Vector_2 delta_left(surfel_intersection, new_point);
        const Point_3 &w_goal_point = Surfel_point::to_point_3(s, u, v, surfel_voxel_index, new_point, xdim, ydim);
        const Point_3 &intersect_w_point = Surfel_point::to_point_3(s, u, v, surfel_voxel_index, surfel_intersection, xdim, ydim);

        //  choose new surfel
        //   get linel on intersection
        std::array<int, 3> linel_start;
        std::array<int, 3> linel_end;
        if (surfel_intersection[0] + surfel_intersection[1] < 1)
            linel_start = s;
        else
            linel_start = helper::nd_vec_add(s, helper::nd_vec_add(u, v)); // s+u+v
        if (surfel_intersection[0] > surfel_intersection[1])
            linel_end = helper::nd_vec_add(s, u); // s+u
        else
            linel_end = helper::nd_vec_add(s, v); // s+v
        //   get surfels around that linel
        std::vector<const Surfel *> adjacent_surfels;
        surface_graph.get_surfels_adjacent_to_linel(
            linel_start[0] + i, linel_start[1] + j, linel_start[2] + k,
            linel_end[0] + i, linel_end[1] + j, linel_end[2] + k,
            std::back_inserter(adjacent_surfels));
        //   take closest to the goal
        const Surfel *new_surfel = nullptr;
        double new_surfel_min_sq_distance = std::numeric_limits<double>::max();
        for (const Surfel *n_surfel : adjacent_surfels)
        {
            if (n_surfel == surfel || n_surfel->timestamp == surfel->timestamp)
                continue;
            if (previous_surfel != nullptr && (n_surfel == previous_surfel || n_surfel->timestamp == previous_surfel->timestamp))
                continue;
            const double &sq_distance = CGAL::squared_distance(w_goal_point, n_surfel->template get_center_point<Point_3>(xdim, ydim));
            if (sq_distance <= new_surfel_min_sq_distance)
            {
                new_surfel = n_surfel;
                new_surfel_min_sq_distance = sq_distance;
            }
        }
        std::array<int, 3> normal = surfel->get_normal();
        if (new_surfel == nullptr)
        {
            // no neighbor on the crossing edge, point is moved along the edge
            // project goal point on surfel and get new surfel (goal point is in the plane defined by surfel)
            const Point_3 surfel_origin(i + s[0], j + s[1], k + s[2]);
            const Vector_3 origin_to_projection(surfel_origin, w_goal_point);
            const Vector_3 goal_project_surfel_vector(
                std::max(std::min(origin_to_projection[0], 1.0), 0.0),
                std::max(std::min(origin_to_projection[1], 1.0), 0.0),
                std::max(std::min(origin_to_projection[2], 1.0), 0.0));
            const Point_3 &project_w_goal_point = surfel_origin + goal_project_surfel_vector;
            // get new surfel
            new_surfel_min_sq_distance = std::numeric_limits<double>::max();
            for (const Surfel *n_surfel : surfel->neighbors)
            {
                if (previous_surfel != nullptr && (n_surfel == previous_surfel || n_surfel->timestamp == previous_surfel->timestamp))
                    continue;
                const double &sq_distance = CGAL::squared_distance(project_w_goal_point, n_surfel->template get_center_point<Point_3>(xdim, ydim));
                if (sq_distance <= new_surfel_min_sq_distance)
                {
                    new_surfel = n_surfel;
                    new_surfel_min_sq_distance = sq_distance;
                }
            }
            // if the distance is stricly greater than half the diag size, choose another surfel (d > sqrt(2)/2 is the same as d^2 > 0.5)
            if (new_surfel == nullptr || CGAL::squared_distance(intersect_w_point, new_surfel->template get_center_point<Point_3>(xdim, ydim)) > 0.5 + epsilon)
            {
                // no neighbor (either delta is too small to reach a surfel corner, or there are no neighbor around that corner)
                // convert goal_project_surfel_vector to surfel point
                const std::uint8_t u_axis = u[0] == 1 ? 0 : u[1] == 1 ? 1 : 2;
                const std::uint8_t v_axis = v[0] == 1 ? 0 : v[1] == 1 ? 1 : 2;
                const Point_2 new_point(goal_project_surfel_vector[u_axis], goal_project_surfel_vector[v_axis]);
                point = new_point;
                return;
            }
        }
        std::array<int, 3> new_surfel_normal = new_surfel->get_normal();
        // transform current point basis using new surfel's normal
        if (new_surfel_normal == u)
            u = helper::nd_vec_mult(normal, -1);
        else if (new_surfel_normal == v)
            v = helper::nd_vec_mult(normal, -1);
        else if (new_surfel_normal == helper::nd_vec_mult(u, -1))
            u = normal;
        else if (new_surfel_normal == helper::nd_vec_mult(v, -1))
            v = normal;
        // get intersection point on new surfel
        const Vector_3 delta_left_w(delta_left[0] * u[0] + delta_left[1] * v[0],
                                    delta_left[0] * u[1] + delta_left[1] * v[1],
                                    delta_left[0] * u[2] + delta_left[1] * v[2]);
        std::array<int, 3> new_s, new_u, new_v;
        new_surfel->get_orientation_vectors(new_s, new_u, new_v);
        const auto &[new_i, new_j, new_k] = helper::voxel_index_to_coord(new_surfel->voxel_index, xdim, ydim);
        //   new_surfel_vector_w vector from surfel lowest point to intersection point
        const Vector_3 new_surfel_vector_w(intersect_w_point[0] - new_i - new_s[0], intersect_w_point[1] - new_j - new_s[1], intersect_w_point[2] - new_k - new_s[2]);
        const std::uint8_t new_u_axis = new_u[0] == 1 ? 0 : new_u[1] == 1 ? 1
                                                                          : 2;
        const std::uint8_t new_v_axis = new_v[0] == 1 ? 0 : new_v[1] == 1 ? 1
                                                                          : 2;
        const Point_2 intesect_new_point(new_surfel_vector_w[new_u_axis], new_surfel_vector_w[new_v_axis]);
        const Vector_2 new_delta_left(delta_left_w[new_u_axis], delta_left_w[new_v_axis]);
        const Surfel *old_surfel = surfel;
        point = intesect_new_point;
        surfel = new_surfel;
        move(new_delta_left, surface_graph, old_surfel); // recursivelly move the delta left. Stops when point + delta is inside a surfel.
    }
};

#endif // SURFEL_POINT_H