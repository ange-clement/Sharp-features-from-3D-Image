#ifndef GET_SMOOTHING_AND_PLANE_INTERSECTION_ENERGY_H
#define GET_SMOOTHING_AND_PLANE_INTERSECTION_ENERGY_H

#include "Get_plane_intersection_energy.h"

namespace Feature_graph_energy
{
    template <typename K, typename Vector>
    struct Get_smoothing_and_plane_intersection_energy
    {
        using Get_plane_intersection_energy = Get_plane_intersection_energy<K, Vector>;
        using Surfel_point = typename Get_plane_intersection_energy::Surfel_point;
        using Vector_2 = typename K::Vector_2;
        using Point_3 = typename K::Point_3;

        Get_plane_intersection_energy &get_plane_intersection_energy;
        const Point_3 &goal_point;
        const double &smooth_factor; // value = i + p*g;
        Get_smoothing_and_plane_intersection_energy(
            Get_plane_intersection_energy &get_plane_intersection_energy,
            const Point_3 &goal_point,
            const double &smooth_factor)
            : get_plane_intersection_energy(get_plane_intersection_energy), goal_point(goal_point),
              smooth_factor(smooth_factor)
        {
        }
        double operator()(const Surfel_point &surfel_point, const double &delta_x, const double &delta_y, const bool &is_corner, const std::size_t &line_index)
        {
            // Compute distance from moved point
            Surfel_point moved_point(surfel_point);
            moved_point.move(Vector_2(delta_x, delta_y), get_plane_intersection_energy.surface_graph);
            const Point_3 &surfel_world_point = moved_point.get_point_3(get_plane_intersection_energy.surface_graph.xdim, get_plane_intersection_energy.surface_graph.ydim);
            const double eval_distance = CGAL::squared_distance(surfel_world_point, goal_point);
            // Compute plane intersection value
            const double eval_intersection = get_plane_intersection_energy(moved_point, 0, 0, is_corner, line_index);
            return eval_intersection + smooth_factor * eval_distance;
        }
    };
} // namespace Feature_graph_energy

#endif