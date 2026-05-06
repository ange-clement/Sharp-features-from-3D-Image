#ifndef POLYLINE_OPTIMIZER_H
#define POLYLINE_OPTIMIZER_H

#include "../polyline/Polylines.h"
#include "../polyline/polyline_snap.h"
#include "../surface/Surface_graph.h"
#include "../helper/image_helper.h"
#include "../helper/cgal_helper.h"
#include "energy/Get_plane_intersection_energy.h"
#include "energy/Get_smoothing_and_plane_intersection_energy.h"
#include "Surfel_point.h"

template <typename K>
struct Optimize_polyline_by_plane_patches
{
    std::string out_filepath = "";
    std::size_t maximum_iteration = 20;
    double min_energy_delta = 1.e-3;
    double start_step_size = 1.0;
    double step_size_factor = 0.9;
    double collapse_distance = 0.5;
    double smooth_factor = 1.0;
    double refine_surfel_distance = 4.0;
    double plane_detection_max_distance = 4.0;

    // returns the sharpest descent vector in the space of the surfel point
    template <typename Vector, typename EnergyFunctor>
    typename K::Vector_2 get_minimizing_gradient(
        const Surfel_point<K, Vector>& surfel_point,
        const double& step_size,
        const bool& is_corner, const std::size_t& line_index,
        EnergyFunctor& get_energy = EnergyFunctor()) const
    {
        return typename K::Vector_2(
            get_energy(surfel_point, -step_size, 0, is_corner, line_index)-get_energy(surfel_point, step_size, 0, is_corner, line_index),
            get_energy(surfel_point, 0, -step_size, is_corner, line_index)-get_energy(surfel_point, 0, step_size, is_corner, line_index)
            );
    }

    void initialize_step_size_factor(
        const double start_step_size = 1.0,
        const double end_step_size = 0.125,
        const std::size_t number_of_iterations = 20)
    {
        step_size_factor = std::pow(end_step_size / start_step_size, 1.0 / static_cast<double>(number_of_iterations));
    }

    template <typename Point, typename Vector>
    Polylines<Point> operator()(
        const Polylines<Point> &input_polylines,
        const Surface_graph<Vector> &surface_graph) const
    {
        using Surfel_point = Surfel_point<K, Vector>;
        const std::size_t &xdim = surface_graph.xdim;
        const std::size_t &ydim = surface_graph.ydim;
        Polylines<Point> polylines = input_polylines; // intiailize as a copy
        std::vector<Point> &points = polylines.points;
        // Initialize 2D points on surface
        std::vector<Surfel_point> points_on_surfels;
        points_on_surfels.reserve(points.size());
        for (const Point &point : points)
        {
            const auto &[projected_surfel, projected_point] = surface_graph.template project_on_surfels<K>(point);
            points_on_surfels.emplace_back(projected_surfel, projected_point, surface_graph);
        }

        std::size_t iteration = 0;
        std::set<std::size_t> processed_corners;
        double check_step_size = start_step_size;
        double previous_system_energy = std::numeric_limits<double>::max();
        double system_energy;
        double system_energy_delta = std::numeric_limits<double>::max();
        while (iteration < maximum_iteration && std::abs(system_energy_delta) > min_energy_delta)
        {
            helper::bBlock("iteration " + std::to_string(iteration));

            helper::bBlock("Gradient descent");
            const double check_step_sq_size = check_step_size * check_step_size;
            system_energy = 0;
            // the structure to compute the energy is recreated at each iterations. It performs the normal sharpening
            Feature_graph_energy::Get_plane_intersection_energy<K, Vector>
                get_intersection_energy(surface_graph, polylines, plane_detection_max_distance);
            processed_corners.clear();
            for (std::size_t line_index = 0, number_of_lines = polylines.lines.size(); line_index < number_of_lines; line_index++)
            {
                const std::vector<std::size_t> &line = polylines.lines[line_index];
                for (std::size_t point_in_line_index = 0, line_size = line.size(); point_in_line_index < line_size; point_in_line_index++)
                {
                    const std::size_t &point_index = line[point_in_line_index];
                    const bool &is_corner = point_in_line_index == 0 || point_in_line_index == line_size - 1;
                    if (is_corner)
                    {
                        const auto &[it, inserted] = processed_corners.insert(point_index);
                        if (!inserted)
                            continue;
                    }
                    Point &point = points[point_index];
                    Surfel_point &surfel_point = points_on_surfels[point_index];
                    Point smooth_point;
                    double point_smooth_factor;
                    if (is_corner)
                    {
                        smooth_point = point;
                        point_smooth_factor = 0.0;
                    }
                    else
                    {
                        const Point &previous_point = points[line[point_in_line_index - 1]];
                        const Point &next_point = points[line[point_in_line_index + 1]];
                        smooth_point = CGAL::midpoint(previous_point, next_point);
                        point_smooth_factor = smooth_factor;
                    }
                    Feature_graph_energy::Get_smoothing_and_plane_intersection_energy<K, Vector>
                        get_smoothing_and_intersection_energy(get_intersection_energy, smooth_point, point_smooth_factor);
                    typename K::Vector_2 gradient_step = get_minimizing_gradient<Vector>(surfel_point, check_step_size, is_corner, line_index, get_smoothing_and_intersection_energy);
                    // clamp the length of the energy gradient to check_step_size
                    double gradient_step_sq_len = gradient_step.squared_length();
                    if (gradient_step_sq_len > check_step_sq_size)
                    {
                        gradient_step = check_step_size / CGAL::sqrt(gradient_step_sq_len) * gradient_step;
                    }
                    surfel_point.move(gradient_step, surface_graph);
                    point = surfel_point.get_point_3(xdim, ydim);
                    // compute the energy after the move to get the current system energy
                    const double &cur_energy = get_smoothing_and_intersection_energy(surfel_point, 0, 0, is_corner, line_index);
                    system_energy += cur_energy;
                }
            }
            system_energy_delta = system_energy - previous_system_energy;
            previous_system_energy = system_energy;
            check_step_size *= step_size_factor;
            helper::eBlock();
#ifdef VERBOSE
            std::cout << system_energy << " ( " << (system_energy_delta >= 0 ? '+' : '-') << std::abs(system_energy_delta) << ")" << std::endl;
#endif
            helper::bLine("collapse");
            Polylines_process::collapse_short_edges(collapse_distance, polylines);
            helper::eLine();
            iteration++;
            helper::eBlock();
        }

        // helper::bLine("snap lines");
        // std::cout << "TODO" << std::endl;
        // Polylines_process::snap_lines<K>(collapse_distance * 0.5, polylines);
        // helper::eLine();

        return Polylines_process::make_polylines_without_floating_points_and_empty_lines(polylines);
    }
};

#endif // POLYLINE_OPTIMIZER_H