#ifndef SHARPNESS_MEASURE_H
#define SHARPNESS_MEASURE_H

#include "helper/dgtal_helper.h" // include before CGAL

#include <CGAL/Image_3.h>

#include "surface/Surface_graph.h"
#include "polyline/Polylines.h"

namespace internal
{
    template <typename Point, typename Vector>
    void add_values_of_pointels_to_graph(const std::vector<Point> &pointels, const std::vector<double> &values, Surface_graph<Vector> &surface_graph)
    {
        using Surfel = typename Surface_graph<Vector>::Surfel;
        std::map<std::size_t, std::list<Surfel>>::iterator surface_it;
        // reset surfels values
        for (auto &[voxel_id, voxel_surfels] : surface_graph.surface)
            for (Surfel &surfel : voxel_surfels)
            {
                surfel.sharpness = 0;
            }
        // for each pointel, find the adjacents surfels and add value
        std::vector<Surfel *> pointel_adjacent_surfels;
        pointel_adjacent_surfels.reserve(12);
        for (std::size_t i = 0, size = pointels.size(); i < size; i++)
        {
            const Point &p = pointels[i];
            const double &v = values[i];
            const std::size_t x = static_cast<std::size_t>(p[0]);
            const std::size_t y = static_cast<std::size_t>(p[1]);
            const std::size_t z = static_cast<std::size_t>(p[2]);
            pointel_adjacent_surfels.clear();
            surface_graph.get_surfels_adjacent_to_point(x, y, z, std::back_inserter(pointel_adjacent_surfels));
            for (Surfel *surfel : pointel_adjacent_surfels)
            {
                surfel->sharpness += v;
            }
        }
        for (auto &[voxel_id, voxel_surfels] : surface_graph.surface)
            for (Surfel &surfel : voxel_surfels)
            {
                surfel.sharpness /= 4;
            }
    }

    template <typename Point, typename InputNormalVector, typename SurfelNormalvector>
    void add_normals_of_surfel_to_graph(const std::vector<Point> &surfels, const std::vector<InputNormalVector> &normals, Surface_graph<SurfelNormalvector> &surface_graph)
    {
        using Surfel = typename Surface_graph<SurfelNormalvector>::Surfel;
        // for each surfel, find the graph's surfel and set normal
        for (std::size_t i = 0, size = surfels.size(); i < size; i++)
        {
            const Point &surfel_pos = surfels[i];
            const InputNormalVector &normal = normals[i];
            Surfel *surfel = surface_graph.get_surfel_from_midpoint(surfel_pos);
            surfel->normal = SurfelNormalvector(normal[0], normal[1], normal[2]);
        }

        // correct null normals
        bool corrected_normal = true;
        while (corrected_normal)
        {
            corrected_normal = false;
            bool has_null_normal = false;
            for (auto &[voxel_id, voxel_surfels] : surface_graph.surface)
                for (Surfel &surfel : voxel_surfels)
                {
                    if (surfel.normal.squared_length() != 0)
                        continue;
                    has_null_normal = true;
                    // found surfel with null normal. Assign to a neighbor with a normal
                    int nb_normal = 0;
                    SurfelNormalvector avg_normal(0, 0, 0);
                    for (const Surfel *n_surfel : surfel.neighbors)
                    {
                        const SurfelNormalvector &n_normal = n_surfel->normal;
                        if (n_normal.squared_length() == 0)
                            continue;
                        nb_normal++;
                        avg_normal += n_normal;
                    }
                    if (nb_normal != 0)
                    {
                        corrected_normal = true;
                        surfel.normal = avg_normal / static_cast<double>(nb_normal);
                    }
                }
            if (has_null_normal && !corrected_normal)
            {
                for (auto &[voxel_id, voxel_surfels] : surface_graph.surface)
                    for (Surfel &surfel : voxel_surfels)
                    {
                        if (surfel.normal.squared_length() != 0)
                            continue;
                        const std::array<int, 3> &image_normal = surfel.get_normal();
                        surfel.normal = SurfelNormalvector(image_normal[0], image_normal[1], image_normal[2]);
                    }
            }
        }
    }

    void make_pointel_value_map(
        const std::vector<helper::DGtal_types::RealPoint> &pointels,
        const std::vector<double> &values,
        std::map<std::array<int, 3>, double> &point_values)
    {
        point_values.clear();
        for (std::size_t i = 0, size = pointels.size(); i < size; i++)
        {
            const helper::DGtal_types::RealPoint &p = pointels[i];
            const double &v = values[i];
            point_values.emplace(std::array<int, 3>({static_cast<int>(p[0]), static_cast<int>(p[1]), static_cast<int>(p[2])}), v);
        }
    }

    void convert_dgtal_coordinates(
        std::vector<helper::DGtal_types::RealPoint> &coordinates)
    {
        for (helper::DGtal_types::RealPoint &coord : coordinates)
            coord = helper::DGtal_types::RealPoint(coord[0] + 0.5, coord[1] + 0.5, coord[2] + 0.5);
    }
} // namespace internal

struct AT_sharpness_measure
{
    double alpha = 0.01;
    double lambda = 0.02;
    double width_start = 2.0;
    double width_end = 0.5;
    double width_divide = 2.0;

    template <typename Vector>
    void operator()(
        const CGAL::Image_3 &image,
        Surface_graph<Vector> &surface_graph,
        std::map<std::array<int, 3>, double> &point_values) const
    {
        // get normals and values
        helper::DGtal_types::Image dgtal_image = helper::cgal_image_to_dgtal_image(image);
        std::vector<helper::DGtal_types::RealPoint> pointels;
        std::vector<double> values; // on pointels
        std::vector<helper::DGtal_types::RealPoint> surfels;
        std::vector<helper::DGtal_types::RealVector> normals; // on surfels
        helper::DGtal_types::Surface_info surface;
        helper::get_surface_from_image(
            /*In */ dgtal_image,
            /*Out*/ surface);
        helper::compute_AT_normals_and_values_from_surface(
            /*In */ surface,
            /*Out*/ surfels, pointels, normals, values, alpha, lambda, width_start, width_end, width_divide);

        // add values to surfel graph
        internal::convert_dgtal_coordinates(pointels);
        internal::convert_dgtal_coordinates(surfels);
        internal::make_pointel_value_map(pointels, values, point_values);
        internal::add_values_of_pointels_to_graph(pointels, values, surface_graph);
        internal::add_normals_of_surfel_to_graph(surfels, normals, surface_graph);
    }
};

struct Curvature_sharpness_measure
{
    double radius = 7.0;
    double grid_step = 1.0;

    template <typename Vector>
    void operator()(
        const CGAL::Image_3 &image,
        Surface_graph<Vector> &surface_graph,
        std::map<std::array<int, 3>, double> &point_values) const
    {
        // get values
        std::vector<double> values;
        helper::DGtal_types::Image dgtal_image = helper::cgal_image_to_dgtal_image(image);
        std::vector<helper::DGtal_types::RealPoint> pointels;
        std::vector<helper::DGtal_types::RealPoint> surfels;
        helper::compute_integral_invariant_values(dgtal_image, surfels, pointels, values, radius, grid_step);

        // normalize values
        const auto &[min_ptr, max_ptr] = std::minmax_element(values.begin(), values.end());
        double min = *min_ptr;
        double max = *max_ptr;
        double inv_range = 1.0 / (max - min);
        std::for_each(values.begin(), values.end(),
                      [&min, &inv_range](double &v)
                      { v = (v - min) * inv_range; });

        // get normals
        std::vector<helper::DGtal_types::RealPoint> normal_surfels;
        std::vector<helper::DGtal_types::RealPoint> normals;
        helper::compute_normals_from_surface(dgtal_image, normal_surfels, normals);

        // add values to surfel graph
        internal::convert_dgtal_coordinates(pointels);
        internal::convert_dgtal_coordinates(normal_surfels);
        internal::make_pointel_value_map(pointels, values, point_values);
        internal::add_values_of_pointels_to_graph(pointels, values, surface_graph);
        internal::add_normals_of_surfel_to_graph(normal_surfels, normals, surface_graph);
    }
};

struct VCM_sharpness_measure
{
    double offset_radius = 5.0;
    double smooth_radius = 4.0;

    template <typename Vector>
    void operator()(
        const CGAL::Image_3 &image,
        Surface_graph<Vector> &surface_graph,
        std::map<std::array<int, 3>, double> &point_values) const
    {
        // get values
        std::vector<double> values;
        helper::DGtal_types::Image dgtal_image = helper::cgal_image_to_dgtal_image(image);
        std::vector<helper::DGtal_types::RealPoint> pointels;
        std::vector<helper::DGtal_types::RealPoint> surfels;
        helper::compute_VCM_values(dgtal_image, pointels, surfels, values, offset_radius, smooth_radius);

        // normalize values
        const auto &[min_ptr, max_ptr] = std::minmax_element(values.begin(), values.end());
        double min = *min_ptr;
        double max = *max_ptr;
        double inv_range = 1.0 / (max - min);
        std::for_each(values.begin(), values.end(),
                      [&min, &inv_range](double &v)
                      { v = (v - min) * inv_range; });

        // get normals
        std::vector<helper::DGtal_types::RealPoint> normal_surfels;
        std::vector<helper::DGtal_types::RealPoint> normals;
        helper::compute_normals_from_surface(dgtal_image, normal_surfels, normals);

        // add values to surfel graph
        internal::convert_dgtal_coordinates(pointels);
        internal::convert_dgtal_coordinates(normal_surfels);
        internal::make_pointel_value_map(pointels, values, point_values);
        internal::add_values_of_pointels_to_graph(pointels, values, surface_graph);
        internal::add_normals_of_surfel_to_graph(normal_surfels, normals, surface_graph);
    }
};

template <typename Point>
struct Feedback_loop_AT_measure
{
    double alpha = 0.01;
    double lambda = 0.02;
    double width_start = 2.0;
    double width_end = 0.5;
    double width_divide = 2.0;
    double refine_surfel_distance = 4.0;

    const Polylines<Point> &feature_graph;

    Feedback_loop_AT_measure(const Polylines<Point> &feature_graph)
        : feature_graph(feature_graph)
    {}

    template <typename Vector>
    void operator()(
        const CGAL::Image_3 &image,
        Surface_graph<Vector> &surface_graph,
        std::map<std::array<int, 3>, double> &point_values) const
    {
        using Surface_graph = Surface_graph<Vector>;
        using Surfel = typename Surface_graph::Surfel;
        // get initial normals
        helper::DGtal_types::Image dgtal_image = helper::cgal_image_to_dgtal_image(image);
        helper::DGtal_types::Surface_info surface;
        helper::get_surface_from_image(
            /*In */ dgtal_image,
            /*Out*/ surface);
        std::vector<helper::DGtal_types::RealPoint> normal_surfels;
        std::vector<helper::DGtal_types::RealPoint> normals; // on surfels
        helper::compute_normals_from_surface(dgtal_image, normal_surfels, normals);
        internal::convert_dgtal_coordinates(normal_surfels);
        internal::add_normals_of_surfel_to_graph(normal_surfels, normals, surface_graph);
        // get normals from normal optimization using the feature graph
        //   get frontier surfels
        std::set<const Surfel *> corner_surfels;
        std::vector<std::set<const Surfel *>> line_frontier_surfels;
        std::unordered_set<const Surfel *> frontier_surfels;
        Surface_process::get_surfel_frontier(surface_graph, feature_graph, corner_surfels, line_frontier_surfels, frontier_surfels);

        {
            std::set<const Surfel*> frontier_surfels_set;
            frontier_surfels_set.insert(frontier_surfels.begin(), frontier_surfels.end());
            Surface_process::write_surface_ply<helper::CGAL_types::K>("tmp/feedback_frontier.ply", surface_graph, Surface_process::internal::Surfels_in_set<Vector>(frontier_surfels_set));
        }

        //   get surfel to sharpen
        std::set<const Surfel*> surfels_to_sharpen;
        Surface_process::get_surfels_near_frontier<Point>(surface_graph, frontier_surfels, refine_surfel_distance, std::inserter(surfels_to_sharpen, surfels_to_sharpen.begin()));

        {
            Surface_process::write_surface_ply<helper::CGAL_types::K>("tmp/feedback_sharpen.ply", surface_graph, Surface_process::internal::Surfels_in_set<Vector>(surfels_to_sharpen));
        }

        //   sharpen surfel normals
        const std::map<const Surfel *, Vector> &sharpened_normals = Surface_process::sharpen_normals(surface_graph, surfels_to_sharpen, frontier_surfels);
        Surface_process::write_surface_graph_infos<helper::CGAL_types::K, Vector>(surface_graph, "tmp/feedback_surface_before");
        Surface_process::apply_normals(sharpened_normals, surface_graph);

        Surface_process::write_surface_graph_infos<helper::CGAL_types::K, Vector>(surface_graph, "tmp/feedback_surface_sharpen");

        //   output sharpened normals in the DGtal surface
        surface.get_surfel_points(normal_surfels);
        Surface_process::extract_normals(surface_graph, normal_surfels, surface.normals);
        // compute using the sharpned normals
        std::vector<helper::DGtal_types::RealPoint> pointels;
        std::vector<double> values; // on pointels
        helper::compute_AT_normals_and_values_from_surface(
            /*In */ surface,
            /*Out*/ normal_surfels, pointels, normals, values, alpha, lambda, width_start, width_end, width_divide);

        // add values to surfel graph
        internal::convert_dgtal_coordinates(pointels);
        internal::convert_dgtal_coordinates(normal_surfels);
        internal::make_pointel_value_map(pointels, values, point_values);

        internal::add_values_of_pointels_to_graph(pointels, values, surface_graph);
        internal::add_normals_of_surfel_to_graph(normal_surfels, normals, surface_graph);
    }
};

#endif // SHARPNESS_MEASURE_H