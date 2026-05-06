#ifndef SHARP_FEATURE_EXTRACTION_H
#define SHARP_FEATURE_EXTRACTION_H

#include <boost/filesystem.hpp>

#include <CGAL/Image_3.h>

#include "polyline/Polylines.h"
#include "surface/Surface_graph.h"
#include "surface/surface_process.h"
#include "helper/cgal_helper.h"

#include "feature_graph/surfel_selector.h"
#include "feature_graph/polyline_extractor.h"
#include "feature_graph/polyline_regularizer.h"
#include "feature_graph/polyline_optimizer.h"

struct TODO
{
};

template <typename SharpnessMeasure, typename SurfaceSelector, typename PolylineExtractor, typename RegularizationFunctor, typename PlacementOptimizer>
Polylines<helper::CGAL_types::Point> detect_sharp_edges_in_image(
    const CGAL::Image_3 &image,
    const SharpnessMeasure &sharpness_measure,
    const SurfaceSelector &surface_selector,
    const PolylineExtractor &polyline_extractor,
    const RegularizationFunctor &regularization_functor,
    const PlacementOptimizer &placement_optimizer,
    const std::string &out_filepath = "")
{
    typedef helper::CGAL_types::K K;
    typedef helper::CGAL_types::Point Point;
    typedef helper::CGAL_types::Vector Vector;
    typedef Surface_graph<Vector> Surface_graph;

    const boost::filesystem::path out_folder = out_filepath;

    helper::bBlock("Create surfel graph");
    Surface_graph surface_graph = Surface_process::make_surface_graph_from_image<Vector>(image, /*pointels_neighbors = */ false);
    helper::eBlock();
    if (!out_filepath.empty())
    {
        Surface_process::write_surface_ply<K, Vector>(boost::filesystem::path(out_folder).append("surface.ply").string(), surface_graph);
    }

    helper::bBlock("Get normals and sharpness");
    std::map<std::array<int, 3>, double> pointel_values;
    sharpness_measure(image, surface_graph, pointel_values);
    helper::eBlock();
    if (!out_folder.empty())
    {
        Surface_process::write_surface_graph_infos<K, Vector>(surface_graph, boost::filesystem::path(out_folder).append("surface").string());
    }

    helper::bBlock("Select surfels");
    std::vector<const Surface_graph::Surfel *> deleted_surfels;
    Surface_graph sharp_surface_graph = surface_selector(surface_graph, std::back_inserter(deleted_surfels));
    helper::eBlock();
    if (!out_folder.empty())
    {
        Surface_process::write_surface_ply<K>(boost::filesystem::path(out_folder).append("selected_surfels.ply").string(), sharp_surface_graph);
        Surface_process::write_surface_graph_infos<K, Vector>(sharp_surface_graph, boost::filesystem::path(out_folder).append("selected_surfels").string());

        std::set<const Surface_graph::Surfel *> surfels_to_display;
        surfels_to_display.insert(deleted_surfels.begin(), deleted_surfels.end());
        Surface_process::write_surface_ply<K>(boost::filesystem::path(out_folder).append("deleted_surfels.ply").string(), surface_graph, Surface_process::internal::Surfels_in_set<Vector>(surfels_to_display));
    }

    helper::bBlock("Make initial feature graph");
    Polylines<Point> sharp_features = polyline_extractor.template operator()<Point, Vector>(sharp_surface_graph, deleted_surfels, pointel_values);
    helper::eBlock();
    if (!out_filepath.empty())
    {
        Polylines_process::write_polylines<Point, Vector>(sharp_features,
                                                             boost::filesystem::path(out_folder).append("sharp_features_initial").string());
    }

    helper::bBlock("Regularize double lines");
    Polylines<Point> sharp_features_regularized = regularization_functor(sharp_features);
    helper::eBlock();
    if (!out_filepath.empty())
    {
        Polylines_process::write_polylines<Point, Vector>(sharp_features_regularized,
                                                             boost::filesystem::path(out_folder).append("sharp_features_regularized").string());
    }

    helper::bBlock("Optimize positions");
    Polylines<Point> sharp_features_optimized = placement_optimizer(sharp_features_regularized, surface_graph);
    helper::eBlock();
    if (!out_filepath.empty())
    {
        Polylines_process::write_polylines<Point, Vector>(sharp_features_optimized,
                                                             boost::filesystem::path(out_folder).append("sharp_features_optimized").string());
    }

    helper::bBlock("Place polylines in world space");
    Polylines<Point> sharp_features_world_space = Polylines_process::make_world_polyline<Vector>(sharp_features_optimized, image);
    helper::eBlock();
    if (!out_filepath.empty())
    {
        Polylines_process::write_polylines<Point, Vector>(sharp_features_world_space,
                                                             boost::filesystem::path(out_folder).append("sharp_features_world_space").string());
    }

    return sharp_features_world_space;
}

template <typename SharpnessMeasure>
Polylines<helper::CGAL_types::Point> detect_sharp_edges_in_image(
    const CGAL::Image_3 &image,
    // measure parameters
    const SharpnessMeasure &sharpness_measure,
    // selection parameters
    const double &selection_threshold = 0.25,
    // regularization parameters
    const double &regularization_distance = 4.0,
    // optimization parameters
    const std::size_t &maximum_iteration = 20,
    const double &start_step_size = 1.0,
    const double &end_step_size = 0.125,
    const double &min_energy_delta = 1.e-3,
    const double &collapse_distance = 0.5,
    const double &smooth_factor = 1.0,
    const double &refine_surfel_distance = 4.0,
    const double &plane_detection_max_distance = 4.0,
    const std::string &out_filepath = "")
{
    Threshold_surfel_selector surface_selector;
    surface_selector.selection_threshold = selection_threshold;
    Sharpness_oriented_thinning_on_element_graph polyline_extractor;
    polyline_extractor.out_filepath = out_filepath;
    Double_lines_regularizer<helper::CGAL_types::Vector> regularization_functor;
    regularization_functor.collapse_double_lines_distance = regularization_distance;
    Optimize_polyline_by_plane_patches<helper::CGAL_types::K> placement_optimizer;
    placement_optimizer.out_filepath = out_filepath;
    placement_optimizer.maximum_iteration = maximum_iteration;
    placement_optimizer.min_energy_delta = min_energy_delta;
    placement_optimizer.collapse_distance = collapse_distance;
    placement_optimizer.smooth_factor = smooth_factor;
    placement_optimizer.refine_surfel_distance = refine_surfel_distance;
    placement_optimizer.plane_detection_max_distance = plane_detection_max_distance;
    placement_optimizer.initialize_step_size_factor(start_step_size, end_step_size, maximum_iteration);
    return detect_sharp_edges_in_image(image, sharpness_measure, surface_selector, polyline_extractor, regularization_functor, placement_optimizer, out_filepath);
}

#endif // SHARP_FEATURE_EXTRACTION_H