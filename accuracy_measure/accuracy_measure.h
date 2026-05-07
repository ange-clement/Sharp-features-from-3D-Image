#ifndef ACCURACY_MEASURE_H
#define ACCURACY_MEASURE_H

#include <optional>

// C3t3
#include <CGAL/Labeled_mesh_domain_3.h>
#include <CGAL/Mesh_domain_with_polyline_features_3.h>
#include <CGAL/Mesh_triangulation_3.h>
#include <CGAL/Mesh_criteria_3.h>
#include <CGAL/Mesh_complex_3_in_triangulation_3.h>
#include <CGAL/make_mesh_3.h>
// C3t3 to Surface mesh
#include <CGAL/facets_in_complex_3_to_triangle_mesh.h>
// AABB tree
#include <CGAL/AABB_tree.h>

#include "../polyline/polyline_process.h"
#include "../helper/cgal_helper.h"
#include "../helper/geometry_helper.h"

namespace Accuracy_measure
{
    namespace internal
    {
        // samples mesh1 and project to mesh2
        double get_surface_mesh_error_and_output(
            const helper::CGAL_types::Mesh& mesh1, const helper::CGAL_types::Mesh& mesh2, const double& nb_points_per_area_unit,
            const std::string &output_samples_ply,
            double &min_value_out, double &max_value_out,
            const std::optional<double>& min_value_in = std::nullopt, const std::optional<double> &max_value_in = std::nullopt)
        {
            typedef helper::CGAL_types::K K;
            typedef helper::CGAL_types::Point Point_3;
            typedef helper::CGAL_types::Mesh Mesh;
            namespace PMP = CGAL::Polygon_mesh_processing;
            std::vector<Point_3> sample_points;
            PMP::sample_triangle_mesh(mesh1, std::back_inserter(sample_points), CGAL::parameters::number_of_points_per_area_unit(nb_points_per_area_unit));

            CGAL::spatial_sort(sample_points.begin(), sample_points.end());

            typedef typename CGAL::GetVertexPointMap<Mesh, CGAL::parameters::Default_named_parameters>::const_type  VPM;
            typedef CGAL::AABB_face_graph_triangle_primitive<Mesh, VPM> Primitive;
            typedef CGAL::AABB_traits_3<K, Primitive> Tree_traits;
            typedef CGAL::AABB_tree<Tree_traits> Tree;

            VPM vpm = get_const_property_map(CGAL::vertex_point, mesh2);

            Tree_traits tgt/*(gt)*/;
            Tree tree(tgt);
            tree.insert(faces(mesh2).first, faces(mesh2).second, mesh2, vpm);

            Point_3 hint = get(vpm, *vertices(mesh2).first);

            double sq_hdist = 0;
            typename K::Compute_squared_distance_3 squared_distance = K().compute_squared_distance_3_object();

            std::vector<double> display_values; // used iff outputSamples is true
            display_values.reserve(sample_points.size());

            for(const Point_3& pt : sample_points)
            {
                hint = tree.closest_point(pt, hint);
                double sq_d = CGAL::squared_distance(hint, pt);
                sq_hdist += sq_d;

                display_values.push_back(CGAL::sqrt(sq_d));
            }

            {
                auto [min_it, max_it] = std::minmax_element(display_values.begin(), display_values.end());
                min_value_out = *min_it;
                max_value_out = *max_it;
                if (min_value_in.has_value() && max_value_in.has_value())
                    helper::internal::write_point_set_values_min_max(output_samples_ply, sample_points, display_values, min_value_in.value(), max_value_in.value(), false);
                else
                    helper::write_point_set_values(output_samples_ply, sample_points, display_values);
            }

            return sq_hdist / static_cast<double>(sample_points.size());
        }

        // samples mesh1 and project to mesh2
        double get_surface_mesh_error(
            const helper::CGAL_types::Mesh& mesh1, const helper::CGAL_types::Mesh& mesh2, const double& nb_points_per_area_unit,
            double &min_value_out, double &max_value_out)
        {
            typedef helper::CGAL_types::K K;
            typedef helper::CGAL_types::Point Point_3;
            typedef helper::CGAL_types::Mesh Mesh;
            namespace PMP = CGAL::Polygon_mesh_processing;
            std::vector<Point_3> sample_points;
            PMP::sample_triangle_mesh(mesh1, std::back_inserter(sample_points), CGAL::parameters::number_of_points_per_area_unit(nb_points_per_area_unit));

            CGAL::spatial_sort(sample_points.begin(), sample_points.end());

            typedef typename CGAL::GetVertexPointMap<Mesh, CGAL::parameters::Default_named_parameters>::const_type  VPM;
            typedef CGAL::AABB_face_graph_triangle_primitive<Mesh, VPM> Primitive;
            typedef CGAL::AABB_traits_3<K, Primitive> Tree_traits;
            typedef CGAL::AABB_tree<Tree_traits> Tree;

            VPM vpm = get_const_property_map(CGAL::vertex_point, mesh2);

            Tree_traits tgt/*(gt)*/;
            Tree tree(tgt);
            tree.insert(faces(mesh2).first, faces(mesh2).second, mesh2, vpm);

            Point_3 hint = get(vpm, *vertices(mesh2).first);

            double sq_hdist = 0;
            typename K::Compute_squared_distance_3 squared_distance = K().compute_squared_distance_3_object();

            for(const Point_3& pt : sample_points)
            {
                hint = tree.closest_point(pt, hint);
                double sq_d = CGAL::squared_distance(hint, pt);
                sq_hdist += sq_d;
            }

            return sq_hdist / static_cast<double>(sample_points.size());
        }

        template <typename Vector, typename Point>
        Point get_random_sample_in_edge(
            const Point& p0, const Point& p1)
        {
            const double sample_position_in_edge = (double)std::rand() / RAND_MAX;
            return p0 + sample_position_in_edge * Vector(p0, p1);
        }

        template <typename Vector, typename Point>
        std::vector<Point> uniform_sample_on_line(
            const Polylines<Point>& poly,
            const double& poly_length,
            const double& nb_sample)
        {
            std::vector<Point> sample_points;
            for (const std::vector<std::size_t>& line : poly.lines)
            {
                if (line.size() == 0)
                    continue;
                for (std::size_t element_id = 0, nb_p = line.size()-1; element_id < nb_p; element_id++)
                {
                    const Point& p0 = poly.points[line[element_id]];
                    const Point& p1 = poly.points[line[element_id+1]];
                    const double& edge_length = CGAL::sqrt(CGAL::squared_distance(p0, p1));
                    double number_of_samples_this_edge = nb_sample * edge_length / poly_length;
                    while (number_of_samples_this_edge > 1.0)
                    {
                        sample_points.push_back(get_random_sample_in_edge<Vector>(p0, p1));
                        number_of_samples_this_edge -= 1.0;
                    }
                    if (number_of_samples_this_edge > 0.0)
                    {
                        const double random_check = (double)std::rand() / RAND_MAX;
                        if (random_check < number_of_samples_this_edge)
                            sample_points.push_back(get_random_sample_in_edge<Vector>(p0, p1));
                    }
                }
            }
            return sample_points;
        }

    } /* namespace internal */

    template <typename Vector, typename Point>
    std::size_t number_of_lines_in_first_not_in_second(
        const Polylines<Point>& poly_first,
        const Polylines<Point>& poly_second,
        const double& distance)
    {
        std::size_t number_of_lines_not_in_second = 0;
        double sq_dist = distance*distance;
        for (const std::vector<std::size_t>& line_second : poly_second.lines)
        {
            for (const std::size_t& point_second_id : line_second)
            {
                const Point& point_second = poly_second.points[point_second_id];
                const Point& projected_point = Polylines_process::project_point_to_polylines<Vector>(poly_first, point_second);
                const double& projection_sq_len = CGAL::squared_distance(point_second, projected_point);
                if (projection_sq_len > sq_dist)
                {
                    number_of_lines_not_in_second++;
                    break;
                }
            }
        }
        return number_of_lines_not_in_second;
    }

    helper::CGAL_types::Mesh make_mesh_from_image(
        const CGAL::Image_3 &image,
        const Polylines<helper::CGAL_types::Point> &polylines = Polylines<helper::CGAL_types::Point>(),
        const double &sizing = 0.005)
    {
        typedef CGAL::Labeled_mesh_domain_3<helper::CGAL_types::K> Image_domain;
        typedef CGAL::Mesh_domain_with_polyline_features_3<Image_domain> Mesh_domain;
    #ifdef CGAL_CONCURRENT_MESH_3
        typedef CGAL::Parallel_tag Concurrency_tag;
    #else
        typedef CGAL::Sequential_tag Concurrency_tag;
    #endif
        typedef CGAL::Mesh_triangulation_3<Mesh_domain,CGAL::Default,Concurrency_tag>::type Tr;
        typedef CGAL::Mesh_criteria_3<Tr> Mesh_criteria;
        typedef CGAL::Mesh_complex_3_in_triangulation_3<Tr> C3t3;

        typedef helper::CGAL_types::Point Point_3;
        typedef C3t3::Vertex_handle Vertex_handle;
        typedef C3t3::Cell_handle Cell_handle;
        typedef C3t3::Facet Facet;
        typedef helper::CGAL_types::Mesh Mesh;
        typedef Mesh::Vertex_index Vertex_index;
        namespace PMP = CGAL::Polygon_mesh_processing;

        std::vector<std::vector<helper::CGAL_types::Point>> &polylines_for_mesh = Polylines_process::make_polyline_for_mesh_3(polylines);

        Mesh_domain domain = Mesh_domain::create_labeled_image_mesh_domain(image, CGAL::parameters::input_features(polylines_for_mesh));
        CGAL::Bbox_3 bbox = domain.bbox();
        const double diag = CGAL::sqrt(CGAL::square(bbox.xmax() - bbox.xmin()) +
                                       CGAL::square(bbox.ymax() - bbox.ymin()) +
                                       CGAL::square(bbox.zmax() - bbox.zmin()));
        const double default_sizing = diag * sizing;
        Mesh_criteria criteria(
            CGAL::parameters::edge_size(default_sizing)
                .edge_min_size(0.25 * default_sizing)
                // .edge_distance(default_approx)
                .facet_angle(25.)
                .facet_size(default_sizing)
                // .facet_distance(default_approx)
                // .facet_min_size(0.1 * default_sizing)
                // .cell_radius_edge_ratio(3.0)
                // .cell_size(default_sizing)
                // .cell_min_size(0.1 * default_sizing)
            );
        std::cout << "make mesh with polyline with " << polylines.points.size() << " elements" << std::endl;
        C3t3 c3t3 = CGAL::make_mesh_3<C3t3>(domain, criteria);
        // make mesh from c3t3
        Mesh sm;
        CGAL::facets_in_complex_3_to_triangle_mesh(c3t3, sm);
        return sm;
    }

    double get_mesh_error_and_output(
        const helper::CGAL_types::Mesh& mesh_to_compare,
        const helper::CGAL_types::Mesh& mesh_target,
        const std::string &output_samples,
        double &min_value_out, double &max_value_out,
        const std::optional<double>& min_value_in = std::nullopt, const std::optional<double> &max_value_in = std::nullopt,
        const double &number_of_samples = 10000)
    {
        const CGAL::Bbox_3& gt_bbox = helper::get_bbox_from_mesh(mesh_target);
        const double diag = CGAL::sqrt(gt_bbox.x_span()*gt_bbox.x_span() + gt_bbox.y_span()*gt_bbox.y_span() + gt_bbox.z_span()*gt_bbox.z_span());
        const double nb_points_per_area_unit = number_of_samples / diag;

        double mesh_error_1 = internal::get_surface_mesh_error_and_output(
            mesh_to_compare, mesh_target, nb_points_per_area_unit, output_samples + "_from_compare_to_target.ply",
            min_value_out, max_value_out, min_value_in, max_value_in);
        double mesh_error_2 = internal::get_surface_mesh_error_and_output(
            mesh_to_compare, mesh_target, nb_points_per_area_unit, output_samples + "_from_target_to_compare.ply",
            min_value_out, max_value_out, min_value_in, max_value_in);
        return mesh_error_1 + mesh_error_2;
    }

    double get_mesh_error(
        const helper::CGAL_types::Mesh& mesh_to_compare,
        const helper::CGAL_types::Mesh& mesh_target,
        double &min_value_out, double &max_value_out,
        const double &number_of_samples = 100000)
    {
        const CGAL::Bbox_3& gt_bbox = helper::get_bbox_from_mesh(mesh_target);
        const double diag = CGAL::sqrt(gt_bbox.x_span()*gt_bbox.x_span() + gt_bbox.y_span()*gt_bbox.y_span() + gt_bbox.z_span()*gt_bbox.z_span());
        const double nb_points_per_area_unit = number_of_samples / diag;

        double mesh_error_1 = internal::get_surface_mesh_error(
            mesh_to_compare, mesh_target, nb_points_per_area_unit,
            min_value_out, max_value_out);
        double mesh_error_2 = internal::get_surface_mesh_error(
            mesh_to_compare, mesh_target, nb_points_per_area_unit,
            min_value_out, max_value_out);
        return mesh_error_1 + mesh_error_2;
    }

    template <typename Vector, typename Point>
    double get_polyline_error(
        const Polylines<Point>& gt_poly,
        const Polylines<Point>& method_poly,
        double &min_error, double &max_error)
    {
        // MSE : sum(d(method, g_poly)^2) / N but asymetric, so do it on both sides
        // error : sum(d(method, gt_poly)^2) / N1 + sum(d(gt, method_poly)^2) / N2
        // d(method, gt_poly) is the distance of a method point to the gt_poly
        // make random uniform sample
        std::srand(0);
        const double& gt_len = Polylines_process::get_polyline_length(gt_poly);
        const double& method_len = Polylines_process::get_polyline_length(method_poly);
        const double total_len = gt_len + method_len;
        const std::size_t total_nb_points = gt_poly.points.size() + method_poly.points.size();
        std::vector<Point> gt_sample = internal::uniform_sample_on_line<Vector>(gt_poly, total_len, total_nb_points);
        std::vector<Point> method_sample = internal::uniform_sample_on_line<Vector>(method_poly, total_len, total_nb_points);
        double error = 0;
        std::size_t n = 0;
        min_error = std::numeric_limits<double>::max();
        max_error = std::numeric_limits<double>::lowest();
        for (const Point& gt_point : gt_sample)
        {
            const Point& projected_point = Polylines_process::project_point_to_polylines<Vector>(method_poly, gt_point);
            const double& projection_sq_len = CGAL::squared_distance(gt_point, projected_point);
            if (projection_sq_len < min_error)
                min_error = projection_sq_len;
            if (projection_sq_len > max_error)
                max_error = projection_sq_len;
            error += projection_sq_len;
            n++;
        }
        for (const Point& method_point : method_sample)
        {
            const Point& projected_point = Polylines_process::project_point_to_polylines<Vector>(gt_poly, method_point);
            const double& projection_sq_len = CGAL::squared_distance(method_point, projected_point);
            if (projection_sq_len < min_error)
                min_error = projection_sq_len;
            if (projection_sq_len > max_error)
                max_error = projection_sq_len;
            error += projection_sq_len;
            n++;
        }
        error /= static_cast<double>(n);
        return error;
    }

} /* namespace Accuracy_measure */

#endif