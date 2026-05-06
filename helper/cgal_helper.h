#ifndef CGAL_HELPER_H
#define CGAL_HELPER_H

// Kernel
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
// Surface mesh
#include <CGAL/Surface_mesh.h>
#include <CGAL/polygon_mesh_processing.h>
#include <CGAL/Polygon_mesh_processing/detect_features.h>
#include <CGAL/Polygon_mesh_processing/IO/polygon_mesh_io.h>
// KD tree
#include <CGAL/Search_traits_3.h>
#include <CGAL/Search_traits_adapter.h>
#include <CGAL/Orthogonal_k_neighbor_search.h>

#include "../polyline/Polylines.h"
#include "../polyline/polyline_process.h"

namespace helper
{
    namespace CGAL_types
    {
        typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
        typedef K::Point_3                                          Point;
        typedef K::Vector_3                                         Vector;
        typedef CGAL::Surface_mesh<Point>                           Mesh;
        typedef CGAL::_image                                        Image;
    }

    Polylines<CGAL_types::Point> extract_sharp_edges_from_mesh(
        CGAL_types::Mesh &mesh,
        const float &angle_in_degree
    )
    {
        namespace PMP = CGAL::Polygon_mesh_processing;
        using Polylines = Polylines<CGAL_types::Point>;

        // get edges
        boost::property_map<CGAL_types::Mesh, CGAL::edge_is_feature_t>::type edge_is_feature_map
            = get(CGAL::edge_is_feature, mesh);
        PMP::detect_sharp_edges(mesh, angle_in_degree, edge_is_feature_map);
        // make graph from edges
        std::vector<CGAL_types::Point> points;
        std::vector<std::vector<std::size_t>> point_neighbors;
        std::map<boost::graph_traits<CGAL_types::Mesh>::vertex_descriptor, std::size_t> vertex_to_point;
        for(boost::graph_traits<CGAL_types::Mesh>::edge_descriptor e : edges(mesh))
        {
            if(!get(edge_is_feature_map, e))
                continue;
            // edge is sharp : get (or insert) point id and create link
            boost::graph_traits<CGAL_types::Mesh>::vertex_descriptor v0 = mesh.vertex(e, 0);
            boost::graph_traits<CGAL_types::Mesh>::vertex_descriptor v1 = mesh.vertex(e, 1);
            const auto& [v0_index_it, v0_inserted] = vertex_to_point.try_emplace(v0, points.size());
            if (v0_inserted)
            {
                points.push_back(mesh.point(v0));
                point_neighbors.resize(point_neighbors.size()+1);
            }
            const auto& [v1_index_it, v1_inserted] = vertex_to_point.try_emplace(v1, points.size());
            if (v1_inserted)
            {
                points.push_back(mesh.point(v1));
                point_neighbors.resize(point_neighbors.size()+1);
            }
            point_neighbors[v0_index_it->second].push_back(v1_index_it->second);
            point_neighbors[v1_index_it->second].push_back(v0_index_it->second);
        }
        // convert to polylines
        Polylines polylines;
        polylines.points = points;
        polylines.lines = Polylines_process::make_lines_from_graph(point_neighbors);
        return polylines;
    }

    // ray origin is outside the box. Returns null optional if the ray does hit the box exactly on corners (float precision)
    template <typename K>
    std::optional<typename K::Point_2> get_ray_bbox_intersection_2d(
        const typename K::Ray_2& ray, const typename K::Point_2& bbox_min, const typename K::Point_2& bbox_max)
    {
        /*
        *  Find which segment to check based on the position of the ray source
        *
        *     UP     |      |            UP
        *   & LEFT   |  UP bbox_max   & RIGHT
        *  _________  ______X _________
        *            |      |
        *     LEFT   |      |   RIGHT
        *  _________ X______| _________
        *        bbox_min   |
        *     DOWN   | DOWN |      DOWN
        *    & LEFT  |      |    & RIGHT
        *
        */
        using Point_2 = typename K::Point_2;
        using Segment_2 = typename K::Segment_2;

        const Point_2& ray_source = ray.source();
        std::vector<Segment_2> segments_to_check;
        segments_to_check.reserve(2);
        // LEFT or RIGHT
        if (ray_source[0] >= bbox_max[0])
            segments_to_check.push_back(Segment_2(bbox_max, Point_2(bbox_max[0], bbox_min[1]))); // RIGHT
        else if (ray_source[0] <= bbox_min[0])
            segments_to_check.push_back(Segment_2(bbox_min, Point_2(bbox_min[0], bbox_max[1]))); // LEFT
        // UP or DOWN
        if (ray_source[1] >= bbox_max[1])
            segments_to_check.push_back(Segment_2(bbox_max, Point_2(bbox_min[0], bbox_max[1]))); // UP
        else if (ray_source[1] <= bbox_min[1])
            segments_to_check.push_back(Segment_2(bbox_min, Point_2(bbox_max[0], bbox_min[1]))); // DOWN
        std::optional<Point_2> min_inter_point;
        typename K::FT min_dist = std::numeric_limits<typename K::FT>::max();
        for (const Segment_2& segment_to_check : segments_to_check)
        {
            const auto& intersect_result = CGAL::intersection(ray, segment_to_check);
            if (intersect_result)
            {
                if (const Point_2* inter_point = std::get_if<Point_2>(&*intersect_result))
                {
                    if (CGAL::squared_distance(*inter_point, ray_source) < min_dist)
                    {
                        min_inter_point = *inter_point;
                    }
                }
            }
        }
        return min_inter_point;
    }

    template <typename K>
    typename K::Point_3 project_at_plane_intersection(
        const std::vector<std::pair<typename K::Point_3, typename K::Vector_3>>& planes,
        const typename K::Point_3& point,
        const std::size_t nb_plane_average_iteration = 64)
    {
        using FT = typename K::FT;
        using Point_3 = typename K::Point_3;
        using Vector_3 = typename K::Vector_3;
        using Plane_3 = typename K::Plane_3;
        if (planes.size() <= 1)
        {
            // 1 plane : project on it
            const auto& [plane_center, plane_normal] = planes[0];
            return Plane_3(plane_center, plane_normal).projection(point);
        }
        else
        {
            // 2 or more planes : reapeat the average projection on each plane
            std::vector<Plane_3> planes_cgal;
            planes_cgal.reserve(planes.size());
            for (const auto& [plane_center, plane_normal] : planes)
            {
                planes_cgal.emplace_back(plane_center, plane_normal);
            }
            Point_3 cur_point = point;
            // repeatedly move on average of projection
            Vector_3 average_plane_projection_point;
            const FT inv_planes_size = 1.0 / static_cast<FT>(planes.size());
            for (std::size_t i = 0; i < nb_plane_average_iteration; i++)
            {
                average_plane_projection_point = Vector_3(0,0,0);
                for (const Plane_3& plane : planes_cgal)
                {
                    average_plane_projection_point += Vector_3(Point_3(0,0,0), plane.projection(cur_point));
                }
                cur_point = Point_3(0,0,0) + average_plane_projection_point * inv_planes_size;
            }
            return cur_point;
        }
    }

    // The first of tuple in Element_type is the point of type K::Point_3
    template <typename Element_type, typename K>
    struct KD_tree
    {
        typedef Element_type                                    Tree_point;
        typedef CGAL::Search_traits_3<K>                        Traits_base;
        typedef CGAL::Search_traits_adapter<Tree_point,
                                            CGAL::Nth_of_tuple_property_map<0, Tree_point>,
                                            Traits_base>                                          Traits;
        typedef CGAL::Orthogonal_k_neighbor_search<Traits>      K_neighbor_search;
        typedef typename K_neighbor_search::Tree                Tree;
        typedef typename K_neighbor_search::Distance            Distance_func;
        typedef typename K::FT                                  Distance_type;

        typedef typename std::vector<std::pair<const Element_type*, Distance_type>>::const_iterator Iterator_type;

        Tree tree;

        std::vector<Element_type> storage;
        std::vector<std::pair<const Element_type*, Distance_type>> closest_points;
        KD_tree() { }
        KD_tree(const KD_tree& other) { } // do nothing
        void reserve(const std::size_t& number)
        {
            tree.reserve(number);
        }
        void push(const Element_type& element)
        {
            tree.insert(element);
        }
        void push(const Element_type* element)
        {
            tree.insert(*element);
        }
        std::array<Iterator_type, 2> get_closest(const Element_type* point, const std::size_t& knn = 1)
        {
            closest_points.clear();
            storage.clear();
            closest_points.reserve(knn);
            storage.reserve(knn);
            K_neighbor_search search(
                tree, point->template get<0>(), knn,
                /*epsilon=*/K::FT(0.0), /*nearest=*/true, Distance_func(), /*sorted=*/true
                );
            for(typename K_neighbor_search::iterator it = search.begin(); it != search.end(); it++)
            {
                const auto& [tree_point, distance] = *it;
                storage.push_back(tree_point);
                closest_points.emplace_back(&(storage.back()), distance);
            }
            return {closest_points.begin(), closest_points.end()};
        }
        std::array<Iterator_type, 2> get_closest(const typename K::Point_3& point, const std::size_t& knn = 1)
        {
            Tree_point search_point(point, 0);
            return get_closest(&search_point, knn);
        }
    };

    CGAL::Bbox_3 get_bbox_from_mesh(const CGAL_types::Mesh& mesh)
    {
        const auto& first_v = *(mesh.points().begin());
        CGAL::Bbox_3 bbox(first_v.x(), first_v.y(), first_v.z(), first_v.x(), first_v.y(), first_v.z());
        for (const auto& v : mesh.points())
        {
            bbox += CGAL::Bbox_3(v.x(), v.y(), v.z(), v.x(), v.y(), v.z());
        }
        return bbox;
    }

    CGAL::Bbox_3 get_normalized_bbox_with_same_ratio(const CGAL::Bbox_3& bbox)
    {
        const double &xrange = bbox.xmax()-bbox.xmin();
        const double &yrange = bbox.ymax()-bbox.ymin();
        const double &zrange = bbox.zmax()-bbox.zmin();
        double inv_maxrange = 1.0 / std::max(xrange, std::max(yrange, zrange));
        return CGAL::Bbox_3(0,0,0, xrange*inv_maxrange, yrange*inv_maxrange, zrange*inv_maxrange);
    }

    template <typename Iterator>
    void scale_points_to_range(
        const Iterator &begin_iterator,
        const Iterator &end_iterator,
        const CGAL::Bbox_3& poly_input_range,
        const CGAL::Bbox_3& output_range)
    {
        using Point_3 = CGAL_types::Point;
        using Vector_3 = CGAL_types::Vector;
        const Vector_3 input_min_point(poly_input_range.xmin(), poly_input_range.ymin(), poly_input_range.zmin());
        const Vector_3 output_min_point(output_range.xmin(), output_range.ymin(), output_range.zmin());
        const Point_3 range_factor(
            (output_range.xmax() - output_range.xmin()) / (poly_input_range.xmax() - poly_input_range.xmin()),
            (output_range.ymax() - output_range.ymin()) / (poly_input_range.ymax() - poly_input_range.ymin()),
            (output_range.zmax() - output_range.zmin()) / (poly_input_range.zmax() - poly_input_range.zmin())
        );
        for (Iterator iterator = begin_iterator; iterator != end_iterator; ++iterator)
        {
            // point = (point - input_min_point) * range_factor + output_min_point
            *iterator = helper::vec_component_mult(*iterator - input_min_point, range_factor)
                + output_min_point;
        }
    }

    void scale_mesh_to_range(
        CGAL_types::Mesh& mesh,
        const CGAL::Bbox_3& poly_input_range,
        const CGAL::Bbox_3& output_range)
    {
        scale_points_to_range(mesh.points().begin(), mesh.points().end(), poly_input_range, output_range);
    }

    template <typename Point>
    void scale_features_to_range(
        Polylines<Point>& features,
        const CGAL::Bbox_3& poly_input_range,
        const CGAL::Bbox_3& output_range)
    {
        scale_points_to_range(features.points.begin(), features.points.end(), poly_input_range, output_range);
    }

} // namespace helper

#endif // CGAL_HELPER_H