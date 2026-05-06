#ifndef SURFACE_PROCES_H
#define SURFACE_PROCES_H

#include <CGAL/Image_3.h>

#include "surface_graph.h"
#include "../helper/cgal_helper.h"
#include "../helper/image_helper.h"
#include "../element_graph/element_graph_process.h"
#include "../polyline/polyline_process.h"

namespace Surface_process
{

    namespace internal
    {

        // Each surface has four points
        // the points have the same index as voxels (lowest corner)
        template <typename Vector>
        void update_surface_graph(const std::array<int, 3> &u, const std::array<int, 3> &v, const typename Surface_graph<Vector>::Surfel &surfel,
                                  const std::size_t &x, const std::size_t &y, const std::size_t &z, const std::size_t &xdim, const std::size_t &ydim,
                                  std::map<std::size_t, std::set<typename Surface_graph<Vector>::Surfel *>> &points_to_contact_surfaces,
                                  Surface_graph<Vector> &surfaces_graph, const bool &pointel_neighboors = true)
        {
            typedef typename Surface_graph<Vector> Surface_graph;
            typedef typename Surface_graph::Surfel Surfel;

            // get surfel
            Surfel *surfel_ptr = surfaces_graph.add_surfel(surfel);
            std::set<Surfel *> &surface_neighbors = surfel_ptr->neighbors;

            // get points and their neighbors
            std::size_t p0 = helper::coords_to_voxel_index(x, y, z, xdim, ydim);
            std::size_t pu = helper::coords_to_voxel_index(x + u[0], y + u[1], z + u[2], xdim, ydim);
            std::size_t pv = helper::coords_to_voxel_index(x + v[0], y + v[1], z + v[2], xdim, ydim);
            std::size_t puv = helper::coords_to_voxel_index(x + u[0] + v[0], y + u[1] + v[1], z + u[2] + v[2], xdim, ydim);
            std::set<Surfel *> &p0_surface_neighbors = points_to_contact_surfaces[p0];
            std::set<Surfel *> &pu_surface_neighbors = points_to_contact_surfaces[pu];
            std::set<Surfel *> &pv_surface_neighbors = points_to_contact_surfaces[pv];
            std::set<Surfel *> &puv_surface_neighbors = points_to_contact_surfaces[puv];

            // update surfel neighbors
            if (pointel_neighboors)
            {
                // surfels are neighbors on points
                // Union of point neighbors surfels : p0 U pu U pv U puv
                surface_neighbors.insert(p0_surface_neighbors.begin(), p0_surface_neighbors.end());
                surface_neighbors.insert(pu_surface_neighbors.begin(), pu_surface_neighbors.end());
                surface_neighbors.insert(pv_surface_neighbors.begin(), pv_surface_neighbors.end());
                surface_neighbors.insert(puv_surface_neighbors.begin(), puv_surface_neighbors.end());
            }
            else
            {
                // surfels are neighbors on edges
                // Union of Intersections of edge points : (p0 ∩ pu) U (po ∩ pv) U (pu ∩ puv) U (pv ∩ puv)
                std::set_intersection(p0_surface_neighbors.begin(), p0_surface_neighbors.end(), pu_surface_neighbors.begin(), pu_surface_neighbors.end(), std::inserter(surface_neighbors, surface_neighbors.begin()));
                std::set_intersection(p0_surface_neighbors.begin(), p0_surface_neighbors.end(), pv_surface_neighbors.begin(), pv_surface_neighbors.end(), std::inserter(surface_neighbors, surface_neighbors.begin()));
                std::set_intersection(pu_surface_neighbors.begin(), pu_surface_neighbors.end(), puv_surface_neighbors.begin(), puv_surface_neighbors.end(), std::inserter(surface_neighbors, surface_neighbors.begin()));
                std::set_intersection(pv_surface_neighbors.begin(), pv_surface_neighbors.end(), puv_surface_neighbors.begin(), puv_surface_neighbors.end(), std::inserter(surface_neighbors, surface_neighbors.begin()));
            }

            // add surfel to neighbor's neighbors
            for (Surfel *neighbor : surface_neighbors)
            {
                neighbor->neighbors.insert(surfel_ptr);
            }

            // update points neighbors
            p0_surface_neighbors.insert(surfel_ptr);
            pu_surface_neighbors.insert(surfel_ptr);
            pv_surface_neighbors.insert(surfel_ptr);
            puv_surface_neighbors.insert(surfel_ptr);
        }

        template <typename Vector>
        struct Make_surface_graph_from_image_struct
        {
            template <typename Word>
            Surface_graph<Vector> operator()(const CGAL::Image_3 &image, const bool &pointel_neighboors = true) const
            {
                using Surface_graph = Surface_graph<Vector>;
                using Surfel = typename Surface_graph::Surfel;

                const helper::CGAL_types::Image *im = image.image();
                const std::size_t &xdim = im->xdim;
                const std::size_t &ydim = im->ydim;
                const std::size_t &zdim = im->zdim;

                Surface_graph surfaces_graph(xdim, ydim, zdim);

                std::map<std::size_t, std::set<Surfel *>> points_to_contact_surfaces;

                for (std::size_t z = 0; z < zdim; z++)
                    for (std::size_t y = 0; y < ydim; y++)
                        for (std::size_t x = 0; x < xdim; x++)
                        {
                            const std::size_t voxel_index = helper::coords_to_voxel_index(x, y, z, xdim, ydim);
                            const Word &w = CGAL::IMAGEIO::static_evaluate<Word>(im, voxel_index) == 0 ? 0 : 1;
                            const Word &wx = helper::safe_get_image_value<Word>(im, x + 1, y, z, xdim, ydim, zdim, static_cast<Word>(0)) == 0 ? 0 : 1;
                            const Word &wy = helper::safe_get_image_value<Word>(im, x, y + 1, z, xdim, ydim, zdim, static_cast<Word>(0)) == 0 ? 0 : 1;
                            const Word &wz = helper::safe_get_image_value<Word>(im, x, y, z + 1, xdim, ydim, zdim, static_cast<Word>(0)) == 0 ? 0 : 1;

                            if (w != wx)
                            {
                                internal::update_surface_graph({0, 1, 0}, {0, 0, 1}, Surfel(voxel_index, Surfel::X, (w < wx)),
                                                               x + 1, y, z, xdim, ydim, points_to_contact_surfaces, surfaces_graph, pointel_neighboors);
                            }
                            if (w != wy)
                            {
                                internal::update_surface_graph({1, 0, 0}, {0, 0, 1}, Surfel(voxel_index, Surfel::Y, (w < wy)),
                                                               x, y + 1, z, xdim, ydim, points_to_contact_surfaces, surfaces_graph, pointel_neighboors);
                            }
                            if (w != wz)
                            {
                                internal::update_surface_graph({1, 0, 0}, {0, 1, 0}, Surfel(voxel_index, Surfel::Z, (w < wz)),
                                                               x, y, z + 1, xdim, ydim, points_to_contact_surfaces, surfaces_graph, pointel_neighboors);
                            }
                        }
                return surfaces_graph;
            }
        };

        struct All_surfels
        {
            template <typename Vector>
            bool operator()(const Surface_graph<Vector> &graph, const typename Surface_graph<Vector>::Surfel &surfel) const
            {
                return true;
            }
        };

        template <typename Vector>
        struct Surfels_in_set
        {
            using Surface_graph = typename Surface_graph<Vector>;
            using Surfel = typename Surface_graph::Surfel;
            const std::set<const Surfel *> &surfels_to_display;

            Surfels_in_set(const std::set<const Surfel *> &surfels_to_display)
                : surfels_to_display(surfels_to_display)
            {
            }

            bool operator()(const Surface_graph &graph, const Surfel &surfel) const
            {
                return surfels_to_display.count(&surfel) > 0;
            }
        };

        template <typename K, typename Vector, typename SurfelPredicate = All_surfels>
        std::vector<typename K::Point_3> graph_to_mesh(
            const Surface_graph<Vector> &graph,
            const SurfelPredicate &surface_predicate = SurfelPredicate())
        {
            using Surface_graph = Surface_graph<Vector>;
            using Surfel = typename Surface_graph::Surfel;
            using Point_3 = typename K::Point_3;

            std::vector<Point_3> points;
            const std::size_t &xdim = graph.xdim;
            const std::size_t &ydim = graph.ydim;
            for (const auto &[voxel_id, voxel_surfels] : graph.surface)
                for (const Surfel &surfel : voxel_surfels)
                {
                    if (!surface_predicate(graph, surfel))
                        continue;
                    std::array<std::array<int, 3>, 4> surfel_points = surfel.get_pointels(xdim, ydim);
                    for (const std::array<int, 3> &surfel_point : surfel_points)
                    {
                        points.push_back(Point_3(surfel_point[0], surfel_point[1], surfel_point[2]));
                    }
                }
            return points;
        }

        template <typename Point>
        void write_surface_mesh_ply(const std::vector<Point> &points, const std::string &filename)
        {
            std::ofstream output_mesh(filename, std::ios_base::binary);
            if (output_mesh.bad())
            {
                output_mesh.close();
                return;
            }
            output_mesh << "ply\nformat ascii 1.0\n"
                        << "element vertex " << points.size() << "\n"
                        << "property double x\nproperty double y\nproperty double z\n"
                        << "property int id\n"
                        << "element face " << points.size() / 4 << "\n"
                        << "property list uchar int vertex_indices\n"
                        << "end_header\n";
            for (std::size_t i = 0, size = points.size(); i < size; i++)
            {
                const Point &p = points[i];
                output_mesh << p[0] << " " << p[1] << " " << p[2] << " "
                            << i << "\n";
            }
            for (std::size_t i = 0, size = points.size(); i < size; i += 4)
            {
                output_mesh << "4 ";
                for (int j = 0; j < 4; j++)
                {
                    output_mesh << i + j << " ";
                }
                output_mesh << "\n";
            }
            output_mesh.close();
        }

        template <typename K, typename Vector>
        void write_surface_sharpness_ply(
            const Surface_graph<Vector> &graph, const std::string &filename,
            const std::optional<double> &min_value = std::nullopt,
            const std::optional<double> &max_value = std::nullopt)
        {
            using Point = typename K::Point_3;
            using Surfel = typename Surface_graph<Vector>::Surfel;
            std::vector<Point> mesh = graph_to_mesh<K, Vector>(graph);
            std::vector<double> surface_values;
            for (const auto &[voxel_id, voxel_surfels] : graph.surface)
                for (const Surfel &surfel : voxel_surfels)
                {
                    surface_values.push_back(1.0 - surfel.sharpness);
                }
            double min, max;
            if (!min_value || !max_value)
            {
                const auto &[computed_min_it, computed_max_it] = std::minmax_element(surface_values.begin(), surface_values.end());
                min = min_value.value_or(*computed_min_it);
                max = max_value.value_or(*computed_max_it);
            }
            else
            {
                min = min_value.value();
                max = max_value.value();
            }
            helper::write_mesh_with_heat_map_as_ply<Point>(filename, mesh, surface_values, true, true, min, max);
        }

        template <typename K, typename Vector>
        void write_surface_normals_ply(const Surface_graph<Vector> &graph,
                                       const std::string &filename)
        {
            using Point = typename K::Point_3;
            using Surfel = typename Surface_graph<Vector>::Surfel;
            std::vector<Point> mesh = graph_to_mesh<K, Vector>(graph);
            std::vector<double> colors;
            for (const auto &[voxel_id, voxel_surfels] : graph.surface)
                for (const Surfel &surfel : voxel_surfels)
                {
                    Vector n = surfel.normal;
                    helper::normalize_if_possible(n);
                    colors.push_back(n[0] < 0 ? -n[0] * 0.5 : n[0]);
                    colors.push_back(n[1] < 0 ? -n[1] * 0.5 : n[1]);
                    colors.push_back(n[2] < 0 ? -n[2] * 0.5 : n[2]);
                }
            helper::write_mesh_with_colors_as_ply<Point>(filename, mesh, colors);
        }

        template <typename Vector>
        struct Surfel_visitor
        {
            typedef typename Surface_graph<Vector> Surface_graph;
            typedef typename Surface_graph::Surfel Surfel;

            typedef std::size_t Size;
            typedef const Surfel* Vertex;

            Size degree(const Vertex& v) const { return v->neighbors.size(); }
            Size bestCapacity() const { return 12; }

            template <typename OutputIterator>
            void writeNeighbors(OutputIterator out, const Vertex& v) const
            {
                const std::set<Surfel*>& edge_neighbors = v->neighbors;
                for (const Surfel* v_n : edge_neighbors)
                    *out++ = v_n;
            }

            template <typename OutputIterator, typename VertexPredicate>
            void writeNeighbors(OutputIterator out, const Vertex& v, const VertexPredicate& p) const
            {
                const std::set<Surfel*>& edge_neighbors = v->neighbors;
                for (const Surfel* v_n : edge_neighbors)
                    if (p(v, v_n))
                        *out++ = v_n;
            }
        };


        template <typename Point, typename Vector>
        struct Surfel_plane_getter
        {
            using Surface_graph = Surface_graph<Vector>;
            using Surfel = typename Surface_graph::Surfel;
            const std::size_t &xdim, ydim;
            Surfel_plane_getter(const std::size_t &xdim, const std::size_t &ydim)
                : xdim(xdim), ydim(ydim)
            {
            }
            Point get_point(const Surfel *surfel) const
            {
                return surfel->template get_center_point<Point>(xdim, ydim);
            }
            Vector get_normal(const Surfel *surfel) const
            {
                return surfel->normal;
            }
            std::pair<Point, Vector> operator()(const Surfel *surfel) const
            {
                return std::make_pair(get_point(surfel), get_normal(surfel));
            }
        };

        // Checks if the squared distance from the starting point to the surfel's center point is less than sq_distance
        template <typename Point, typename Vector>
        struct Distance_predicate
        {
            using Surface_graph = Surface_graph<Vector>;
            using Surfel = typename Surface_graph::Surfel;
            typedef const Surfel *Vertex;

            const double sq_distance;
            const Point from_position;
            const Surfel_plane_getter<Point, Vector> &surfel_plane_getter; // used to get the surfel center point

            Distance_predicate(const double &sq_distance, const Point &position, const Surfel_plane_getter<Point, Vector> &surfel_plane_getter)
                : sq_distance(sq_distance), from_position(position), surfel_plane_getter(surfel_plane_getter)
            {
            }
            bool operator()(const Surfel *v, const Surfel *vn) const
            {
                return CGAL::squared_distance(from_position, surfel_plane_getter.get_point(vn)) < sq_distance;
            }
        };

        // Checks if the surfel is not in selected and that the squared distance from the starting point to the surfel's center point is less than sq_distance
        template <typename Point, typename Vector>
        struct Distance_without_selected_surfels_predicate : public Distance_predicate<Point, Vector>
        {
            typedef Distance_predicate<Point, Vector> Base;
            using Surfel = typename Base::Surfel;
            typedef const Surfel *Vertex;

            const std::unordered_set<const Surfel *> &wall_surfels; // acts as a wall for the search

            Distance_without_selected_surfels_predicate(const double &sq_distance, const Point &position, const Surfel_plane_getter<Point, Vector> &surfel_plane_getter,
                                                        const std::unordered_set<const Surfel *> &wall_surfels)
                : Base(sq_distance, position, surfel_plane_getter), wall_surfels(wall_surfels)
            {
            }
            bool operator()(const Surfel *v, const Surfel *vn) const
            {
                return wall_surfels.count(v) == 0 && wall_surfels.count(vn) == 0 && Base::operator()(v, vn);
            }
        };

    } // namespace internal

    template <typename Vector>
    Surface_graph<Vector>
    make_surface_graph_from_image(const CGAL::Image_3 &image, const bool &pointel_neighboors = true)
    {
        const internal::Make_surface_graph_from_image_struct<Vector> functor;
        CGAL_IMAGE_IO_CASE(image.image(), return functor.template operator()<Word>(image, pointel_neighboors));
        return Surface_graph<Vector>();
    }

    template <typename K, typename Vector, typename SurfelPredicate = internal::All_surfels>
    void write_surface_ply(const std::string &filename,
                           const Surface_graph<Vector> &graph,
                           const SurfelPredicate &surfel_selector = SurfelPredicate())
    {
        using Point_3 = typename K::Point_3;
        std::vector<Point_3> mesh = internal::graph_to_mesh<K, Vector>(graph, surfel_selector);
        internal::write_surface_mesh_ply<Point_3>(mesh, filename);
    }

    template <typename K, typename Vector>
    void write_surface_graph_infos(const Surface_graph<Vector> &surface_graph,
                                   const std::string &out_folder)
    {
        using Point = typename K::Point_3;
        internal::write_surface_sharpness_ply<K, Vector>(surface_graph, out_folder + "_values.ply", 0.0, 1.0);
        internal::write_surface_normals_ply<K, Vector>(surface_graph, out_folder + "_normals.ply");
    }

    template <typename Vector, typename SurfelIterator>
    Surface_graph<Vector> make_sub_graph_from_surfels(const Surface_graph<Vector> &graph, SurfelIterator begin, SurfelIterator end)
    {
        typedef typename Surface_graph<Vector> Surface_graph;
        typedef typename Surface_graph::Surfel Surfel;
        Surface_graph subgraph(graph.xdim, graph.ydim, graph.zdim);
        std::map<const Surfel *, Surfel *> graph_surfel_to_new_graph_surfel;
        for (; begin != end; ++begin)
        {
            const Surfel *surfel = *begin;
            Surfel *newSurfel = subgraph.add_surfel(Surfel(surfel->voxel_index, surfel->orientation, surfel->inverse_direction));
            newSurfel->normal = surfel->normal;
            newSurfel->sharpness = surfel->sharpness;
            graph_surfel_to_new_graph_surfel[surfel] = newSurfel;
        }
        for (; begin != end; ++begin)
        {
            const Surfel *surfel = *begin;
            const std::set<Surfel *> &neighbors = surfel->neighbors;
            Surfel *newSurfel = graph_surfel_to_new_graph_surfel.at(surfel);
            std::set<Surfel *> &newSurfel_neighbors = newSurfel->neighbors;
            for (const Surfel *point_neighbor : neighbors)
            {
                if (graph_surfel_to_new_graph_surfel.find(point_neighbor) != graph_surfel_to_new_graph_surfel.end())
                    newSurfel_neighbors.insert(graph_surfel_to_new_graph_surfel.at(point_neighbor));
            }
        }
        return subgraph;
    }

    template <typename Vector, typename Predicate, typename OutputIterator>
    void get_surfel_rings_within_predicate(
        const Surface_graph<Vector>& surface_graph,
        const typename Surface_graph<Vector>::Surfel* surfel,
        const Predicate& predicate,
        OutputIterator out)
    {
        using Surface_graph = Surface_graph<Vector>;
        using Surfel = typename Surface_graph::Surfel;
        using Surfel_visitor = internal::Surfel_visitor<Vector>;
        Surfel_visitor surfel_visitor;
        DGtal::BreadthFirstVisitor<Surfel_visitor, std::set<const Surfel*> > visitor
            ( surfel_visitor, surfel );
        while ( ! visitor.finished() )
        {
            const auto& [explored_surfel, topological_distance] = visitor.current();
            *out++ = explored_surfel;
            visitor.expand(predicate);
        }
    }

    template <typename Point, typename Vector, typename OutputIterator>
    void get_close_surfels_on_surface(
        const typename Surface_graph<Vector>::Surfel *surfel, const Surface_graph<Vector> &surface_graph,
        const double &distance,
        OutputIterator out)
    {
        using Surface_graph = Surface_graph<Vector>;
        using Surfel = typename Surface_graph::Surfel;
        const double sq_distance = distance*distance;
        internal::Surfel_plane_getter<Point, Vector> surfel_plane_getter(surface_graph.xdim, surface_graph.ydim);
        const internal::Distance_predicate<Point, Vector> distance_predicate(
            sq_distance,
            surfel_plane_getter.get_point(surfel), // distance from the center of the surfel
            surfel_plane_getter);
        Surface_process::get_surfel_rings_within_predicate(surface_graph, surfel, distance_predicate, out);
    }

    template <typename Vector, typename Point>
    void get_surfel_frontier(
        const Surface_graph<Vector> &surface_graph,
        const Polylines<Point> &polylines,
        std::set<const typename Surface_graph<Vector>::Surfel *> &corner_surfels,
        std::vector<std::set<const typename Surface_graph<Vector>::Surfel *>> &line_frontier_surfels,
        std::unordered_set<const typename Surface_graph<Vector>::Surfel *> &frontier_surfels,
        const double &step_size = 1.0,
        const double &min_advancement_epsilon = 0.1) // garantees termination
    {
        using Surface_graph = Surface_graph<Vector>;
        using Surfel = typename Surface_graph::Surfel;

        corner_surfels.clear();
        line_frontier_surfels.clear();
        line_frontier_surfels.resize(polylines.lines.size());
        frontier_surfels.clear();
        frontier_surfels.reserve(polylines.lines.size());
        const std::vector<Point> &points = polylines.points;
        Surface_process::internal::Surfel_plane_getter<Point, Vector> get_plane(surface_graph.xdim, surface_graph.ydim);
        for (std::size_t line_index = 0, number_of_lines = polylines.lines.size(); line_index < number_of_lines; line_index++)
        {
            const std::vector<std::size_t> &line = polylines.lines[line_index];
            std::set<const Surfel *> &cur_line_frontier_surfels = line_frontier_surfels[line_index];
            std::size_t current_line_point = 0;
            double current_edge_factor = 0;
            // start from the front corner and find surfels while traversing the line to the back corner
            const Surfel *current_surfel = surface_graph.template get_closest_surfel<Point>(points[line.front()]);
            corner_surfels.insert(current_surfel);
            while (current_line_point < line.size() - 1)
            {
                cur_line_frontier_surfels.insert(current_surfel);
                // get aim point at current_line_point, current_edge_factor + step_size
                const auto &[aim_line_point, aim_edge_factor] = Polylines_process::get_line_point_at_geodesic_distance<Vector>(polylines, line_index, current_line_point, current_edge_factor, step_size);
                const Point &aim_point = Polylines_process::get_line_point_with_edge_factor<Vector>(polylines, line_index, aim_line_point, aim_edge_factor);
                const Surfel *closest_surfel = Surface_process::get_surfel_neighbor_closest_to_point(surface_graph, current_surfel, aim_point);
                // travel to closest surfel
                current_surfel = closest_surfel;
                // update current_line_point
                if (aim_line_point == line.size() - 1)
                {
                    // aim is at the end of the line, end loop
                    current_line_point = aim_line_point;
                    current_edge_factor = 0.0;
                }
                else
                {
                    // setup next iteration : project the current line point in the range [current + epsilon ; aim + step_size]
                    const auto &[min_projection_line_point, min_projection_edge_factor] = Polylines_process::get_line_point_at_geodesic_distance<Vector>(polylines, line_index, current_line_point, current_edge_factor, min_advancement_epsilon);
                    const auto &[max_projection_line_point, max_projection_edge_factor] = Polylines_process::get_line_point_at_geodesic_distance<Vector>(polylines, line_index, aim_line_point, aim_edge_factor, step_size);
                    const auto &[proj_line_point, proj_edge_factor] = Polylines_process::project_point_to_line_point_between_start_end<Vector>(polylines, line_index,
                                                                                                                                                min_projection_line_point, min_projection_edge_factor,
                                                                                                                                                max_projection_line_point, max_projection_edge_factor,
                                                                                                                                                get_plane.get_point(closest_surfel));
                    current_line_point = proj_line_point;
                    current_edge_factor = proj_edge_factor;
                }
            }
            cur_line_frontier_surfels.insert(current_surfel);
            frontier_surfels.insert(cur_line_frontier_surfels.begin(), cur_line_frontier_surfels.end());
            corner_surfels.insert(current_surfel);
        }
    }

    template <typename Point, typename Vector, typename OutputIterator>
    void get_surfels_near_frontier(
        const Surface_graph<Vector> &surface_graph,
        const std::unordered_set<const typename Surface_graph<Vector>::Surfel*> &frontier_surfels,
        const double &distance,
        OutputIterator out)
    {
        using Surface_graph = Surface_graph<Vector>;
        using Surfel = typename Surface_graph::Surfel;
        internal::Surfel_plane_getter<Point, Vector> surfel_plane_getter(surface_graph.xdim, surface_graph.ydim);
        for (const Surfel *surfel : frontier_surfels)
            get_close_surfels_on_surface<Point>(surfel, surface_graph, distance, out);
    }

    // outputs the normals of the surface graph following the surfel order given by surfels
    //     i.e. normals[i] correspond to the surfel at surfels[i]
    template <typename SurfelOrderPoint, typename OutputNormal, typename Vector>
    void extract_normals(const Surface_graph<Vector>& surface_graph, const std::vector<SurfelOrderPoint>& surfels, std::vector<OutputNormal>& normals)
    {
        using Surface_graph = Surface_graph<Vector>;
        using Surfel = typename Surface_graph::Surfel;
        std::map<std::array<int, 3>, std::size_t> doubled_int_point_to_index;
        for (std::size_t surfel_index = 0, nb_surfels = surfels.size(); surfel_index < nb_surfels; surfel_index++)
        {
            const SurfelOrderPoint& surfel_center = surfels[surfel_index];
            std::array<int, 3> int_point = {static_cast<int>(surfel_center[0]*2+1.5),
                                            static_cast<int>(surfel_center[1]*2+1.5),
                                            static_cast<int>(surfel_center[2]*2+1.5)};
            doubled_int_point_to_index.emplace(int_point, surfel_index);
        }
        normals.reserve(surfels.size());
        internal::Surfel_plane_getter<SurfelOrderPoint, Vector> surfel_plane_getter(surface_graph.xdim, surface_graph.ydim);
        for (const auto& [voxel_id, voxel_surfels] : surface_graph.surface)
            for (const Surfel& surfel : voxel_surfels)
            {
                const SurfelOrderPoint& surfel_center = surfel_plane_getter.get_point(&surfel);
                const Vector& surfel_normal = surfel.normal;
                std::array<int, 3> int_point = {static_cast<int>(surfel_center[0]*2+0.5),
                                                static_cast<int>(surfel_center[1]*2+0.5),
                                                static_cast<int>(surfel_center[2]*2+0.5)};
                const auto& find_it = doubled_int_point_to_index.find(int_point);
                if (find_it == doubled_int_point_to_index.end())
                {
                    std::cout << "WARNING : a surfel has no point in surface, ignored" << std::endl;
                    continue;
                }
                const std::size_t& surfel_index = find_it->second;
                normals[surfel_index] = OutputNormal(surfel_normal[0], surfel_normal[1], surfel_normal[2]);
            }
    }

    template <typename Vector>
    void apply_normals(
        const std::map<const typename Surface_graph<Vector>::Surfel *, Vector> &sharpened_normals,
        Surface_graph<Vector>& surface_graph)
    {
        using Surface_graph = Surface_graph<Vector>;
        using Surfel = typename Surface_graph::Surfel;
        for (auto& [surfel, normal] : sharpened_normals)
        {
            const_cast<Surfel*>(surfel)->normal = normal;
        }
    }

    // TODO : solve only for the surfels that changed group during the optimization process
    template <typename Vector>
    std::map<const typename Surface_graph<Vector>::Surfel*, Vector> sharpen_normals(
        const Surface_graph<Vector>& surface_graph,
        const std::set<const typename Surface_graph<Vector>::Surfel*>& surfels_to_sharpen,
        const std::unordered_set<const typename Surface_graph<Vector>::Surfel*>& frontier_surfels,
        const std::size_t &maximum_number_of_iteration = 1000,
        const double& minimum_value_change = 1.e-3)
    {
        using Surface_graph = Surface_graph<Vector>;
        using Surfel = typename Surface_graph::Surfel;
        using Poisson_graph = Element_graph<Vector>;
        using Poisson_graph_element = typename Poisson_graph::Graph_element;

        helper::bLine("build graph");
        std::map<const Surfel*, Graph_poisson_surfel_element<Vector>*> surfel_to_element;
        std::map<const Poisson_graph_element*, const Surfel*> element_to_surfel;
        std::shared_ptr<Poisson_graph> poisson_graph_ptr = Element_graph_process::build_poisson_element_graph_from_surfels(
            surface_graph, surfels_to_sharpen, frontier_surfels,
            surfel_to_element, element_to_surfel);
        Poisson_graph &poisson_graph = *poisson_graph_ptr;
        helper::eLine();

        helper::bLine("solve normals");
        std::vector<std::optional<Vector>> element_values =
            Element_graph_process::propagate_values_via_average_iterations(poisson_graph, maximum_number_of_iteration, minimum_value_change);
        helper::eLine();

        // output with surfels
        std::map<const Surfel*, Vector> surfel_values;
        for (const auto& [surfel, surfel_element] : surfel_to_element)
        {
            const std::size_t &element_index = surfel_element->get_index();
            const std::optional<Vector> &value = element_values[element_index];
            if (!value.has_value())
                continue;
            Vector final_vector = value.value();
            helper::normalize_if_possible(final_vector); // the solver don't normalize values
            surfel_values.emplace(surfel, final_vector);
        }

        return surfel_values;
    }

    template <typename Point_3, typename Vector>
    const typename Surface_graph<Vector>::Surfel*
    get_surfel_neighbor_closest_to_point(
        const Surface_graph<Vector> &surface_graph,
        const typename Surface_graph<Vector>::Surfel* surfel,
        const Point_3 &aim_point)
    {
        using Surface_graph = Surface_graph<Vector>;
        using Surfel = typename Surface_graph::Surfel;
        const Surfel *closest_surfel = nullptr;
        double min_sq_distance = std::numeric_limits<double>::max();
        for (const Surfel *surfel : surfel->neighbors)
        {
            const Point_3 &surfel_point = surfel->typename get_center_point<Point_3>(surface_graph.xdim, surface_graph.ydim);
            double sq_distance = CGAL::squared_distance(surfel_point, aim_point);
            if (sq_distance < min_sq_distance)
            {
                min_sq_distance = sq_distance;
                closest_surfel = surfel;
            }
        }
        return closest_surfel;
    }

} // namespace Surface_process

#endif // SURFACE_PROCES_H
