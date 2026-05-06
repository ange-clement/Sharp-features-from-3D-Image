#ifndef GET_PLANE_INTERSECTION_ENERGY
#define GET_PLANE_INTERSECTION_ENERGY

#include "../../surface/Surface_graph.h"
#include "../../surface/surface_process.h"
#include "../../helper/cgal_helper.h"
#include "../Surfel_point.h"



#include "../../helper/geometry_helper.h"

namespace Feature_graph_energy
{
    namespace internal
    {

        template <typename Surfel>
        bool surfel_group_contains_surfel(const std::vector<std::set<const Surfel *>> &surfel_groups, const Surfel *surfel)
        {
            for (const std::set<const Surfel *> &other_group : surfel_groups)
            {
                if (other_group.count(surfel) > 0)
                    return true;
            }
            return false;
        }

        // returns true if the two groups share a common surfel
        //   assumes that the groups are sorted
        template <typename Surfel>
        bool surfel_group_intersects_group(const std::set<const Surfel *> &group, const std::set<const Surfel *> &other_group)
        {
            using it = typename std::set<const Surfel *>::const_iterator;
            it group_it = group.begin(), group_end = group.end();
            it other_group_it = other_group.begin(), other_group_end = other_group.end();
            while (group_it != group_end && other_group_it != other_group_end)
            {
                if (*group_it < *other_group_it)
                    ++group_it;
                else if (*other_group_it < *group_it)
                    ++other_group_it;
                else
                    return true;
            }
            return false;
        }

        // returns true if the group list share a common surfel with the group
        template <typename Surfel>
        bool surfel_group_list_intersects_group(const std::vector<std::set<const Surfel *>> &surfel_groups, const std::set<const Surfel *> &group)
        {
            // Assumes that the sets are sorted
            for (const std::set<const Surfel *> &other_group : surfel_groups)
            {
                if (surfel_group_intersects_group(group, other_group))
                    return true;
            }
            return false;
        }

    } // namespace internal

    template <typename K, typename Vector>
    struct Get_plane_intersection_energy
    {
        using Surface_graph = Surface_graph<Vector>;
        using Surfel = typename Surface_graph::Surfel;

        using Point_3 = typename K::Point_3;
        using Vector_3 = typename K::Vector_3;

        using Plane = std::pair<Point_3, Vector_3>;
        using Frontier_surfel_and_planes = std::pair<const Surfel *, std::vector<Plane>>;
        using KD_tree = helper::KD_tree<boost::tuple<Point_3, std::size_t>, K>;
        using KD_tree_iterator_type = typename KD_tree::Iterator_type;

        using Surfel_point = Surfel_point<K, Vector>;

        KD_tree corners_surfel_frontier_kd_tree;            // all corners
        std::vector<KD_tree> lines_surfel_frontier_kd_tree; // only line surfels
        std::vector<Frontier_surfel_and_planes> surfel_and_planes;
        const Surface_graph &surface_graph;

        void get_surfel_groups_and_close_surfels(
            const Surfel *surfel,
            const std::unordered_set<const Surfel *> &frontier_surfels,
            const double &plane_detection_distance,
            const Surface_graph &surface_graph,
            std::vector<std::set<const Surfel *>> &surfel_groups,
            std::set<const Surfel *> &close_surfels) const
        {
            using Surfel_group = std::set<const Surfel *>;
            Surface_process::internal::Surfel_plane_getter<Point_3, Vector> surfel_plane_getter(surface_graph.xdim, surface_graph.ydim);
            const double distance_wall = plane_detection_distance + 0.5; // avoids floating point error
            const double sq_distance_wall = distance_wall * distance_wall;
            const Surface_process::internal::Distance_without_selected_surfels_predicate<Point_3, Vector> distance_wall_predicate(
                sq_distance_wall,
                surfel_plane_getter.get_point(surfel),
                surfel_plane_getter, frontier_surfels); // allows to get all surfels on surface near surfel without crossing the supporting surfels
            std::vector<const Surfel *> close_surfels_on_surface;
            Surface_process::get_close_surfels_on_surface<Point_3>(surfel, surface_graph, plane_detection_distance, std::back_inserter(close_surfels_on_surface));
            for (const Surfel *close_surfel : close_surfels_on_surface)
            {
                close_surfels.insert(close_surfel);
                // do not start a search if the starting surfel is already in a group
                if (internal::surfel_group_contains_surfel(surfel_groups, close_surfel))
                    continue;
                Surfel_group group;
                Surface_process::get_surfel_rings_within_predicate(surface_graph, close_surfel, distance_wall_predicate, std::inserter(group, group.begin()));
                if (group.size() > 1)
                    surfel_groups.push_back(group);
            }
        }

        void solve_incorrect_line_surfels(
            std::map<const Surfel *, std::vector<std::set<const Surfel *>>> &cur_line_surfel_groups,
            std::list<const Surfel *> &incorrect_surfels) const
        {
            // Solve line surfels that have more than 2 planes:
            //  remove groups that do not intersect a neighbor line surfel with 2 planes
            typename std::list<const Surfel *>::iterator incorrect_surfels_it;
            bool solved_at_least_one_surfel = true;
            while (solved_at_least_one_surfel)
            {
                incorrect_surfels_it = incorrect_surfels.begin();
                solved_at_least_one_surfel = false;
                while (incorrect_surfels_it != incorrect_surfels.end())
                {
                    bool solved_current_surfel = false;
                    const Surfel *surfel = *incorrect_surfels_it;
                    // for each neighbor, erase the unwanted groups of this surfel until there are only two groups left
                    for (const Surfel *neighbor : surfel->neighbors)
                    {
                        const auto &neighbor_surfel_groups_iterator = cur_line_surfel_groups.find(neighbor);
                        if (neighbor_surfel_groups_iterator == cur_line_surfel_groups.end()) // skips neighbor surfel does not belong to the line
                            continue;
                        // erases all groups that do not intersect with the neighbor's group
                        std::vector<std::set<const Surfel *>> &surfel_groups = cur_line_surfel_groups.at(surfel);
                        const std::vector<std::set<const Surfel *>> &neighbor_groups = neighbor_surfel_groups_iterator->second;
                        typename std::vector<std::set<const Surfel *>>::iterator group_iterator = surfel_groups.begin();
                        while (group_iterator != surfel_groups.end() && surfel_groups.size() > 2)
                        {
                            if (!internal::surfel_group_list_intersects_group(neighbor_groups, *group_iterator))
                            {
                                group_iterator = surfel_groups.erase(group_iterator);
                            }
                            else
                                group_iterator++;
                        }
                        if (surfel_groups.size() <= 2) // we sucessfully cleaned the groups of this surfel
                        {
                            solved_at_least_one_surfel = true;
                            solved_current_surfel = true;
                            break;
                        }
                    }
                    if (solved_current_surfel)
                        incorrect_surfels_it = incorrect_surfels.erase(incorrect_surfels_it);
                    else
                        incorrect_surfels_it++;
                }
            }
            if (!incorrect_surfels.empty())
            {
                std::cout << "WARNING: Could not solve the groups of some surfels, they are treated as corners" << std::endl;
            }
        }

        std::vector<std::map<const Surfel *, std::vector<Plane>>> make_surfel_planes(
            const std::vector<std::set<const Surfel *>> &line_frontier_surfels,
            const std::unordered_set<const Surfel *> &frontier_surfels,
            const std::set<const Surfel *> &corner_surfels,
            const Surface_graph &surface_graph,
            const double &plane_detection_distance) const
        {
            using Surfel_group = std::set<const Surfel *>;
            using Line_surfel_to_groups = std::map<const Surfel *, std::vector<Surfel_group>>;

            // Make groups via breath first search, and use neighbors in ambiguous cases
            helper::bLine("Get surfel patches");
            std::vector<Line_surfel_to_groups> all_surfel_groups; // for each line, map the line's surfel to its groups
            all_surfel_groups.reserve(line_frontier_surfels.size());
            std::set<const Surfel *> surfels_too_close_to_feature; // contains all close surfels, used for the normal rectification
            for (const std::set<const Surfel *> &cur_line_surfels : line_frontier_surfels)
            {
                all_surfel_groups.resize(all_surfel_groups.size() + 1);
                Line_surfel_to_groups &cur_line_surfel_groups = all_surfel_groups.back();
                std::list<const Surfel *> incorrect_surfels; // line surfels with more than two planes
                for (const Surfel *line_surfel : cur_line_surfels)
                {
                    std::vector<Surfel_group> &surfel_groups = cur_line_surfel_groups[line_surfel];
                    get_surfel_groups_and_close_surfels(line_surfel, frontier_surfels, plane_detection_distance, surface_graph, surfel_groups, surfels_too_close_to_feature);
                    // store surfels that are not a corner with more than two planes
                    if (corner_surfels.count(line_surfel) == 0 && surfel_groups.size() > 2)
                        incorrect_surfels.push_back(line_surfel);
                }
                solve_incorrect_line_surfels(cur_line_surfel_groups, incorrect_surfels);
            }
            helper::eLine();

            // Solve normals on surfels too close
            helper::bBlock("Sharpen normals");
            const std::map<const Surfel *, Vector> &sharpened_normals = Surface_process::sharpen_normals(surface_graph, surfels_too_close_to_feature, frontier_surfels);
            helper::eBlock();

            // Make planes
            std::vector<std::map<const Surfel *, std::vector<Plane>>> all_surfel_planes;
            all_surfel_planes.reserve(line_frontier_surfels.size());
            Surface_process::internal::Surfel_plane_getter<Point_3, Vector> get_plane(surface_graph.xdim, surface_graph.ydim);
            for (const auto &cur_line_surfel_groups : all_surfel_groups)
            {
                all_surfel_planes.resize(all_surfel_planes.size() + 1);
                std::map<const Surfel *, std::vector<Plane>> &cur_line_planes = all_surfel_planes.back();
                for (const auto &[surfel, surfel_groups] : cur_line_surfel_groups)
                {
                    std::vector<Plane> &pointel_planes = cur_line_planes[surfel];
                    for (const std::set<const Surfel *> &group : surfel_groups)
                    {
                        Vector_3 group_centroid(0, 0, 0);
                        Vector_3 group_normal(0, 0, 0);
                        std::size_t nb_surfel_taken = group.size();
                        // Create a plane as the average point and normal of surfels that are not too close of the feature
                        for (const Surfel *surfel : group)
                        {
                            const auto &[surfel_center, surfel_normal] = get_plane(surfel);
                            Vector_3 surfel_normal_for_plane;
                            // use new normal for close surfels
                            const auto &sharp_normal_iterator = sharpened_normals.find(surfel);
                            if (sharp_normal_iterator != sharpened_normals.end())
                            {
                                surfel_normal_for_plane = sharp_normal_iterator->second;
                            }
                            else
                            {
                                surfel_normal_for_plane = surfel_normal;
                            }
                            group_centroid += Vector_3(Point_3(0, 0, 0), surfel_center);
                            group_normal += surfel_normal_for_plane;
                        }
                        const double nb_surfel_taken_double = static_cast<double>(nb_surfel_taken);
                        const Point_3 plane_center = Point_3(0, 0, 0) + (group_centroid / nb_surfel_taken_double);
                        const Vector_3 plane_normal = group_normal / nb_surfel_taken_double;
                        pointel_planes.emplace_back(plane_center, plane_normal);
                    }
                }
            }
            return all_surfel_planes;
        }

        template <typename Point>
        Get_plane_intersection_energy(
            const Surface_graph &surface_graph, const Polylines<Point> &polylines,
            double plane_detection_distance = 4.0)
            : surface_graph(surface_graph)
        {
            using Polylines = Polylines<Point>;
            helper::bLine("Find surfels supporting polylines");
            std::set<const Surfel *> corner_surfels;
            std::vector<std::set<const Surfel *>> line_frontier_surfels;
            std::unordered_set<const Surfel *> frontier_surfels;
            Surface_process::get_surfel_frontier(surface_graph, polylines, corner_surfels, line_frontier_surfels, frontier_surfels);
            helper::eLine();
            helper::bBlock("Make planes and initialize Kd_trees");
            surfel_and_planes.reserve(frontier_surfels.size());
            Surface_process::internal::Surfel_plane_getter<Point_3, Vector> get_plane(surface_graph.xdim, surface_graph.ydim);
            const double sq_distance = plane_detection_distance * plane_detection_distance;
            std::vector<std::map<const Surfel *, std::vector<Plane>>> all_surfel_planes =
                make_surfel_planes(line_frontier_surfels, frontier_surfels, corner_surfels,
                                   surface_graph, plane_detection_distance);
            std::set<const Surfel *> processed_surfels;
            lines_surfel_frontier_kd_tree.resize(line_frontier_surfels.size());
            for (std::size_t line_index = 0, number_of_lines = line_frontier_surfels.size(); line_index < number_of_lines; line_index++)
            {
                const std::set<const Surfel *> &cur_line_surfels = line_frontier_surfels[line_index];
                const std::map<const Surfel *, std::vector<Plane>> &cur_line_surfel_planes = all_surfel_planes[line_index];
                KD_tree &cur_line_kd_tree = lines_surfel_frontier_kd_tree[line_index];
                for (const Surfel *surfel : cur_line_surfels)
                {
                    const bool is_surfel_processed = !processed_surfels.insert(surfel).second;
                    if (is_surfel_processed)
                        continue;
                    const bool is_corner = corner_surfels.count(surfel) > 0;
                    const std::vector<Plane> &surfel_planes = cur_line_surfel_planes.at(surfel);
                    const auto &surfel_tree_point = boost::make_tuple(get_plane.get_point(surfel), surfel_and_planes.size());
                    if (is_corner || surfel_planes.size() > 2)
                    {
                        corners_surfel_frontier_kd_tree.push(surfel_tree_point);
                    }
                    else
                    {
                        cur_line_kd_tree.push(surfel_tree_point);
                    }
                    surfel_and_planes.emplace_back(surfel, surfel_planes);
                }
            }
            helper::eBlock();
        }

        KD_tree &get_kd_tree_of_point(const bool &is_corner, const std::size_t &line_index)
        {
            return (is_corner) ? corners_surfel_frontier_kd_tree : lines_surfel_frontier_kd_tree[line_index];
        }

        double operator()(const Surfel_point &surfel_point, const double &delta_x, const double &delta_y, const bool &is_corner, const std::size_t &line_index)
        {
            constexpr std::size_t point_maximum_number_of_incident_surfel = 12;
            constexpr double maximum_squared_distance = 0.6; // d > sqrt(2)/2 is the same as d^2 > 0.5, then add epsilon
            // 1. move point on surface
            Surfel_point moved_point(surfel_point);
            moved_point.move(typename K::Vector_2(delta_x, delta_y), surface_graph);
            // 2. get closest surfels and planes to the moved point
            const Point_3 &surfel_world_point = moved_point.get_point_3(surface_graph.xdim, surface_graph.ydim);
            auto [it_begin, it_end] = get_kd_tree_of_point(is_corner, line_index).get_closest(surfel_world_point, point_maximum_number_of_incident_surfel);
            // if the tree is empty, look into the other tree
            if (it_begin == it_end)
            {
                const auto &[other_tree_it_begin, other_tree_it_end] = get_kd_tree_of_point(!is_corner, line_index).get_closest(surfel_world_point, point_maximum_number_of_incident_surfel);
                it_begin = other_tree_it_begin;
                it_end = other_tree_it_end;
                // if both tree are empty, return max value as error
                if (it_begin == it_end)
                {
                    return std::numeric_limits<double>::max();
                }
            }
            // 3. average planes projection on surfels adjacent to the point
            //  usefull in cases where the point is close to surfel edges and corners
            Vector_3 avg_plane_intersection(0, 0, 0);
            std::size_t nb_adjacent_surfels = 0;
            for (KD_tree_iterator_type it = it_begin; it != it_end; it++)
            {
                const auto &[tree_point, squared_distance] = *it;
                if (squared_distance > maximum_squared_distance)
                    continue;
                const std::size_t &point_id = tree_point->template get<1>();
                const auto &nearest_surfel_and_planes = surfel_and_planes[point_id];
                const Point_3 plane_intersection = helper::project_at_plane_intersection<K>(nearest_surfel_and_planes.second, surfel_world_point);
                avg_plane_intersection += Vector_3(Point_3(0, 0, 0), plane_intersection);
                nb_adjacent_surfels++;
            }
            if (nb_adjacent_surfels == 0)
            {
                // no surfel close enough, taking the closest
                const auto &[tree_point, distance] = *it_begin;
                nb_adjacent_surfels++;
                const std::size_t &point_id = tree_point->template get<1>();
                const auto &nearest_surfel_and_planes = surfel_and_planes[point_id];
                const Point_3 plane_intersection = helper::project_at_plane_intersection<K>(nearest_surfel_and_planes.second, surfel_world_point);
                avg_plane_intersection += Vector_3(Point_3(0, 0, 0), plane_intersection);
            }
            const Point_3 plane_intersection = Point_3(0, 0, 0) + avg_plane_intersection / static_cast<double>(nb_adjacent_surfels);
            return CGAL::sqrt(CGAL::squared_distance(surfel_world_point, plane_intersection));
        }
    };
} // namespace Feature_graph_energy

#endif // GET_PLANE_INTERSECTION_ENERGY