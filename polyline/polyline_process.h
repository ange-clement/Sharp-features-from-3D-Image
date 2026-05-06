#ifndef POLYLINE_PROCESS_H
#define POLYLINE_PROCESS_H

#include "Polylines.h"
#include "../helper/geometry_helper.h"
#include "../element_graph/Element_graph.h"

namespace Polylines_process
{
    namespace internal
    {
        void get_one_line_from_graph_recursive(
            const std::vector<std::vector<std::size_t>> &graph_vector,
            const std::size_t current, const std::size_t previous,
            std::vector<bool> &processed_elements, std::vector<std::size_t> &line)
        {
            const std::vector<std::size_t> &current_neighbors = graph_vector[current];
            if (current_neighbors.size() != 2 || processed_elements[current])
            {
                line.clear();
                line.push_back(current);
            }
            else
            {
                processed_elements[current] = true;
                // the next neighbor is the one different from the previous (faster than looking in processed elements)
                std::size_t next = current_neighbors[0];
                if (next == previous)
                    next = current_neighbors[1];
                get_one_line_from_graph_recursive(graph_vector, next, current, processed_elements, line);
                line.push_back(current);
            }
        }

        template <typename Point>
        void write_as_text(const Polylines<Point> &polylines, const std::string &filename)
        {
            std::ofstream polylines_file(filename);
            for (const std::vector<std::size_t> &poly : polylines.lines)
            {
                polylines_file << poly.size();
                for (const std::size_t &p_id : poly)
                {
                    const Point &p = polylines.points[p_id];
                    polylines_file << " " << p[0] << " " << p[1] << " " << p[2];
                }
                polylines_file << "\n";
            }
            polylines_file.close();
        }

        template <typename Point, typename Vector>
        void make_mesh_from_polylines(
            const Polylines<Point> &polylines,
            const double &edge_size, const double &sphere_size,
            std::vector<Point> &points,
            std::vector<std::array<std::size_t, 3>> &triangles,
            std::vector<std::size_t> &triangle_color_index,
            const std::size_t &sphere_subidvision = 2,
            const std::size_t &cylinder_number_of_faces = 4)
        {
            std::unordered_set<std::size_t> processed_corner;
            for (std::size_t poly_id = 0, poly_size = polylines.lines.size(); poly_id < poly_size; ++poly_id)
            {
                const std::vector<std::size_t> &poly = polylines.lines[poly_id];
                if (poly.size() == 0)
                    continue;
                // corners
                std::vector<std::size_t> poly_corner;
                if (poly.size() == 1)
                    poly_corner.push_back(poly.front());
                else if (poly.size() > 1)
                {
                    const std::size_t &first_corner = poly.front();
                    const std::size_t &second_corner = poly.back();
                    if (first_corner != second_corner) // not a loop
                    {
                        poly_corner.push_back(first_corner);
                        poly_corner.push_back(second_corner);
                    }
                }
                for (const std::size_t &corner : poly_corner)
                {
                    const auto &[it, inserted] = processed_corner.insert(corner);
                    if (!inserted)
                        continue;
                    const Point &corner_point = polylines.points[corner];
                    helper::add_mesh_sphere_from_point<Point, Vector>(corner_point, points, triangles, sphere_size, sphere_subidvision);
                    triangle_color_index.resize(triangles.size(), 0);
                }
                // edges
                if (poly.size() < 2)
                    continue;
                helper::add_mesh_cylinder_from_polyline<Point, Vector>(poly, polylines.points, edge_size, cylinder_number_of_faces, points, triangles);
                triangle_color_index.resize(triangles.size(), poly_id + 1);
            }
        }

        // returns a vector that, for each corner, gives the adjacent lines
        template <typename Point>
        std::vector<std::vector<std::size_t>> get_corner_to_line_adjacency(const Polylines<Point> &polylines)
        {
            std::vector<std::vector<std::size_t>> corner_to_lines;
            corner_to_lines.resize(polylines.points.size());
            for (std::size_t line_index = 0, number_of_lines = polylines.lines.size(); line_index < number_of_lines; line_index++)
            {
                const std::vector<std::size_t> &line = polylines.lines[line_index];
                if (line.empty())
                    continue;
                corner_to_lines[line.front()].push_back(line_index);
                if (line.size() > 1 && line.front() != line.back()) // the line has more than one element and is not a loop
                    corner_to_lines[line.back()].push_back(line_index);
            }
            return corner_to_lines;
        }

        template <typename Point>
        std::array<Point, 2> get_polylines_bbox(const Polylines<Point> &polylines)
        {
            if (polylines.points.empty())
                return {Point(0, 0, 0), Point(0, 0, 0)};
            std::array<double, 3> min_coords{polylines.points[0][0], polylines.points[0][1], polylines.points[0][2]};
            std::array<double, 3> max_coords{polylines.points[0][0], polylines.points[0][1], polylines.points[0][2]};
            for (const Point &p : polylines.points)
            {
                for (int i = 0; i < 3; i++)
                {
                    const double &point_coord = p[i];
                    if (point_coord < min_coords[i])
                        min_coords[i] = point_coord;
                    if (point_coord > max_coords[i])
                        max_coords[i] = point_coord;
                }
            }
            return {Point(min_coords[0], min_coords[1], min_coords[2]), Point(max_coords[0], max_coords[1], max_coords[2])};
        }

        void remove_line_from_corner_lines(
            const std::size_t &line_index,
            const std::size_t &corner_index,
            std::vector<std::vector<std::size_t>> &corner_to_lines)
        {
            std::vector<std::size_t> &lines_of_corner = corner_to_lines[corner_index];
            for (std::vector<std::size_t>::iterator corner_line_iterator = lines_of_corner.begin(); corner_line_iterator != lines_of_corner.end(); ++corner_line_iterator)
            {
                if (*corner_line_iterator == line_index)
                {
                    lines_of_corner.erase(corner_line_iterator);
                    return;
                }
            }
        }


    } // namespace internal

    std::vector<std::vector<std::size_t>>
    make_lines_from_graph(
        const std::vector<std::vector<std::size_t>> &graph_vector)
    {
        std::vector<std::vector<std::size_t>> lines;
        // merge one_ring
        std::vector<bool> processed_elements(graph_vector.size());
        for (std::size_t point_index = 0, size = graph_vector.size(); point_index < size; point_index++)
        {
            const std::vector<std::size_t> &neighbors = graph_vector[point_index];
            // skip non corner
            if (neighbors.size() == 2)
                continue;
            // set corner as processed
            processed_elements[point_index] = true;
            // create and fill a new line for each neighbor
            for (const std::size_t neigh : neighbors)
            {
                // skip processed neighbor
                if (processed_elements[neigh])
                    continue;
                lines.resize(lines.size() + 1);
                std::vector<std::size_t> &line = lines.back();
                internal::get_one_line_from_graph_recursive(graph_vector, neigh, point_index, processed_elements, line);
                line.push_back(point_index);
            }
        }
        // handle ring (loop) cases :
        //   process each unprocessed vertex as a corner
        for (std::size_t point_index = 0, size = graph_vector.size(); point_index < size; point_index++)
        {
            // insert or skip processed element
            if (processed_elements[point_index])
                continue;
            processed_elements[point_index] = true;
            // create a new line from the first neighbor
            const std::vector<std::size_t> &neighbors = graph_vector[point_index];
            lines.resize(lines.size() + 1);
            std::vector<std::size_t> &line = lines.back();
            internal::get_one_line_from_graph_recursive(graph_vector, neighbors[0], point_index, processed_elements, line);
            line.push_back(point_index);
        }

        return lines;
    }

    template <typename Point>
    Polylines<Point> make_polylines_without_floating_points_and_empty_lines(const Polylines<Point> &polylines)
    {
        Polylines<Point> new_polylines;
        std::vector<std::optional<std::size_t>> old_point_index_to_new_point_index(polylines.points.size(), std::nullopt);
        for (const std::vector<std::size_t> &line : polylines.lines)
        {
            if (line.size() == 0)
                continue;
            new_polylines.lines.resize(new_polylines.lines.size() + 1);
            std::vector<std::size_t> &new_line = new_polylines.lines.back();
            new_line.reserve(line.size());
            for (const std::size_t &point_index : line)
            {
                std::optional<std::size_t> &new_point_index = old_point_index_to_new_point_index[point_index];
                if (!new_point_index.has_value())
                {
                    new_point_index = new_polylines.points.size();
                    new_polylines.points.push_back(polylines.points[point_index]);
                }
                new_line.push_back(new_point_index.value());
            }
        }
        return new_polylines;
    }

    template <typename Point>
    double line_length(const Polylines<Point> &polylines, const std::size_t &line_index, const double &stop_length_computation_after_squared_length = std::numeric_limits<double>::max())
    {
        const std::vector<std::size_t> &line = polylines.lines[line_index];
        double line_length = 0;
        for (std::size_t p_id = 1, max_p = line.size(); p_id < max_p; p_id++)
        {
            line_length += CGAL::sqrt(CGAL::squared_distance(polylines.points[line[p_id - 1]], polylines.points[line[p_id]]));
            if (line_length * line_length > stop_length_computation_after_squared_length)
            {
                return line_length;
            }
        }
        return line_length;
    }

    template <typename Point, typename Vector>
    void write_polylines(
        const Polylines<Point> &polylines,
        const std::string &filename,
        const bool &write_as_lines_mesh = true,
        const bool &write_as_dark_lilnes_mesh = true,
        const bool &write_as_light_lilnes_mesh = true,
        const bool &write_as_text = true,
        const double &edge_size = 0.002,
        const double &sphere_size = 0.005)
    {
        if (write_as_lines_mesh || write_as_dark_lilnes_mesh || write_as_light_lilnes_mesh)
        {
            std::array<Point, 2> polylines_bbox = internal::get_polylines_bbox(polylines);
            const double diagonal_size = CGAL::sqrt(CGAL::squared_distance(polylines_bbox[0], polylines_bbox[1]));
            std::vector<Point> points;
            std::vector<std::array<std::size_t, 3>> triangles;
            std::vector<std::size_t> triangle_color_index;
            internal::make_mesh_from_polylines<Point, Vector>(
                /*In  */ polylines, edge_size * diagonal_size, sphere_size * diagonal_size,
                /*Out */ points, triangles, triangle_color_index);
            if (write_as_lines_mesh)
            {
                helper::write_mesh_as_ply(filename + "_color.ply", points, triangles, triangle_color_index, true);
            }
            if (write_as_dark_lilnes_mesh)
            {
                helper::write_mesh_as_ply(filename + "_dark.ply", points, triangles, triangle_color_index, false);
            }
            if (write_as_light_lilnes_mesh)
            {
                helper::write_mesh_as_ply(filename + "_light.ply", points, triangles, triangle_color_index, false, true);
            }
        }
        if (write_as_text)
        {
            internal::write_as_text(polylines, filename + ".polylines.txt");
        }
    }

    // returns the linearly interpolated point in the line
    template <typename Vector, typename Point>
    Point get_line_point_with_edge_factor(
        const Polylines<Point> &polylines,
        const std::size_t &line_index,
        const std::size_t &line_point_index, const double &edge_factor)
    {
        const std::vector<std::size_t> &line = polylines.lines[line_index];
        const Point edge_start = polylines.points[line[line_point_index]];
        if (edge_factor == 0 || line_point_index >= line.size() - 1)
            return edge_start;
        const Point edge_end = polylines.points[line[line_point_index + 1]];
        return edge_start + edge_factor * Vector(edge_start, edge_end);
    }

    // returns a line point: the index of the point and the factor to the next point
    template <typename Vector, typename Point>
    std::pair<std::size_t, double> get_line_point_at_geodesic_distance(
        const Polylines<Point> &polylines,
        const std::size_t &line_index,
        std::size_t line_point_index, double edge_factor, double step_size)
    {
        const double epsilon = 1.e-3;
        const std::vector<std::size_t> &line = polylines.lines[line_index];
        if (line_point_index >= line.size() - 1) // do not move if already on the last point
            return {line_point_index, 0.0};
        const Point edge_start = polylines.points[line[line_point_index]];
        const Point edge_end = polylines.points[line[line_point_index + 1]];
        const Point start_point = edge_start + edge_factor * Vector(edge_start, edge_end);
        Point previous_point = start_point;
        while (line_point_index < line.size() - 1)
        {
            const Point edge_end = polylines.points[line[line_point_index + 1]];
            const double &distance_to_next_point = CGAL::sqrt(CGAL::squared_distance(previous_point, edge_end));
            if (step_size < distance_to_next_point - epsilon) // step_size strictly under distance_to_next_point
            {
                // do not update line point, update edge value
                if (edge_factor == 0)
                    return {line_point_index, step_size / distance_to_next_point};
                else
                {
                    double new_edge_factor = edge_factor + (step_size / distance_to_next_point) * (1.0 - edge_factor);
                    return {line_point_index, new_edge_factor};
                }
            }
            else if (step_size < distance_to_next_point + epsilon) // step_size equals (in range + or - epsilon) distance_to_next_point
            {
                // start of next point
                return {line_point_index + 1, 0.0};
            }
            // step_size is strictly over distance_to_next_point
            previous_point = edge_end;
            edge_factor = 0;
            line_point_index++;
            step_size -= distance_to_next_point;
        }
        // step size is greater than the whole line, returns the last point
        return {line_point_index, 0.0};
    }

    template <typename Vector, typename Point>
    std::pair<std::size_t, double> project_point_to_line_point_between_start_end(
        const Polylines<Point> &polylines,
        const std::size_t &line_index,
        const std::size_t &line_point_start, const double &edge_factor_start,
        const std::size_t &line_point_end, const double &edge_factor_end,
        const Point &point_to_project)
    {
        const std::vector<std::size_t> &line = polylines.lines[line_index];
        if (line_point_start >= line.size() - 1) // do not move if already on the last point
            return {line_point_start, 0.0};
        Point previous_point = get_line_point_with_edge_factor<Vector>(polylines, line_index, line_point_start, edge_factor_start);
        // project the point to the edges in the range and take the closest projection
        std::size_t closest_line_point = line_point_start;
        Point closest_projection_point = previous_point;
        double closest_sq_distance = std::numeric_limits<double>::max();
        for (std::size_t current_line_point = line_point_start; current_line_point <= line_point_end; current_line_point++)
        {
            // project on the segment from the current point to the start of the next line segment (or the end of the projection range)
            const bool is_on_same_edge = current_line_point == line_point_end;
            const Point &edge_end = is_on_same_edge ? get_line_point_with_edge_factor<Vector>(polylines, line_index, line_point_end, edge_factor_end)
                                                    : polylines.points[line[current_line_point + 1]];
            const Point &projected_point = helper::get_point_segment_projection<Point, Vector>(point_to_project, previous_point, edge_end);
            const double &segment_sq_distance = CGAL::squared_distance(projected_point, point_to_project);
            if (segment_sq_distance < closest_sq_distance)
            {
                closest_sq_distance = segment_sq_distance;
                closest_line_point = current_line_point;
                closest_projection_point = projected_point;
            }
            previous_point = edge_end;
        }
        // convert the point to the edge factor
        double closest_edge_factor;
        const Point &line_segment_start = polylines.points[line[closest_line_point]];
        const double &segment_sq_length = CGAL::squared_distance(line_segment_start, polylines.points[line[closest_line_point + 1]]);
        if (segment_sq_length == 0)
            closest_edge_factor = 0.0;
        else
            closest_edge_factor = CGAL::sqrt(CGAL::squared_distance(line_segment_start, closest_projection_point) / segment_sq_length);
        if (closest_edge_factor >= 1.0)
        {
            closest_line_point++;
            closest_edge_factor = 0.0;
        }
        return {closest_line_point, closest_edge_factor};
    }

    template <typename Vector, typename Point>
    Point project_point_to_polylines(
        const Polylines<Point> &polylines,
        const Point &point_to_project)
    {
        Point projected_point = point_to_project;
        double min_point_sq_distance = std::numeric_limits<double>::max();
        for (const std::vector<std::size_t> &line : polylines.lines)
        {
            for (std::size_t element_id = 0, max_element_id = line.size() - 1; element_id < max_element_id; element_id++)
            {
                const Point &element_segment_begin = polylines.points[line[element_id]];
                const Point &element_segment_end = polylines.points[line[element_id + 1]];
                const Point &edge_projection = helper::get_point_segment_projection<Point, Vector>(point_to_project, element_segment_begin, element_segment_end);
                const double &edge_projection_sq_len = CGAL::squared_distance(point_to_project, edge_projection);
                if (edge_projection_sq_len < min_point_sq_distance)
                {
                    min_point_sq_distance = edge_projection_sq_len;
                    projected_point = edge_projection;
                }
            }
        }
        return projected_point;
    }

    // removes the polyline edges with a length less than collapse_under_distance
    //   The points are not removed from the polylines's point vector.
    template <typename Point>
    void collapse_short_edges(
        const double &collapse_under_distance,
        Polylines<Point> &polylines)
    {
        const std::vector<Point> &points = polylines.points;
        std::vector<std::vector<std::size_t>> &lines = polylines.lines;
        std::vector<std::size_t>::iterator line_iterator;
        const double collapse_under_sq_distance = collapse_under_distance * collapse_under_distance;
        for (std::vector<std::size_t> &line : lines)
        {
            if (line.size() <= 2)
                continue;
            line_iterator = line.begin();
            Point previous_point = points[*line_iterator];
            line_iterator++;
            while (line_iterator != std::prev(line.end()))
            {
                const Point &curent_point = points[*line_iterator];
                const double &sq_distance = CGAL::squared_distance(previous_point, curent_point);
                if (sq_distance < collapse_under_sq_distance)
                {
                    line_iterator = line.erase(line_iterator);
                }
                else
                {
                    previous_point = curent_point;
                    line_iterator++;
                }
            }
            // for the last element, if the edge is too short then remove the element before last
            const Point &curent_point = points[*line_iterator];
            const double &sq_distance = CGAL::squared_distance(previous_point, curent_point);
            if (sq_distance < collapse_under_sq_distance)
                line.erase(std::prev(line_iterator)); // line_iterator invalidated
        }
    }

    template <typename Point>
    std::vector<std::vector<Point>> make_polyline_for_mesh_3(const Polylines<Point> &polylines)
    {
        std::vector<std::vector<Point>> point_format_polyline = polylines.make_point_format();

        // call polyline to protect : merge identical points and split where the angle is acute
        std::vector<std::vector<Point>> polylines_copy = point_format_polyline;
        // polylines_copy.resize(point_format_polyline.size());
        // for (std::size_t i = 0, size = point_format_polyline.size(); i < size; i++)
        // {
        //     const std::vector<Point>& line = point_format_polyline[i];
        //     std::vector<Point>& line_copy = polylines_copy[i];
        //     line_copy.insert(line_copy.begin(), line.begin(), line.end());
        // }
        point_format_polyline.clear();
        CGAL::polylines_to_protect(point_format_polyline,
                                polylines_copy.begin(),
                                polylines_copy.end());
        return point_format_polyline;
    }

    template <typename Vector, typename Point, typename Image>
    Polylines<Point> make_world_polyline(
        const Polylines<Point> &polylines,
        const Image &image)
    {
        Polylines<Point> polylines_world = polylines.clone();
        Vector offset(-0.5, -0.5, -0.5);
        std::set<std::size_t> processed_corners;
        for (const std::vector<std::size_t> &line : polylines_world.lines)
        {
            if (line.size() == 0)
                continue;
            {
                const auto& [it, inserted] = processed_corners.insert(line.front());
                if (inserted)
                {
                    polylines_world.points[line.front()] = helper::image_pos_to_world_pos<Point>(polylines_world.points[line.front()] + offset, image);
                }
            }
            for (std::size_t l_p_id = 1, size = line.size()-1; l_p_id < size; l_p_id++)
            {
                const std::size_t &p_id = line[l_p_id];
                Point& point = polylines_world.points[p_id];
                point = helper::image_pos_to_world_pos<Point>(point + offset, image);
            }
            {
                const auto& [it, inserted] = processed_corners.insert(line.back());
                if (inserted)
                {
                    polylines_world.points[line.back()] = helper::image_pos_to_world_pos<Point>(polylines_world.points[line.back()] + offset, image);
                }
            }
        }
        return polylines_world;
    }

    template <typename Point>
    double get_polyline_length(
        const Polylines<Point>& poly
    )
    {
        double length = 0;
        for (const std::vector<std::size_t>& line : poly.lines)
        {
            if (line.size() == 0)
                continue;
            for (std::size_t element_id = 0, nb_p = line.size()-1; element_id < nb_p; element_id++)
            {
                const Point& p0 = poly.points[line[element_id]];
                const Point& p1 = poly.points[line[element_id+1]];
                length += CGAL::sqrt(CGAL::squared_distance(p0, p1));
            }
        }
        return length;
    }

} // namespace Polylines_process

#endif // POLYLINE_PROCESS_H