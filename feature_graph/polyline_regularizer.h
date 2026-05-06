#ifndef POLYLINE_REGULIZER_H
#define POLYLINE_REGULIZER_H

#undef NDEBUG
#include <assert.h>
#define NDEBUG

#include "../helper/geometry_helper.h"
#include "../polyline/Polylines.h"
#include "../polyline/polyline_process.h"
#include "../element_graph/graph_element.h"
#include "../element_graph/element_graph_process.h"

template <typename Vector>
struct Double_lines_regularizer
{
    std::string out_filepath = "";
    double collapse_double_lines_distance;

    enum Line_status
    {
        NOT_DOUBLE,
        ISTHMUS,
        SIMPLE,
        OTHER
    };

    // the double line condition:
    //   true iff the maximum projection squared distance from line A to line B is less than collapse_double_lines_sq_distance
    //   this criterion is not the same from line A to B as for line B to A
    template <typename Point>
    bool is_double_line_against_other_line(
        const std::size_t &line_to_check_index,
        const std::size_t &other_line_index,
        const Polylines<Point> &polylines,
        const double &collapse_double_lines_sq_distance) const
    {
        using Graph_element = Element_graph<Point>::Graph_element;
        const std::vector<std::size_t> &line_to_check = polylines.lines[line_to_check_index];
        const std::vector<std::size_t> &other_line = polylines.lines[other_line_index];
        if (line_to_check.empty() || other_line.empty())
            return false;
        if (other_line.size() == 1)
        {
            // the other line is a point: search for a point in the line to check whose distance is larger than the criterion
            const Point &other_line_point = polylines.points[other_line.front()];
            for (const std::size_t &point_index : line_to_check)
            {
                const Point &point = polylines.points[point_index];
                double sq_distance = CGAL::squared_distance(point, other_line_point);
                if (sq_distance > collapse_double_lines_sq_distance)
                {
                    // we found a point whose projection length is greater than collapse_under_max_distance
                    return false;
                }
            }
            // every point is under the threshold distance so we have a double line
            return true;
        }
        // the other line is composed of segments. Project each points on each segments
        for (const std::size_t &point_index : line_to_check)
        {
            // get the minimum projection distance from the point in line to segments of other_line
            const Point &point = polylines.points[point_index];
            double min_point_sq_distance = std::numeric_limits<double>::max();
            for (std::size_t other_element_id = 0, max_other_element_id = other_line.size() - 1; other_element_id < max_other_element_id; other_element_id++)
            {
                const Point &other_element_segment_begin = polylines.points[other_line[other_element_id]];
                const Point &other_element_segment_end = polylines.points[other_line[other_element_id + 1]];
                const Point &edge_projection = helper::get_point_segment_projection<Point, Vector>(point, other_element_segment_begin, other_element_segment_end);
                const double &edge_projection_sq_len = CGAL::squared_distance(point, edge_projection);
                if (edge_projection_sq_len < min_point_sq_distance)
                {
                    min_point_sq_distance = edge_projection_sq_len;
                }
            }
            if (min_point_sq_distance > collapse_double_lines_sq_distance)
            {
                // we found a point whose projection length is greater than collapse_under_max_distance
                return false;
            }
        }
        // every point is under the threshold distance so we have a double line
        return true;
    }

    // checks is the line is close to another line (double line criterion):
    //   true iff there is another line where the maximum projection squared distance to that line is less than (or equals) collapse_double_lines_sq_distance
    template <typename Point>
    bool should_collapse_line(
        const std::size_t &line_to_check_index,
        const Polylines<Point> &polylines,
        const double &collapse_double_lines_sq_distance) const
    {
        // empty lines are already deleted
        if (polylines.lines[line_to_check_index].empty())
            return false;
        // small lines are collapsed
        const double &line_length = Polylines_process::line_length(polylines, line_to_check_index, collapse_double_lines_sq_distance);
        if (line_length*line_length <= collapse_double_lines_sq_distance)
            return true;
        // test for every lines that is not the line to check: is this line double beacuse of the other line ?
        for (std::size_t other_line_index = 0; other_line_index < line_to_check_index; other_line_index++)
        {
            if (is_double_line_against_other_line<Point>(line_to_check_index, other_line_index, polylines, collapse_double_lines_sq_distance))
                return true;
        }
        for (std::size_t other_line_index = line_to_check_index + 1, nb_lines = polylines.lines.size(); other_line_index < nb_lines; other_line_index++)
        {
            if (is_double_line_against_other_line<Point>(line_to_check_index, other_line_index, polylines, collapse_double_lines_sq_distance))
                return true;
        }
        return false;
    }

    // computes the lines that are double lines because they are close to the input line
    template <typename Point>
    std::vector<std::size_t> get_double_lines_from_line(
        const std::size_t &line_to_check_index,
        const Polylines<Point> &polylines,
        const double &collapse_double_lines_sq_distance) const
    {
        // empty lines cannot make another line double
        if (polylines.lines[line_to_check_index].empty())
            return {};
        std::vector<std::size_t> close_lines;
        // test for every lines that is not the line to check: is the other line double beacuse of this line ?
        for (std::size_t other_line_index = 0; other_line_index < line_to_check_index; other_line_index++)
        {
            if (is_double_line_against_other_line<Point>(other_line_index, line_to_check_index, polylines, collapse_double_lines_sq_distance))
                close_lines.push_back(other_line_index);
        }
        for (std::size_t other_line_index = line_to_check_index + 1, nb_lines = polylines.lines.size(); other_line_index < nb_lines; other_line_index++)
        {
            if (is_double_line_against_other_line<Point>(other_line_index, line_to_check_index, polylines, collapse_double_lines_sq_distance))
                close_lines.push_back(other_line_index);
        }
        return close_lines;
    }

    // get the status of a line: is it double? if so, what type of double line is it?
    template <typename Point>
    Line_status compute_line_status(
        const std::size_t &line_index,
        const Polylines<Point> &polylines,
        const std::vector<std::vector<std::size_t>> &corner_to_lines,
        const Element_graph<Point> &dual_graph,
        const std::vector<typename Element_graph<Point>::Graph_element *> &line_index_to_dual_element,
        const double &collapse_double_lines_sq_distance) const
    {
        using Element_graph = Element_graph<Point>;
        using Graph_element = Element_graph::Graph_element;
        const bool &collapse_line =
            should_collapse_line<Point>(line_index, polylines, collapse_double_lines_sq_distance);
        if (!collapse_line)
            return Line_status::NOT_DOUBLE;
        // get type of line
        const std::vector<std::size_t> &line = polylines.lines[line_index];
        //   check isthmus: there is a corner with only one element (floating loops are considered isthmuses)
        if (line.size() < 2 || corner_to_lines[line.front()].size() == 1 || corner_to_lines[line.back()].size() == 1)
        {
            return Line_status::ISTHMUS;
        }
        //   check simple
        const Graph_element *dual_element = line_index_to_dual_element[line_index];
        std::shared_ptr<Element_graph> neighbors_subgraph =
            Element_graph_process::internal::make_sub_graph_from_elements(dual_graph, dual_element->neighbors);
        if (Element_graph_process::internal::is_graph_connected(*neighbors_subgraph))
            return Line_status::SIMPLE;
        return Line_status::OTHER;
    }

    // Updates the corner lines, merges the lines in the polylines, and merges in the dual graph
    template <typename Point>
    void merge_lines_following_order(
        const std::size_t &line_container_index,
        const bool &should_reverse_line_container,
        const std::size_t &line_to_push_index,
        const bool &should_reverse_line_to_push,
        Polylines<Point> &polylines,
        std::vector<std::vector<std::size_t>> &corner_to_lines,
        const std::vector<typename Element_graph<Point>::Graph_element *> &line_index_to_dual_element,
        std::set<std::size_t> &lines_updated_by_merge,
        std::set<std::size_t> &lines_deleted_by_merge,
        std::vector<Line_status> &lines_status) const
    {
        using Graph_element = Element_graph<Point>::Graph_element;

        std::vector<std::size_t> &line_container = polylines.lines[line_container_index];
        std::vector<std::size_t> &line_to_push = polylines.lines[line_to_push_index];
        // reverse if needed
        if (should_reverse_line_container)
            std::reverse(line_container.begin(), line_container.end());
        if (should_reverse_line_to_push)
            std::reverse(line_to_push.begin(), line_to_push.end());
        const std::size_t &common_corner = line_container.back();

        assert((void("Trying to merge lines that are not consecutives"), common_corner == line_to_push.front()));

        // clear corner (that is not a corner anymore)
        corner_to_lines[common_corner].clear();
        const std::size_t &other_corner_to_push = line_to_push.back();
        // update the other corner of the line to push : remove the cleared line and add the merged line (if it is not a loop)
        Polylines_process::internal::remove_line_from_corner_lines(line_to_push_index, other_corner_to_push, corner_to_lines);
        if (other_corner_to_push != line_container.front())
            corner_to_lines[other_corner_to_push].push_back(line_container_index);
        // merge the lines (do not duplicate the corner)
        line_container.insert(line_container.end(), std::next(line_to_push.begin()), line_to_push.end());
        line_to_push.clear();
        // merge in the dual graph: add all neighbors of the line to push to the line container, and disconnect the line to push
        Graph_element *element_container = line_index_to_dual_element[line_container_index];
        Graph_element *element_to_push = line_index_to_dual_element[line_to_push_index];
        std::set<Graph_element *> &element_container_neighbors = element_container->neighbors;
        for (auto &element_to_push_neighbor : element_to_push->neighbors)
        {
            std::set<Graph_element *> &neighbors_of_element_to_push_neighbor = element_to_push_neighbor->neighbors;
            neighbors_of_element_to_push_neighbor.erase(element_to_push);
            if (element_to_push_neighbor != element_container)
            {
                element_container_neighbors.insert(element_to_push_neighbor);
                neighbors_of_element_to_push_neighbor.insert(element_container);
            }
        }
        element_to_push->neighbors.clear();
        // mark the merged line
        lines_updated_by_merge.insert(line_container_index);
        lines_deleted_by_merge.insert(line_to_push_index);
        // update line status
        lines_status[line_to_push_index] = Line_status::NOT_DOUBLE;
    }

    // merges the lines from the common corner
    template <typename Point>
    void merge_lines_on_corner(
        const std::size_t &line_index_A,
        const std::size_t &line_index_B,
        const std::size_t &corner_index,
        Polylines<Point> &polylines,
        std::vector<std::vector<std::size_t>> &corner_to_lines,
        const std::vector<typename Element_graph<Point>::Graph_element *> &line_index_to_dual_element,
        std::set<std::size_t> &lines_updated_by_merge,
        std::set<std::size_t> &lines_deleted_by_merge,
        std::vector<Line_status> &lines_status) const
    {
        // do not merge a line into itself (possible input because of loops)
        if (line_index_A == line_index_B)
            return;
        std::vector<std::size_t> &line_A = polylines.lines[line_index_A];
        std::vector<std::size_t> &line_B = polylines.lines[line_index_B];
        // do not merge if a line is a loop
        if (line_A.front() == line_A.back() || line_B.front() == line_B.back())
            return;

        assert((void("Trying to merge line with one or less element"), line_A.size() > 1 && line_B.size() > 1));
        assert((void("Trying to merge lines that are not adjacent"), (line_A.front() == corner_index || line_A.back() == corner_index) && (line_B.front() == corner_index || line_B.back() == corner_index)));

        if (line_A.front() == corner_index)
        {
            if (line_A.front() == line_B.front())
            {
                // ( line_A.back() ) ----- ( line_A.front() [corner_index] line_B.front() ) ----- ( line_B.back() )
                // take the reverse of line_A and push line_B
                merge_lines_following_order(line_index_A, true, line_index_B, false, polylines, corner_to_lines, line_index_to_dual_element, lines_updated_by_merge, lines_deleted_by_merge, lines_status);
            }
            else if (line_A.front() == line_B.back())
            {
                // ( line_A.back() ) ----- ( line_A.front() [corner_index] line_B.back() ) ----- ( line_B.front() )
                // take line_B and push line_A
                merge_lines_following_order(line_index_B, false, line_index_A, false, polylines, corner_to_lines, line_index_to_dual_element, lines_updated_by_merge, lines_deleted_by_merge, lines_status);
            }
        }
        else // line_A.back() == corner_index
        {
            if (line_A.back() == line_B.front())
            {
                // ( line_A.front() ) ----- ( line_A.back() [corner_index] line_B.front() ) ----- ( line_B.back() )
                // take line_A and push line_B
                merge_lines_following_order(line_index_A, false, line_index_B, false, polylines, corner_to_lines, line_index_to_dual_element, lines_updated_by_merge, lines_deleted_by_merge, lines_status);
            }
            else if (line_A.back() == line_B.back())
            {
                // ( line_A.front() ) ----- ( line_A.back() [corner_index] line_B.back() ) ----- ( line_B.front() )
                // take line_A and push the reverse of line_B
                merge_lines_following_order(line_index_A, false, line_index_B, true, polylines, corner_to_lines, line_index_to_dual_element, lines_updated_by_merge, lines_deleted_by_merge, lines_status);
            }
        }
    }

    template <typename Point>
    void remove_line(
        const std::size_t &line_index,
        Polylines<Point> &polylines,
        std::vector<std::vector<std::size_t>> &corner_to_lines,
        const std::vector<typename Element_graph<Point>::Graph_element *> &line_index_to_dual_element,
        std::set<std::size_t> &lines_updated_by_merge,
        std::set<std::size_t> &lines_deleted_by_merge,
        std::vector<Line_status> &lines_status) const
    {
        using Element_graph = Element_graph<Point>;
        using Graph_element = Element_graph::Graph_element;

        std::vector<std::size_t> &line_to_delete = polylines.lines[line_index];
        if (line_to_delete.size() < 2)
        {
            line_to_delete.clear();
            return;
        }
        const std::size_t &front_corner = line_to_delete.front();
        const std::size_t &back_corner = line_to_delete.back();
        // update corner's adjacent lines
        Polylines_process::internal::remove_line_from_corner_lines(line_index, front_corner, corner_to_lines);
        Polylines_process::internal::remove_line_from_corner_lines(line_index, back_corner, corner_to_lines);
        // merge the neighbors lines if the corner disapear (there is only two adjacent lines)
        const std::vector<std::size_t> &front_corner_neighbor_lines = corner_to_lines[front_corner];
        const std::vector<std::size_t> &back_corner_neighbor_lines = corner_to_lines[back_corner];
        if (front_corner_neighbor_lines.size() == 2)
            merge_lines_on_corner(front_corner_neighbor_lines.front(), front_corner_neighbor_lines.back(), front_corner, polylines, corner_to_lines, line_index_to_dual_element, lines_updated_by_merge, lines_deleted_by_merge, lines_status);
        if (back_corner_neighbor_lines.size() == 2)
            merge_lines_on_corner(back_corner_neighbor_lines.front(), back_corner_neighbor_lines.back(), back_corner, polylines, corner_to_lines, line_index_to_dual_element, lines_updated_by_merge, lines_deleted_by_merge, lines_status);
        // collapse line
        line_to_delete.clear();
        // disconnect in dual_graph
        Graph_element *dual_element = line_index_to_dual_element[line_index];
        for (Graph_element *dual_element_neighbor : dual_element->neighbors)
        {
            dual_element_neighbor->neighbors.erase(dual_element);
        }
        dual_element->neighbors.clear();
        // update line status
        lines_status[line_index] = Line_status::NOT_DOUBLE;
    }

    // returns adjacent lines from a line without duplicates
    template <typename Point>
    std::vector<std::size_t> get_line_adjacent_lines(
        const std::size_t &line_index,
        const Polylines<Point> &polylines,
        const std::vector<std::vector<std::size_t>> &corner_to_lines) const
    {
        std::vector<std::size_t> adjacent_lines;
        const std::vector<std::size_t> &line = polylines.lines[line_index];
        const std::vector<std::size_t> &front_corner_adjacent_lines = corner_to_lines[line.front()];
        const std::vector<std::size_t> &back_corner_adjacent_lines = corner_to_lines[line.back()];
        adjacent_lines.reserve(front_corner_adjacent_lines.size() + back_corner_adjacent_lines.size());
        adjacent_lines.insert(adjacent_lines.end(), front_corner_adjacent_lines.begin(), front_corner_adjacent_lines.end());
        adjacent_lines.insert(adjacent_lines.end(), back_corner_adjacent_lines.begin(), back_corner_adjacent_lines.end());
        std::sort(adjacent_lines.begin(), adjacent_lines.end());
        auto last = std::unique(adjacent_lines.begin(), adjacent_lines.end());
        adjacent_lines.erase(last, adjacent_lines.end());
        return adjacent_lines;
    }

    template <typename Point>
    void collapse_line_to_middle_element(
        const std::size_t &line_index,
        Polylines<Point> &polylines,
        std::vector<std::vector<std::size_t>> &corner_to_lines,
        const std::vector<typename Element_graph<Point>::Graph_element *> &line_index_to_dual_element,
        std::vector<Line_status> &lines_status) const
    {
        using Element_graph = Element_graph<Point>;
        using Graph_element = Element_graph::Graph_element;

        std::vector<std::size_t> &line_to_delete = polylines.lines[line_index];
        if (line_to_delete.size() < 2)
            return;
        const std::size_t &middle_element_line_index = line_to_delete.size() / 2;
        const std::size_t &middle_element = line_to_delete[middle_element_line_index];
        std::vector<std::size_t> &middle_corner_neighbor_lines = corner_to_lines[middle_element];
        if (line_to_delete.size() > 2) // clear only if the middle element is not a corner
            middle_corner_neighbor_lines.clear();
        // update corners' adjacent lines
        const std::size_t &front_corner = line_to_delete.front();
        const std::size_t &back_corner = line_to_delete.back();
        Polylines_process::internal::remove_line_from_corner_lines(line_index, front_corner, corner_to_lines);
        Polylines_process::internal::remove_line_from_corner_lines(line_index, back_corner, corner_to_lines);
        //   connect the new corner to the adjacent lines.
        middle_corner_neighbor_lines = get_line_adjacent_lines(line_index, polylines, corner_to_lines);
        //   clear old corners
        if (middle_element != front_corner)
            corner_to_lines[front_corner].clear();
        if (middle_element != back_corner)
            corner_to_lines[back_corner].clear();
        // update the adjacent lines
        for (const std::size_t &adjacent_line_index : middle_corner_neighbor_lines)
        {
            std::vector<std::size_t> &adjacent_line = polylines.lines[adjacent_line_index];
            if (adjacent_line.front() == front_corner || adjacent_line.front() == back_corner)
                adjacent_line.front() = middle_element;
            if (adjacent_line.back() == front_corner || adjacent_line.back() == back_corner)
                adjacent_line.back() = middle_element;
        }
        // collapse line
        line_to_delete.clear();
        // in dual_graph, disconnect the line and connect neighbors between them
        Graph_element *dual_element = line_index_to_dual_element[line_index];
        std::set<Graph_element *> &dual_element_neighbors = dual_element->neighbors;
        for (Graph_element *dual_element_neighbor : dual_element_neighbors)
        {
            std::set<Graph_element *> &dual_element_neighbor_neighbors = dual_element_neighbor->neighbors;
            dual_element_neighbor_neighbors.erase(dual_element);
            for (Graph_element *other_dual_element_neighbor : dual_element_neighbors)
            {
                if (other_dual_element_neighbor != dual_element_neighbor)
                {
                    dual_element_neighbor_neighbors.insert(other_dual_element_neighbor);
                }
            }
        }
        dual_element_neighbors.clear();
        // update line status
        lines_status[line_index] = Line_status::NOT_DOUBLE;
    }

    template <typename Point>
    Polylines<Point> operator()(
        const Polylines<Point> &input_polylines) const
    {
        using Polylines = Polylines<Point>;
        using Element_graph = Element_graph<Point>;
        using Graph_element = Element_graph::Graph_element;
        using Element_iterator = Element_graph::Element_iterator;
        const double epsilon = 1.e-3; // Use a bigger radius to detect lines to update
        const double collapse_double_lines_sq_distance = collapse_double_lines_distance * collapse_double_lines_distance;

        Polylines polylines = input_polylines; // intiailize as a copy
        // For each corner, save the adjacent lines. Updated when polylines updates
        std::vector<std::vector<std::size_t>> corner_to_lines = Polylines_process::internal::get_corner_to_line_adjacency(polylines);

        std::vector<std::size_t> lines_isthmus_marked_deletion;
        std::vector<std::size_t> lines_simple_marked_deletion;
        std::vector<std::size_t> lines_other_marked_deletion;

        // Make the dual graph of lines to check the simple property on lines. The graph is updated on line deletion
        std::vector<Graph_element *> line_index_to_dual_element;
        std::shared_ptr<Element_graph> dual_graph_ptr = Element_graph_process::build_dual_graph_from_polylines(polylines, line_index_to_dual_element);
        Element_graph &dual_graph = *dual_graph_ptr; // updated on line deletion

        // Initialize lines_status to know the type of line (not double, isthmus, simple or other)
        std::vector<Line_status> lines_status;
        lines_status.resize(polylines.lines.size());
        for (std::size_t line_index = 0, nb_lines = polylines.lines.size(); line_index < nb_lines; line_index++)
        {
            lines_status[line_index] = compute_line_status<Point>(line_index, polylines, corner_to_lines, dual_graph, line_index_to_dual_element, collapse_double_lines_sq_distance);
        }

        bool collapsed = true;
        std::size_t iteration_i = 0;
        std::set<std::size_t> lines_to_delete;
        std::set<std::size_t> lines_to_update_status;
        std::set<std::size_t> lines_updated_by_merge;
        std::set<std::size_t> lines_deleted_by_merge;
        std::set<std::size_t> adjacent_lines; // lines adjacent to a deleted line

        while (collapsed)
        {
            helper::bLine("collapse iteration " + std::to_string(iteration_i));
            iteration_i++;

            lines_isthmus_marked_deletion.clear();
            lines_simple_marked_deletion.clear();
            lines_other_marked_deletion.clear();

            // Detect and classify lines to delete
            for (std::size_t line_index = 0, nb_lines = polylines.lines.size(); line_index < nb_lines; line_index++)
            {
                const Graph_element *dual_element = line_index_to_dual_element[line_index];
                const Line_status &line_status = lines_status[line_index];
                if (line_status == Line_status::NOT_DOUBLE)
                    continue;
                else if (line_status == Line_status::ISTHMUS) // always fill isthmus lines
                    lines_isthmus_marked_deletion.push_back(line_index);
                else if (lines_isthmus_marked_deletion.empty() && line_status == Line_status::SIMPLE) // fill simple lines if no isthmus lines have been found yet (to lower unused memory cost)
                    lines_simple_marked_deletion.push_back(line_index);
                else if (lines_isthmus_marked_deletion.empty() && lines_simple_marked_deletion.empty() && line_status == Line_status::OTHER) // fill other lines if no isthmus nor simple lines have been found yet (to lower unused memory cost)
                    lines_other_marked_deletion.push_back(line_index);
            }

            // Select the lines to remove this iteration (select in order : isthmus, simple or other lines)
            lines_to_delete.clear();
            Line_status line_deletion_status;
            if (!lines_isthmus_marked_deletion.empty())
            {
                line_deletion_status = Line_status::ISTHMUS;
                lines_to_delete.insert(lines_isthmus_marked_deletion.begin(), lines_isthmus_marked_deletion.end());
#ifdef VERBOSE
                std::cout << " remove " << lines_to_delete.size() << " isthmus lines";
#endif
            }
            else if (!lines_simple_marked_deletion.empty())
            {
                line_deletion_status = Line_status::SIMPLE;
                // delete only the longest line
                std::size_t longest_line;
                std::size_t max_length = 0;
                for (const std::size_t &line_index : lines_simple_marked_deletion)
                {
                    const std::vector<std::size_t> &line = polylines.lines[line_index];
                    if (line.size() > max_length)
                    {
                        max_length = line.size();
                        longest_line = line_index;
                    }
                }
                lines_to_delete.insert(longest_line);
#ifdef VERBOSE
                std::cout << " remove the longest simple line ";
#endif
            }
            else if (!lines_other_marked_deletion.empty())
            {
                line_deletion_status = Line_status::OTHER;
                // delete only one line
                lines_to_delete.insert(lines_other_marked_deletion.front());
#ifdef VERBOSE
                std::cout << " remove one other line ";
#endif
            }
            else
            {
#ifdef VERBOSE
                std::cout << " there are no double lines ";
#endif
                // no more lines to delete
                collapsed = false;
                helper::eLine();
                break;
            }

            {
                Polylines display_lines_to_delete;
                display_lines_to_delete.points = polylines.points;
                for (const std::size_t line_index : lines_to_delete)
                {
                    const std::vector<std::size_t>& line = polylines.lines[line_index];
                    display_lines_to_delete.lines.push_back(line);
                }
                Polylines_process::write_polylines<Point, Vector>(display_lines_to_delete, "tmp/lines_to_delete_iter_" + std::to_string(iteration_i) + "_not_double", false, false, false, true);
            }

            // Start deletion of lines
            lines_to_update_status.clear();
            lines_updated_by_merge.clear();
            lines_deleted_by_merge.clear();
            adjacent_lines.clear();
            // Before collapsing the lines, collect the indices of the adjacent lines, and the lines that are double because of a deleted line.
            // they will be used to determine which lines to update
            for (const std::size_t &line_index : lines_to_delete)
            {
                const std::vector<std::size_t> &deleted_line_adjacent_lines = get_line_adjacent_lines(line_index, polylines, corner_to_lines);
                adjacent_lines.insert(deleted_line_adjacent_lines.begin(), deleted_line_adjacent_lines.end());
                lines_to_update_status.insert(deleted_line_adjacent_lines.begin(), deleted_line_adjacent_lines.end());
                // Conservative check: add an epsilon to the distance
                const std::vector<std::size_t> &near_lines = get_double_lines_from_line<Point>(line_index, polylines, collapse_double_lines_sq_distance + epsilon);
                lines_to_update_status.insert(near_lines.begin(), near_lines.end());
            }
            // Delete or collapse lines. Updates polylines, corner_to_lines, dual_graph, and lines_status.
            // marks the lines updated or deleted by the merge operation
            if (line_deletion_status == Line_status::ISTHMUS)
            {
                for (const std::size_t &line_index : lines_to_delete)
                {
                    if (lines_updated_by_merge.count(line_index) == 0) // do not remove the line if it has been updated by a merge operation
                        remove_line(line_index, polylines, corner_to_lines, line_index_to_dual_element, lines_updated_by_merge, lines_deleted_by_merge, lines_status);
                }
            }
            else if (line_deletion_status == Line_status::SIMPLE)
            {
                remove_line(*(lines_to_delete.begin()), polylines, corner_to_lines, line_index_to_dual_element, lines_updated_by_merge, lines_deleted_by_merge, lines_status);
            }
            else // OTHER
            {
                collapse_line_to_middle_element(*(lines_to_delete.begin()), polylines, corner_to_lines, line_index_to_dual_element, lines_status);
            }

            // Update the status of:
            //   - lines adjacent to a deleted line (already collected)
            //   - lines that are double because of a deleted line (already collected)
            //   - lines that are double because of an adjacent line
            for (const std::size_t &adjacent_line_index : adjacent_lines)
            {
                // double lines from deleted lines are already collected
                if (lines_to_delete.count(adjacent_line_index) > 0)
                    continue;
                // Conservative check: add an epsilon to the distance
                const std::vector<std::size_t> &near_lines = get_double_lines_from_line<Point>(adjacent_line_index, polylines, collapse_double_lines_sq_distance + epsilon);
                lines_to_update_status.insert(near_lines.begin(), near_lines.end());
            }
            // Update lines_status.
            for (const std::size_t &line_to_update : lines_to_update_status)
            {
                // lines_to_update_status may contain deleted of merged lines.
                if (lines_to_delete.count(line_to_update) > 0 || lines_deleted_by_merge.count(line_to_update) > 0)
                    continue;
                lines_status[line_to_update] = compute_line_status<Point>(line_to_update, polylines, corner_to_lines, dual_graph, line_index_to_dual_element, collapse_double_lines_sq_distance);
            }

            helper::eLine();
        }

        return Polylines_process::make_polylines_without_floating_points_and_empty_lines(polylines);
    }
};

#endif // POLYLINE_REGULIZER_H