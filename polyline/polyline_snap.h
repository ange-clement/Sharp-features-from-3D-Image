#ifndef POLYLINE_SNAP_H
#define POLYLINE_SNAP_H

#include "Polylines.h"
#include "polyline_process.h"
#include "../feature_graph/polyline_regularizer.h"

#include <CGAL/Distance_3/Segment_3_Segment_3.h>

namespace Polylines_process
{
    namespace internal
    {
        // the point is not removed from the polyline
        template <typename Point>
        void snap_corners(
            const std::size_t corner_container_index,
            const std::size_t corner_to_push_index,
            Polylines<Point> &polylines,
            std::vector<std::vector<std::size_t>> &corner_to_lines)
        {
            std::cout << "  Snap corner " << corner_to_push_index << " into " << corner_container_index << std::endl;
            std::vector<std::size_t> &corner_container_incident_lines = corner_to_lines[corner_container_index];
            std::vector<std::size_t> &corner_to_push_incident_lines = corner_to_lines[corner_to_push_index];
            std::cout << "  got incident lines" << std::endl;
            for (const std::size_t &corner_to_push_incident_line_index : corner_to_push_incident_lines)
            {
                std::vector<std::size_t> &incident_line = polylines.lines[corner_to_push_incident_line_index];
                if (incident_line.size() < 2)
                {
                    incident_line.clear();
                    continue;
                }
                std::cout << "  incident_line " << corner_to_push_incident_line_index << " was " << std::flush;
                for (const std::size_t &p_id : incident_line)
                    std::cout << " " << p_id << std::flush;
                std::cout << std::endl;
                if (incident_line.front() == corner_to_push_index)
                    incident_line.front() = corner_container_index;
                else
                    std::cout << "no, " << incident_line.front() << " != " << corner_to_push_index << std::endl;
                if (incident_line.back() == corner_to_push_index)
                    incident_line.back() = corner_container_index;
                else
                    std::cout << "no, " << incident_line.back() << " != " << corner_to_push_index << std::endl;
                std::cout << "  incident_line " << corner_to_push_incident_line_index << " is " << std::flush;
                for (const std::size_t &p_id : incident_line)
                    std::cout << " " << p_id << std::flush;
                std::cout << std::endl;
                if (incident_line.size() == 2 && incident_line.front() == incident_line.back())
                {
                    remove_line_from_corner_lines(corner_to_push_incident_line_index, corner_container_index, corner_to_lines);
                    incident_line.clear();
                    continue;
                }
                corner_container_incident_lines.push_back(corner_to_push_incident_line_index);
            }
            corner_to_push_incident_lines.clear();
            // remove duplicates
            std::sort(corner_container_incident_lines.begin(), corner_container_incident_lines.end());
            auto last = std::unique(corner_container_incident_lines.begin(), corner_container_incident_lines.end());
            corner_container_incident_lines.erase(last, corner_container_incident_lines.end());
        }

        template <typename Point>
        bool check_and_snap_corners(
            const Point &corner_container_point,
            const std::size_t &corner_container_index,
            const std::size_t &corner_to_push_index,
            const double &collapse_under_distance_sq,
            const std::vector<bool> &processed_corners,
            Polylines<Point> &polylines,
            std::vector<std::vector<std::size_t>> &corner_to_lines)
        {
            if (processed_corners[corner_to_push_index])
                return false;
            if (CGAL::squared_distance(corner_container_point, polylines.points[corner_to_push_index]) < collapse_under_distance_sq)
            {
                snap_corners(corner_container_index, corner_to_push_index, polylines, corner_to_lines);
                return true;
            }
            return false;
        }

        // check and snap to corners in the line after the current line
        template <typename Point>
        void snap_corner_to_corners_after(
            const std::size_t &line_index,
            const std::size_t corner_to_check,
            const double &collapse_under_distance_sq,
            Polylines<Point> &polylines,
            std::vector<std::vector<std::size_t>> &corner_to_lines,
            std::vector<bool> &processed_corners)
        {
            if (processed_corners[corner_to_check])
                return;
            processed_corners[corner_to_check] = true;
            // find close corner
            const Point &corner_point = polylines.points[corner_to_check];
            for (std::size_t other_line_index = line_index + 1, number_of_lines = polylines.lines.size(); other_line_index < number_of_lines; other_line_index++)
            {
                const std::vector<std::size_t> &other_line = polylines.lines[other_line_index];
                if (other_line.size() == 0)
                    continue;
                const std::size_t &other_line_front_corner = other_line.front();
                check_and_snap_corners(corner_point, corner_to_check, other_line_front_corner, collapse_under_distance_sq, processed_corners, polylines, corner_to_lines);
                if (other_line.size() == 1)
                    continue;
                const std::size_t &other_line_back_corner = other_line.back();
                check_and_snap_corners(corner_point, corner_to_check, other_line_back_corner, collapse_under_distance_sq, processed_corners, polylines, corner_to_lines);
            }
        }

        template <typename Point>
        void snap_corner_to_corners(
            const double &collapse_under_distance,
            Polylines<Point> &polylines,
            std::vector<std::vector<std::size_t>> &corner_to_lines)
        {
            // TODO Optimization: KD tree of corners
            double collapse_under_distance_sq = collapse_under_distance * collapse_under_distance;
            std::vector<bool> processed_corners(polylines.points.size(), false);
            // iterate on corner, and snap on corners after it
            for (std::size_t line_index = 0, number_of_lines = polylines.lines.size(); line_index < number_of_lines; line_index++)
            {
                const std::vector<std::size_t> &line = polylines.lines[line_index];
                if (line.size() == 0)
                    continue;
                // snap the front corner to corners on lines after
                const std::size_t &front_corner = line.front();
                snap_corner_to_corners_after(line_index, front_corner, collapse_under_distance_sq, polylines, corner_to_lines, processed_corners);
                if (line.size() <= 1)
                    continue;
                // snap the back corner to the front corner
                const std::size_t &back_corner = line.back();
                if (front_corner == back_corner)
                    continue;
                const Point &front_corner_point = polylines.points[front_corner];
                if (!check_and_snap_corners(front_corner_point, front_corner, back_corner, collapse_under_distance_sq, processed_corners, polylines, corner_to_lines))
                    snap_corner_to_corners_after(line_index, back_corner, collapse_under_distance_sq, polylines, corner_to_lines, processed_corners);
            }
        }

        // May invalidate the line iterators
        template <typename K, typename Point>
        void snap_current_corner_to_edge(
            const std::size_t &line_index, // can be -1 if the corner does not belong to a line
            const std::size_t &corner_to_snap,
            const std::size_t &other_line_index, // can be the same as line_index
            const std::size_t &edge_start_index,
            Polylines<Point> &polylines,
            std::vector<std::vector<std::size_t>> &corner_to_lines)
        {
            std::cout << "new line. Corner to snap is " << corner_to_snap << std::endl;
            std::cout << "edge to snap is on line " << other_line_index << " with index " << edge_start_index << " [" << polylines.lines[other_line_index][edge_start_index] << " " << polylines.lines[other_line_index][edge_start_index+1] << "]" << std::endl;
            std::cout << "other line size : " << polylines.lines[other_line_index].size() << std::endl;
            // split the compared line at the current edge and join on the corner
            //   create a line with the corner and the rest of the compared line
            std::size_t new_line_index = polylines.lines.size();
            std::cout << "new line index is " << new_line_index << std::endl;
            polylines.lines.resize(new_line_index + 1); // may invalidate the line iterators
            std::cout << "number of lines is " << polylines.lines.size() << std::endl;
            std::vector<std::size_t> &new_line = polylines.lines[new_line_index];
            std::vector<std::size_t> &other_line = polylines.lines[other_line_index]; // get the other line after adding the new line
            std::cout << new_line_index << " add " << corner_to_snap << std::endl;
            new_line.push_back(corner_to_snap);
            for (std::size_t id = edge_start_index + 1; id < other_line.size(); id++)
                std::cout << new_line_index << " add " << other_line[id] << std::endl;
            new_line.insert(new_line.end(), other_line.begin() + edge_start_index + 1, other_line.end());
            std::cout << "new line is " << std::flush;
            for (const std::size_t &p_id : new_line)
                std::cout << " " << p_id << std::flush;
            std::cout << " end" << std::endl;
            std::cout << "other line" << std::endl;
            //   on the compared line, remove the rest of the line and add the corner
            other_line.resize(edge_start_index + 2);
            other_line.back() = corner_to_snap;
            std::cout << "corner adjacency" << std::endl;
            // update the corner adjacency
            std::vector<std::size_t> &corner_adjacency = corner_to_lines[corner_to_snap];
            if (line_index != other_line_index) // snap the corner on a different line
                corner_adjacency.push_back(other_line_index);
            corner_adjacency.push_back(new_line_index);
            std::cout << "back corner adjacency" << std::endl;
            //   update the other corner because the line has changed index
            const std::size_t &back_corner = new_line.back();
            if (corner_to_snap == back_corner) // if a loop is created, the adjacency is already valid
                return;
            std::vector<std::size_t> &back_corner_adjacency = corner_to_lines[back_corner];
            //   replace the index, or add if the other line was a loop
            if (back_corner == other_line.front())
                back_corner_adjacency.push_back(new_line_index);
            else
            {
                for (std::size_t &adjacent_line : back_corner_adjacency)
                {
                    if (adjacent_line == other_line_index)
                    {
                        adjacent_line = new_line_index;
                        break; // there are no repeats, so there is maximum one element to change
                    }
                }
            }
        }

        template <typename K, typename Point>
        void snap_current_corner_to_line(
            const std::size_t &line_index,
            const std::size_t &corner_to_check,
            const Point &corner_point,
            const std::size_t &other_line_index, // can be the same as line_index
            const double &collapse_under_distance_sq,
            Polylines<Point> &polylines,
            std::vector<std::vector<std::size_t>> &corner_to_lines)
        {
            using Segment_3 = typename K::Segment_3;
            std::vector<std::size_t> &other_line = polylines.lines[other_line_index];
            if (other_line.size() < 2)
                return;
            std::size_t start_edge_start_index = (other_line.front() == corner_to_check) ? 1 : 0;
            std::size_t end_edge_start_index = (other_line.back() == corner_to_check) ? other_line.size() - 2 : other_line.size() - 1;
            for (std::size_t edge_start_index = start_edge_start_index; edge_start_index < end_edge_start_index; edge_start_index++)
            {
                const std::size_t &edge_point_0_index = other_line[edge_start_index];
                const std::size_t &edge_point_1_index = other_line[edge_start_index + 1];
                Segment_3 segment(polylines.points[edge_point_0_index], polylines.points[edge_point_1_index]);
                if (CGAL::squared_distance(corner_point, segment) < collapse_under_distance_sq)
                {
                    snap_current_corner_to_edge<K>(line_index, corner_to_check, other_line_index, edge_start_index, polylines, corner_to_lines); // line pointers may be invalidated
                    return; // the snap will cut other_line, so there are no remaining edges to compare
                }
            }
        }

        template <typename K, typename Point>
        void snap_current_corner_to_all_edges(
            const std::size_t &line_index,
            const std::size_t &corner_to_check,
            const double &collapse_under_distance_sq,
            Polylines<Point> &polylines,
            std::vector<std::vector<std::size_t>> &corner_to_lines,
            std::vector<bool> &processed_corners)
        {
            if (processed_corners[corner_to_check])
                return;
            processed_corners[corner_to_check] = true;
            const Point &corner_point = polylines.points[corner_to_check];
            for (std::size_t other_line_index = 0; other_line_index < polylines.lines.size(); other_line_index++) // the number of line can increase
                snap_current_corner_to_line<K>(line_index, corner_to_check, corner_point, other_line_index, collapse_under_distance_sq, polylines, corner_to_lines);  // line pointers may be invalidated
        }

        template <typename K, typename Point>
        void snap_corner_to_edges(
            const double &collapse_under_distance,
            Polylines<Point> &polylines,
            std::vector<std::vector<std::size_t>> &corner_to_lines)
        {
            // TODO Optimization: KD tree of lines so that far lines are not checked
            using Segment_3 = typename K::Segment_3;

            double collapse_under_distance_sq = collapse_under_distance * collapse_under_distance;
            std::vector<bool> processed_corners(polylines.points.size(), false);

            // iterate on corners, and snap edges
            for (std::size_t line_index = 0; line_index < polylines.lines.size(); line_index++) // the number of line can increase
            {
                {
                    const std::vector<std::size_t> &line = polylines.lines[line_index];
                    if (line.size() == 0)
                        continue;
                    snap_current_corner_to_all_edges<K>(line_index, line.front(), collapse_under_distance_sq, polylines, corner_to_lines, processed_corners);  // line pointers may be invalidated
                }
                {
                    const std::vector<std::size_t> &line = polylines.lines[line_index]; // need to fetch the line again as it may be invalid
                    if (line.size() == 1)
                        continue;
                    snap_current_corner_to_all_edges<K>(line_index, line.back(), collapse_under_distance_sq, polylines, corner_to_lines, processed_corners);  // line pointers may be invalidated
                }
            }
        }

        template <typename K, typename Point>
        void check_and_snap_corner_adjacent_point(
            const std::size_t &corner_adjacent_point_in_line_index,
            const std::size_t &corner_point_index,
            const Point &corner_point,
            const double &collapse_under_distance_sq,
            std::vector<std::size_t> &incident_line,
            Polylines<Point> &polylines,
            std::vector<std::vector<std::size_t>> &corner_to_lines)
        {
            const Point &corner_adjacent_point = polylines.points[incident_line[corner_adjacent_point_in_line_index]];
            if (CGAL::squared_distance(corner_adjacent_point, corner_point) >= collapse_under_distance_sq)
                return;
            std::cout << "snap adjacent point too close : " << incident_line[corner_adjacent_point_in_line_index] << std::endl;
            std::cout << "incident_line was " << std::flush;
            for (const std::size_t &p_id : incident_line)
                std::cout << " " << p_id << std::flush;
            std::cout << std::endl;
            if (incident_line.size() == 2)
                snap_corners(corner_point_index, incident_line[corner_adjacent_point_in_line_index], polylines, corner_to_lines);
            else
                incident_line.erase(incident_line.begin() + corner_adjacent_point_in_line_index);
            std::cout << "incident_line is " << std::flush;
            for (const std::size_t &p_id : incident_line)
                std::cout << " " << p_id << std::flush;
            std::cout << std::endl;
        }

        template <typename K, typename Point>
        void snap_current_edge_to_edge(
            const std::size_t &line_index,
            const std::size_t &edge_start_index,
            const typename K::Segment_3 &edge_segment,
            const std::size_t &other_line_index, // can be the same as line_index
            const std::size_t &other_edge_start_index,
            const typename K::Segment_3 &other_edge_segment,
            const CGAL::Distance_3::internal::Segment_3_Segment_3_Result<K> &closest_point_result,
            const double &collapse_under_distance_sq,
            const K& k,
            Polylines<Point> &polylines,
            std::vector<std::vector<std::size_t>> &corner_to_lines)
        {
            std::cout << "        snap line " << line_index << " edge at " << edge_start_index << " [" << polylines.lines[line_index][edge_start_index] << " " << polylines.lines[line_index][edge_start_index+1] << "]"
                      << " to line " << other_line_index << " edge at " << other_edge_start_index << " [" << polylines.lines[other_line_index][other_edge_start_index] << " " << polylines.lines[other_line_index][other_edge_start_index+1] << "]" << std::endl;
            {
                std::cout << "BEFORE SNAP STATUS :" << std::endl;
                for (std::size_t l_id = 0; l_id < polylines.lines.size(); l_id++)
                {
                    std::cout << l_id << " :" << std::flush;
                    const std::vector<std::size_t> &line = polylines.lines[l_id];
                    for (const std::size_t &p_id : line)
                        std::cout << " " << p_id;
                    std::cout << std::endl;
                }
                for (std::size_t p_id = 0; p_id < corner_to_lines.size(); p_id++)
                {
                    if (corner_to_lines[p_id].empty())
                        continue;
                    std::cout << p_id << " -->" << std::flush;
                    for (const std::size_t& adj_l_id : corner_to_lines[p_id])
                        std::cout << " " << adj_l_id << std::flush;
                    std::cout << std::endl;
                }
            }
            // Create new point
            std::cout << "new point" << std::endl;
            const Point closest_point_on_edge = edge_segment.source() + closest_point_result.x*edge_segment.to_vector();
            const Point closest_point_on_other_edge = other_edge_segment.source() + closest_point_result.y*other_edge_segment.to_vector();
            const Point &new_point = CGAL::midpoint(closest_point_on_edge, closest_point_on_other_edge);
            const std::size_t &new_point_index = polylines.points.size();
            polylines.points.push_back(new_point);
            corner_to_lines.resize(new_point_index + 1);
            // Snap edges to the new point (snap first to the other edge because it is after in the case where the lines are the same)
            std::cout << "snap other line" << std::endl;
            snap_current_corner_to_edge<K>(static_cast<std::size_t>(-1), new_point_index, other_line_index, other_edge_start_index, polylines, corner_to_lines); // line pointers may be invalidated
            std::cout << "snap line" << std::endl;
            snap_current_corner_to_edge<K>(other_line_index, new_point_index, line_index, edge_start_index, polylines, corner_to_lines); // line pointers may be invalidated
            {
                std::cout << "after snap new corner " << new_point_index << " :" << std::endl;
                for (const std::size_t &incident_line_index : corner_to_lines[new_point_index])
                {
                    std::cout << "line " << incident_line_index << " :" << std::flush;
                    for (const std::size_t &p_id : polylines.lines[incident_line_index])
                        std::cout << " " << p_id << std::flush;
                    std::cout << std::endl;
                }
            }
            // Snap the new corner adjacent point that are too close
            for (const std::size_t &incident_line_index : corner_to_lines[new_point_index])
            {
                std::vector<std::size_t> &incident_line = polylines.lines[incident_line_index];
                if (incident_line.empty())
                    continue;
                if (incident_line.front() == new_point_index)
                    check_and_snap_corner_adjacent_point<K>(1, new_point_index, new_point, collapse_under_distance_sq, incident_line, polylines, corner_to_lines);
                if (incident_line.size() > 1 && incident_line.back() == new_point_index)
                    check_and_snap_corner_adjacent_point<K>(incident_line.size() - 2, new_point_index, new_point, collapse_under_distance_sq, incident_line, polylines, corner_to_lines);
            }

            {
                std::cout << "after snap too close to " << new_point_index << " :" << std::endl;
                for (const std::size_t &incident_line_index : corner_to_lines[new_point_index])
                {
                    std::cout << "line " << incident_line_index << " :" << std::flush;
                    for (const std::size_t &p_id : polylines.lines[incident_line_index])
                        std::cout << " " << p_id << std::flush;
                    std::cout << std::endl;
                }
            }
            {
                std::cout << "AFTER SNAP STATUS :" << std::endl;
                for (std::size_t l_id = 0; l_id < polylines.lines.size(); l_id++)
                {
                    std::cout << l_id << " :" << std::flush;
                    const std::vector<std::size_t> &line = polylines.lines[l_id];
                    for (const std::size_t &p_id : line)
                        std::cout << " " << p_id;
                    std::cout << std::endl;
                }
                for (std::size_t p_id = 0; p_id < corner_to_lines.size(); p_id++)
                {
                    if (corner_to_lines[p_id].empty())
                        continue;
                    std::cout << p_id << " -->" << std::flush;
                    for (const std::size_t& adj_l_id : corner_to_lines[p_id])
                        std::cout << " " << adj_l_id << std::flush;
                    std::cout << std::endl;
                }
            }
        }

        template <typename K, typename Point>
        void snap_current_edge_to_line(
            const std::size_t &line_index,
            const std::size_t &edge_start_index,
            const std::size_t &other_line_index, // can be the same as line_index
            const double &collapse_under_distance_sq,
            const K& k,
            Polylines<Point> &polylines,
            std::vector<std::vector<std::size_t>> &corner_to_lines)
        {
            using Segment_3 = typename K::Segment_3;
            using Segment_3_Segment_3_Result = CGAL::Distance_3::internal::Segment_3_Segment_3_Result<K>;

            // std::cout << "      check other line  " << other_line_index << std::endl;

            // snap the start corner edge to the other edges of the line
            std::vector<std::size_t> &line = polylines.lines[line_index];
            if (line.size() < 2)
                return;
            std::vector<std::size_t> &other_line = polylines.lines[other_line_index];
            if (other_line.size() < 2)
                return;
            const std::size_t& edge_start_point = line[edge_start_index];
            const std::size_t& edge_end_point = line[edge_start_index + 1];
            const Segment_3 edge_segment(polylines.points[edge_start_point], polylines.points[edge_end_point]);
            std::size_t other_edge_start_index = (line_index == other_line_index) ? edge_start_index + 2 : 0; // on the same line, start on the second next edge
            for (; other_edge_start_index < other_line.size() - 1; other_edge_start_index++) // the size of the line can decrease in a snap, avoid possible underflow
            {
                const std::size_t& other_edge_start_point = other_line[other_edge_start_index];
                const std::size_t& other_edge_end_point = other_line[other_edge_start_index + 1];
                if (other_edge_start_point == edge_start_point || other_edge_start_point == edge_end_point
                    || other_edge_end_point == edge_start_point || other_edge_end_point == edge_end_point) // does the two edges share a point ?
                    continue;
                const Segment_3 other_edge_segment(polylines.points[other_line[other_edge_start_index]], polylines.points[other_line[other_edge_start_index+1]]);
                Segment_3_Segment_3_Result closest_point_result = CGAL::Distance_3::internal::squared_distance(edge_segment, other_edge_segment, k);
                if (closest_point_result.squared_distance >= collapse_under_distance_sq)
                    continue;
                // the two edges are too close : collapse
                snap_current_edge_to_edge(
                    line_index, edge_start_index, edge_segment,
                    other_line_index, other_edge_start_index, other_edge_segment,
                    closest_point_result, collapse_under_distance_sq, k,
                    polylines, corner_to_lines); // line pointers may be invalidated
                return; // the snap will cut the line, so there are no remaining edges to compare
            }
        }

        template <typename K, typename Point>
        void snap_current_edge_to_edges_after(
            const std::size_t &line_index,
            const std::size_t &edge_start_index,
            const double &collapse_under_distance_sq,
            const K& k,
            Polylines<Point> &polylines,
            std::vector<std::vector<std::size_t>> &corner_to_lines)
        {
            std::cout << "    check edge " << edge_start_index << std::endl;
            // snap the current edge to edges that are after
            for (std::size_t other_line_index = line_index; other_line_index < polylines.lines.size(); other_line_index++) // the number of line can increase
                snap_current_edge_to_line(line_index, edge_start_index, other_line_index, collapse_under_distance_sq, k, polylines, corner_to_lines); // line pointers may be invalidated
        }

        template <typename K, typename Point>
        void snap_edge_to_edges(
            const double &collapse_under_distance,
            Polylines<Point> &polylines,
            std::vector<std::vector<std::size_t>> &corner_to_lines)
        {
            // TODO Optimization: KD tree of lines so that far lines are not checked
            using Segment_3 = typename K::Segment_3;

            double collapse_under_distance_sq = collapse_under_distance * collapse_under_distance;

            const K k; // for Segment_3 to Segment_3 distance

            // snap edges (with different source and target points) that are too close
            for (std::size_t line_index = 0; line_index < polylines.lines.size(); line_index++) // the number of line can increase
            {
                std::cout << "  check line " << line_index << std::endl;
                if (polylines.lines[line_index].size() < 2)
                    continue;
                for (std::size_t edge_start_index = 0; polylines.lines[line_index].size() >= 1 && edge_start_index < polylines.lines[line_index].size()-1; edge_start_index++) // the size of the line can decrease in a snap, avoid possible underflow
                    snap_current_edge_to_edges_after(line_index, edge_start_index, collapse_under_distance_sq, k, polylines, corner_to_lines);  // line pointers may be invalidated
            }
        }
    } // namespace internal

    // collapse polyline elements with short distance
    // corners are protected
    // may leave floating points and empty lines
    template <typename K, typename Point>
    void snap_lines(
        const double &collapse_under_distance,
        Polylines<Point> &polylines)
    {
        {
            for (std::size_t l_id = 0; l_id < polylines.lines.size(); l_id++)
            {
                std::cout << l_id << " :" << std::flush;
                const std::vector<std::size_t> &line = polylines.lines[l_id];
                for (const std::size_t &p_id : line)
                    std::cout << " " << p_id;
                std::cout << std::endl;
            }
        }

        Polylines_process::write_polylines<Point, typename K::Vector_3>(polylines, "tmp/snap_init");

        std::cout << "corner adjacency" << std::endl;
        std::vector<std::vector<std::size_t>> corner_to_lines = internal::get_corner_to_line_adjacency(polylines);

        {
            for (std::size_t l_id = 0; l_id < polylines.lines.size(); l_id++)
            {
                std::cout << l_id << " :" << std::flush;
                const std::vector<std::size_t> &line = polylines.lines[l_id];
                for (const std::size_t &p_id : line)
                    std::cout << " " << p_id;
                std::cout << std::endl;
            }
            for (std::size_t p_id = 0; p_id < corner_to_lines.size(); p_id++)
            {
                if (corner_to_lines[p_id].empty())
                    continue;
                std::cout << p_id << " -->" << std::flush;
                for (const std::size_t& adj_l_id : corner_to_lines[p_id])
                    std::cout << " " << adj_l_id << std::flush;
                std::cout << std::endl;
            }
        }

        std::cout << "snap corner to corner" << std::endl;
        internal::snap_corner_to_corners(collapse_under_distance, polylines, corner_to_lines);
        Polylines_process::write_polylines<Point, typename K::Vector_3>(polylines, "tmp/snap_corner_corner");

        {
            for (std::size_t l_id = 0; l_id < polylines.lines.size(); l_id++)
            {
                std::cout << l_id << " :" << std::flush;
                const std::vector<std::size_t> &line = polylines.lines[l_id];
                for (const std::size_t &p_id : line)
                    std::cout << " " << p_id;
                std::cout << std::endl;
            }
            for (std::size_t p_id = 0; p_id < corner_to_lines.size(); p_id++)
            {
                if (corner_to_lines[p_id].empty())
                    continue;
                std::cout << p_id << " -->" << std::flush;
                for (const std::size_t& adj_l_id : corner_to_lines[p_id])
                    std::cout << " " << adj_l_id << std::flush;
                std::cout << std::endl;
            }
        }

        std::cout << "snap corner to edges" << std::endl;
        internal::snap_corner_to_edges<K>(collapse_under_distance, polylines, corner_to_lines);
        Polylines_process::write_polylines<Point, typename K::Vector_3>(polylines, "tmp/snap_corner_edges");

        {
            for (std::size_t l_id = 0; l_id < polylines.lines.size(); l_id++)
            {
                std::cout << l_id << " :" << std::flush;
                const std::vector<std::size_t> &line = polylines.lines[l_id];
                for (const std::size_t &p_id : line)
                    std::cout << " " << p_id;
                std::cout << std::endl;
            }
            for (std::size_t p_id = 0; p_id < corner_to_lines.size(); p_id++)
            {
                if (corner_to_lines[p_id].empty())
                    continue;
                std::cout << p_id << " -->" << std::flush;
                for (const std::size_t& adj_l_id : corner_to_lines[p_id])
                    std::cout << " " << adj_l_id << std::flush;
                std::cout << std::endl;
            }
        }

        std::cout << "snap edges to edges" << std::endl;
        internal::snap_edge_to_edges<K>(collapse_under_distance, polylines, corner_to_lines);
        Polylines_process::write_polylines<Point, typename K::Vector_3>(polylines, "tmp/snap_edges_edges");

        {
            for (std::size_t l_id = 0; l_id < polylines.lines.size(); l_id++)
            {
                std::cout << l_id << " :" << std::flush;
                const std::vector<std::size_t> &line = polylines.lines[l_id];
                for (const std::size_t &p_id : line)
                    std::cout << " " << p_id;
                std::cout << std::endl;
            }
            for (std::size_t p_id = 0; p_id < corner_to_lines.size(); p_id++)
            {
                if (corner_to_lines[p_id].empty())
                    continue;
                std::cout << p_id << " -->" << std::flush;
                for (const std::size_t& adj_l_id : corner_to_lines[p_id])
                    std::cout << " " << adj_l_id << std::flush;
                std::cout << std::endl;
            }
        }

        polylines = Polylines_process::make_polylines_without_floating_points_and_empty_lines(polylines);

        // collapse double lines because snap can create simple lines
        Double_lines_regularizer<typename K::Vector_3> regularization_functor;
        regularization_functor.collapse_double_lines_distance = collapse_under_distance;
        polylines = regularization_functor(polylines);
    }
} // namespace Polylines_process

#endif