
#ifndef POLYLINES_H
#define POLYLINES_H

#include <vector>

#include "../helper/mesh_generator.h"

template <typename Point>
struct Polylines
{
    std::vector<Point> points;
    std::vector<std::vector<std::size_t>> lines;

    void clear()
    {
        points.clear();
        lines.clear();
    }

    // assumes corners points are exactly in the same coords
    void insert_from_point_format(const std::vector<std::vector<Point>> &point_poly)
    {
        lines.reserve(point_poly.size());
        std::map<Point, std::size_t> inserted_points_ids;
        for (const std::vector<Point> &point_line : point_poly)
        {
            lines.resize(lines.size() + 1);
            std::vector<std::size_t> &line = lines.back();
            line.reserve(point_line.size());
            for (const Point &p : point_line)
            {
                const auto &[iterator, inserted] = inserted_points_ids.emplace(p, points.size());
                if (inserted)
                {
                    points.push_back(p);
                }
                line.push_back(iterator->second);
            }
        }
    }

    std::vector<std::vector<Point>> make_point_format() const
    {
        std::vector<std::vector<Point>> point_polyline;
        point_polyline.reserve(lines.size());
        for (const std::vector<std::size_t> &line : lines)
        {
            std::vector<Point> point_line;
            point_line.reserve(line.size());
            for (const std::size_t &p_id : line)
            {
                const Point &p = points[p_id];
                point_line.push_back(p);
            }
            point_polyline.push_back(point_line);
        }
        return point_polyline;
    }

    Polylines<Point> clone() const
    {
        Polylines<Point> clone;
        clone.points = points;
        clone.lines = lines;
        return clone;
    }
};

#endif // POLYLINES_H