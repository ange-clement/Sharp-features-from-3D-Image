#ifndef GEMETRY_HELPER_H
#define GEMETRY_HELPER_H

#include <CGAL/number_utils.h>
#include <CGAL/Kernel/global_functions_3.h>

#include "mesh_generator.h" // get_rgb_color_heat_map

namespace helper
{

    namespace internal
    {
        template <typename Point, typename Value>
        void write_point_set_values_min_max(
            const std::string& filename,
            const std::vector<Point>& points,
            const std::vector<Value>& values,
            const Value& min, const Value& max,
            bool invert_values = false)
        {
            Value inv_range = 1.0 / (max - min);
            std::ofstream output_point_set(filename, std::ios_base::binary);
            output_point_set << "ply\nformat ascii 1.0\n"
                            << "element vertex " << points.size() << "\n"
                            << "property double x\nproperty double y\nproperty double z\n"
                            << "property double nx\nproperty double ny\nproperty double nz\n"
                            << "property uchar red\nproperty uchar green\nproperty uchar blue\n"
                            << "end_header\n";
            for (std::size_t i = 0, size = points.size(); i < size; i++)
            {
                const Point& p = points[i];
                const Value& c = values[i];
                Value value = (c - min) * inv_range;
                if (invert_values)
                    value = 1.0 - value;
                std::array<Value, 3> rgbColor = internal::get_rgb_color_heat_map(value);
                const unsigned int& c0 = static_cast<unsigned int>(255 * rgbColor[0]);
                const unsigned int& c1 = static_cast<unsigned int>(255 * rgbColor[1]);
                const unsigned int& c2 = static_cast<unsigned int>(255 * rgbColor[2]);
                output_point_set << p[0] << " " << p[1] << " " << p[2] << " "
                                << "1 0 0 "
                                << c0 << " " << c1 << " " << c2 << "\n";
            }
            output_point_set.close();
        }
    }

    template <typename Vector>
    bool normalize_if_possible(Vector &vec)
    {
        double sqLen = vec[0] * vec[0] + vec[1] * vec[1] + vec[2] * vec[2];
        if (sqLen == 0)
            return false;
        vec /= CGAL::sqrt(sqLen);
        return true;
    }

    template <typename Point, typename Vector>
    Point get_ray_plane_intersection(
        const Point& ray_origin, const Vector& ray_direction,
        const Point& plane_center, const Vector& plane_normal
        )
    {
        float cosAngle = CGAL::scalar_product(ray_direction, plane_normal);
        if (cosAngle == 0.0)
            return ray_origin;
        return ray_origin + (CGAL::scalar_product(Vector(ray_origin, plane_center), plane_normal) / cosAngle) * ray_direction;
    }

    template <typename Point, typename Vector>
    Point get_point_segment_projection(const Point& point_to_project, const Point& segment_begin, const Point& segment_end)
    {
        Vector segment_vector(segment_begin, segment_end);
        auto segment_length = segment_vector.squared_length();
        if (segment_length == 0)
            return segment_begin;
        segment_length /= CGAL::sqrt(segment_length);
        segment_vector /= segment_length;
        const auto& projection_length = CGAL::scalar_product(segment_vector, Vector(segment_begin, point_to_project));
        if (projection_length <= 0)
            return segment_begin;
        if (projection_length >= segment_length)
            return segment_end;
        return segment_begin + projection_length * segment_vector;
    }

    template <typename T, std::size_t dim>
    std::array<T, dim> nd_vec_add(const std::array<T, dim>& a, const std::array<T, dim>& b)
    {
        std::array<T, dim> res;
        for (std::size_t i = 0; i < dim; i++)
        {
            res[i] = a[i] + b[i];
        }
        return res;
    }

    template <typename T, std::size_t dim>
    std::array<T, dim> nd_vec_mult(const std::array<T, dim>& p, const T& factor)
    {
        std::array<T, dim> res;
        for (std::size_t i = 0; i < dim; i++)
            res[i] = p[i] * factor;
        return res;
    }

    template <typename P>
    P vec_component_mult(const P& p1, const P& p2)
    {
        return P(p1[0]*p2[0], p1[1]*p2[1], p1[2]*p2[2]);
    }

    template <typename Point, typename Value>
    void write_point_set_values(
        const std::string& filename,
        const std::vector<Point>& points,
        const std::vector<Value>& values,
        bool invert_values = false,
        bool normalize_range = true)
    {
        Value min, max;
        if (normalize_range)
        {
            auto [min_it, max_it] = std::minmax_element(values.begin(), values.end());
            min = *min_it;
            max = *max_it;
        }
        else
        {
            min = 0;
            max = 1;
        }
        internal::write_point_set_values_min_max(filename, points, values, min, max, invert_values);
    }

    template <typename Point, typename Image>
    Point image_pos_to_world_pos(const Point &p, const Image &image)
    {
        return Point(p[0] * image.vx() + image.tx(), p[1] * image.vy() + image.ty(), p[2] * image.vz() + image.tz());
    }
} // namespace helper

#endif // GEMETRY_HELPER_H