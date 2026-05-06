#ifndef MESH_GENERATOR_H
#define MESH_GENERATOR_H

#include <vector>
#include <array>

#include <CGAL/number_utils.h>
#include <CGAL/Kernel/global_functions_3.h>
#include <CGAL/global_functions_distance_3.h>

#include "geometry_helper.h"

namespace helper
{

    namespace internal
    {
        // https://stackoverflow.com/questions/3018313/algorithm-to-convert-rgb-to-hsv-and-hsv-to-rgb-in-range-0-255-for-both
        std::array<double, 3> hsv2rgb(const std::array<double, 3> &in)
        {
            double                hh, p, q, t, ff;
            long                  i;
            std::array<double, 3> out;

            if(in[1] <= 0.0) {       // < is bogus, just shuts up warnings
                out[0] = in[2];
                out[1] = in[2];
                out[2] = in[2];
                return out;
            }
            hh = in[0] * 360;
            if(hh >= 360.0) hh = 0.0;
            hh /= 60.0;
            i = (long)hh;
            ff = hh - i;
            p = in[2] * (1.0 - in[1]);
            q = in[2] * (1.0 - (in[1] * ff));
            t = in[2] * (1.0 - (in[1] * (1.0 - ff)));

            switch(i) {
            case 0:
                out[0] = in[2];
                out[1] = t;
                out[2] = p;
                break;
            case 1:
                out[0] = q;
                out[1] = in[2];
                out[2] = p;
                break;
            case 2:
                out[0] = p;
                out[1] = in[2];
                out[2] = t;
                break;

            case 3:
                out[0] = p;
                out[1] = q;
                out[2] = in[2];
                break;
            case 4:
                out[0] = t;
                out[1] = p;
                out[2] = in[2];
                break;
            case 5:
            default:
                out[0] = in[2];
                out[1] = p;
                out[2] = q;
                break;
            }
            return out;
        }

        std::array<double, 3> get_rgb_color_heat_map(const double& value)
        {
            if (value >= 0)
                return hsv2rgb({0.2 - value * 0.2, 0.75 + value * 0.25, 0.5 + value * 0.5});
            else
                return hsv2rgb({0.2 - value * 0.2, 0.75 - value * 0.25, 0.5 - value * 0.5});
        }
    } // namespace internal


    template <typename Point, typename Vector>
    void add_mesh_sphere_from_point(
        const Point& point,
        std::vector<Point> & mesh_points,
        std::vector<std::array<std::size_t, 3>> & mesh_faces,
        float sphereSize, int subdivision
        )
    {
        // 12 initial points and 20 initial faces, for each subdivision:
        // the number of triangles is multiplied by 4 (nt+1 = nt * 4, nt_0 = 20)
        // to the points, add 3 points per triangle (np+1 = np + 3 * nt, np_0 = 12, nt_0 = 20)
        std::size_t number_of_triangles = 20 * pow(4, subdivision);
        std::size_t number_of_point = 12 + 20 * (pow(4, subdivision) - 1);

        mesh_points.reserve(mesh_points.size() + number_of_point);
        mesh_faces.reserve(mesh_faces.size() + number_of_triangles);

        // Icosahedron
        std::vector<Vector> icosahedron_points = {
            Vector( 0,         0.525731,  0.850651),
            Vector( 0.850651,  0,         0.525731),
            Vector( 0.525731,  0.850651,  0),
            Vector( 0,        -0.525731,  0.850651),
            Vector( 0.525731, -0.850651,  0),
            Vector( 0,        -0.525731, -0.850651),
            Vector(-0.850651,  0,        -0.525731),
            Vector( 0,         0.525731, -0.850651),
            Vector(-0.525731,  0.850651,  0),
            Vector( 0.850651,  0,        -0.525731),
            Vector(-0.850651,  0,         0.525731),
            Vector(-0.525731, -0.850651,  0)
        };
        icosahedron_points.reserve(number_of_point);

        std::vector<std::size_t> icosahedron_triangles = {
            0,1,2,
            3,1,0,
            3,4,1,
            5,6,7,
            8,0,2,
            8,2,7,
            8,7,6,
            9,1,4,
            9,5,7,
            9,2,1,
            9,7,2,
            9,4,5,
            10,3,0,
            10,0,8,
            10,8,6,
            11,4,3,
            11,3,10,
            11,6,5,
            11,5,4,
            11,10,6,
        };
        icosahedron_triangles.reserve(number_of_triangles);

        // Subdivide
        std::vector<std::size_t> subdivided_triangles;
        subdivided_triangles.reserve(number_of_triangles);

        std::vector<std::size_t>* source_triangles_pt = &icosahedron_triangles;
        std::vector<std::size_t>* final_triangles_pt = &subdivided_triangles;

        for (int div = 0; div < subdivision; div++)
        {
            const std::vector<std::size_t>& source_triangles = *source_triangles_pt;
            std::vector<std::size_t>& final_triangles = *final_triangles_pt;
            final_triangles.clear();

            // subdivide each source triangle into final triangle
            for (std::size_t t = 0, size = source_triangles.size(); t < size; t += 3)
            {
                const std::size_t& t0 = source_triangles[t];
                const std::size_t& t1 = source_triangles[t+1];
                const std::size_t& t2 = source_triangles[t+2];
                const Vector& pt0 = icosahedron_points[t0];
                const Vector& pt1 = icosahedron_points[t1];
                const Vector& pt2 = icosahedron_points[t2];

                const std::size_t& point_index = icosahedron_points.size();

                Vector m01 ((pt0[0] + pt1[0]) * 0.5, (pt0[1] + pt1[1]) * 0.5, (pt0[2] + pt1[2]) * 0.5);
                Vector m12 ((pt1[0] + pt2[0]) * 0.5, (pt1[1] + pt2[1]) * 0.5, (pt1[2] + pt2[2]) * 0.5);
                Vector m20 ((pt2[0] + pt0[0]) * 0.5, (pt2[1] + pt0[1]) * 0.5, (pt2[2] + pt0[2]) * 0.5);

                m01 *= 1.0 / CGAL::sqrt(CGAL::squared_length(m01));
                m12 *= 1.0 / CGAL::sqrt(CGAL::squared_length(m12));
                m20 *= 1.0 / CGAL::sqrt(CGAL::squared_length(m20));

                icosahedron_points.push_back(m01); // point_index
                icosahedron_points.push_back(m12); // point_index+1
                icosahedron_points.push_back(m20); // point_index+2

                final_triangles.insert(final_triangles.end(), {
                                                                t0, point_index, point_index+2,
                                                                point_index, t1, point_index+1,
                                                                point_index+2, point_index+1, t2,
                                                                point_index, point_index+1, point_index+2
                                                            });
            }

            std::swap(source_triangles_pt, final_triangles_pt);
        }

        // make mesh
        const std::vector<std::size_t>& final_triangles = *source_triangles_pt;
        std::size_t point_index = mesh_points.size();
        for (const Vector v : icosahedron_points)
        {
            mesh_points.push_back(point + sphereSize * v);
        }

        for (std::size_t t = 0, size = final_triangles.size(); t < size; t += 3)
        {
            mesh_faces.push_back({point_index + final_triangles[t  ],
                                point_index + final_triangles[t+1],
                                point_index + final_triangles[t+2]});
        }
    }

    template <typename Point, typename Vector>
    void add_mesh_cylinder_from_polyline(
        const std::vector<std::size_t> & polyline,
        const std::vector<Point> &polyline_points,
        float sizeFactor, int resCylinder,
        std::vector<Point> &points,
        std::vector<std::array<std::size_t, 3>> &faces)
    {
        std::vector<Point> cylinder_shape; // a circle
        cylinder_shape.reserve(resCylinder);
        float angleIncr = 2.0 * CGAL_PI / (float) resCylinder;
        for (std::size_t i = 0; i < resCylinder; i++){
            float angle = i * angleIncr;
            float c = sizeFactor * cos(angle);
            float s = sizeFactor * sin(angle);
            cylinder_shape.push_back(Point(c, s, 0));
        }
        std::size_t shape_size = cylinder_shape.size();

        points.reserve(points.size() + shape_size *  polyline.size());
        std::size_t number_of_added_faces = shape_size * (polyline.size()-1);
        faces.reserve( faces.size() + number_of_added_faces * 3);
        // construct basis on a plane defined by curVec :
        //      Y is the plane normal, X,Z are vectors in the plane
        Vector vecY;//
        std::size_t start_index = 0;
        double vecY_sqlen = 0;
        while (vecY_sqlen == 0 && start_index < polyline.size()-1)
        {
            vecY = Vector(polyline_points[polyline[start_index]], polyline_points[polyline[start_index+1]]);
            vecY_sqlen = vecY.squared_length();
            start_index++;
        }
        if (start_index == polyline.size())
        {
            std::cout << "ERROR: mesh cylinder: invalid polyline, all length are 0 !" << std::endl;
            return;
        }
        start_index--;
        vecY /= CGAL::sqrt(vecY_sqlen);
        Vector vecZ (-vecY[1], vecY[0], 0);
        if (vecZ[0] == 0 && vecZ[1] == 0)
            vecZ = Vector(1,0,0);
        else
            vecZ /= CGAL::sqrt(vecZ.squared_length());
        Vector vecX = CGAL::cross_product(vecZ, vecY);
        for (const Point& p : cylinder_shape)
        {
            const Point in_plane_point = polyline_points[polyline[0]] + p[0] * vecX + p[1] * vecZ;
            points.push_back(in_plane_point);
        }

        std::vector<Point> shape;
        shape.reserve(cylinder_shape.size());
        for (std::size_t p_id = 1, size = polyline.size(); p_id < size; ++p_id)
        {
            // polyline :
            //
            //       _(i)--curVec-->(i+1)
            //       /\
            //    prevVec
            //     /
            //    /
            // (i-1)
            //
            // The current cylinder is constructed over curVec.
            // construct vectors
            const Point& p_cur  = polyline_points[polyline[p_id]];
            const Point& p_prev = polyline_points[polyline[p_id-1]];
            Vector prevVec = Vector(p_prev, p_cur);
            helper::normalize_if_possible(prevVec);
            Vector curVec;
            if (p_id < size-1)
            {
                Point p_next;
                p_next = polyline_points[polyline[p_id+1]];
                curVec = Vector(p_cur, p_next);
                helper::normalize_if_possible(curVec);
            }
            else
                curVec = prevVec;
            // construct plane
            Vector plane_normal;
            if (CGAL::scalar_product(prevVec, curVec) < -0.707) // disable small angles   -0.707 less than 45°
                plane_normal = prevVec - curVec;
            else
                plane_normal = prevVec + curVec;
            helper::normalize_if_possible(plane_normal);
            // get shape
            shape.clear();
            shape.insert(shape.end(), points.end() - shape_size, points.end());
            // Construct cylinder
            std::size_t start_index = points.size();
            for (const Point& p : shape)
            {
                points.push_back(helper::get_ray_plane_intersection(p, prevVec, p_cur, plane_normal));
            }
            for (std::size_t f_id = 1; f_id < shape_size; f_id++)
            {
                const std::size_t& index = start_index + f_id;
                faces.push_back({index-shape_size-1,
                                index-shape_size,
                                index});
                faces.push_back({index-shape_size-1,
                                index,
                                index-1});
            }
            const std::size_t& index = points.size()-1;
            faces.push_back({index-shape_size,
                            start_index-shape_size,
                            start_index});
            faces.push_back({index-shape_size,
                            start_index,
                            index});
        }
    }

    template <typename Point>
    void write_mesh_as_ply(
        const std::string& filename,
        const std::vector<Point> & points,
        const std::vector<std::array<std::size_t, 3>> & triangles,
        const std::vector<std::size_t> & triangles_values,
        const bool& multicolor = true,
        const bool& white_mode = false)
    {
        std::map<std::size_t, std::array<double, 3>> color_map;
        std::srand(1); // consistant colors
        for (const std::size_t& value : triangles_values)
        {
            if (color_map.count(value) != 0)
                continue;
            if (multicolor)
            {
                color_map[value] = internal::hsv2rgb({std::rand() /(double) RAND_MAX,                 // [0.0;1.0]
                                                      0.5 + (std::rand() /(double)(RAND_MAX*4)),      // [0.5;0.75]
                                                      0.5 + (std::rand() /(double)(RAND_MAX*4))});    // [0.5;0.75]
            }
            else
            {
                if (white_mode)
                {
                    // white mode : around HSV(235, 72, 255)
                    color_map[value] = internal::hsv2rgb({0.6 + (std::rand() /(double)(RAND_MAX*10)),     // [0.6;0.7]
                                                          0.2 + (std::rand() /(double)(RAND_MAX*5)),      // [0.2;0.4]
                                                          0.8 + (std::rand() /(double)(RAND_MAX*5))});    // [0.8;1.0]
                }
                else
                {
                    color_map[value] = {0.3, 0.3, 0.3};
                }
            }
        }
        // color_map[0] is the color of corners
        if (!multicolor && white_mode)
            color_map[0] = {0.8, 0.6, 0.6};
        else
            color_map[0] = {0.5, 0.2, 0.2};
        std::ofstream output_mesh_value(filename, std::ios_base::binary);
        output_mesh_value << "ply\nformat ascii 1.0\n"
                        << "element vertex " << points.size() << "\n"
                        << "property double x\nproperty double y\nproperty double z\n"
                        << "property int id\n"
                        << "element face " << triangles.size() << "\n"
                        << "property list uchar int vertex_indices\n"
                        << "property uchar red\nproperty uchar green\nproperty uchar blue\n"
                        << "end_header\n";
        for (std::size_t i = 0, size = points.size(); i < size; i++)
        {
            const Point& p = points[i];
            output_mesh_value << p[0] << " " << p[1] << " " << p[2] << " "
                            << i << "\n";
        }
        for (std::size_t i = 0; i < triangles.size(); i++)
        {
            const std::array<std::size_t, 3>& triangle = triangles[i];
            const std::size_t& triangle_value = triangles_values[i];
            output_mesh_value << 3 << " ";
            for (int j = 0; j < 3; j++)
            {
                output_mesh_value << triangle[j] << " ";
            }
            std::array<double, 3> rgbColor = color_map.at(triangle_value);
            const unsigned int& c0 = static_cast<unsigned int>(255 * std::min(1.0, std::max(0.0, rgbColor[0])));
            const unsigned int& c1 = static_cast<unsigned int>(255 * std::min(1.0, std::max(0.0, rgbColor[1])));
            const unsigned int& c2 = static_cast<unsigned int>(255 * std::min(1.0, std::max(0.0, rgbColor[2])));
            output_mesh_value << c0 << " " << c1 << " " << c2 << "\n";
        }
        output_mesh_value.close();
    }

    template <typename Point>
    void write_mesh_with_heat_map_as_ply(
        const std::string& filename,
        const std::vector<Point> points, const std::vector<double>& values,
        const bool hsvColors = true, const bool normalise_negatives = true,
        const std::optional<double>& min_value = std::nullopt,
        const std::optional<double>& max_value = std::nullopt)
    {
        double min, max;
        if (!min_value || !max_value)
        {
            const auto& [computed_min_it, computed_max_it] = std::minmax_element(values.begin(), values.end());
            min = min_value.value_or(*computed_min_it);
            max = max_value.value_or(*computed_max_it);
        }
        else
        {
            min = min_value.value();
            max = max_value.value();
        }
        if (!normalise_negatives)
        {
            max = std::max(max, std::abs(min));
            min = 0;
        }
        double inv_range = 1.0 / (max - min);

        std::ofstream output_mesh_value(filename, std::ios_base::binary);
        output_mesh_value << "ply\nformat ascii 1.0\n"
                        << "element vertex " << points.size() << "\n"
                        << "property double x\nproperty double y\nproperty double z\n"
                        << "property int id\n"
                        << "element face " << points.size() / 4 << "\n"
                        << "property list uchar int vertex_indices\n"
                        << "property uchar red\nproperty uchar green\nproperty uchar blue\n"
                        << "end_header\n";
        for (std::size_t i = 0, size = points.size(); i < size; i++)
        {
            const Point& p = points[i];
            output_mesh_value << p[0] << " " << p[1] << " " << p[2] << " "
                            << i << "\n";
        }
        for (std::size_t i = 0, size = points.size()/4; i < size; i++)
        {
            output_mesh_value << "4 ";
            for (int j = 0; j < 4; j++)
            {
                output_mesh_value << i*4+j << " ";
            }

            double value = (values[i] - min) * inv_range;
            std::array<double, 3> rgbColor;
            if (hsvColors)
                rgbColor = internal::get_rgb_color_heat_map(value);
            else
                rgbColor = {value, value, value};
            const unsigned int& c0 = static_cast<unsigned int>(255 * rgbColor[0]);
            const unsigned int& c1 = static_cast<unsigned int>(255 * rgbColor[1]);
            const unsigned int& c2 = static_cast<unsigned int>(255 * rgbColor[2]);
            output_mesh_value << c0 << " " << c1 << " " << c2 << "\n";
        }
        output_mesh_value.close();
    }

    template <typename Point>
    void write_mesh_with_colors_as_ply(
        const std::string& filename,
        const std::vector<Point>& points,
        const std::vector<double>& colors)
    {
        std::ofstream output_mesh_value(filename, std::ios_base::binary);
        output_mesh_value << "ply\nformat ascii 1.0\n"
                        << "element vertex " << points.size() << "\n"
                        << "property double x\nproperty double y\nproperty double z\n"
                        << "property int id\n"
                        << "element face " << points.size() / 4 << "\n"
                        << "property list uchar int vertex_indices\n"
                        << "property uchar red\nproperty uchar green\nproperty uchar blue\n"
                        << "end_header\n";
        for (std::size_t i = 0, size = points.size(); i < size; i++)
        {
            const Point& p = points[i];
            output_mesh_value << p[0] << " " << p[1] << " " << p[2] << " "
                            << i << "\n";
        }
        for (std::size_t i = 0, size = points.size()/4; i < size; i++)
        {
            output_mesh_value << "4 ";
            for (int j = 0; j < 4; j++)
            {
                output_mesh_value << i*4+j << " ";
            }
            const unsigned int& c0 = static_cast<unsigned int>(255 * colors[i*3]  );
            const unsigned int& c1 = static_cast<unsigned int>(255 * colors[i*3+1]);
            const unsigned int& c2 = static_cast<unsigned int>(255 * colors[i*3+2]);
            output_mesh_value << c0 << " " << c1 << " " << c2 << "\n";
        }
        output_mesh_value.close();
    }


} // namespace helper

#endif // MESH_GENERATOR_H