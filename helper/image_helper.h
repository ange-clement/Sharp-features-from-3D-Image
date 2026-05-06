#ifndef IMAGE_HELPER_H
#define IMAGE_HELPER_H

#include "cgal_helper.h"
#include "dgtal_helper_types.h"

#include <CGAL/Image_3.h>

namespace helper
{

    namespace internal
    {
        struct Is_label_condition
        {
            const std::size_t label;
            Is_label_condition(const std::size_t& label) : label(label) { }
            template <typename Word>
            bool operator()(const Word& w) const
            {
                return w == label;
            }
        };

        struct Is_not_label_condition
        {
            const std::size_t label;
            Is_not_label_condition(const std::size_t& label) : label(label) { }
            template <typename Word>
            bool operator()(const Word& w) const
            {
                return w != label;
            }
        };

        template <typename Word, typename DomainCondition>
        std::array<std::size_t, 6> get_voxels_with_condition_bbox_known_word(const CGAL::Image_3 &image, const DomainCondition &cond = DomainCondition())
        {
            std::size_t xmin = std::numeric_limits<std::size_t>::max();
            std::size_t ymin = std::numeric_limits<std::size_t>::max();
            std::size_t zmin = std::numeric_limits<std::size_t>::max();
            std::size_t xmax = 0;
            std::size_t ymax = 0;
            std::size_t zmax = 0;

            const CGAL_types::Image *im = image.image();

            const std::size_t &xdim = im->xdim;
            const std::size_t &ydim = im->ydim;
            const std::size_t &zdim = im->zdim;

            for (std::size_t z = 0; z < zdim; z++)
                for (std::size_t y = 0; y < ydim; y++)
                    for (std::size_t x = 0; x < xdim; x++)
                    {
                        const Word &w = CGAL::IMAGEIO::static_evaluate<Word>(im, x, y, z);
                        if (cond(w))
                        {
                            if (x < xmin)
                                xmin = x;
                            if (y < ymin)
                                ymin = y;
                            if (z < zmin)
                                zmin = z;

                            if (x > xmax)
                                xmax = x;
                            if (y > ymax)
                                ymax = y;
                            if (z > zmax)
                                zmax = z;
                        }
                    }
            return {xmin, ymin, zmin, xmax, ymax, zmax};
        }

        template <typename DomainCondition>
        struct Get_voxels_with_condition_bbox_struct
        {
            template <typename Word>
            std::array<std::size_t, 6> operator()(CGAL::Image_3 &image, const DomainCondition &cond = DomainCondition()) const
            {
                return get_voxels_with_condition_bbox_known_word<Word, DomainCondition>(image, cond);
            }
        };

        template <typename Word>
        CGAL::Image_3 crop_image_known_word(const CGAL::Image_3& image, const std::array<int, 6>& image_space_bbox, const Word& outside_value = 0)
        {
            auto [xmin, ymin, zmin, xmax, ymax, zmax] = image_space_bbox;

            const std::size_t xrange = xmax - xmin;
            const std::size_t yrange = ymax - ymin;
            const std::size_t zrange = zmax - zmin;

            const CGAL_types::Image* src_im = image.image();

            CGAL_types::Image* res_im = CGAL::_createImage(xrange, yrange, zrange, 1,
                                        src_im->vx, src_im->vy, src_im->vz, 1,
                                        src_im->wordKind, src_im->sign);

            res_im->tx = static_cast<float>(src_im->tx + xmin * src_im->vx);
            res_im->ty = static_cast<float>(src_im->ty + ymin * src_im->vy);
            res_im->tz = static_cast<float>(src_im->tz + zmin * src_im->vz);

            const std::size_t& xdim = src_im->xdim;
            const std::size_t& ydim = src_im->ydim;
            const std::size_t& zdim = src_im->zdim;

            for (std::size_t z = 0; z < zrange; z++)
                for (std::size_t y = 0; y < yrange; y++)
                    for (std::size_t x = 0; x < xrange; x++)
                    {
                        const std::size_t src_x_coord = x + xmin;
                        const std::size_t src_y_coord = y + ymin;
                        const std::size_t src_z_coord = z + zmin;
                        if (src_x_coord >= xdim || src_y_coord >= ydim || src_z_coord >= zdim)
                            CGAL::IMAGEIO::static_evaluate<Word>(res_im, x, y, z) = outside_value;
                        else
                        {
                            const Word& w = CGAL::IMAGEIO::static_evaluate<Word>(src_im, x + xmin, y + ymin, z + zmin);
                            CGAL::IMAGEIO::static_evaluate<Word>(res_im, x, y, z) = w;
                        }
                    }

            return CGAL::Image_3(res_im);
        }

        // does not copy the values
        CGAL_types::Image *create_image_from_image_dimensions(const CGAL_types::Image * src)
        {
            CGAL_types::Image * res = CGAL::_createImage(src->xdim, src->ydim, src->zdim, src->vdim,
                                        src->vx, src->vy, src->vz, src->wdim,
                                        src->wordKind, src->sign);

            res->cx = src->cx;
            res->cy = src->cy;
            res->cz = src->cz;

            res->tx = src->tx;
            res->ty = src->ty;
            res->tz = src->tz;

            res->rx = src->rx;
            res->ry = src->ry;
            res->rz = src->rz;

            return res;
        }

        // Removes floating voxels: groups that have less voxels than voxel_count_threshold
        //      Also remove groups with background value
        template <typename Word>
        void remove_floating_voxels_word(
            const CGAL::Image_3& src_image, const Word &background,
            CGAL_types::Image* &res_im,
            const unsigned int voxel_count_threshold = 27)
        {
            using CGAL::IMAGEIO::static_evaluate;
            const CGAL_types::Image* src_im = src_image.image();
            res_im = create_image_from_image_dimensions(src_im);
            std::fill_n(static_cast<Word*>(res_im->data), std::size_t(res_im->xdim)*std::size_t(res_im->ydim)*std::size_t(res_im->zdim), static_cast<Word>(0));

            std::vector<bool> processed_voxels(src_image.size(), false);
            std::vector<std::size_t> component_coords;
            component_coords.reserve(voxel_count_threshold);

            const unsigned int& xdim = static_cast<unsigned int>(src_im->xdim);
            const unsigned int& ydim = static_cast<unsigned int>(src_im->ydim);
            const unsigned int& zdim = static_cast<unsigned int>(src_im->zdim);
            // look for unprocessed voxel of a new object
            for(unsigned int z=0; z < zdim; z++){
                for(unsigned int y=0; y < ydim; y++){
                    for(unsigned int x=0; x < xdim; x++){
                        std::size_t point_index = helper::coords_to_voxel_index(x, y, z, xdim, ydim);
                        if (processed_voxels[point_index])
                            continue;
                        // new voxel : save component coordinates and copy if its voxel count is over the threshold
                        Word src_value = static_evaluate<Word>(src_im, point_index);
                        component_coords.clear();
                        // search
                        DGtal_types::ImageVoxelVisitor voxel_visitor(src_image);
                        const DGtal_types::CGAL_IsLabelPredicate is_label_predicate(src_im, src_value);
                        DGtal::DepthFirstVisitor<DGtal_types::ImageVoxelVisitor, std::set<DGtal_types::ImageVoxelVisitor::Vertex> > visitor
                            ( voxel_visitor, {x,y,z} );
                        while ( ! visitor.finished() && component_coords.size() < voxel_count_threshold )
                        {
                            const auto& [visited_voxel, topological_distance] = visitor.current();
                            std::size_t visited_point_index = helper::coords_to_voxel_index(visited_voxel[0], visited_voxel[1], visited_voxel[2], xdim, ydim);
                            component_coords.push_back(visited_point_index);
                            processed_voxels[visited_point_index] = true;
                            visitor.expand(is_label_predicate);
                        }
                        if (component_coords.size() < voxel_count_threshold) // component is noise
                        {
                            // background group that is too small : replace with a neighbor value (choose from a previous adjacent voxel)
                            // set the value for replacement
                            std::array<unsigned int, 3> previous_voxel = {x, y, z};
                            if (previous_voxel[0] != 0)
                            {
                                previous_voxel[0]--;
                                src_value = static_evaluate<Word>(src_im, previous_voxel[0], previous_voxel[1], previous_voxel[2]);
                            }
                            else if (previous_voxel[1] != 0)
                            {
                                previous_voxel[1]--;
                                src_value = static_evaluate<Word>(src_im, previous_voxel[0], previous_voxel[1], previous_voxel[2]);
                            }
                            else if (previous_voxel[2] != 0)
                            {
                                previous_voxel[2]--;
                                src_value = static_evaluate<Word>(src_im, previous_voxel[0], previous_voxel[1], previous_voxel[2]);
                            }
                            else
                            {
                                src_value = background;
                            }
                        }
                        if (src_value != background)
                        {
                            // copy the processed voxels then copy the rest of the component
                            for (const std::size_t& component_visited_voxel : component_coords)
                            {
                                static_evaluate<Word>(res_im, component_visited_voxel) = src_value;
                            }
                        }
                    }
                }
            }
        }
    } // namespace internal


    template <typename T>
    inline T
    coords_to_voxel_index(const T &i, const T &j, const T &k, const T &xdim, const T &ydim)
    {
        return static_cast<T>((k * ydim + j) * xdim + i);
    }

    inline std::array<std::size_t, 3>
    voxel_index_to_coord(const std::size_t &voxel_id, const std::size_t &xdim, const std::size_t &ydim)
    {
        std::size_t xydim = xdim * ydim;
        return {voxel_id % xdim, (voxel_id % xydim) / xdim, voxel_id / xydim};
    }

    template <typename Word>
    const Word& safe_get_image_value(const CGAL_types::Image* im,
                                    const long i, const long j, const long k,
                                    const long& xdim, const long& ydim, const long& zdim,
                                    const Word& background)
    {
        if (i < 0 || i >= xdim
            || j < 0 || j >= ydim
            || k < 0 || k >= zdim)
            return background;
        return CGAL::IMAGEIO::static_evaluate<Word>(im, i, j, k);
    }

    template <typename Word>
    const Word& safe_get_image_value(const CGAL_types::Image* im,
                                    const long i, const long j, const long k,
                                    const Word& background)
    {
        return safe_get_image_value<Word>(im, i, j, k, im->xdim, im->ydim, im->zdim, background);
    }

    template <typename DomainCondition>
    std::array<std::size_t, 6> get_voxels_with_condition_bbox(CGAL::Image_3 &image, const DomainCondition &condition = DomainCondition())
    {
        internal::Get_voxels_with_condition_bbox_struct<DomainCondition> op;
        CGAL_IMAGE_IO_CASE(image.image(), return op.template operator()<Word>(image, condition));
        return {0, 0, 0, 0, 0, 0};
    }

    CGAL::Image_3 crop_image(CGAL::Image_3& image, const std::array<int, 6>& image_space_bbox, const double& outside_value = 0)
    {
        CGAL_IMAGE_IO_CASE(image.image(), return internal::crop_image_known_word<Word>(image, image_space_bbox, static_cast<Word>(outside_value)));
        return CGAL::Image_3();
    }

    void write_image(CGAL::Image_3& image, const std::string &filename)
    {
        _writeImage(image.image(), filename.c_str());
    }

    CGAL::Image_3 remove_floating_voxels(const CGAL::Image_3& src_image, const unsigned int voxel_count_threshold = 27)
    {
        CGAL_types::Image* res_im = nullptr;
        CGAL_IMAGE_IO_CASE(src_image.image(), internal::remove_floating_voxels_word<Word>(src_image, CGAL::IMAGEIO::static_evaluate<Word>(src_image.image(), 0), res_im, voxel_count_threshold));
        return CGAL::Image_3(res_im);
    }
} // namespace helper

#endif // IMAGE_HELPER_H