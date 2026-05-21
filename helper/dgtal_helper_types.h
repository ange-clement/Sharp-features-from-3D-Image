#ifndef DGTAL_HELPER_TYPES_H
#define DGTAL_HELPER_TYPES_H

#include <DGtal/base/Common.h>
#include <DGtal/helpers/Shortcuts.h>
#include <DGtal/helpers/ShortcutsGeometry.h>
#include <DGtal/helpers/StdDefs.h>
#include <DGtal/topology/DigitalSurface.h>
#include <DGtal/images/ImageSelector.h>
// #include <DGtal/shapes/SurfaceMesh.h>
#include <DGtal/kernel/domains/HyperRectDomain.h>

#include "cgal_helper.h"

namespace helper
{

    struct DGtal_types
    {
        typedef DGtal::Z3i::Space Space;
        typedef DGtal::Z3i::KSpace KSpace;
        typedef DGtal::HyperRectDomain<Space> Domain;

        // typedef KSpace::Cell Cell;
        typedef DGtal::Z3i::RealPoint RealPoint;
        typedef DGtal::Z3i::RealVector RealVector;

        typedef DGtal::ImageSelector<Domain, unsigned char>::Type Image;

        // typedef KSpace::SurfelSet SurfelSet;
        // typedef DGtal::SetOfSurfels<KSpace, SurfelSet> DigitalSurfaceContainer_Set;
        // typedef DGtal::DigitalSurface<DigitalSurfaceContainer_Set> DigitalSurface_Set;

        // typedef DGtal::SurfaceMesh<RealPoint, RealVector> SurfaceMesh;

        typedef DGtal::Shortcuts<KSpace> SH;
        typedef DGtal::ShortcutsGeometry<KSpace> SHG;

        struct Surface_info
        {
            KSpace space;
            DGtal::CountedPtr<SH::DigitalSurface> surface;
            SH::SurfelRange surfels;
            SH::PointelRange pointels;
            SH::RealVectors normals;
            DGtal::CanonicSCellEmbedder<KSpace> embedder;
            DGtal::CanonicCellEmbedder<KSpace> uembedder;

            void update_normals(const std::vector<RealVector> &new_normals)
            {
                std::transform(new_normals.cbegin(), new_normals.cend(), normals.begin(),
                               [](const RealVector &v)
                               { return v; });
            }

            void get_surfel_points(std::vector<DGtal_types::RealPoint> &surfel_points) const
            {
                surfel_points.resize(surfels.size());
                std::transform(surfels.cbegin(), surfels.cend(), surfel_points.begin(),
                            [this](const DGtal_types::SH::SCell &c)
                            { return embedder(c); });
            }
        };

        struct ImageVoxelVisitor
        {
            typedef std::size_t Size;
            typedef std::array<unsigned int, 3> Vertex;

            const std::array<std::size_t, 3> image_dimensions;

            ImageVoxelVisitor(const Image &image)
                : image_dimensions({static_cast<std::size_t>(image.extent()[0]),
                                    static_cast<std::size_t>(image.extent()[1]),
                                    static_cast<std::size_t>(image.extent()[2])})
            {
            }

            ImageVoxelVisitor(const CGAL::Image_3 &image)
                : image_dimensions({static_cast<std::size_t>(image.xdim()),
                                    static_cast<std::size_t>(image.ydim()),
                                    static_cast<std::size_t>(image.zdim())})
            {
            }

            Size degree(const Vertex &v) const
            {
                std::size_t degree = 6;
                if (v[0] == 0 || v[0] == image_dimensions[0] - 1)
                    degree--;
                if (v[1] == 1 || v[1] == image_dimensions[1] - 1)
                    degree--;
                if (v[2] == 0 || v[2] == image_dimensions[2] - 1)
                    degree--;
                return degree;
            }
            Size bestCapacity() const { return 6; }

            template <typename OutputIterator>
            void writeNeighbors(OutputIterator out, const Vertex &v) const
            {
                if (v[2] != 0)
                    *out++ = {v[0], v[1], v[2] - 1};
                if (v[2] != image_dimensions[2] - 1)
                    *out++ = {v[0], v[1], v[2] + 1};
                if (v[1] != 0)
                    *out++ = {v[0], v[1] - 1, v[2]};
                if (v[1] != image_dimensions[1] - 1)
                    *out++ = {v[0], v[1] + 1, v[2]};
                if (v[0] != 0)
                    *out++ = {v[0] - 1, v[1], v[2]};
                if (v[0] != image_dimensions[0] - 1)
                    *out++ = {v[0] + 1, v[1], v[2]};
            }

            template <typename OutputIterator, typename VertexPredicate>
            void insert_if_predicate(OutputIterator out, const Vertex &v, const Vertex &vn, const VertexPredicate &p) const
            {
                if (p(v, vn))
                    *out++ = vn;
            }

            template <typename OutputIterator, typename VertexPredicate>
            void writeNeighbors(OutputIterator out, const Vertex &v, const VertexPredicate &p) const
            {
                if (v[2] != 0)
                    insert_if_predicate(out, v, {v[0], v[1], v[2] - 1}, p);
                if (v[2] != image_dimensions[2] - 1)
                    insert_if_predicate(out, v, {v[0], v[1], v[2] + 1}, p);
                if (v[1] != 0)
                    insert_if_predicate(out, v, {v[0], v[1] - 1, v[2]}, p);
                if (v[1] != image_dimensions[1] - 1)
                    insert_if_predicate(out, v, {v[0], v[1] + 1, v[2]}, p);
                if (v[0] != 0)
                    insert_if_predicate(out, v, {v[0] - 1, v[1], v[2]}, p);
                if (v[0] != image_dimensions[0] - 1)
                    insert_if_predicate(out, v, {v[0] + 1, v[1], v[2]}, p);
            }
        };

        struct IsLabelPredicate
        {
            const DGtal_types::Image &image;
            const unsigned char &value;
            IsLabelPredicate(const DGtal_types::Image &image, const unsigned char &value)
                : image(image), value(value)
            {
            }
            bool operator()(const std::array<unsigned int, 3> &/*v*/, const std::array<unsigned int, 3> &v_n) const
            {
                return image(DGtal_types::Image::Point(v_n[0], v_n[1], v_n[2])) == value;
            }
        };

        template <typename Word>
        struct CGAL_IsLabelPredicate
        {
            const CGAL_types::Image*& image;
            const Word& value;
            CGAL_IsLabelPredicate(const CGAL_types::Image*& image, const Word& value)
                : image(image), value(value)
            { }
            bool operator()(const std::array<unsigned int, 3>& /*v*/, const std::array<unsigned int, 3>& v_n) const
            {
                return CGAL::IMAGEIO::static_evaluate<Word>(image, v_n[0], v_n[1], v_n[2]) == value;
            }
        };
    };
}

#endif // DGTAL_HELPER_TYPES_H