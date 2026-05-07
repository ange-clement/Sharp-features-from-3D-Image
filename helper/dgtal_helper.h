// This file must be included before CGAL files that uses Eigen.
#ifndef DGTAL_HELPER_H
#define DGTAL_HELPER_H

#include <DGtal/images/SimpleThresholdForegroundPredicate.h>

#include <DGtal/base/Common.h>
#include <DGtal/helpers/Shortcuts.h>
#include <DGtal/helpers/ShortcutsGeometry.h>
#include <DGtal/helpers/StdDefs.h>
#include <DGtal/topology/helpers/Surfaces.h>
#include <DGtal/topology/DigitalSurface.h>
#include <DGtal/images/ImageSelector.h>
#include <DGtal/shapes/SurfaceMesh.h>
// This file must be included before CGAL files that uses Eigen.
#include <DGtal/dec/ATSolver2D.h>
#include <DGtal/dec/DiscreteExteriorCalculusFactory.h>
#include <DGtal/topology/MetricAdjacency.h>
#include <DGtal/shapes/MeshVoxelizer.h>
#include <DGtal/kernel/domains/HyperRectDomain.h>
#include <DGtal/io/readers/MeshReader.h>
#include <DGtal/images/ImageContainerBySTLVector.h>
#include <DGtal/topology/ImplicitDigitalSurface.h>
#include <DGtal/base/ConceptUtils.h>
#include <DGtal/geometry/volumes/KanungoNoise.h>
// #include <DGtal/images/IntervalForegroundPredicate.h>
// #include <DGtal/base/Clock.h>
// // #include <DGtal/base/CowPtr.h>

#include <CGAL/Image_3.h>

#include "dgtal_helper_types.h"
#include "timing_helper.h"
#include "image_helper.h"

namespace helper
{
    namespace internal
    {
        std::vector<DGtal_types::Image>
        split_in_objects_dgtal_image(const DGtal_types::Image &src_image, const unsigned char &background = 0)
        {
            typename DGtal_types::Image::Domain domain = src_image.domain();
            typename DGtal_types::Image::Vector extent = src_image.extent();
            std::size_t w = extent[0];
            std::size_t h = extent[1];
            std::size_t e = extent[2];

            std::vector<DGtal_types::Image> images;
            std::set<std::size_t> processed_voxels;

            // look for unprocessed voxel of a new object
            for (unsigned int z = 0; z < e; z++)
                for (unsigned int y = 0; y < h; y++)
                    for (unsigned int x = 0; x < w; x++)
                    {
                        DGtal_types::Image::Point pt(x, y, z);
                        const unsigned char &src_value = src_image(pt);
                        if (src_value == background)
                            continue;
                        std::size_t point_index = z * w * h + y * w + x;
                        const auto &[it, inserted] = processed_voxels.emplace(point_index);
                        if (!inserted)
                            continue;
                        // new voxel : create image and fill with this object
                        images.push_back(DGtal_types::Image(domain));
                        DGtal_types::Image &image = images.back();
                        // search
                        DGtal_types::ImageVoxelVisitor voxel_visitor(src_image);
                        const DGtal_types::IsLabelPredicate is_label_predicate(src_image, src_value);
                        DGtal::BreadthFirstVisitor<DGtal_types::ImageVoxelVisitor, std::set<DGtal_types::ImageVoxelVisitor::Vertex>> visitor(voxel_visitor, {x, y, z});
                        while (!visitor.finished())
                        {
                            const auto &[visited_voxel, topological_distance] = visitor.current();
                            std::size_t visited_point_index = visited_voxel[2] * w * h + visited_voxel[1] * w + visited_voxel[0];
                            image.setValue(DGtal_types::Image::Vertex(visited_voxel[0], visited_voxel[1], visited_voxel[2]), src_value);
                            processed_voxels.insert(visited_point_index);
                            visitor.expand(is_label_predicate);
                        }
                    }

            return images;
        }

        void add_pointel_values_from_surfel_value(
            const DGtal_types::RealPoint &surfel_point,
            const double &surfel_value,
            std::map<std::array<int, 3>, std::pair<double, std::size_t>> &point_value_and_number)
        {
            int i2 = static_cast<int>(surfel_point[0] * 2);
            int j2 = static_cast<int>(surfel_point[1] * 2);
            int k2 = static_cast<int>(surfel_point[2] * 2);
            std::array<int, 3> start;
            std::array<int, 3> axisU;
            std::array<int, 3> axisV;
            if ((i2 & 1) != 0) // surfel normal is the x axis
            {
                start = {i2 / 2 + 1, j2 / 2, k2 / 2};
                axisU = {0, 1, 0};
                axisV = {0, 0, 1};
            }
            else if ((j2 & 1) != 0) // surfel normal is the y axis
            {
                start = {i2 / 2, j2 / 2 + 1, k2 / 2};
                axisU = {0, 0, 1};
                axisV = {1, 0, 0};
            }
            else if ((k2 & 1) != 0) // surfel normal is the z axis
            {
                start = {i2 / 2, j2 / 2, k2 / 2 + 1};
                axisU = {1, 0, 0};
                axisV = {0, 1, 0};
            }
            for (int ic = 0; ic < 2; ic++)
                for (int jc = 0; jc < 2; jc++)
                {
                    std::array<int, 3> coord = {
                        start[0] + ic * axisU[0] + jc * axisV[0],
                        start[1] + ic * axisU[1] + jc * axisV[1],
                        start[2] + ic * axisU[2] + jc * axisV[2]};
                    auto [pt_value_it, inserted] = point_value_and_number.emplace(coord, std::make_pair(0, 0));
                    pt_value_it->second.first += surfel_value;
                    pt_value_it->second.second++;
                }
        }

        void make_average_pointel_values(
            const std::map<std::array<int, 3>, std::pair<double, std::size_t>> &point_value_and_number,
            std::vector<DGtal_types::RealPoint> &valuePoints,
            std::vector<double> &values)
        {
            valuePoints.reserve(valuePoints.size() + point_value_and_number.size());
            values.reserve(values.size() + point_value_and_number.size());
            for (const auto &[intPoint, value_and_number] : point_value_and_number)
            {
                valuePoints.push_back(DGtal_types::RealPoint(intPoint[0] - .5, intPoint[1] - .5, intPoint[2] - .5));
                values.push_back(value_and_number.first / (double)value_and_number.second);
            }
        }

        // values are per valuePoints, do not normalize values
        // add info to result
        void compute_integral_invariant_values_single_object_image(
            const DGtal_types::Image &src_image,
            std::vector<DGtal_types::RealPoint> &surfaces,
            std::vector<DGtal_types::RealPoint> &valuePoints,
            std::vector<double> &values,
            double radius, double grid_step = 1.0, float threshold = 0.5)
        {
            bBlock("transform image -> build digital surface");

            typedef DGtal::functors::SimpleThresholdForegroundPredicate<DGtal_types::Image> ImagePredicate;
            typedef DGtal::LightImplicitDigitalSurface<DGtal::Z3i::KSpace, ImagePredicate> MyLightImplicitDigitalSurface;
            typedef DGtal::DigitalSurface<MyLightImplicitDigitalSurface> MyDigitalSurface;

            ImagePredicate predicate = ImagePredicate(src_image, static_cast<DGtal_types::Image::Value>(threshold));

            DGtal::Z3i::Domain domain = src_image.domain();

            DGtal::Z3i::KSpace KSpaceShape;

            bool space_ok = KSpaceShape.init(domain.lowerBound(), domain.upperBound(), true);
            if (!space_ok)
            {
                DGtal::trace.error() << "Error in the Khalimsky space construction." << std::endl;
                eBlock();
                return;
            }

            DGtal::SurfelAdjacency<DGtal::Z3i::KSpace::dimension> SAdj(true);
            DGtal::Z3i::KSpace::Surfel bel = DGtal::Surfaces<DGtal::Z3i::KSpace>::findABel(KSpaceShape, predicate, 100000);
            MyLightImplicitDigitalSurface LightImplDigSurf(KSpaceShape, predicate, SAdj, bel);
            MyDigitalSurface digSurf(LightImplDigSurf);
            eBlock();

            bBlock("Integral Invariant measure");
            typedef DGtal::functors::IIPrincipalCurvatures3DFunctor<DGtal::Z3i::Space> MyIICurvatureFunctor;
            typedef DGtal::IntegralInvariantCovarianceEstimator<DGtal::Z3i::KSpace, ImagePredicate, MyIICurvatureFunctor> MyIICurvatureEstimator;

            typedef MyIICurvatureFunctor::Value Value;

            typedef DGtal::DepthFirstVisitor<MyDigitalSurface> Visitor;
            typedef DGtal::GraphVisitorRange<Visitor> VisitorRange;
            typedef typename VisitorRange::ConstIterator SurfelConstIterator;

            VisitorRange range(new Visitor(digSurf, *digSurf.begin()));
            SurfelConstIterator abegin = range.begin();
            SurfelConstIterator aend = range.end();

            MyIICurvatureFunctor curvatureFunctor;    // Functor used to convert volume -> curvature
            curvatureFunctor.init(grid_step, radius); // Initialisation for a grid step and a given Euclidean radius of convolution kernel

            MyIICurvatureEstimator curvatureEstimator(curvatureFunctor);
            curvatureEstimator.attach(KSpaceShape, predicate); // Setting a KSpace and a predicate on the object to evaluate
            curvatureEstimator.setParams(radius / grid_step);  // Setting the digital radius of the convolution kernel
            curvatureEstimator.init(grid_step, abegin, aend);  // Initialisation for a given h, and a range of surfels

            std::vector<Value> results;
            std::back_insert_iterator<std::vector<Value>> resultsIt(results); // output iterator for results of Integral Invariant curvature computation
            curvatureEstimator.eval(abegin, aend, resultsIt);                 // Computation
            eBlock();

            bBlock("output curvature values");
            std::vector<double> KMeanAbs_surfel_values;
            KMeanAbs_surfel_values.resize(results.size());
            for (std::size_t i = 0, size = results.size(); i < size; i++)
            {
                const auto &[k1, k2] = results[i];
                KMeanAbs_surfel_values[i] = -0.5 * (abs(k1) + abs(k2));
            }

            DGtal::CanonicSCellEmbedder<DGtal::Z3i::KSpace> embedder = DGtal_types::SH::getSCellEmbedder(KSpaceShape);

            std::map<std::array<int, 3>, std::pair<double, std::size_t>> point_value_and_number;
            std::size_t surfel_id = 0;

            surfaces.reserve(surfaces.size() + digSurf.size());

            VisitorRange range2(new Visitor(digSurf, *digSurf.begin()));
            aend = range2.end();
            // set surfel values and add them on pointels
            for (abegin = range2.begin(); abegin != aend; ++abegin)
            {
                const MyDigitalSurface::SCell &surfel = *abegin;
                DGtal_types::RealPoint surfel_point = embedder(surfel);
                const double &surfel_value = KMeanAbs_surfel_values[surfel_id];
                surfaces.push_back(surfel_point);
                add_pointel_values_from_surfel_value(surfel_point, surfel_value, point_value_and_number);
                surfel_id++;
            }

            make_average_pointel_values(point_value_and_number, valuePoints, values);
            eBlock();
        }

        // values are per valuePoints
        // does not normalize values
        void compute_VCM_values_single_object_image(
            const DGtal_types::Image &src_image,
            std::vector<DGtal_types::RealPoint> &surfaces,
            std::vector<DGtal_types::RealPoint> &valuePoints,
            std::vector<double> &values,
            double offset_radius = 5.0, double smooth_radius = 4.0)
        {
            // typedef DGtal::ImageSelector<DGtal_types::Domain, bool>::Type BinaryImage;
            typedef DGtal::functors::IntervalForegroundPredicate<DGtal_types::Image> ThresholdedImage;
            typedef DGtal::ImplicitDigitalSurface<DGtal_types::KSpace, ThresholdedImage> DigitalSurfaceContainer;

            typedef DGtal_types::KSpace::Surfel Surfel;
            typedef DGtal::DigitalSurface<DigitalSurfaceContainer> DigitalSurface;

            typedef DGtal::ExactPredicateLpSeparableMetric<DGtal_types::Space, 2> Metric;        // L2-metric type
            typedef DGtal::functors::HatPointFunction<DGtal::Z3i::Point, double> KernelFunction; // chi function type
            typedef DGtal::VoronoiCovarianceMeasureOnDigitalSurface<DigitalSurfaceContainer, Metric,
                                                                    KernelFunction>
                VCMOnSurface;
            typedef VCMOnSurface::Surfel2Normals::const_iterator S2NConstIterator;

            int thresholdMin = 0;
            int thresholdMax = 255;
            const double R = offset_radius;
            const double r = smooth_radius;
            const double trivial_r = 3;

            DGtal_types::KSpace ks;
            // Reads the volume
            bBlock("Loading image into memory and build digital surface.");
            ThresholdedImage thresholdedImage(src_image, thresholdMin, thresholdMax);
            ks.init(src_image.domain().lowerBound(),
                    src_image.domain().upperBound(), true);
            DGtal::SurfelAdjacency<DGtal_types::KSpace::dimension> surfAdj(true); // interior in all directions.
            Surfel bel = DGtal::Surfaces<DGtal_types::KSpace>::findABel(ks, thresholdedImage, 10000);
            DigitalSurfaceContainer *container =
                new DigitalSurfaceContainer(ks, thresholdedImage, surfAdj, bel, false);
            DigitalSurface surface(container); // acquired
            eBlock();

            bBlock("Initializing VCM");
            DGtal::Surfel2PointEmbedding embType = DGtal::Pointels; // Could be Pointels|InnerSpel|OuterSpel;
            Metric l2;                                              // Euclidean L2 metric
            KernelFunction chi(1.0, r);                             // hat function with support of radius r
            VCMOnSurface vcm_surface(surface, embType, R, r,
                                     chi, trivial_r, l2, true /* verbose */);
            eBlock();

            std::size_t index_first_added_surfel = surfaces.size();
            surfaces.reserve(surfaces.size() + surface.size());

            bBlock("Computing eigenvalues");
            DGtal_types::RealVector lambda; // eigenvalues of chi-vcm
            std::vector<double> surfel_values;
            for (S2NConstIterator it = vcm_surface.mapSurfel2Normals().begin(),
                                  itE = vcm_surface.mapSurfel2Normals().end();
                 it != itE; ++it)
            {
                Surfel s = it->first;
                DGtal::Z3i::Point kp = ks.sKCoords(s);
                DGtal_types::RealPoint rp(0.5 * (double)kp[0], 0.5 * (double)kp[1], 0.5 * (double)kp[2]);
                DGtal_types::RealVector n = it->second.vcmNormal;
                vcm_surface.getChiVCMEigenvalues(lambda, s);
                double ratio = -(lambda[0] + lambda[1]) / (lambda[0] + lambda[1] + lambda[2]);
                surfaces.push_back(rp);
                surfel_values.push_back(ratio);
            }
            eBlock();

            bBlock("Output VCM values");
            std::map<std::array<int, 3>, std::pair<double, std::size_t>> point_value_and_number;
            for (std::size_t surface_id = 0, surface_nb = surface.size(); surface_id < surface_nb; ++surface_id)
            {
                const DGtal_types::RealPoint &surfel_point = surfaces[index_first_added_surfel + surface_id];
                const double &surfel_value = surfel_values[surface_id];

                add_pointel_values_from_surfel_value(surfel_point, surfel_value, point_value_and_number);
            }

            make_average_pointel_values(point_value_and_number, valuePoints, values);
            eBlock();
        }

        // From DGtalTools
        void fillInterior(
            DGtal::ImageContainerBySTLVector<DGtal::Z3i::Domain, unsigned char> &image,
            const unsigned char &fillValue)
        {
            using namespace DGtal;

            // Flag image
            ImageContainerBySTLVector<Z3i::Domain, bool> imageFlag(image.domain());
            for (auto &p : imageFlag.domain())
            {
                if (image(p) != 0)
                    imageFlag.setValue(p, true);
            }

            std::stack<Z3i::Point> pstack;
            pstack.push(*(image.domain().begin()));
            FATAL_ERROR_MSG(image(pstack.top()) == 0, "Starting point of the domain must be equal to zero.");

            // 6-Pencil for the exterior propagation
            std::vector<Z3i::Point> pencil6 = {{1, 0, 0}, {0, 1, 0}, {0, 0, 1}, {-1, 0, 0}, {0, -1, 0}, {0, 0, -1}};

            while (!pstack.empty())
            {
                Z3i::Point p = pstack.top();
                pstack.pop();
                imageFlag.setValue(p, true);

                for (auto &delta : pencil6)
                    if ((image.domain().isInside(p + delta)) &&
                        (imageFlag(p + delta) == false))
                        pstack.push(p + delta);
            }

            for (auto &p : image.domain())
                if ((image(p) == 0) && (!imageFlag(p)))
                    image.setValue(p, fillValue);
        }

        // From DgtalTools
        template <unsigned int SEP>
        DGtal::ImageContainerBySTLVector<DGtal::Z3i::Domain, unsigned char>
        voxelizeAndExport(const std::string &inputFilename,
                          const unsigned int resolution,
                          const unsigned int margin,
                          const unsigned char fillVal,
                          std::pair<DGtal::Z3i::RealPoint, DGtal::Z3i::RealPoint> &meshBbox)
        {
            using namespace DGtal;

            using Domain = Z3i::Domain;
            using PointR3 = Z3i::RealPoint;
            using PointZ3 = Z3i::Point;

            Mesh<PointR3> inputMesh;
            inputMesh << inputFilename.c_str();
            meshBbox = inputMesh.getBoundingBox();

            const double smax = (meshBbox.second - meshBbox.first).max();
            if (resolution != 0)
            {
                const double factor = resolution / smax;
                const PointR3 translate = -meshBbox.first;
                for (auto it = inputMesh.vertexBegin(), itend = inputMesh.vertexEnd();
                     it != itend; ++it)
                {
                    // scale + translation
                    *it += translate;
                    *it *= factor;
                }
            }
            Domain aDomain(PointZ3().diagonal(-(int)margin), PointZ3().diagonal(resolution + margin));
            if (resolution == 0)
            {
                aDomain = Domain(meshBbox.first, meshBbox.second);
            }

            // Digitization step
            Z3i::DigitalSet mySet(aDomain);
            MeshVoxelizer<Z3i::DigitalSet, SEP> aVoxelizer;
            aVoxelizer.voxelize(mySet, inputMesh, 1.0);

            // Export the digital set to a vol file
            ImageContainerBySTLVector<Domain, unsigned char> image(aDomain);
            for (auto &p : mySet)
                image.setValue(p, fillVal);

            fillInterior(image, fillVal);

            return image;
        }

        template <unsigned int SEP>
        DGtal::ImageContainerBySTLVector<DGtal::Z3i::Domain, unsigned char>
        voxelizeAndExport(const std::string &inputFilename,
                          const unsigned int resolution,
                          const unsigned int margin,
                          const unsigned char fillVal)
        {
            std::pair<DGtal::Z3i::RealPoint, DGtal::Z3i::RealPoint> meshBbox;
            return voxelizeAndExport(inputFilename, resolution, margin, fillVal, meshBbox);
        }

        struct DGtal_is_not_background_image_predicate
        {
            typedef DGtal_types::Image::Point Point;
            typedef DGtal_types::Image::Value Value;
            const DGtal_types::Image& image;
            const Value label;
            DGtal_is_not_background_image_predicate(const DGtal_types::Image& image)
            : image(image), label(image(Point(0,0,0)))
            { }

            bool operator()(const Point& p) const
            {
                return image(p) != label;
            }
        };

    }; // namespace internal

    DGtal_types::Image create_dgtal_image(const int& xdim, const int& ydim, const int& zdim)
    {
        DGtal_types::Image::Point firstPoint = DGtal_types::Image::Point::zero;
        DGtal_types::Image::Point lastPoint;
        lastPoint[0] = xdim - 1;
        lastPoint[1] = ydim - 1;
        lastPoint[2] = zdim - 1;
        DGtal_types::Domain domain(firstPoint, lastPoint);
        return DGtal_types::Image(domain);
    }

    DGtal_types::Image cgal_image_to_dgtal_image(const CGAL::Image_3 &src_image)
    {
        std::size_t w = src_image.xdim();
        std::size_t h = src_image.ydim();
        std::size_t e = src_image.zdim();

        DGtal_types::Image image = create_dgtal_image(w, h, e);
        const CGAL_types::Image *im = src_image.image();

        for (unsigned int z = 0; z < e; z++)
        {
            for (unsigned int y = 0; y < h; y++)
            {
                for (unsigned int x = 0; x < w; x++)
                {
                    typename DGtal_types::Image::Point pt;
                    pt[0] = x;
                    pt[1] = y;
                    pt[2] = z;
                    image.setValue(pt, CGAL::IMAGEIO::static_evaluate<unsigned char>(im, x, y, z));
                }
            }
        }

        return image;
    }

    void get_surface_from_image(const DGtal_types::Image &src_image,
                                DGtal_types::Surface_info &surface_info,
                                const bool &compute_normals = true)
    {
        typedef DGtal_types::SH::BinaryImage BinaryImage;
        bBlock("transform image -> build digital surface -> estimate II normals.");
        DGtal::Parameters params = DGtal_types::SH::defaultParameters() | DGtal_types::SHG::defaultParameters();
        DGtal_types::Domain shapeDomain = src_image.domain();
        DGtal::CountedPtr<BinaryImage> bimage(new BinaryImage(shapeDomain));
        BinaryImage::iterator it = bimage->begin();
        for (std::size_t i = 0, size = src_image.size(); i < size; i++)
        {
            *(it++) = src_image.at(i);
        }

        DGtal_types::KSpace K = DGtal_types::SH::getKSpace(bimage, params);
        surface_info.surface = DGtal_types::SH::makeDigitalSurface(bimage, K, params);
        surface_info.surfels = DGtal_types::SH::getSurfelRange(surface_info.surface, params);
        surface_info.pointels = DGtal_types::SH::getCellRange(surface_info.surface, 0);
        if (compute_normals)
            surface_info.normals = DGtal_types::SHG::getIINormalVectors(bimage, surface_info.surfels, params);
        surface_info.uembedder = DGtal_types::SH::getCellEmbedder(K);
        surface_info.embedder = DGtal_types::SH::getSCellEmbedder(K);
        eBlock();
    }

    // normals are per surfaces
    // values are per valuePoints
    void compute_AT_normals_and_values_from_surface(
        const DGtal_types::Surface_info &surface_info,
        std::vector<DGtal_types::RealPoint> &surfaces,
        std::vector<DGtal_types::RealPoint> &valuePoints,
        std::vector<DGtal_types::RealVector> &normals,
        std::vector<double> &values,
        const double alpha = 0.05,
        const double lambda = 0.02,
        const double start_width = 2.0,
        const double end_width = 0.25,
        const double width_divide = 2.0)
    {
        bBlock("Creating AT solver for digital surface");
        typedef DGtal::DiscreteExteriorCalculusFactory<DGtal::EigenLinearAlgebraBackend> CalculusFactory;
        const auto calculus = CalculusFactory::createFromNSCells<2>(surface_info.surfels.begin(), surface_info.surfels.end());
        DGtal::ATSolver2D<DGtal_types::KSpace> at_solver(calculus, 1);
        at_solver.initInputVectorFieldU2(surface_info.normals, surface_info.surfels.cbegin(), surface_info.surfels.cend());
        // TODO check normal inversions
        at_solver.setUp(alpha, lambda);
        at_solver.solveGammaConvergence(start_width, end_width, width_divide);
        eBlock();

        bBlock("output AT normals and values");
        normals.resize(surface_info.normals.size());
        at_solver.getOutputVectorFieldU2(normals, surface_info.surfels.cbegin(), surface_info.surfels.cend());
        surfaces.resize(surface_info.surfels.size());
        std::transform(surface_info.surfels.cbegin(), surface_info.surfels.cend(), surfaces.begin(),
                       [surface_info](const DGtal_types::SH::SCell &c)
                       { return surface_info.embedder(c); });

        DGtal_types::SH::Scalars features(surface_info.pointels.size());
        at_solver.getOutputScalarFieldV0(features, surface_info.pointels.cbegin(), surface_info.pointels.cend(),
                                         at_solver.Maximum);

        valuePoints.resize(surface_info.pointels.size());
        std::transform(surface_info.pointels.cbegin(), surface_info.pointels.cend(), valuePoints.begin(),
                       [surface_info](const DGtal_types::SH::Cell &v)
                       { return surface_info.uembedder(v); });
        values.resize(features.size());
        std::transform(features.cbegin(), features.cend(), values.begin(),
                       [](const double &v) { return v; });
        eBlock();
    }

    void compute_integral_invariant_values(
        const DGtal_types::Image &src_image,
        std::vector<DGtal_types::RealPoint> &surfaces,
        std::vector<DGtal_types::RealPoint> &valuePoints,
        std::vector<double> &values,
        double radius = 7.0, double grid_step = 1.0, float threshold = 0.5)
    {
        bBlock("split image");
        std::vector<DGtal_types::Image> separated_images = internal::split_in_objects_dgtal_image(src_image, src_image(DGtal_types::Image::Point(0, 0, 0)));
        message(std::to_string(separated_images.size()) + " images");
        eBlock();
        for (const DGtal_types::Image &image : separated_images)
        {
            internal::compute_integral_invariant_values_single_object_image(
                /*In */ image,
                /*Out */ surfaces, valuePoints, values,
                /*Optional In*/ radius, grid_step, threshold);
        }
    }

    void compute_normals_from_surface(
        const DGtal_types::Image &src_image,
        std::vector<DGtal_types::RealPoint> &surfaces,
        std::vector<DGtal_types::RealVector> &normals)
    {
        DGtal::Parameters params = DGtal_types::SH::defaultParameters() | DGtal_types::SHG::defaultParameters();
        DGtal_types::Domain shapeDomain = src_image.domain();
        DGtal::CountedPtr<DGtal_types::SH::BinaryImage> bimage(new DGtal_types::SH::BinaryImage(shapeDomain));
        DGtal_types::SH::BinaryImage::iterator it = bimage->begin();
        for (std::size_t i = 0, size = src_image.size(); i < size; i++)
        {
            *(it++) = src_image.at(i);
        }
        DGtal_types::KSpace K = DGtal_types::SH::getKSpace(bimage, params);
        DGtal::CanonicSCellEmbedder<DGtal_types::KSpace> embedder = DGtal_types::SH::getSCellEmbedder(K);
        DGtal::CountedPtr<DGtal_types::SH::DigitalSurface> surface = DGtal_types::SH::makeDigitalSurface(bimage, K, params);
        DGtal_types::SH::SurfelRange dg_surfaces = DGtal_types::SH::getSurfelRange(surface, params);
        normals = DGtal_types::SHG::getIINormalVectors(bimage, dg_surfaces, params);
        surfaces.resize(dg_surfaces.size());
        std::transform(dg_surfaces.cbegin(), dg_surfaces.cend(), surfaces.begin(),
                       [embedder](const DGtal_types::SH::SCell &c)
                       { return embedder(c); });
    }

    void compute_VCM_values(
        const DGtal_types::Image &src_image,
        std::vector<DGtal_types::RealPoint> &surfaces,
        std::vector<DGtal_types::RealPoint> &valuePoints,
        std::vector<double> &values,
        double offset_radius, double smooth_radius)
    {
        bBlock("split image");
        std::vector<DGtal_types::Image> separated_images = internal::split_in_objects_dgtal_image(src_image, src_image(DGtal_types::Image::Point(0, 0, 0)));
        message(std::to_string(separated_images.size()) + " images");
        eBlock();
        for (const DGtal_types::Image &image : separated_images)
        {
            internal::compute_VCM_values_single_object_image(
                /*In */ image,
                /*Out */ surfaces, valuePoints, values,
                /*Optional In*/ offset_radius, smooth_radius);
        }
    }

    template <typename DGTalImage>
    CGAL::Image_3 dgtal_image_to_cgal_image(const DGTalImage &src_image)
    {
        typename DGTalImage::Domain domain = src_image.domain();
        typename DGTalImage::Vector extent = src_image.extent();
        std::size_t xdim = extent[0];
        std::size_t ydim = extent[1];
        std::size_t zdim = extent[2];

        CGAL_types::Image *im = CGAL::_createImage(xdim, ydim, zdim, 1,
                                  1, 1, 1, 1,
                                  CGAL::WK_FIXED, CGAL::SGN_UNSIGNED);
        std::fill_n(static_cast<unsigned char *>(im->data),
                    std::size_t(im->xdim) * std::size_t(im->ydim) * std::size_t(im->zdim),
                    static_cast<unsigned char>(0));

        for (unsigned int z = 0; z < zdim; z++)
        {
            for (unsigned int y = 0; y < ydim; y++)
            {
                for (unsigned int x = 0; x < xdim; x++)
                {
                    typename DGTalImage::Point pt(x, y, z);
                    pt += domain.lowerBound();
                    CGAL::IMAGEIO::static_evaluate<unsigned char>(im, x, y, z) = static_cast<unsigned char>(src_image(pt));
                }
            }
        }

        return CGAL::Image_3(im);
    }

    CGAL::Image_3 make_cgal_image_from_mesh(
        const std::string &input_mesh,
        const std::size_t &resolution,
        const std::size_t margin = 1,
        const unsigned char fill_value = 128)
    {
        // mesh2vol
        std::pair<DGtal_types::RealPoint, DGtal_types::RealPoint> meshBbox;
        auto STLImage = internal::voxelizeAndExport<26>(input_mesh, resolution, margin, fill_value, meshBbox);
        // make cgal
        CGAL::Image_3 cgal_image = dgtal_image_to_cgal_image(STLImage);
        // apply mesh bbox
        std::array<std::size_t, 6> cgal_image_bbox = get_voxels_with_condition_bbox(cgal_image, internal::Is_not_label_condition(0));
        CGAL_types::Image *im = cgal_image.image();
        im->tx = meshBbox.first[0];
        im->ty = meshBbox.first[1];
        im->tz = meshBbox.first[2];
        DGtal::Z3i::RealPoint meshDiag = meshBbox.second - meshBbox.first;
        im->vx = meshDiag[0] / static_cast<float>(cgal_image_bbox[3] - cgal_image_bbox[0]);
        im->vy = meshDiag[1] / static_cast<float>(cgal_image_bbox[4] - cgal_image_bbox[1]);
        im->vz = meshDiag[2] / static_cast<float>(cgal_image_bbox[5] - cgal_image_bbox[2]);
        // cut image to its bbox
        std::array<int, 6> cgal_image_bbox_with_margin = {
            static_cast<int>(cgal_image_bbox[0]) - static_cast<int>(2 * margin),
            static_cast<int>(cgal_image_bbox[1]) - static_cast<int>(2 * margin),
            static_cast<int>(cgal_image_bbox[2]) - static_cast<int>(2 * margin),
            static_cast<int>(cgal_image_bbox[3]) + static_cast<int>(2 * margin),
            static_cast<int>(cgal_image_bbox[4]) + static_cast<int>(2 * margin),
            static_cast<int>(cgal_image_bbox[5]) + static_cast<int>(2 * margin)};
        return crop_image(cgal_image, cgal_image_bbox_with_margin, 0);
    }

    // Add noise then filter the image
    //      Noise: a voxel is inverted with the probability of pow(a, distance_to_border) (0<a<1)
    //      Filter: floating groups of voxels that have less than voxel_count_threshold voxels are replaced
    CGAL::Image_3 make_noise_in_image(
        const CGAL::Image_3 &image,
        const double &noise_amount,
        const unsigned int voxel_count_threshold = 27)
    {
        helper::bBlock("add noise with DGtal");
        DGtal_types::Image dgtal_image = cgal_image_to_dgtal_image(image);
        auto params = DGtal_types::SH::defaultParameters();
        internal::DGtal_is_not_background_image_predicate image_predicate(dgtal_image);
        typedef DGtal::KanungoNoise<internal::DGtal_is_not_background_image_predicate, DGtal_types::Domain> KanungoPredicate;
        const DGtal_types::Domain shapeDomain = dgtal_image.domain();
        KanungoPredicate noisy_dshape(image_predicate, shapeDomain, noise_amount);
        std::transform(shapeDomain.begin(), shapeDomain.end(),
                        dgtal_image.begin(),
                        [&noisy_dshape](const DGtal_types::Image::Point &p)
                        { return noisy_dshape(p); });
        helper::eBlock();
        helper::bBlock("filter noise to remove floating voxels");
        CGAL::Image_3 cgal_noise_image = dgtal_image_to_cgal_image(dgtal_image);
        CGAL::Image_3 result_image = remove_floating_voxels(cgal_noise_image, voxel_count_threshold);
        // put the image in the same world space
        result_image.tx() = image.tx();
        result_image.ty() = image.ty();
        result_image.tz() = image.tz();
        result_image.image()->vx = image.vx();
        result_image.image()->vy = image.vy();
        result_image.image()->vz = image.vz();

        helper::eBlock();
        return result_image;
    }

} // namespace helper

#endif // DGTAL_HELPER_H
