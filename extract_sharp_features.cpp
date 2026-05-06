#define VERBOSE

#include <boost/filesystem.hpp>
#include <boost/program_options.hpp>

#include <CGAL/assertions.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <CGAL/Image_3.h>
#include <CGAL/Polygon_mesh_processing/IO/polygon_mesh_io.h>
#include <CGAL/Surface_mesh/IO/PLY.h>

#include "sharpness_measure.h"
#include "polyline/Polylines.h"
#include "helper/cgal_helper.h"
#include "sharp_feature_extraction.h"
#include "accuracy_measure/accuracy_measure.h"

// #include "Surface/Surface_process.h"
// #include "sharp_edge_extraction.h"

struct Input_parameters
{
    std::vector<std::string> filenames;
    std::string filepath = "./";
    std::string output_filepath = "./";
    int mesh_image_resolution = 512;
    double extract_line_angle = 30.0;
    bool add_noise = false;
    double noise_amount;

    boost::program_options::options_description make_option_descriptor()
    {
        namespace po = boost::program_options;
        po::options_description desc("Input options");
        auto op = desc.add_options();
        op("input,i", po::value<std::vector<std::string>>(&filenames), "the input files, either images (.inr) or meshes (.off)");
        op("extract_line_angle", po::value<double>(&extract_line_angle)->default_value(extract_line_angle), "if the input are meshes, specify the angle used for sharp line extraction used as ground truth lines in the accuracy measurement");
        // op("add_noise", po::value<bool>(&add_noise)->default_value(add_noise), "true if noise should be added to the image");
        op("noise_amount", po::value<double>(&noise_amount), "add Kanungo noise to the image");
        op("filepath", po::value<std::string>(&filepath)->default_value(filepath), "set the input files location");
        op("output_filepath", po::value<std::string>(&output_filepath), "set the output files location. One folder is created for each input");
        op("mesh_image_resolution", po::value<int>(&mesh_image_resolution)->default_value(mesh_image_resolution), "if the input is a mesh, specify the resolution of the image created");
        return desc;
    }

    void process_arguments(const boost::program_options::variables_map &vm)
    {
        if (vm.count("noise_amount"))
        {
            add_noise = true;
        }
    }
};

struct Sharpness_ambrosio_tortorelli_parameters
{
    bool active = true;
    double alpha = 0.05;
    double lambda = 0.02;
    double width_start = 2.0;
    double width_end = 0.5;
    double width_divide = 2.0;

    boost::program_options::options_description make_option_descriptor()
    {
        namespace po = boost::program_options;
        po::options_description desc("Amrosio-Tortorelli sharpness options");
        auto op = desc.add_options();
        op("ambrosio_tortorelli,a", po::value<bool>(&active)->default_value(active), "is the Ambrosio-Tortorelli measure active");
        op("at_alpha", po::value<double>(&alpha)->default_value(alpha), "set alpha, the normal field approximation factor");
        op("at_lambda", po::value<double>(&lambda)->default_value(lambda), "set lambda, the discontinuity edges length factor");
        op("at_width", po::value<double>(&width_end)->default_value(width_end), "set epsilon, the final width of the discontinuity edges");
        op("at_width_start", po::value<double>(&width_start)->default_value(width_start), "set the starting width of the discontinuity edges");
        op("at_width_divide", po::value<double>(&width_divide)->default_value(width_divide), "set the dividing factor for the width of the discontinuity edges");
        return desc;
    }

    void process_arguments(const boost::program_options::variables_map & /*vm*/) {}
};

struct Sharpness_curvature_parameters
{
    bool active = false;
    double radius = 7.0;
    double grid_step = 1.0;

    boost::program_options::options_description make_option_descriptor()
    {
        namespace po = boost::program_options;
        po::options_description desc("Integral invariant curvature sharpness options");
        auto op = desc.add_options();
        op("curvature,c", po::value<bool>(&active)->default_value(active), "is the Integral invariant curvature measure active");
        op("c_radius", po::value<double>(&radius)->default_value(radius), "set the radius of the measure kernel");
        op("c_grid_step", po::value<double>(&grid_step)->default_value(grid_step), "set the grid step");
        return desc;
    }

    void process_arguments(const boost::program_options::variables_map & /*vm*/) {}
};

struct Sharpness_voronoi_covariance_measure_parameters
{
    bool active = false;
    double offset_radius = 5.0;
    double smooth_radius = 4.0;

    boost::program_options::options_description make_option_descriptor()
    {
        namespace po = boost::program_options;
        po::options_description desc("Voronoi covariance measure sharpness options");
        auto op = desc.add_options();
        op("voronoi_covariance,v", po::value<bool>(&active)->default_value(active), "is the Voronoi covariance measure active");
        op("vcm_offset_radius", po::value<double>(&offset_radius)->default_value(offset_radius), "set the offset radius");
        op("vcm_smooth_radius", po::value<double>(&smooth_radius)->default_value(smooth_radius), "set the smooth radius");
        return desc;
    }

    void process_arguments(const boost::program_options::variables_map & /*vm*/) {}
};

struct Sharpness_feedback_ambrosio_tortorelli_parameters
{
    bool active = false;
    double alpha = 0.05;
    double lambda = 0.02;
    double width_start = 2.0;
    double width_end = 0.5;
    double width_divide = 2.0;
    double refine_surfel_distance = 4.0;

    boost::program_options::options_description make_option_descriptor()
    {
        namespace po = boost::program_options;
        po::options_description desc("Feedback with Amrosio-Tortorelli sharpness options");
        auto op = desc.add_options();
        op("feedback,f", po::value<bool>(&active)->default_value(active), "is the feedback measure active");
        op("feedback_at_alpha", po::value<double>(&alpha)->default_value(alpha), "set alpha, the normal field approximation factor");
        op("feedback_at_lambda", po::value<double>(&lambda)->default_value(lambda), "set lambda, the discontinuity edges length factor");
        op("feedback_at_width", po::value<double>(&width_end)->default_value(width_end), "set epsilon, the final width of the discontinuity edges");
        op("feedback_at_width_start", po::value<double>(&width_start)->default_value(width_start), "set the starting width of the discontinuity edges");
        op("feedback_at_width_divide", po::value<double>(&width_divide)->default_value(width_divide), "set the dividing factor for the width of the discontinuity edges");
        op("feedback_at_refine_surfel_distance", po::value<double>(&refine_surfel_distance)->default_value(refine_surfel_distance), "set the distance to the lines where the normal of surfels are refined");
        return desc;
    }

    void process_arguments(const boost::program_options::variables_map & /*vm*/) {}
};

struct Optimization_parameters
{
    std::size_t maximum_iteration = 20;
    double start_step_size = 1.0;
    double end_step_size = 0.125;
    double min_energy_delta = 1.e-3;
    double collapse_distance = 0.5;
    double smooth_factor = 1.0;
    double refine_surfel_distance = 4.0;
    double plane_detection_max_distance = 4.0;

    boost::program_options::options_description make_option_descriptor()
    {
        namespace po = boost::program_options;
        po::options_description desc("Optimization options");
        auto op = desc.add_options();
        op("iteration", po::value<std::size_t>(&maximum_iteration)->default_value(maximum_iteration), "set the maximal number of iteration");
        op("start_step_size", po::value<double>(&start_step_size)->default_value(start_step_size), "set the starting step size of the gradient descent");
        op("end_step_size", po::value<double>(&end_step_size)->default_value(end_step_size), "set the ending step size of the gradient descent");
        op("min_energy_delta", po::value<double>(&min_energy_delta)->default_value(min_energy_delta), "set the minimum energy delta to stop the iterations");
        op("collapse_distance", po::value<double>(&collapse_distance)->default_value(collapse_distance), "set collapse distance in the optimization");
        op("refine_surfel_distance", po::value<double>(&refine_surfel_distance)->default_value(refine_surfel_distance), "set the distance to the lines where the normal of surfels are refined");
        op("smooth_factor", po::value<double>(&smooth_factor)->default_value(smooth_factor), "set the smoothing factor of the energy");
        op("plane_detection_max_distance", po::value<double>(&plane_detection_max_distance)->default_value(plane_detection_max_distance), "set the search maximum distance when creating patches of surfel for plane detection");
        return desc;
    }

    void process_arguments(const boost::program_options::variables_map & /*vm*/) {}
};

struct Mesh_generation_parameters
{
    bool active = true;
    double triangle_size = 0.02;
    double edge_size = 0.02;
    double edge_min_size = 0.001;

    boost::program_options::options_description make_option_descriptor()
    {
        namespace po = boost::program_options;
        po::options_description desc("Mesh generation options");
        auto op = desc.add_options();
        op("mesh", po::value<bool>(&active)->default_value(active), "if true, two tetrahedral meshes will be generated with and without the sharp edges");
        op("triangle_size", po::value<double>(&triangle_size)->default_value(triangle_size), "set the triangle size as a factor of the diagonal of the axis aligned bounding box");
        op("edge_size", po::value<double>(&edge_size)->default_value(edge_size), "set the sharp edges size as a factor of the diagonal of the axis aligned bounding box");
        op("edge_min_size", po::value<double>(&edge_min_size)->default_value(edge_min_size), "set the minimum size of the sharp edges");
        return desc;
    }

    void process_arguments(const boost::program_options::variables_map & /*vm*/) {}
};

struct Program_parameters
{
    Input_parameters input_params;

    Sharpness_ambrosio_tortorelli_parameters at_params;
    Sharpness_curvature_parameters curv_params;
    Sharpness_voronoi_covariance_measure_parameters vcm_params;
    Sharpness_feedback_ambrosio_tortorelli_parameters feedback_at_params;

    double selection_threshold = 0.25;

    double regularization_distance = 4.0;

    Optimization_parameters op_params;

    Mesh_generation_parameters mesh_params;

    bool should_measure_line_accuracy = true;
    bool should_measure_mesh_accuracy = true;

    boost::program_options::options_description make_option_descriptor()
    {
        namespace po = boost::program_options;
        po::options_description desc("Program options");
        auto op = desc.add_options();
        desc.add(input_params.make_option_descriptor());
        desc.add(at_params.make_option_descriptor());
        desc.add(curv_params.make_option_descriptor());
        desc.add(vcm_params.make_option_descriptor());
        desc.add(feedback_at_params.make_option_descriptor());
        desc.add(op_params.make_option_descriptor());
        desc.add(mesh_params.make_option_descriptor());
        op("help,h", "produce help message");
        op("selection_threshold", po::value<double>(&selection_threshold)->default_value(selection_threshold), "set the threshold for the selection step");
        op("regularization_distance", po::value<double>(&regularization_distance)->default_value(regularization_distance), "set the distance to delete lines in the regularization step, proportional to the maximum voxel edge length");
        op("measure_line_accuracy", po::value<bool>(&should_measure_line_accuracy)->default_value(should_measure_line_accuracy), "should the line accuracy be measured (only for meshes as input)");
        op("measure_mesh_accuracy", po::value<bool>(&should_measure_mesh_accuracy)->default_value(should_measure_mesh_accuracy), "should the mesh accuracy be measured (only for meshes as input)");
        return desc;
    }

    void process_arguments(const boost::program_options::variables_map &vm)
    {
        input_params.process_arguments(vm);
        at_params.process_arguments(vm);
        curv_params.process_arguments(vm);
        vcm_params.process_arguments(vm);
        op_params.process_arguments(vm);
        mesh_params.process_arguments(vm);
    }
};

bool setup_commands(int &argc, char *argv[], Program_parameters &parameters)
{
    namespace po = boost::program_options;

    try
    {
        // Declare the supported options.
        po::options_description desc = parameters.make_option_descriptor();

        po::positional_options_description pos_desc;
        pos_desc.add("input", -1);

        po::variables_map vm;
        po::store(po::command_line_parser(argc, argv).options(desc).positional(pos_desc).run(), vm);
        po::notify(vm);

        parameters.process_arguments(vm);

        if (vm.count("help"))
        {
            std::cout << desc << std::endl;
        }
    }
    catch (std::exception &e)
    {
        std::cout << e.what() << "\n";
        return false;
    }
    return true;
}

// holds ground truth informations for the accuracy measures
struct Accuracy_measure_structures
{
    Polylines<helper::CGAL_types::Point> lines;
    helper::CGAL_types::Mesh mesh;
};

void display_line_accuracy(
    const CGAL::Image_3 &image,
    const Accuracy_measure_structures &ground_truth_structures,
    const Polylines<helper::CGAL_types::Point> &feature_graph,
    const double voxel_distance_too_far = 2.0)
{
    // get bbox
    const CGAL::Bbox_3& bbox = helper::get_bbox_from_mesh(ground_truth_structures.mesh);
    const double distance_too_far = voxel_distance_too_far / (std::max(std::max(image.xdim(), image.ydim()), image.zdim()));
    CGAL::Bbox_3 normalized_bbox = helper::get_normalized_bbox_with_same_ratio(bbox);
    // gt poly
    Polylines<helper::CGAL_types::Point> gt_poly = ground_truth_structures.lines.clone();
    helper::scale_features_to_range(gt_poly, bbox, normalized_bbox);
    // method poly
    Polylines<helper::CGAL_types::Point> method_poly = feature_graph.clone(); // copy poly
    helper::scale_features_to_range(method_poly, bbox, normalized_bbox);
    std::size_t nb_line_not_in_gt = Accuracy_measure::number_of_lines_in_first_not_in_second<helper::CGAL_types::Vector>(gt_poly, method_poly, distance_too_far);
    std::size_t nb_line_not_in_method = Accuracy_measure::number_of_lines_in_first_not_in_second<helper::CGAL_types::Vector>(method_poly, gt_poly, distance_too_far);
    // compare
    double min_error, max_error;
    double error = Accuracy_measure::get_polyline_error<helper::CGAL_types::Vector>(gt_poly, method_poly, min_error, max_error);
    // output min max and average distances
    std::cout << "feature graph error to mesh sharp lines :" << std::endl;
    std::cout << "      min error = " << min_error << std::endl;
    std::cout << "      max error = " << max_error << std::endl;
    std::cout << "      avg error = " << error << std::endl;
    std::cout << "      number of lines in feature graph but not in mesh sharp lines = " << nb_line_not_in_gt << std::endl;
    std::cout << "      number of lines in mesh sharp lines but not in feature graph = " << nb_line_not_in_method << std::endl;
}

void display_mesh_accuracy(
    const CGAL::Image_3 &image,
    const std::string &out_filepath,
    const Accuracy_measure_structures &ground_truth_structures,
    const Polylines<helper::CGAL_types::Point> &feature_graph)
{
    helper::bBlock("measure mesh accuracy");
    const boost::filesystem::path out_folder = out_filepath;
    helper::CGAL_types::Mesh gt_mesh = ground_truth_structures.mesh; // copy mesh
    const CGAL::Bbox_3& gt_bbox = helper::get_bbox_from_mesh(ground_truth_structures.mesh);
    CGAL::Bbox_3 normalized_bbox = helper::get_normalized_bbox_with_same_ratio(gt_bbox);
    helper::scale_mesh_to_range(gt_mesh, gt_bbox, normalized_bbox); // normalize the bbox of the ground truth mesh to have comparable data
    if (!out_filepath.empty())
    {
        std::ofstream output(boost::filesystem::path(out_folder).append("mesh_ground_truth.ply").string());
        CGAL::IO::write_PLY(output, gt_mesh);
    }
    helper::bBlock("make meshes");
    helper::CGAL_types::Mesh mesh_without_features = Accuracy_measure::make_mesh_from_image(image);
    helper::scale_mesh_to_range(mesh_without_features, gt_bbox, normalized_bbox); // apply the same scaling as the ground truth mesh
    if (!out_filepath.empty())
    {
        std::ofstream output(boost::filesystem::path(out_folder).append("mesh_without_features.ply").string());
        CGAL::IO::write_PLY(output, mesh_without_features);
    }
    helper::CGAL_types::Mesh mesh_with_features = Accuracy_measure::make_mesh_from_image(image, feature_graph);
    helper::scale_mesh_to_range(mesh_with_features, gt_bbox, normalized_bbox); // apply the same scaling as the ground truth mesh
    if (!out_filepath.empty())
    {
        std::ofstream output(boost::filesystem::path(out_folder).append("mesh_with_features.ply").string());
        CGAL::IO::write_PLY(output, mesh_with_features);
    }
    helper::eBlock();
    helper::bBlock("measure accuracy");
    double min_value_without_features;
    double max_value_without_features;
    double mesh_error_without_features;
    double min_value_with_features;
    double max_value_with_features;
    double mesh_error_with_features;
    if (!out_filepath.empty())
    {
        mesh_error_without_features = Accuracy_measure::get_mesh_error_and_output(mesh_without_features, gt_mesh,
            boost::filesystem::path(out_folder).append("mesh_error_without_features").string(),
            min_value_without_features, max_value_without_features);
        mesh_error_with_features = Accuracy_measure::get_mesh_error_and_output(mesh_with_features, gt_mesh,
            boost::filesystem::path(out_folder).append("mesh_error_with_features").string(),
            min_value_with_features, max_value_with_features,
            min_value_without_features, max_value_without_features); // output with the same color scale
    }
    else
    {
        mesh_error_without_features = Accuracy_measure::get_mesh_error(mesh_without_features, gt_mesh,
            min_value_without_features, max_value_without_features);
        mesh_error_with_features = Accuracy_measure::get_mesh_error(mesh_with_features, gt_mesh,
            min_value_with_features, max_value_with_features);
    }
    helper::eBlock();
    // output min max and average distance of sampled points
    std::cout << "mesh error without features :" << std::endl;
    std::cout << "      min error = " << min_value_without_features << std::endl;
    std::cout << "      max error = " << max_value_without_features << std::endl;
    std::cout << "      avg error = " << mesh_error_without_features << std::endl;
    std::cout << "mesh error with features : " << std::endl;
    std::cout << "      min error = " << min_value_with_features << std::endl;
    std::cout << "      max error = " << max_value_with_features << std::endl;
    std::cout << "      avg error = " << mesh_error_with_features << std::endl;
    helper::eBlock();
}

template <typename SharpnessMeasure>
Polylines<helper::CGAL_types::Point> process_image_with_measure(
    const CGAL::Image_3 &image,
    const std::string &out_filepath,
    const Program_parameters &parameters,
    const Accuracy_measure_structures &ground_truth_structures,
    const SharpnessMeasure &sharpness_measure = SharpnessMeasure())
{
    // extract features
    helper::bBlock("Extract sharp feature from the image");
    Polylines<helper::CGAL_types::Point> feature_graph =
        detect_sharp_edges_in_image(image, sharpness_measure,
                                parameters.selection_threshold,
                                parameters.regularization_distance,
                                parameters.op_params.maximum_iteration,
                                parameters.op_params.start_step_size,
                                parameters.op_params.end_step_size,
                                parameters.op_params.min_energy_delta,
                                parameters.op_params.collapse_distance,
                                parameters.op_params.smooth_factor,
                                parameters.op_params.refine_surfel_distance,
                                parameters.op_params.plane_detection_max_distance,
                                out_filepath);
    helper::eBlock();
    if (parameters.should_measure_line_accuracy)
    {
        display_line_accuracy(image, ground_truth_structures, feature_graph);
    }
    if (parameters.should_measure_mesh_accuracy)
    {
        display_mesh_accuracy(image, out_filepath, ground_truth_structures, feature_graph);
    }
    return feature_graph;
}

template <typename SharpnessMeasure>
void process_using_measure(
    const CGAL::Image_3 &image,
    const std::string &out_filepath,
    const Program_parameters &parameters,
    const Accuracy_measure_structures &ground_truth_structures,
    const SharpnessMeasure &sharpness_measure = SharpnessMeasure())
{
    const boost::filesystem::path out_folder = out_filepath;
#ifdef VERBOSE
    helper::clear_trace_ostream();
    helper::add_trace_ostream(std::cout);
    std::ofstream trace_fd(boost::filesystem::path(out_folder).append("trace.txt").string());
    helper::add_trace_ostream(trace_fd, false);
#endif
    Polylines<helper::CGAL_types::Point> feature_graph = process_image_with_measure(image, out_filepath, parameters, ground_truth_structures, sharpness_measure);

    // feedback loop
    if (parameters.feedback_at_params.active)
    {
        helper::bBlock("Using the feedback loop with Ambrosio-Tortorelli measure");
        const boost::filesystem::path feedback_at_folder = boost::filesystem::path(out_folder).append("feedback_at");
        boost::filesystem::create_directory(feedback_at_folder);
        Feedback_loop_AT_measure feedback_at_sharpness(feature_graph);
        feedback_at_sharpness.alpha = parameters.feedback_at_params.alpha;
        feedback_at_sharpness.lambda = parameters.feedback_at_params.lambda;
        feedback_at_sharpness.width_start = parameters.feedback_at_params.width_start;
        feedback_at_sharpness.width_end = parameters.feedback_at_params.width_end;
        feedback_at_sharpness.width_divide = parameters.feedback_at_params.width_divide;
        feedback_at_sharpness.refine_surfel_distance = parameters.feedback_at_params.refine_surfel_distance;
        process_image_with_measure(image, feedback_at_folder.string(), parameters, ground_truth_structures, feedback_at_sharpness);
        helper::eBlock();
    }

#ifdef VERBOSE
    helper::clear_trace_ostream();
    helper::add_trace_ostream(std::cout);
#endif
}

void extract_sharp_edges(
    const CGAL::Image_3 &image,
    const std::string &out_filepath,
    const Program_parameters &parameters,
    const Accuracy_measure_structures &ground_truth_structures)
{
    const boost::filesystem::path out_folder = out_filepath;
    if (parameters.at_params.active)
    {
        helper::bBlock("Using the Ambrosio-Tortorelli measure");
        const boost::filesystem::path at_folder = boost::filesystem::path(out_folder).append("at");
        boost::filesystem::create_directory(at_folder);
        AT_sharpness_measure at_sharpness;
        at_sharpness.alpha = parameters.at_params.alpha;
        at_sharpness.lambda = parameters.at_params.lambda;
        at_sharpness.width_start = parameters.at_params.width_start;
        at_sharpness.width_end = parameters.at_params.width_end;
        at_sharpness.width_divide = parameters.at_params.width_divide;
        process_using_measure(image, at_folder.string(), parameters, ground_truth_structures, at_sharpness);
        helper::eBlock();
    }
    if (parameters.curv_params.active)
    {
        helper::bBlock("Using the integral invariant measure");
        const boost::filesystem::path cuvature_folder = boost::filesystem::path(out_folder).append("curvature");
        boost::filesystem::create_directory(cuvature_folder);
        Curvature_sharpness_measure c_sharpness;
        c_sharpness.radius = parameters.curv_params.radius;
        c_sharpness.grid_step = parameters.curv_params.grid_step;
        process_using_measure(image, cuvature_folder.string(), parameters, ground_truth_structures, c_sharpness);
        helper::eBlock();
    }
    if (parameters.vcm_params.active)
    {
        helper::bBlock("Using the Voronoi curvature measure");
        const boost::filesystem::path vcm_folder = boost::filesystem::path(out_folder).append("vcm");
        boost::filesystem::create_directory(vcm_folder);
        VCM_sharpness_measure v_sharpness;
        v_sharpness.offset_radius = parameters.vcm_params.offset_radius;
        v_sharpness.smooth_radius = parameters.vcm_params.smooth_radius;
        process_using_measure(image, vcm_folder.string(), parameters, ground_truth_structures, v_sharpness);
        helper::eBlock();
    }
}

void sharp_edge_process_from_file(const std::string &filename, const std::string &out_filepath, const Program_parameters &parameters)
{
    const boost::filesystem::path out_folder = out_filepath;
    CGAL::Image_3 image;
    Accuracy_measure_structures ground_truth_structures;
    // identify input : try to read the image
    helper::bLine("Read image " + filename);
    bool image_read = image.read(filename);
    helper::eLine();
    // if not an image, read the mesh and make the image
    if (!image_read)
    {
        helper::message("Input file is not an image: " + filename);
        helper::bBlock("Constructing an image from the mesh, resolution: " + std::to_string(parameters.input_params.mesh_image_resolution));
        image = helper::make_cgal_image_from_mesh(filename, parameters.input_params.mesh_image_resolution);
        helper::eBlock();
        {
            helper::bBlock("Write image");
            helper::write_image(image, boost::filesystem::path(out_folder).append("image.inr").string());
            helper::eBlock();
        }
        if (parameters.should_measure_mesh_accuracy || parameters.should_measure_line_accuracy)
        {
            helper::bBlock("Read mesh with cgal");
            CGAL::Polygon_mesh_processing::IO::read_polygon_mesh(filename, ground_truth_structures.mesh);
            helper::eBlock();
            if (parameters.should_measure_line_accuracy)
            {
                helper::bBlock("Extract ground truth lines, angle threshold: " + std::to_string(parameters.input_params.extract_line_angle));
                ground_truth_structures.lines = helper::extract_sharp_edges_from_mesh(ground_truth_structures.mesh, parameters.input_params.extract_line_angle);
                helper::eBlock();
                {
                    helper::bBlock("Write ground truth lines");
                    Polylines_process::write_polylines<helper::CGAL_types::Point, helper::CGAL_types::Vector>(ground_truth_structures.lines, boost::filesystem::path(out_folder).append("groud_truth_lines").string());
                    helper::eBlock();
                }
            }
        }
    }

    if (parameters.input_params.add_noise)
    {
        helper::bBlock("Adding noise, value: " + std::to_string(parameters.input_params.noise_amount));
        image = helper::make_noise_in_image(image, parameters.input_params.noise_amount);
        helper::eBlock();
        {
            helper::bBlock("Write image with noise");
            helper::write_image(image, boost::filesystem::path(out_folder).append("image_noise.inr").string());
            helper::eBlock();
        }
    }

    extract_sharp_edges(image, out_filepath, parameters, ground_truth_structures);
}

int main(int argc, char *argv[])
{
    Program_parameters parameters;
    if (!setup_commands(argc, argv, parameters))
        return EXIT_FAILURE;

    for (const std::string &filename_str : parameters.input_params.filenames)
    {
        const boost::filesystem::path user_filename = filename_str;
        const boost::filesystem::path user_filepath;
        boost::filesystem::path input_filepath;
        if (user_filename.is_absolute())
            input_filepath = user_filename;
        else
            input_filepath = boost::filesystem::path(parameters.input_params.filepath).append(user_filename);
        const boost::filesystem::path output_filepath = boost::filesystem::path(parameters.input_params.output_filepath).append(user_filename.stem());
        helper::bBlock("Extract sharp edges from " + input_filepath.string() + ", outputing to " + output_filepath.string());
        boost::filesystem::create_directory(output_filepath);
        sharp_edge_process_from_file(input_filepath.string(), output_filepath.string(), parameters);
#ifdef VERBOSE
        helper::clear_trace_ostream();
        helper::add_trace_ostream(std::cout);
#endif
        helper::eBlock();
    }

    return EXIT_SUCCESS;
}