#ifndef POLYLINE_EXTRACTOR_H
#define POLYLINE_EXTRACTOR_H

#include "../polyline/Polylines.h"
#include "../surface/Surface_graph.h"
#include "../element_graph/Element_graph.h"
#include "../element_graph/Element_graph_process.h"
#include "../helper/image_helper.h"

struct Sharpness_oriented_thinning_on_element_graph
{
    std::string out_filepath = "";

    template <typename Point, typename Vector>
    Polylines<Point> operator()(
        const Surface_graph<Vector> &sharp_surface_graph,
        const std::vector<const typename Surface_graph<Vector>::Surfel *> &deleted_surfels,
        const std::map<std::array<int, 3>, double> &pointel_values) const
    {
        using Polylines = Polylines<Point>;
        using Element_graph = Element_graph<Point>;
        using Graph_element = Element_graph::Graph_element;

        const boost::filesystem::path out_folder = out_filepath;

        helper::bBlock("Build element graph");
        std::shared_ptr<Element_graph> element_graph_ptr =
            Element_graph_process::build_element_graph_from_surfels<Point, Vector>(sharp_surface_graph, deleted_surfels, pointel_values);
        Element_graph &element_graph = *element_graph_ptr;
        helper::eBlock();
        if (!out_filepath.empty())
        {
            Element_graph_process::write_graph_with_infos(
                boost::filesystem::path(out_folder).append("element_graph").string(), element_graph);
        }

        helper::bBlock("Thin graph");
        Element_graph_process::thin_element_graph_with_sharpness(element_graph, Element_graph_process::internal::Invert_value_priority());
        Polylines thinned_graph = Element_graph_process::create_polyline_from_graph<Point>(element_graph);
        helper::eBlock();

        return thinned_graph;
    }
};

#endif // POLYLINE_EXTRACTOR_H