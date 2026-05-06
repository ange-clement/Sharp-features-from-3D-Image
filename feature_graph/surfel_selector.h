#ifndef SURFEL_SELECTOR_H
#define SURFEL_SELECTOR_H

#include "../surface/Surface_graph.h"
#include "../surface/surface_process.h"

struct Threshold_surfel_selector
{
    double selection_threshold;

    template <typename Vector, typename OutputIterator>
    Surface_graph<Vector> operator()(
        const Surface_graph<Vector>& surface_graph,
        OutputIterator &deleted_surfel_output) const
    {
        typedef Surface_graph<Vector>    Surface_graph;
        typedef Surface_graph::Surfel    Surfel;

        std::set<const Surfel*> sharp_surfels;
        const float over_selection_threshold = 1.0 - selection_threshold;
        for (auto& [voxel_id, voxel_surfels] : surface_graph.surface)
            for (const Surfel& surfel : voxel_surfels)
            {
                if (surfel.sharpness < over_selection_threshold)
                {
                    sharp_surfels.insert(&surfel);
                }
                else
                {
                    *deleted_surfel_output++ = &surfel;
                }
            }
        return Surface_process::make_sub_graph_from_surfels(surface_graph, sharp_surfels.cbegin(), sharp_surfels.cend());
    }
};

#endif // SURFEL_SELECTOR_H