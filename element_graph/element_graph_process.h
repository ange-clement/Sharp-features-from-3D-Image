#ifndef ELEMENT_GRAPH_PROCESS_H
#define ELEMENT_GRAPH_PROCESS_H

#include "Element_graph.h"
#include "graph_element.h"
#include "../helper/geometry_helper.h"
#include "../surface/Surface_graph.h"
#include "../polyline/Polylines.h"
#include "../polyline/polyline_process.h"

namespace Element_graph_process
{
    namespace internal
    {
        template <typename Point>
        struct Element_graph_visitor
        {
            using Element_graph = Element_graph<Point>;
            using Graph_element = typename Element_graph::Graph_element;

            typedef std::size_t Size;
            typedef Graph_element* Vertex;

            Size degree(const Vertex& v) const { return v->neighbors.size(); }
            Size bestCapacity() const { return 8; }

            template <typename OutputIterator>
            void writeNeighbors(OutputIterator out, const Vertex& v) const
            {
                const std::set<Graph_element*>& edge_neighbors = v->neighbors;
                for (const auto& v_n : edge_neighbors)
                    *out++ = v_n;
            }

            template <typename OutputIterator, typename VertexPredicate>
            void writeNeighbors(OutputIterator out, const Vertex& v, const VertexPredicate& p) const
            {
                const std::set<Graph_element*>& edge_neighbors = v->neighbors;
                for (const auto& v_n : edge_neighbors)
                    if (p(v))
                        *out++ = v_n;
            }
        };

        template <typename Point>
        void write_graph_as_polylines(const std::string &filename, const Element_graph<Point> &graph)
        {
            using Graph_element = typename Element_graph<Point>::Graph_element;
            std::ofstream elements_stream(filename);
            for (const std::unique_ptr<Graph_element> &element : graph.elements)
            {
                const Point &coord = element->point;
                for (const auto &n_element: element->neighbors)
                {
                    const Point &n_coord = n_element->point;
                    elements_stream << "2 "
                                    << coord[0] << " " << coord[1] << " " << coord[2] << " "
                                    << n_coord[0] << " " << n_coord[1] << " " << n_coord[2] << "\n";
                }
            }
            elements_stream.close();
        }

        template <typename Point>
        void write_graph_as_point_set_with_infos(const std::string &filename, const Element_graph<Point> &graph)
        {
            using Graph_element = typename Element_graph<Point>::Graph_element;
            std::vector<Point> points;
            std::vector<double> values;
            points.reserve(graph.elements.size());
            values.reserve(graph.elements.size());
            for (const std::unique_ptr<Graph_element>& pointel : graph.elements)
            {
                points.push_back(pointel->point);
                values.push_back(pointel->value);
            }
            helper::write_point_set_values(filename, points, values, /*invert_values=*/true, /*normalize_range=*/false);
        }

        template <typename Point>
        bool is_graph_connected ( const Element_graph<Point>& graph )
        {
            using Element_graph = Element_graph<Point>;
            using Graph_element = typename Element_graph::Graph_element;
            if (graph.size() <= 1)
                return true;
            Element_graph_visitor<Point> graph_visitor;

            DGtal::BreadthFirstVisitor<Element_graph_visitor<Point>, std::set<Graph_element*> > visitor
                ( graph_visitor, graph.get_first_element() );
            while ( ! visitor.finished() )
            {
                visitor.expand();
            }
            return visitor.markedVertices().size() == graph.size();
        }

        // template <typename Point>
        // bool is_graph_connected ( const Element_graph<Point>& graph )
        // {
        //     Graph& nc_graph = const_cast<Graph&>(graph);
        //     return is_connected_print(nc_graph);
        // }

        template <typename Point>
        std::shared_ptr<Element_graph<Point>> make_sub_graph_from_elements(
                const Element_graph<Point>& graph,
                const std::set<typename Element_graph<Point>::Graph_element*>& subgraph_elements)
        {
            using Element_graph = Element_graph<Point>;
            using Graph_element = typename Element_graph::Graph_element;
            std::shared_ptr<Element_graph> subgraph(new Element_graph());
            std::map<Graph_element*, Graph_element*> graph_element_to_new_graph_element;
            // Add all elements to the new graph and create map
            for (auto& element : subgraph_elements)
            {
                Graph_element* new_element = subgraph->add_element(element->construct_new());
                new_element->point = element->point;
                new_element->value = element->value;
                graph_element_to_new_graph_element[element] = new_element;
            }

            // Add all connection of existing element using the map
            for (auto& element : subgraph_elements)
            {
                const std::set<Graph_element*>& neighbors = element->neighbors;
                Graph_element* new_element = graph_element_to_new_graph_element.at(element);
                std::set<Graph_element*>& new_element_neighbors = new_element->neighbors;
                for (auto& element_neighbor : neighbors)
                {
                    auto element_convert_iterator = graph_element_to_new_graph_element.find(element_neighbor);
                    if ( element_convert_iterator != graph_element_to_new_graph_element.end() )
                        new_element_neighbors.insert( element_convert_iterator->second );
                }
            }

            return subgraph;
        }

        // true iif :
        //    not an isthmus pointel : more than 1 neighbor or not a pointel
        //    is simple :
        //      do not create hole : is next to a hole (a neighbor element has been deleted)
        //      do not merge hole : the subgraph of neighbors is connected
        template <typename Point>
        bool is_element_removable( const Element_graph<Point>& graph, typename Element_graph<Point>::Graph_element* element )
        {
            using Element_graph = Element_graph<Point>;
            using Graph_element = typename Element_graph::Graph_element;
            if (element->is_removable.value_or(false))
            {
                // cached simple property
                return true;
            }
            // Get direct neighbors without the point
            const std::set<Graph_element*>& neighbors = element->neighbors;
            if (element->is_protected())  // with isthmus protection
            {
                element->is_removable = false;
                return false;
            }
            // Check if a hole is created
            if (!element->is_next_to_hole)
            {
                element->is_removable = false;
                return false;
            }
            // Check if holes are merged
            std::shared_ptr<Element_graph> neighbors_subgraph = make_sub_graph_from_elements(graph, neighbors);
            if (!is_graph_connected(*neighbors_subgraph))
            {
                // store result in cache
                element->is_removable = false;
                return false;
            }
            // All check passed, set chache
            element->is_removable = true;
            return true;
        }

        struct No_priority
        {
            template <typename Point>
            double operator()(const Element_graph<Point>& graph, const typename Element_graph<Point>::Graph_element& element) const
            {
                return 0;
            }
        };
        struct Invert_value_priority
        {
            template <typename Point>
            double operator()(const Element_graph<Point>& graph, const typename Element_graph<Point>::Graph_element& element) const
            {
                return 1.0 - element.value;
            }
        };
        struct Value_priority
        {
            template <typename Point>
            double operator()(const Element_graph<Point>& graph, const typename Element_graph<Point>::Graph_element& element) const
            {
                return element.value;
            }
        };

        // delete the first removable element from the graph
        // start search from "priority_begin"
        template <typename Point, typename PriorityValue = No_priority>
        bool thinning_remove_one_element( Element_graph<Point>& graph,
                                          std::map<double, std::list<typename Element_graph<Point>::Element_iterator>>& priority_map,
                                 typename std::map<double, std::list<typename Element_graph<Point>::Element_iterator>>::iterator& priority_begin,
                                          const PriorityValue& priority = PriorityValue() )
        {
            using Element_graph = Element_graph<Point>;
            using Graph_element = typename Element_graph::Graph_element;
            using Element_iterator = typename Element_graph::Element_iterator;
            using Prioriry_iterator = typename std::map<double, std::list<Element_iterator>>::iterator;

            typename std::list<Element_iterator>::iterator value_element_it;
            typename std::list<Element_iterator>::iterator value_element_end;

            for (Prioriry_iterator sorted_element_end = priority_map.end(); priority_begin != sorted_element_end; priority_begin++)
            {
                std::list<Element_iterator>& value_elements = priority_begin->second;

                for (value_element_it = value_elements.begin(), value_element_end = value_elements.end();
                    value_element_it != value_element_end; ++value_element_it)
                {
                    Element_iterator& element_it = *value_element_it;
                    Graph_element* element = element_it->get();
                    if (is_element_removable(graph, element))
                    {
                        // update the start of search (priority_begin) to the lowest neighbor
                        const double& element_value = priority(graph, *element);
                        double lowest_neighbor_value = element_value;
                        for ( auto& element_neigh : element->neighbors )
                        {
                            const double& element_value = priority(graph, *element_neigh);
                            if (element_value < lowest_neighbor_value)
                            {
                                lowest_neighbor_value = element_value;
                            }
                        }
                        // erase from graph
                        graph.remove_element(element_it);
                        // erase from priority map
                        value_elements.erase(value_element_it);
                        // if (value_elements.size() == 0)
                        //     priority_map.erase(priority_begin);
                        if (lowest_neighbor_value != element_value)
                        {
                            priority_begin = priority_map.find(lowest_neighbor_value);
                        }
                        return true;
                    }
                }
            }

            return false;
        }

        template <typename Point>
        std::vector<std::set<typename Element_graph<Point>::Graph_element*>> get_connex_groups(const Element_graph<Point>& graph)
        {
            using Element_graph = Element_graph<Point>;
            using Graph_element = typename Element_graph::Graph_element;
            using Element_graph_visitor = Element_graph_visitor<Point>;
            std::vector<std::set<Graph_element*>> connex_groups;
            std::set<Graph_element*> processed_elements;
            Element_graph_visitor graph_visitor;
            for (auto& graph_element_ptr : graph.elements)
            {
                auto& [it, inserted] = processed_elements.insert(graph_element_ptr.get());
                if (!inserted)
                    continue;
                DGtal::BreadthFirstVisitor<Element_graph_visitor, std::set<Graph_element*> > visitor
                    ( graph_visitor, graph_element_ptr.get() );
                while ( ! visitor.finished() )
                {
                    visitor.expand();
                }
                connex_groups.push_back(visitor.visitedVertices());
                processed_elements.insert(connex_groups.back().begin(), connex_groups.back().end());
            }
            return connex_groups;
        }

        template <typename Point>
        struct Graph_poisson_surfel_element_frontier_checker
        {
            using Element_graph = Element_graph<Point>;
            using Graph_element = typename Element_graph::Graph_element;
            bool operator()(const Graph_element* graph_element)
            {
                return static_cast<const Graph_poisson_surfel_element<Point>*>(graph_element)->is_frontier();
            }
        };

        // returns a vector of the elements and their topological distance to the frontier
        //   sorted by the topological distance
        template <typename Point, typename IsFrontierChecker = Graph_poisson_surfel_element_frontier_checker<Point>>
        std::vector<std::pair<typename Element_graph<Point>::Graph_element*, std::size_t>>
        sort_elements_by_topological_distance_to_frontier(
            const std::set<typename Element_graph<Point>::Graph_element*>& connex_group,
            const IsFrontierChecker& is_frontier = IsFrontierChecker())
        {
            using Element_graph = Element_graph<Point>;
            using Graph_element = typename Element_graph::Graph_element;
            using Graph_poisson_surfel_element = Graph_poisson_surfel_element<Point>;
            std::set<Graph_element*> processed_elements;
            std::vector<std::pair<Graph_element*, std::size_t>> elements_and_topological_distance;
            // init elements with topological distance 0 (on the frontier)
            for (Graph_element* element : connex_group)
            {
                Graph_poisson_surfel_element* surfel_element = static_cast<Graph_poisson_surfel_element*>(element);
                if (surfel_element->is_frontier())
                {
                    processed_elements.insert(element);
                    elements_and_topological_distance.emplace_back(element, 0);
                }
            }
            // process elements with topological_distance >= 0
            std::size_t start_index = 0;
            std::size_t current_topological_distance = 1;
            while (start_index < elements_and_topological_distance.size())
            {
                const std::size_t end_index = elements_and_topological_distance.size(); // the [start, end[ range contains elements with the previous topological distance
                for (std::size_t sorted_element_index = start_index; sorted_element_index < end_index; sorted_element_index++)
                {
                    Graph_element* element = elements_and_topological_distance[sorted_element_index].first;
                    for (Graph_element* neighbor_element : element->neighbors)
                    {
                        const bool &inserted = processed_elements.insert(neighbor_element).second;
                        if (inserted)
                            elements_and_topological_distance.emplace_back(neighbor_element, current_topological_distance);
                    }
                }
                start_index = end_index; // start of element with the current topological_distance
                current_topological_distance++;
            }
            return elements_and_topological_distance;
        }

        template <typename Vector>
        struct Poisson_surfel_element_value_extractor
        {
            static const int value_dim = 3;
            using Value = std::array<double, value_dim>;

            std::optional<Value> operator()(const Graph_poisson_surfel_element<Vector>* e) const
            {
                if (!e->is_constraint())
                    return std::nullopt;
                std::array<double, 3> element_value = {static_cast<double>(e->normal()[0]),
                                                       static_cast<double>(e->normal()[1]),
                                                       static_cast<double>(e->normal()[2])};
                return element_value;
            }
            Vector construct_value(const std::array<double, value_dim>& values) const
            {
                return Vector(values[0], values[1], values[2]);
            }
        };

    } // namespace internal

    template <typename Point>
    void write_graph_with_infos(const std::string &filename, const Element_graph<Point> &graph)
    {
        internal::write_graph_as_polylines(filename + ".polylines.txt", graph);
        internal::write_graph_as_point_set_with_infos(filename + ".ply", graph);
    }

    template <typename Point, typename Vector>
    std::shared_ptr<Element_graph<Point>>
    build_element_graph_from_surfels(
        const Surface_graph<Vector>& sharp_surface_graph,
        const std::vector<const typename Surface_graph<Vector>::Surfel*> &deleted_surfels,
        const std::map<std::array<int, 3>, double>& pointel_values)
    {
        using Surface_graph = Surface_graph<Vector>;
        using Surfel = typename Surface_graph::Surfel;
        using Int_point = std::array<int, 3>;

        using Element_graph = Element_graph<Point>;
        using Graph_element = typename Element_graph::Graph_element;

        std::shared_ptr<Element_graph> element_graph(new Element_graph());

        const std::map<std::size_t, std::list<Surfel>>& surface = sharp_surface_graph.surface;

        std::map<Int_point, Graph_element*> pointel_to_element;
        std::map<Int_point, Graph_element*> linel_to_element;
        const std::size_t& xdim = sharp_surface_graph.xdim;
        const std::size_t& ydim = sharp_surface_graph.ydim;
        Int_point s, u, v;
        // for each voxel
        for (const auto& [voxel_id, surfels] : surface)
        {
            const auto& [x, y, z] = helper::voxel_index_to_coord(voxel_id, xdim, ydim);
            // for each surfel
            for (const Surfel& surfel : surfels)
            {
                surfel.get_orientation_vectors(s, u, v); // a pointel (in local space) and two orthogonal directions
                const Int_point p = {static_cast<int>(x + s[0]), static_cast<int>(y + s[1]), static_cast<int>(z + s[2])};
                // for each elements
                /*  6 - 5 - 4
                *   | \ | / |
                *   7 - 8 - 3
                *   | / | \ |
                *   0 - 1 - 2
                */
                std::array<Point, 9> surfel_elements = {
                    Point(p[0],                   p[1],                   p[2]),
                    Point(p[0]+u[0]*0.5,          p[1]+u[1]*0.5,          p[2]+u[2]*0.5),
                    Point(p[0]+u[0],              p[1]+u[1],              p[2]+u[2]),
                    Point(p[0]+u[0]+v[0]*0.5,     p[1]+u[1]+v[1]*0.5,     p[2]+u[2]+v[2]*0.5),
                    Point(p[0]+u[0]+v[0],         p[1]+u[1]+v[1],         p[2]+u[2]+v[2]),
                    Point(p[0]+u[0]*0.5+v[0],     p[1]+u[1]*0.5+v[1],     p[2]+u[2]*0.5+v[2]),
                    Point(p[0]+v[0],              p[1]+v[1],              p[2]+v[2]),
                    Point(p[0]+v[0]*0.5,          p[1]+v[1]*0.5,          p[2]+v[2]*0.5),
                    Point(p[0]+u[0]*0.5+v[0]*0.5, p[1]+u[1]*0.5+v[1]*0.5, p[2]+u[2]*0.5+v[2]*0.5)
                };
                // get (or add) pointels
                std::array<Graph_element*, 9> elements;
                for (int i = 0; i < 8; i+=2)
                {
                    const Point& coord = surfel_elements[i];
                    Int_point key = {static_cast<int>(coord[0]), static_cast<int>(coord[1]), static_cast<int>(coord[2])};
                    typename std::map<Int_point, Graph_element*>::iterator pointel_it
                        = pointel_to_element.find(key);
                    if (pointel_it == pointel_to_element.end())
                    {
                        auto it_and_inserted = pointel_to_element.emplace(key, element_graph->add_pointel_element());
                        Graph_element* new_element = it_and_inserted.first->second;
                        new_element->point = coord;
                        new_element->value = pointel_values.at(key);
                        new_element->is_next_to_hole = false;
                        pointel_it = it_and_inserted.first;
                    }
                    elements[i] = pointel_it->second;
                }
                // get (or add) linels
                for (int i = 1; i < 8; i+=2)
                {
                    const Point& coord = surfel_elements[i];
                    Int_point key = {static_cast<int>(coord[0]*2), static_cast<int>(coord[1]*2), static_cast<int>(coord[2]*2)};
                    typename std::map<Int_point, Graph_element*>::iterator linel_it
                        = linel_to_element.find(key);
                    if (linel_it == linel_to_element.end())
                    {
                        auto it_and_inserted = linel_to_element.emplace(key, element_graph->add_linel_element());
                        Graph_element* new_element = it_and_inserted.first->second;
                        new_element->point = coord;
                        // linel element value is the average of the two adjacent pointel elements' values (ignore surfel element)
                        new_element->value = 0.5*(elements[i-1]->value + elements[(i+1)&7]->value);
                        new_element->is_next_to_hole = false;
                        linel_it = it_and_inserted.first;
                    }
                    elements[i] = linel_it->second;
                }
                // add surfel
                Graph_element* new_element = element_graph->add_surfel_element();
                new_element->point = surfel_elements[8];
                // surfel element value is the average of the four adjacent pointel elements' values (ignore linel element)
                new_element->value = 0.25*
                                                (elements[0]->value + elements[2]->value +
                                                elements[4]->value + elements[6]->value);
                new_element->is_next_to_hole = false;
                elements[8] = new_element;
                // make the connections
                Graph_element* element_surfel = elements[8];
                auto& element_surfel_neigh = element_surfel->neighbors;
                for (int i = 0; i < 8; i++)
                {
                    Graph_element* element_0 = elements[i];
                    Graph_element* element_1 = elements[(i+1)&7];
                    // linel
                    element_0->neighbors.insert(element_1);
                    element_1->neighbors.insert(element_0);
                    // surfel
                    element_0->neighbors.insert(element_surfel);
                    element_surfel->neighbors.insert(element_0);
                }
            }
        }
        // update "next_to_hole" property of elements for all element adjacent to deleted_surfels
        for (const Surfel* surfel : deleted_surfels)
        {
            const auto& [x, y, z] = helper::voxel_index_to_coord(surfel->voxel_index, xdim, ydim);
            surfel->get_orientation_vectors(s, u, v);
            const Int_point p = {static_cast<int>(x + s[0]), static_cast<int>(y + s[1]), static_cast<int>(z + s[2])};
            // for each elements
            /*  6 - 5 - 4
            *   | \ | / |
            *   7 - X - 3
            *   | / | \ |
            *   0 - 1 - 2
            */
            std::array<Int_point, 8> surfel_elements = {
                 p,
                {p[0]*2+u[0],          p[1]*2+u[1],          p[2]*2+u[2]},
                {p[0]  +u[0],          p[1]  +u[1],          p[2]  +u[2]},
                {(p[0] +u[0])*2+v[0], (p[1]  +u[1])*2+v[1], (p[2]  +u[2])*2+v[2]},
                {p[0]  +u[0]   +v[0],  p[1]  +u[1]   +v[1],  p[2]  +u[2]   +v[2]},
                {(p[0] +v[0])*2+u[0], (p[1]  +v[1])*2+u[1], (p[2]  +v[2])*2+u[2]}, // {(p[0]+v[0])*2+v[0], (p[1]+v[1])*2+v[1], (p[2]+v[2])*2+v[2]},
                {p[0]  +v[0],          p[1]  +v[1],          p[2]  +v[2]},
                {p[0]*2+v[0],          p[1]*2+v[1],          p[2]*2+v[2]}
            };
            // update pointels if present
            for (int i = 0; i < 8; i+=2)
            {
                const Int_point& coord = surfel_elements[i];
                typename std::map<Int_point, Graph_element*>::iterator pointel_it
                    = pointel_to_element.find(coord);
                if (pointel_it == pointel_to_element.end())
                    continue;
                pointel_it->second->is_next_to_hole = true;
            }
            // update linels if present
            for (int i = 1; i < 8; i+=2)
            {
                const Int_point& coord = surfel_elements[i];
                typename std::map<Int_point, Graph_element*>::iterator linel_it
                    = linel_to_element.find(coord);
                if (linel_it == linel_to_element.end())
                    continue;
                linel_it->second->is_next_to_hole = true;
            }
        }

        return element_graph;
    }

    // adds only the poisson surfel elements
    template <typename Vector>
    std::shared_ptr<Element_graph<Vector>>
    build_poisson_element_graph_from_surfels(
        const Surface_graph<Vector>& surface_graph,
        const std::set<const typename Surface_graph<Vector>::Surfel*> &surfels_to_update,
        const std::unordered_set<const typename Surface_graph<Vector>::Surfel*> &frontier_surfels,
        std::map<const typename Surface_graph<Vector>::Surfel*, Graph_poisson_surfel_element<Vector>*> &surfel_to_element,
        std::map<const typename Element_graph<Vector>::Graph_element*, const typename Surface_graph<Vector>::Surfel*> &element_to_surfel)
    {
        using Surface_graph = Surface_graph<Vector>;
        using Surfel = typename Surface_graph::Surfel;

        using Element_graph = Element_graph<Vector>;
        using Graph_element = typename Element_graph::Graph_element;
        using Graph_poisson_surfel_element = Graph_poisson_surfel_element<Vector>;

        std::shared_ptr<Element_graph> element_graph(new Element_graph());

        const std::map<std::size_t, std::list<Surfel>>& surface = surface_graph.surface;

        surfel_to_element.clear();
        element_to_surfel.clear();

        // build nodes and values
        for (const auto& [voxel_id, voxel_surfels] : surface_graph.surface)
        {
            for (const Surfel& surfel : voxel_surfels)
            {
                Graph_poisson_surfel_element* poisson_surfel_element = static_cast<Graph_poisson_surfel_element*>(element_graph->add_element(new Graph_poisson_surfel_element()));
                poisson_surfel_element->normal() = surfel.normal;
                poisson_surfel_element->is_constraint() = (surfels_to_update.count(&surfel) == 0);
                poisson_surfel_element->is_frontier() = (frontier_surfels.count(&surfel) > 0);
                surfel_to_element.emplace(&surfel, poisson_surfel_element);
                element_to_surfel.emplace(poisson_surfel_element, &surfel);
            }
        }
        // build topology
        for (auto& [surfel, element] : surfel_to_element)
        {
            Graph_poisson_surfel_element* surfel_element = static_cast<Graph_poisson_surfel_element*>(element);
            if (surfel_element->is_frontier())
            {
                // Frontier case : connect to the non frontier neighbor with the most similar normal. If all the neighbor are frontiers, ignore the frontier property.
                const Vector& element_normal = surfel_element->normal();
                Graph_poisson_surfel_element* closest_normal_neighbor = nullptr;
                Graph_poisson_surfel_element* closest_frontier_normal_neighbor = nullptr;
                double closest_normal_dot = -2;
                double closest_frontier_normal_dot = -2;
                for (const Surfel* neighbor_surfel : surfel->neighbors)
                {
                    Graph_poisson_surfel_element* neighbor_element = surfel_to_element.at(neighbor_surfel);
                    const Vector& neighbor_normal = neighbor_element->normal();
                    const double& dot = CGAL::scalar_product(element_normal, neighbor_normal);
                    if (neighbor_element->is_frontier())
                    {
                        if (dot > closest_frontier_normal_dot)
                        {
                            closest_frontier_normal_dot = dot;
                            closest_frontier_normal_neighbor = neighbor_element;
                        }
                    }
                    else
                    {
                        if (dot > closest_normal_dot)
                        {
                            closest_normal_dot = dot;
                            closest_normal_neighbor = neighbor_element;
                        }
                    }
                }
                if (closest_normal_neighbor == nullptr)
                    closest_normal_neighbor = closest_frontier_normal_neighbor;
                if (closest_normal_neighbor != nullptr)
                {
                    surfel_element->neighbors.insert(closest_normal_neighbor);
                    closest_normal_neighbor->neighbors.insert(surfel_element);
                }
            }
            else
            {
                // Default case : connect to all non-frontier element
                for (const Surfel* neighbor_surfel : surfel->neighbors)
                {
                    Graph_poisson_surfel_element* neighbor_element = surfel_to_element.at(neighbor_surfel);
                    if (neighbor_element->is_frontier())
                        continue;
                    surfel_element->neighbors.insert(neighbor_element);
                    neighbor_element->neighbors.insert(surfel_element);
                }
            }
        }

        // handle cases where no nodes are contraints : set the furthest element from frontier as constraint
        const std::vector<std::set<Graph_element*>> connex_groups = internal::get_connex_groups(*element_graph);
        for (const std::set<Graph_element*>& group : connex_groups)
        {
            if (group.size() <= 1) // filter out isolated element
                continue;
            bool has_constraint_element = false;
            for (const Graph_element* element : group)
            {
                const Graph_poisson_surfel_element* surfel_element = static_cast<const Graph_poisson_surfel_element*>(element);
                if (surfel_element->is_constraint())
                {
                    has_constraint_element = true;
                    break;
                }
            }
            if (has_constraint_element)
                continue;
            // there are no constraint elements: set element with maximum topological distance as constraint
            std::vector<std::pair<Graph_element*, std::size_t>> elements_and_topological_distance = internal::sort_elements_by_topological_distance_to_frontier<Vector>(group);
            Graph_poisson_surfel_element* max_surfel_element;
            max_surfel_element = static_cast<Graph_poisson_surfel_element*>(elements_and_topological_distance.back().first);
            max_surfel_element->is_constraint() = true;
        }

        return element_graph;
    }

    // ElementValueExtractor allows to retrieve a value from an element. It contains
    //      int value_dim :
    //          The dimension of the extracted value.
    //      std::optional<std::array<double, value_dim>> operator()(const Graph_poisson_surfel_element* e) :
    //          Returns the value of the constraint element,
    //          or an empty optional if the element is not a constraint
    //      Point construct_point(const std::array<double, value_dim>& values) :
    //          construct a Point object from an array of double
    // This method return the values of non constrained elements
    template <typename Point, typename ElementValueExtractor = internal::Poisson_surfel_element_value_extractor<Point>>
    std::vector<std::optional<Point>>
    propagate_values_via_average_iterations(
        const Element_graph<Point>& graph, const std::size_t& maximum_number_of_iteration, const double& minimum_value_change,
        const ElementValueExtractor& extract_value = ElementValueExtractor())
    {
        using Element_graph = Element_graph<Point>;
        using Graph_element = typename Element_graph::Graph_element;
        using Graph_poisson_surfel_element = Graph_poisson_surfel_element<Point>;
        const int& value_dim = ElementValueExtractor::value_dim;
        using Value = std::array<double, value_dim>;
        // Construct graph of elements with value to compute and initialize values
        std::vector<std::optional<Value>> element_values;
        element_values.resize(graph.elements.size());
        std::vector<const Graph_poisson_surfel_element*> elements_to_fill;
        for (const auto& element_ptr : graph.elements)
        {
            const Graph_poisson_surfel_element* element = static_cast<const Graph_poisson_surfel_element*>(element_ptr.get());
            const std::optional<Value>& element_value = extract_value(element);
            const std::size_t &element_index = element->get_index();
            element_values[element_index] = element_value;
            if (!element_value.has_value())
                elements_to_fill.push_back(element);
        }

        // solve : each element with a value to fill is set as the average of its neighbors values
        // stops when all elements have a value and
        //      either a certain number of loop is reached or the biggest value change is under min_value_change
        const std::size_t &number_of_values_to_fill = elements_to_fill.size();
        std::size_t previous_number_of_filled_values = 0;
        std::size_t number_of_filled_values;
        for (std::size_t iteration = 0; iteration < maximum_number_of_iteration;) // iteration incremented when all values are diffused at least once.
        {
            number_of_filled_values = 0;
            double max_value_change = 0;
            for (std::size_t value_to_fill_index = 0; value_to_fill_index < number_of_values_to_fill; value_to_fill_index++)
            {
                const Graph_poisson_surfel_element* element_to_fill = elements_to_fill[value_to_fill_index];
                std::optional<Value> &value_to_fill = element_values[element_to_fill->get_index()];
                Value average_value_of_neighbors;
                average_value_of_neighbors.fill(0);
                std::size_t number_of_values_in_average = 0;
                for (const Graph_element* neighbor : element_to_fill->neighbors)
                {
                    const Graph_poisson_surfel_element* neighor_surfel_element =  static_cast<const Graph_poisson_surfel_element*>(neighbor);
                    const std::optional<Value> &neighbor_value_opt = element_values[neighor_surfel_element->get_index()];
                    if (neighbor_value_opt.has_value())
                    {
                        const Value& neighbor_value = neighbor_value_opt.value();
                        for (int d = 0; d < value_dim; d++)
                            average_value_of_neighbors[d] += neighbor_value[d];
                        number_of_values_in_average++;
                    }
                }
                if (number_of_values_in_average == 0)
                    continue;
                double inv_number_of_values_in_average = 1.0 / static_cast<double>(number_of_values_in_average);
                for (int d = 0; d < value_dim; d++)
                    average_value_of_neighbors[d] *= inv_number_of_values_in_average;
                // update max_value_change
                if (!value_to_fill.has_value())
                    max_value_change = std::numeric_limits<double>::max();
                else
                {
                    const Value& previous_value = value_to_fill.value();
                    for (int d = 0; d < value_dim; d++)
                        max_value_change = std::max(max_value_change, std::abs(previous_value[d] - average_value_of_neighbors[d]));
                }
                // update value
                value_to_fill = average_value_of_neighbors;
                number_of_filled_values++;
            }
            if (number_of_filled_values == previous_number_of_filled_values)
            {
                iteration++; // start the iteration count when all values have been diffused at least once.
            }
            previous_number_of_filled_values = number_of_filled_values;
            if (max_value_change < minimum_value_change)
                break;
        }

        // output
        std::vector<std::optional<Point>> element_points;
        const Value null_value = {0,0,0};
        element_points.resize(element_values.size());
        std::transform(element_values.begin(), element_values.end(), element_points.begin(),
                [&null_value, extract_value] (const std::optional<Value>& value_opt) -> std::optional<Point> 
                {
                    if (value_opt.has_value())
                        return extract_value.construct_value(value_opt.value());
                    else
                        return std::nullopt;
                });

        return element_points;
    }

    template <typename Point, typename PriorityValue = internal::No_priority>
    void thin_element_graph_with_sharpness(Element_graph<Point>& graph, const PriorityValue& priority = PriorityValue())
    {
        using Element_graph = Element_graph<Point>;
        using Graph_element = typename Element_graph::Graph_element;
        using Elements_it = typename Element_graph::Element_iterator;

        bool deleted_simple = false;
        int iteration_number = 0;
        Element_graph::Element_container& graph_elements = graph.elements;

        // sort surfels by mapping values to list of surfels. Also strores iterators for future deletion
        std::map<double, std::list<Elements_it>> sorted_elements;
        for (Elements_it element_it = graph_elements.begin(); element_it != graph_elements.end(); ++element_it)
        {
            sorted_elements[priority(graph, **element_it)].push_back(element_it);
        }
        typename std::map<double, std::list<Elements_it>>::iterator priority_begin = sorted_elements.begin();

#ifdef VERBOSE
        std::chrono::seconds start_ms = std::chrono::duration_cast< std::chrono::seconds >(std::chrono::system_clock::now().time_since_epoch());
        std::chrono::seconds cur_sec(0);
        std::size_t max_number_of_iterations = graph.elements.size();
#endif

        do
        {
            deleted_simple = internal::thinning_remove_one_element( graph, sorted_elements, priority_begin, priority );
            iteration_number++;

#ifdef VERBOSE
            std::chrono::seconds cur_ms = std::chrono::duration_cast< std::chrono::seconds >(std::chrono::system_clock::now().time_since_epoch());
            if (cur_ms - start_ms > cur_sec)
            {
                cur_sec = cur_ms - start_ms;
                std::cout << "\riter : " << iteration_number << " / " << max_number_of_iterations << " (" << graph.elements.size() << " left, "<< (1.0 - (double)graph.elements.size()/(double)max_number_of_iterations)*100 << "%)" << std::flush;
            }
#endif
        }
        while (deleted_simple);
#ifdef VERBOSE
        std::cout << "\r" << std::flush;
#endif
        helper::message("thining done in "+std::to_string(iteration_number)+" iterations");
    }

    template <typename Point>
    Polylines<Point> create_polyline_from_graph(const Element_graph<Point>& graph)
    {
        using Element_graph = Element_graph<Point>;
        using Graph_element = typename Element_graph::Graph_element;

        Polylines<Point> polylines;
        // fill points (same order)
        polylines.points.reserve(graph.elements.size());
        std::map<const Graph_element*, std::size_t> element_to_index;
        for (const std::unique_ptr<Graph_element> &element_ptr : graph.elements)
        {
            const Graph_element* element = element_ptr.get();
            element_to_index.emplace(element, polylines.points.size());
            polylines.points.push_back(element->point);
        }
        // convert graph to index graph
        polylines.lines.resize(graph.elements.size());
        for (const std::unique_ptr<Graph_element> &element_ptr : graph.elements)
        {
            const Graph_element* element = element_ptr.get();
            const std::size_t& element_index = element_to_index.at(element);
            std::vector<std::size_t>& line = polylines.lines[element_index];
            for (const Graph_element* neighbor_element : element->neighbors)
            {
                line.push_back(element_to_index.at(neighbor_element)); // safe only if the neighbor is in the element list of graph
            }
        }
        // order into lines
        polylines.lines = Polylines_process::make_lines_from_graph(polylines.lines);
        return polylines;
    }

    // constructed elements are Graph_pointel_element
    template <typename Point>
    std::shared_ptr<Element_graph<Point>>
    build_element_graph_from_polylines(
        const Polylines<Point>& polylines,
        std::vector<typename Element_graph<Point>::Graph_element*> &point_index_to_graph_element)
    {
        using Element_graph = Element_graph<Point>;
        using Graph_element = typename Element_graph::Graph_element;

        std::shared_ptr<Element_graph> element_graph(new Element_graph());
        point_index_to_graph_element.resize(polylines.points.size());
        // create one element per point
        for (std::size_t point_index = 0, nb_points = polylines.points.size(); point_index < nb_points; point_index++)
        {
            Graph_element* new_element = element_graph->add_pointel_element();
            new_element->point = polylines.points[point_index];
            point_index_to_graph_element[point_index] = new_element;
        }
        // connect elements
        for (const std::vector<std::size_t> &line : polylines.lines)
        {
            for (std::size_t line_element = 1, line_size = line.size(); line_element < line_size; line_element++)
            {
                const std::size_t &current_point_index = line[line_element];
                const std::size_t &previous_point_index = line[line_element-1];
                Graph_element* current_element = point_index_to_graph_element[current_point_index];
                Graph_element* previous_element = point_index_to_graph_element[previous_point_index];
                current_element->neighbors.insert(previous_element);
                previous_element->neighbors.insert(current_element);
            }
        }
        return element_graph;
    }

    //  // constructed elements are Graph_dual_line_element
    template <typename Point>
    std::shared_ptr<Element_graph<Point>>
    build_dual_graph_from_polylines(
        const Polylines<Point>& polylines,
        std::vector<typename Element_graph<Point>::Graph_element*> &line_index_to_dual_element)
    {
        using Element_graph = Element_graph<Point>;
        using Graph_element = typename Element_graph::Graph_element;

        std::shared_ptr<Element_graph> element_graph(new Element_graph());
        line_index_to_dual_element.resize(polylines.lines.size());
        std::map<std::size_t, std::vector<std::size_t>> corner_point_adjacent_lines;
        // create one element per line
        for (std::size_t line_index = 0, nb_lines = polylines.lines.size(); line_index < nb_lines; line_index++)
        {
            const std::vector<std::size_t>& line = polylines.lines[line_index];
            if (line.empty())
                continue;
            corner_point_adjacent_lines[line.front()].push_back(line_index);
            if (line.size() > 1 && line.back() != line.front()) // the line has more than one element and is not a loop
                corner_point_adjacent_lines[line.back()].push_back(line_index);
            Graph_element* new_element = element_graph->add_element(new Graph_element());
            line_index_to_dual_element[line_index] = new_element;
            // for display
            new_element->point = polylines.points[line[line.size()/2]];
            new_element->value = 0;
        }
        // connect elements
        for (std::size_t line_index = 0, nb_lines = polylines.lines.size(); line_index < nb_lines; line_index++)
        {
            const std::vector<std::size_t>& line = polylines.lines[line_index];
            if (line.empty())
                continue;
            Graph_element* dual_element = line_index_to_dual_element[line_index];
            const std::vector<std::size_t>& front_corner_neighbors = corner_point_adjacent_lines.at(line.front());
            std::set<Graph_element*> &dual_element_neihbors = dual_element->neighbors;
            for (const std::size_t& neighbor_line : front_corner_neighbors)
            {
                if (neighbor_line != line_index)
                    dual_element_neihbors.insert(line_index_to_dual_element[neighbor_line]);
            }
            if (line.size() > 1 && line.back() != line.front()) // the line has more than one element and is not a loop
            {
                const std::vector<std::size_t>& back_corner_neighbors = corner_point_adjacent_lines.at(line.back());
                for (const std::size_t& neighbor_line : back_corner_neighbors)
                {
                    if (neighbor_line != line_index)
                        dual_element_neihbors.insert(line_index_to_dual_element[neighbor_line]);
                }
            }
        }
        return element_graph;
    }

    template <typename Point>
    std::shared_ptr<Element_graph<Point>>
    build_dual_graph_from_polylines(
        const Polylines<Point>& polylines)
    {
        std::vector<Graph_dual_line_element<Point>*> line_index_to_dual_element;
        build_dual_graph_from_polylines(polylines, line_index_to_dual_element);
    }

} // namespace Element_graph_process

#endif // ELEMENT_GRAPH_PROCESS_H