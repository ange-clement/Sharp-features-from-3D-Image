#ifndef GRAPH_H
#define GRAPH_H

#include "graph_element.h"

template <typename Point>
struct Element_graph
{
    typedef Graph_element<Point> Graph_element;
    typedef Graph_pointel_element<Point> Graph_pointel_element;
    typedef Graph_linel_element<Point> Graph_linel_element;
    typedef Graph_surfel_element<Point> Graph_surfel_element;

    typedef typename std::list<std::unique_ptr<Graph_element>> Element_container;
    typedef typename Element_container::iterator Element_iterator;

    Element_container elements;
    std::size_t timestamp = 0; // stamp elements, incremented and never decremented

    Element_graph() { }
    ~Element_graph() { } // elements are automatically deleted
    Element_graph( const Element_graph& ) = delete; // non construction-copyable
    Element_graph& operator=( const Element_graph& ) = delete; // non copyable

    Graph_element* add_element(Graph_element* new_element)
    {
        new_element->neighbors.clear();
        new_element->timestamp = timestamp;
        elements.emplace_back(new_element);
        timestamp++;
        return new_element;
    }

    Graph_element* add_pointel_element()
    {
        return add_element(new Graph_pointel_element());
    }
    Graph_element* add_linel_element()
    {
        return add_element(new Graph_linel_element());
    }
    Graph_element* add_surfel_element()
    {
        return add_element(new Graph_surfel_element());
    }

    // remove an element and update the iterator
    void remove_element(typename Element_iterator& element_it)
    {
        // delete element :
        //  delete it from its neighbors's neighbor list
        //  set its neigbors as next to a hole
        Graph_element* element = element_it->get();
        for ( auto& element_neigh : element->neighbors )
        {
            element_neigh->neighbors.erase(element); // &(*element)
            element_neigh->is_removable = std::nullopt;
            element_neigh->is_next_to_hole = true;
        }
        //  then delete the surfel
        element_it = elements.erase(element_it);
    }

    void collapse_element(typename Element_iterator& element_it,
                          Graph_element* to_element)
    {
        // add all neighbors of the element to collapse to the other element
        const Graph_element* element = element_it->get();
        std::set<Graph_element*>& to_element_neighbors = to_element->neighbors;
        for ( auto& element_neigh : element->neighbors )
        {
            if (element_neigh != to_element)
            {
                to_element_neighbors.insert(element_neigh);
                element_neigh->neighbors.insert(to_element);
            }
        }
        remove_element(element_it);
    }

    Graph_element* get_first_element() const
    {
        return elements.front().get();
    }

    std::size_t size() const
    {
        return elements.size();
    }
};


#endif // GRAPH_H