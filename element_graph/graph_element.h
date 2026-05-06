#ifndef GRAPH_ELEMENT_H
#define GRAPH_ELEMENT_H

#include <CGAL/assertions.h> // CGAL_unreachable

template <typename Point>
struct Graph_element
{
    Graph_element() { }

    std::size_t timestamp = static_cast<std::size_t>(-1);

    bool is_next_to_hole = false;
    std::optional<bool> is_removable = std::nullopt;
    std::set<Graph_element*> neighbors;

    Point point;
    double value;

    virtual bool is_protected() const
    {
        CGAL_unreachable();
        return false;
    };

    virtual Graph_element* construct_new() const { return new Graph_element<Point>(); }
};

template <typename Point>
struct Graph_pointel_element : public Graph_element<Point>
{
    using Base = Graph_element<Point>;
    Graph_pointel_element() : Base() { }

    virtual bool is_protected() const override
    {
        return Base::neighbors.size() == 1;
    }

    virtual Base* construct_new() const override { return new Graph_pointel_element<Point>(); }
};

template <typename Point>
struct Graph_linel_element : public Graph_element<Point>
{
    using Base = Graph_element<Point>;
    Graph_linel_element() : Base() { }

    virtual bool is_protected() const override
    {
        return false;
    }

    virtual Base* construct_new() const override { return new Graph_linel_element<Point>(); }
};

template <typename Point>
struct Graph_surfel_element : public Graph_element<Point>
{
    using Base = Graph_element<Point>;
    Graph_surfel_element() : Base() { }

    virtual bool is_protected() const override
    {
        return false;
    }

    virtual Base* construct_new() const override { return new Graph_surfel_element<Point>(); }
};

template <typename Vector>
struct Graph_poisson_surfel_element : public Graph_element<Vector>
{
    using Base = Graph_element<Vector>;

    Graph_poisson_surfel_element() : Base() { }

    virtual bool is_protected() const override
    {
        CGAL_unreachable();
        return false;
    }

    virtual Base* construct_new() const override
    {
        CGAL_unreachable();
        return new Graph_poisson_surfel_element<Vector>();
    }

    const bool& is_frontier() const {
        return Base::is_next_to_hole;
    }
    bool& is_frontier() {
        return Base::is_next_to_hole;
    }

    const bool& is_constraint() const {
        return is_constraint_bool;
    }
    bool& is_constraint() {
        return is_constraint_bool;
    }

    const Vector& normal() const {
        return Base::point;
    }
    Vector& normal() {
        return Base::point;
    }

    // no Graph_poisson_surfel_element is removed
    const std::size_t& get_index() const {
        return Base::timestamp;
    }

private:
    bool is_constraint_bool;
};

#endif // GRAPH_ELEMENT_H