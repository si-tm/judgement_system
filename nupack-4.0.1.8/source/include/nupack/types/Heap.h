#pragma once

#include "../common/Config.h"
#include <algorithm>

namespace nupack {

/******************************************************************************************/

/// Heap that cannot have its size changed once constructed
template <class V, class Comp=less_t>
class StaticHeap {
protected:
    V c;

public:
    using container_type = V;
    using value_type = value_type_of<V>;

    StaticHeap(container_type v) : c(std::move(v)) {}

    template <class ...Ts>
    explicit StaticHeap(Ts &&...ts) : c(fw<Ts>(ts)...) {
        std::make_heap(std::begin(c), std::end(c), Comp());
    };

    auto begin() const {return std::begin(c);}
    auto end() const {return std::end(c);}

    /// Retrieve maximum element and remove it from the heap
    value_type pop() {
        std::pop_heap(std::begin(c), std::end(c), Comp());
        value_type out = std::move(back(c));
        c.pop_back();
        return out;
    }

    container_type const & contents() const {return c;}

    /// Retrieve maximum element
    value_type const & top() const {return *std::begin(c);}

    container_type sorted() && {
        std::sort_heap(std::begin(c), std::end(c), Comp());
        return std::move(c);
    }

    container_type sorted() const & {return StaticHeap(*this).sorted();}

    constexpr auto comparator() const {return Comp();}
    auto size() const {return std::size(c);}
    bool empty() const {return std::empty(c);}
};

NUPACK_DEFINE_TEMPLATE(isStaticHeap, StaticHeap, class, class);

/******************************************************************************************/

/// Heap that can have its size changed once constructed
template <class V, class Comp=less_t>
class Heap : public StaticHeap<V, Comp> {
    using StaticHeap<V, Comp>::c;

public:
    template <class ...Ts>
    void emplace(Ts &&...ts) {
        c.emplace_back(fw<Ts>(ts)...);
        std::push_heap(begin_of(c), end_of(c), Comp());
    };

    template <class B, class E>
    void insert(B b, E e) {
        for (auto it = c.insert(c.end(), b, e); it != end_of(c); ++it)
            std::push_heap(begin_of(c), std::next(it), Comp());
    }

    void clear() {c.clear();}
};

NUPACK_DEFINE_TEMPLATE(isHeap, Heap, class, class);

/******************************************************************************************/

/// Heap that can have its size changed up to a maximum size
template <class V, class Comp=less_t>
class MaxSizeHeap : public StaticHeap<V, Comp> {
    usize m_max;
    using StaticHeap<V, Comp>::c;

public:
    using base_type = StaticHeap<V, Comp>;

    template <class ...Ts>
    explicit MaxSizeHeap(usize m, Ts &&...ts) : StaticHeap<V, Comp>(fw<Ts>(ts)...), m_max(m){};

    template <class ...Ts>
    void emplace_if(Ts &&...ts) {
        if (!m_max) return;
        else if (this->size() < m_max) {
            c.emplace_back(fw<Ts>(ts)...);
            std::push_heap(begin_of(c), end_of(c), Comp());
        } else {
            typename base_type::value_type t{fw<Ts>(ts)...};
            if (Comp()(t, base_type::top())) {
                std::pop_heap(begin_of(c), end_of(c), Comp());
                c.back() = std::move(t);
                std::push_heap(begin_of(c), end_of(c), Comp());
            }
        }
    }

    auto max() const {return m_max;}

    void clear() {c.clear();}
};

NUPACK_DEFINE_TEMPLATE(isMaxSizeHeap, MaxSizeHeap, class, class);

/******************************************************************************************/

}
