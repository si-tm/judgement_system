#pragma once
#include "../math/SIMD.h"
#include "Common.h"

namespace nupack::thermo {

/******************************************************************************************/

template <class V>
decltype(auto) flat_index(V &&v, Index i) noexcept(!Debug) {
    if constexpr(Debug) {
        if (i >= std::size(v)) NUPACK_ERROR("Index too high", i, std::size(v));
        if (i < 0) NUPACK_ERROR("Index too low", i, 0);
    }
    // using T = std::remove_reference_t<decltype(v.begin()[i])>;
    // if constexpr(std::is_const_v<T> && std::is_arithmetic_v<std::decay_t<T>>)
    //     NUPACK_QUICK_ASSERT(std::isfinite(v.begin()[i]), v.begin()[i], i);
    return v.begin()[i];
}

/******************************************************************************************/

template <class V, int N>
decltype(auto) flat_index(V &&v, simd::Chunk<N> i) {return simd::load(i, v.data());}

/******************************************************************************************/

template <class V>
auto flat_index(V &&v, span const &s) noexcept(!Debug);

/******************************************************************************************/

template <class T>
struct Slice : View<T> {
    using View<T>::View;

    T data() const {return this->begin();}

    template <int N>
    decltype(auto) operator()(simd::Chunk<N> const &i) const {return flat_index(*this, i);}

    decltype(auto) operator()(Index i) const {return flat_index(*this, i);}

    auto operator()(span const &i) const {return flat_index(*this, i);}

    void fill(value_type_of<T> const &t) const {::nupack::fill(*this, t);}

    template <class U>
    void set(Index i, U &&u) const {flat_index(*this, i) = fw<U>(u);}
};

/******************************************************************************************/

template <class V>
auto flat_index(V &&v, span const &s) noexcept(!Debug) {
    if constexpr(Debug) {
        if (s.stop() > std::size(v)) NUPACK_ERROR("Span too high", s.start(), s.stop(), std::size(v));
        if (s.start() < 0) NUPACK_ERROR("Span too low", s.start(), s.stop(), 0);
    }
    using Iter = decltype(v.data() + s.start());
    auto r = Slice<Iter>(v.data() + s.start(), v.data() + s.stop());
    // for (auto &&x : r) if (std::is_const_v<std::remove_reference_t<decltype(x)>>) NUPACK_QUICK_ASSERT(std::isfinite(x), x);
    return r;
}

/******************************************************************************************/


template <class R, class T>
auto addressed_element(R const &recursion, T const &t, Index i, Index j);

template <class R, class T>
auto addressed_element(R const &recursion, T const &t, Index i, span j);

template <class R, class T>
auto addressed_element(R const &recursion, T const &t, span i, Index j);

template <class R>
struct Address {
    Index i, j;
    R recursion;
    NUPACK_REFLECT(Address, i, j, recursion);
};

template <class T, class R>
struct AddressedElement : Addressed {
    T value;
    Index i, j;
    R recursion;

    Address<R> address() const {return {i, j, recursion};}

    AddressedElement(R const &r, T t, Index i, Index j) : value(t), i(i), j(j), recursion(r) {}
};

template <class T, class R>
struct AddressedRow {
    T value;
    Index i;
    span j;
    R recursion;

    auto size() const {return value.size();}

    AddressedRow(R const &r, T t, Index i, span j) : value(t), i(i), j(j), recursion(r) {}

    auto operator()(Index k) const {return addressed_element(recursion, value(k), i, j.start() + k);}
};

template <class T, class R>
struct AddressedColumn {
    T value;
    span i;
    Index j;
    R recursion;

    auto operator()(Index k) const {return addressed_element(recursion, value(k), i.start() + k, j);}

    auto size() const {return value.size();}

    AddressedColumn(R const &r, T t, span i, Index j) : value(t), i(i), j(j), recursion(r) {}
};

/******************************************************************************************/

template <class R, class T>
auto addressed_element(R const &recursion, T const &t, Index i, Index j) {return AddressedElement<T, R>(recursion, t, i, j);}

template <class R, class T>
auto addressed_element(R const &recursion, T const &t, Index i, span j) {return AddressedRow<T, R>(recursion, t, i, j);}

template <class R, class T>
auto addressed_element(R const &recursion, T const &t, span i, Index j) {return AddressedColumn<T, R>(recursion, t, i, j);}

/******************************************************************************************/

/// Holds a full rectangle of data indexable (i, j), sliceable on the first index
struct Rectangle : MemberComparable {
    std::array<Index, 2> shape = {0, 0};
    NUPACK_REFLECT(Rectangle, shape);

    void check(Index i, Index j) const {
        NUPACK_QUICK_REQUIRE(i, >=, 0);
        NUPACK_QUICK_REQUIRE(j, >=, 0);
        NUPACK_QUICK_REQUIRE(i, <, shape[0], i, j, shape);
        NUPACK_QUICK_REQUIRE(j, <, shape[1], i, j, shape);
    }
    void check(Index i, span j) const {if (!j.empty()) {check(i, j.start()); check(i, j.stop()-1);}}
    void check(span i, Index j) const {if (!i.empty()) {check(i.start(), j); check(i.stop()-1, j);}}

    template <int N>
    constexpr auto operator()(Index i, simd::Chunk<N> j) const {return i * shape[1] + j;}

    constexpr Index operator()(Index i, Index j) const {check(i, j); return i * shape[1] + j;}
    constexpr span operator()(Index i, span j) const {check(i, j); return i * shape[1] + j;}

    std::size_t resize(std::size_t i, std::size_t j) {shape[0] = i; shape[1] = j; return i * j;}
    std::array<Index, 1> set(Index i, Index j) const {return {(*this)(i, j)};}
};

/// Indexable as (i, j) where i < j. may be sliced on either index
struct SpanBoth : Rectangle {
    constexpr Index operator()(Index i, Index j) const {check(i, j); return i * shape[1] + j;}
    constexpr span operator()(Index i, span j) const {check(i, j); return i * shape[1] + j;}
    constexpr span operator()(span i, Index j) const {check(i, j); return (shape[1] + j) * shape[0] + i;}

    template <int N>
    constexpr auto operator()(simd::Chunk<N> i, Index j) const {return j * shape[0] + i;}
    template <int N>
    constexpr auto operator()(Index i, simd::Chunk<N> j) const {return i * shape[0] + j;}

    std::size_t resize(std::size_t i, std::size_t j) {shape[0] = i; shape[1] = j; return 2 * i * j;}
    std::array<Index, 2> set(Index i, Index j) const {return {i * shape[1] + j, (shape[1] + j) * shape[0] + i};}
};

/// Indexable as (i, j) where i < j. may be sliced on the first index.
struct SpanFirst : Rectangle {
    constexpr Index operator()(Index i, Index j) const {check(i, j); return j * shape[0] + i;}
    constexpr span operator()(span i, Index j) const {check(i, j); return j * shape[0] + i;}

    std::array<Index, 1> set(Index i, Index j) const {return {(*this)(i, j)};}
};

struct SpanSecond : Rectangle {};

/******************************************************************************************/

extern std::atomic<std::size_t> AllocatedBytes;

template <class T>
T* allocate_tensor(std::size_t n) {
    AllocatedBytes.fetch_add(sizeof(T) * n, std::memory_order_relaxed);
    return simd::allocator<T>().allocate(n);
}

template <class T>
void deallocate_tensor(T *t, std::size_t n) noexcept {
    AllocatedBytes.fetch_sub(sizeof(T) * n, std::memory_order_relaxed);
    simd::allocator<T>().deallocate(t, n);
}

template <class T>
class SharedArray : public TotallyOrdered {
    std::shared_ptr<T> ptr;
    std::size_t length = 0;

public:
    static_assert(std::is_scalar_v<T>);
    using value_type = T;


    NUPACK_REFLECT(SharedArray, ptr, length);

    struct Deleter {
        std::size_t n;
        void operator()(T *t) const {deallocate_tensor(t, n);}
    };

    bool operator<(SharedArray const &a) const {
        return ptr != a.ptr && std::lexicographical_compare(begin(), end(), a.begin(), a.end());
    }

    bool operator==(SharedArray const &a) const {
        return ptr == a.ptr || std::equal(begin(), end(), a.begin(), a.end());
    }

    SharedArray() = default;
    SharedArray(SharedArray &&) noexcept = default;
    SharedArray & operator=(SharedArray &&) noexcept = default;

    SharedArray(SharedArray const &a) {
        allocate(a.size());
        std::copy(a.begin(), a.end(), begin());
    }

    SharedArray & operator=(SharedArray const &a) {
        allocate(a.size());
        std::copy(a.begin(), a.end(), begin());
        return *this;
    }

    std::shared_ptr<T const> pointer() const {return ptr;}
    std::shared_ptr<T> pointer() {return ptr;}

    bool allocate(std::size_t n) {
        if (n != length) {
            if (n) ptr.reset(allocate_tensor<T>(n), Deleter{n});
            else ptr.reset();
            length = n;

            if (!Release) {
                if constexpr(std::is_floating_point_v<T>) fill(std::numeric_limits<T>::quiet_NaN());
                if constexpr(std::is_integral_v<T>) fill(std::numeric_limits<T>::max());
            }

            return true;
        } else return false;
    }

    void fill(T const &t) noexcept {std::fill(begin(), end(), t);}

    T *data() {return ptr.get();}
    T const *data() const {return ptr.get();}

    T * begin() {return ptr.get();}
    T * end() {return begin() + size();}

    T const * begin() const {return ptr.get();}
    T const * end() const {return begin() + size();}

    auto size() const {NUPACK_QUICK_ASSERT(!length || (ptr && length), length); return length;}
    explicit operator bool() const {return bool(ptr);}

    void clear() noexcept {ptr.reset(); length = 0;}
};

/******************************************************************************************/

template <class From, class To, class SFINAE=void>
struct ConvertArray {
    template <class Converter>
    static void copy(SharedArray<From> const &m, SharedArray<To> &out, Converter const &c) {
        zip(m, out, c);
    }
};

template <class From, class To, class Converter>
void convert_array(SharedArray<From> const &m, SharedArray<To> &out, Converter const &c) {
    if constexpr(std::is_same_v<From, To>) {
        out = m;
    } else {
        out.allocate(m.size());
        ConvertArray<From, To>::copy(m, out, c);
    }
}

/******************************************************************************************/

struct MatrixBase {};

/******************************************************************************************/

template <class T, class Layout=Rectangle>
struct Matrix : MatrixBase, MemberComparable {
    using value_type = T;
    SharedArray<T> storage;
    Layout layout;

    NUPACK_REFLECT(Matrix, storage, layout);

    Matrix() = default;
    Matrix(std::size_t i, std::size_t j) {allocate(i, j);}

    auto shape() const {return layout.shape;}

    template <class U, class Converter>
    void assign_and_clear(Matrix<U, Layout> &&m, Converter const &c) {layout = m.layout; convert_array(m.storage, storage, c); m.clear();}
    
    // template <class U>
    // explicit Matrix(Matrix<U, Layout> m) {assign_and_clear(std::move(m));}

    template <class I, class J>
    decltype(auto) operator()(I const &i, J const &j) {
        return flat_index(storage, layout(i, j));
    }

    template <class I, class J>
    decltype(auto) operator()(I const &i, J const &j) const {
        if constexpr(Debug) try {return flat_index(storage, layout(i, j));} catch (...) {BEEP(*this, i, j); throw;}
        else return flat_index(storage, layout(i, j));
    }

    template <class R, class I, class J>
    auto addressed(R const &r, I const &i, J const &j) const {return addressed_element(r, (*this)(i, j), i, j);}

    bool resize(std::size_t i, std::size_t j) {return storage.allocate(layout.resize(i, j));}

    void fill(T const &t) {storage.fill(t);}

    void allocate(std::size_t i, std::size_t j) {resize(i, j);}

    auto corner(bool i, bool j) const {return (*this)(i ? layout.shape[0] - 1u : 0u, j ? layout.shape[1] - 1u : 0u);}

    friend std::ostream & operator<<(std::ostream &os, Matrix const &m) {
        if (!all_of(m.shape())) return os << "[]";
        std::ios::fmtflags f(os.flags());
        os << std::setprecision(6) << "[\n";
        for (auto i : range(m.shape()[0])) {
            if (i) os << ",\n";
            os << "    [";
            dump_os(os, std::setw(11), const_cast<Matrix&>(m)(i, 0));
            for (auto j : range(1, m.shape()[1])) dump_os(os << ", ", std::setw(11), const_cast<Matrix&>(m)(i, j));
            os << ']';
        }
        os << "\n]";
        os.flags(f);
        return os;
    }

    auto begin() {return storage.begin();}
    auto end() {return storage.end();}
    auto begin() const {return storage.begin();}
    auto end() const {return storage.end();}
    void clear() noexcept {storage.clear(); layout = Layout();}
};

/******************************************************************************************/

template <class T>
struct Column {
    SharedArray<T> storage;
    NUPACK_REFLECT(Column, storage);

    template <class I>
    decltype(auto) operator()(I const &i) const noexcept {return flat_index(storage, i);}

    template <class I>
    decltype(auto) operator()(I const &i) noexcept {return flat_index(storage, i);}

    void resize(std::size_t n) {storage.allocate(n);}

    void fill(T const &t) {std::fill(storage.begin(), storage.end(), t);}

    auto begin() {return storage.begin();}
    auto end() {return storage.end();}
    auto begin() const {return storage.begin();}
    auto end() const {return storage.end();}
    void clear() noexcept {storage.clear();}
};

/******************************************************************************************/

}

namespace nupack::io {
    template <class T> struct PrintAsContainer<thermo::Slice<T>> : PrintAsList {};
}

namespace nupack::memory {

template <class T>
struct custom<thermo::SharedArray<T>> : True {};

template <class T>
struct impl<thermo::SharedArray<T>> {
    std::size_t operator()(thermo::SharedArray<T> const &m) const {
        return sizeof(m) + sizeof(T) * m.size();
    }
    void release(thermo::SharedArray<T> &m) const {m.clear();}
};

/******************************************************************************************/

}
