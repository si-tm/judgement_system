#pragma once
#include <nupack/thermo/Compute.h>
#include <nupack/thermo/Block.h>
#include <nupack/thermo/CachedModel.h>

namespace nupack::thermo {

template <class T>
struct MemoryFactor<NoStacking, T> {
    static constexpr std::size_t strand_memory(Ignore, Ignore) {return 0;}
    static constexpr std::size_t complex_memory(Complex const &x, Ignore) {return 22 * sizeof(T) * len(x.back()) * len(x.front());}
};

template <class T>
struct MemoryFactor<AllDangles, T> {
    static constexpr std::size_t strand_memory(Sequence const &s, Alphabet const &a) {
        auto const n = a.length();
        return len(s) * sizeof(T) * (3 + 2 * n * n);
    }
    static constexpr std::size_t complex_memory(Complex const &x, Ignore) {return 20 * sizeof(T) * len(x.back()) * len(x.front());}
};

template <class T>
struct MemoryFactor<MinDangles, T> : MemoryFactor<AllDangles, T> {};

template <class T>
struct MemoryFactor<Stacking, T> {
    static constexpr std::size_t strand_memory(Sequence const &s, Alphabet const &a) {
        auto const n = a.length();
        return sizeof(T) * len(s) * (2 * cube(n) + sq(n) + n);
    }
    static constexpr std::size_t complex_memory(Complex const &x, Ignore) {return 29 * sizeof(T) * len(x.back()) * len(x.front());}
};

/******************************************************************************************/

template <class Model, class Ensemble, int Complexity>
struct PFStrategy : MemoryFactor<Ensemble, value_type_of<Model>> {
    using T = typename Model::value_type;
    using ensemble = Ensemble;
    using model = CachedModel<PF<T>, Model>;
    using block = Block<T, Ensemble>;
    using strand = typename block::StrandData;
    using forward = ScalarAlgebra<PF<T>>;
    using backward = BackwardAlgebra<PF<T>, true>;
    using cache = PF<T>;
    using complexity = size_constant<Complexity>;
    using is_partition_function = True;
    using backup = False;
};

template <class Model, class Ensemble, int Complexity>
struct MFEStrategy : MemoryFactor<Ensemble, value_type_of<Model>> {
    using T = typename Model::value_type;
    using ensemble = Ensemble;
    using model = CachedModel<MFE<T>, Model>;
    using block = Block<T, Ensemble>;
    using strand = typename block::StrandData;
    using forward = ScalarAlgebra<MFE<T>>;
    using backward = BackwardAlgebra<MFE<T>, false>;
    using cache = MFE<T>;
    using complexity = size_constant<Complexity>;
    using is_partition_function = False;
    using backup = False;
};

template <class Model1, class Model2, class Ensemble, int Complexity>
struct DynamicStrategy : MemoryFactor<Ensemble, value_type_of<Model1>> {
    using ensemble = Ensemble;
    using complexity = size_constant<Complexity>;
    using is_partition_function = True;

    using T = typename Model1::value_type;
    using block = Block<T, Ensemble>;
    using model = CachedModel<PF<T>, Model1>;
    using forward = ScalarAlgebra<PF<T>>;
    using backward = BackwardAlgebra<PF<T>, true>;
    using strand = typename block::StrandData;

    struct backup : True {
        using T = typename Model2::value_type;
        using block = Block<Big<T>, Ensemble>;
        using model = CachedModel<PF<T>, Model2>;
        using forward = BigAlgebra<PF<T>>;
        using backward = BackwardAlgebra<PF<Big<T>>, true>;
        using strand = typename block::StrandData;
    };

    struct converter {
        void operator()(value_type_of<Model1> const &t, value_type_of<Model2> &u) const {u = t;}
        void operator()(value_type_of<Model1> const &t, value_type_of<Model2> &m, Exponent<value_type_of<Model2>> &e) const {
            if (t == 0) {m = 0; e = 0;}
            else std::tie(m, e) = simd::ifrexp(t);
        }
    };
};

/******************************************************************************************/

}