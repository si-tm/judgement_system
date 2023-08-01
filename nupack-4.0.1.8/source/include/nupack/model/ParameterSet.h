#pragma once
#include "../standard/Array.h" // for multi_array
#include "../common/Constants.h" // for temperature
#include "../types/IO.h" // for goto_line_after
#include "../types/Alphabet.h"
#include "../types/Matrix.h"
#include "../standard/Ptr.h"
#include "../reflect/Serialize.h"
#include "ParameterStorage.h"


namespace nupack {

/******************************************************************************************/

struct ParameterFile : MemberOrdered {

    string path;
    NUPACK_REFLECT(ParameterFile, path);

    /* Take a parameter file name and load in the parameters
    Priority ranking:
        1) the specified path
        2) "NUPACKHOME" environment variable
        3) then the NUPACK parameters folder (set by CMake)
    */
    ParameterFile(string name="rna");

    json open() const;
};

NUPACK_DEFINE_TYPE(is_parameter_file, ParameterFile);

/******************************************************************************************/

template <class T>
struct ParameterArray : TotallyOrdered {
    using value_type = T;
    std::shared_ptr<T> array; // should be T[] but seems screwed up to compile on my clang
    std::size_t m_size = 0;

    static auto allocate(std::size_t n) {
        return n ? std::shared_ptr<T>(new T[n], [](T *t) {delete[] t;}) : std::shared_ptr<T>();
    }

    ParameterArray deep_copy() const {
        auto out = *this;
        if (array) {
            out.array = allocate(out.size());
            std::copy(begin(), end(), out.begin());
        }
        return out;
    }

    NUPACK_REFLECT(ParameterArray, array, m_size);

    std::size_t size() const {return m_size;}
    T const *begin() const {return array.get();}
    T const *end() const {return array.get() + size();}

    T *begin() {return array.get();}
    T *end() {return array.get() + size();}

    bool operator==(ParameterArray const &p) const {
        return array == p.array || std::equal(begin(), end(), p.begin(), p.end());
    }

    bool operator<(ParameterArray const &p) const {
        if (array == p.array || !array || !p.array) return array < p.array;
        return std::lexicographical_compare(begin(), end(), p.begin(), p.end());
    }

    ParameterArray() = default;
    
    explicit ParameterArray(std::size_t n) : array(allocate(n)), m_size(n) {}
    // explicit ParameterArray(std::shared_ptr<T> ptr) : array(std::move(ptr)) {}

    template <class U>
    explicit ParameterArray(ParameterArray<U> const &p) : ParameterArray(p.size()) {
        std::copy(p.begin(), p.end(), begin());
    }

    template <class I>
    auto span(I const &i, std::uint32_t n) const {return view(begin() + i.start(n), begin() + i.stop(n));}

    template <class I>
    auto span(I const &i, std::uint32_t n) {return view(begin() + i.start(n), begin() + i.stop(n));}

    auto save_repr() const {return vec<T>(begin(), end());}
    void load_repr(vec<T> const &v) {
        *this = ParameterArray(v.size());
        std::copy(v.begin(), v.end(), begin());
    }
};

/******************************************************************************************/

struct ParameterInfo {
    using is_member_ordered = True;

    ParameterFile file;
    string kind = "dG";
    real loop_bias = 0;
    real temperature = DefaultTemperature;

    NUPACK_REFLECT(ParameterInfo, temperature, loop_bias, kind, file);
};

/******************************************************************************************/

struct ParameterIndex {
    // Fixed size arrays
    auto multi_base()        const {return ParameterDescriptor<0>();}
    auto multi_init()        const {return ParameterDescriptor<1>();}
    auto multi_pair()        const {return ParameterDescriptor<2>();}
    auto log_loop_penalty()  const {return ParameterDescriptor<3>();}
    auto join_penalty()      const {return ParameterDescriptor<4>();}
    auto ninio()             const {return ParameterDescriptor<5, 5>();}
    auto interior_size()     const {return ParameterDescriptor<10, 30>();}
    auto bulge_size()        const {return ParameterDescriptor<40, 30>();}
    auto hairpin_size()      const {return ParameterDescriptor<70, 30>();} // up to 100

    std::uint32_t alphabet_length;
    bool is_condensed;

    // Variable size arrays
    ParameterGrid<4> stack; // up to 356
    ParameterGrid<3> dangle5; // 420
    ParameterGrid<3> dangle3; // 484
    ParameterGrid<2> terminal_penalty; // 500
    ParameterGrid<4> coaxial_stack; // 756
    ParameterGrid<4> interior_mismatch; // 1012
    ParameterGrid<2> interior_mismatch_1; // 1028
    ParameterGrid<4> terminal_mismatch; // 1284
    ParameterGrid<4> hairpin_mismatch; // 1540
    ParameterGrid<6> hairpin_tetra; // 5636
    ParameterGrid<5> hairpin_tri; // 6660
    // These are optional!
    ParameterGrid<6> interior_1_1; // 10756
    ParameterGrid<7> interior_1_2; // 27140
    ParameterGrid<8> interior_2_2; // 92676

    NUPACK_REFLECT(ParameterIndex, alphabet_length, is_condensed, interior_mismatch, 
        interior_mismatch_1, terminal_mismatch, stack, coaxial_stack, hairpin_tetra, 
        hairpin_tri, hairpin_mismatch, dangle5, dangle3, terminal_penalty, 
        interior_1_1, interior_1_2, interior_2_2);

    ParameterIndex() = default;
    void set_length(std::uint32_t n) noexcept;

    std::uint32_t calculate_size(std::uint32_t n) const noexcept;
};

/******************************************************************************************/

struct ParameterBase : ParameterIndex, TotallyOrdered {
    using base_type = ParameterIndex;
    
    BasePairing pairing;
    Alphabet alphabet;

    ParameterInfo info;
    std::string material;

    bool operator==(ParameterBase const &p) const {return info == p.info;}
    bool operator<(ParameterBase const &p) const {return info < p.info;}

    NUPACK_EXTEND_REFLECT(ParameterBase, ParameterIndex, info, alphabet, material, pairing);

    friend std::ostream & operator<<(std::ostream &os, ParameterBase const &p) {
        return os << "ParameterSet(" << p.info.file.path << ", " << p.info.kind << ", " << p.info.temperature << " K)";
    }
};

/******************************************************************************************/

BasePairing load_pairing(Alphabet const &, json const &j);

ParameterArray<real> load_parameter_data(Alphabet const &, ParameterIndex &, json const &j);
json save_parameter_data(Alphabet const &a, ParameterIndex const &i, ParameterArray<real> const &p);

ParameterSet<real> load_parameter_set(ParameterInfo i);

/******************************************************************************************/

template <class T>
struct ParameterSet : ParameterBase {
    using base_type = ParameterBase;
    using value_type = T;
    ParameterArray<T> data;
    ParameterArray<if_t<std::is_same_v<T, real>, float, real>> conversion_cache;
    
    ParameterSet() = default;

    template <class U>
    explicit ParameterSet(ParameterSet<U> o) : base_type(std::move(static_cast<base_type &>(o))), 
        data(std::move(o.conversion_cache)), conversion_cache(std::move(o.data)) {
        static_assert(!is_same<T, U>, "Should use normal copy constructor");
    }

    template <class F>
    void modify(F &&f) {
        data = data.deep_copy();
        conversion_cache = conversion_cache.deep_copy();
        f(*this);
        ParameterSet<if_t<std::is_same_v<T, real>, float, real>> other(*this);
        f(other);
    }

    NUPACK_EXTEND_REFLECT(ParameterSet, ParameterBase, data, conversion_cache);

    bool valid() const noexcept {return bool(data.array);}

    void load_repr(ParameterInfo const &i) {*this = ParameterSet(load_parameter_set(i));}
    ParameterInfo save_repr() const {return info;}

    template <class I, class ...Is>
    T operator()(I const &i, Is...is) const {
        auto const index = i.index(this->alphabet_length, is...);
        NUPACK_QUICK_ASSERT(valid());
        NUPACK_QUICK_REQUIRE(index, <, data.size());
        return data.begin()[index];
    }

};

/******************************************************************************************/

// Return a view of all parameter values of a given type
template <class Parameters, class Index>
auto parameter_view(Parameters &&p, Index const &i) {
    return view(p.data, i.start(p.alphabet_length), i.stop(p.alphabet_length));
} 

/******************************************************************************************/

}
