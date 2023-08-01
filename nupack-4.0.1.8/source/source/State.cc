/** \file State.cc
 * @brief Defines any non-templated non-inline functions for State classes
 */
#include <nupack/state/State.h>
#include <nupack/model/Model.h>

namespace nupack {

/******************************************************************************************/

template <>
real64 Model<real64>::structure_energy(Complex const &x, PairList p, bool d) const {return ::nupack::structure_energy(x, std::move(p), *this, d);}

template <>
real32 Model<real32>::structure_energy(Complex const &x, PairList p, bool d) const {return ::nupack::structure_energy(x, std::move(p), *this, d);}

/******************************************************************************************/

StateBase::StateBase(System y, PairList p) : sys(std::move(y)), pairs(std::move(p)) {
    if (!sys.data) NUPACK_ERROR("Empty System pointer");
    auto const n = len(sys.sequence());
    if (pairs.empty()) pairs = PairList(n);
    else if (len(pairs) == sys.n_bases())
        pairs = pairs.with_null_bases(vec<iseq>{
            indirect_view(sys.strands(), [](auto const &s) {return len(s) - 2u;})});
    else if (len(pairs) != n)
        NUPACK_ERROR("pair list length doesn't match sequence length", len(pairs), n);
    // for_pseudoknots(pairs, [](auto i, auto j, auto k, auto l) { // doesn't work if strands out of order.
    //     NUPACK_ERROR("structure contains at least one pseudoknot", i, j, k, l);
    // });
}

/******************************************************************************************/

template <class X, class V, class F>
bool for_each_aligned_pairs(std::size_t s, X &data, V const &self, V const &other, F &&f) {
    if (s == len(self)) return f();
    for (auto i : indices<iseq>(other)) // for each possible matching index
        if (self[i] == other[s] && !contains(key_view(data), i)) {
            data.push_back({i, iseq(sum(view(self, 0, i), len)), iseq(len(other[s]))});
            if (!for_each_aligned_pairs(s+1, data, self, other, f)) return false;
            data.pop_back();
        }
    return true;
}

vec<PairList> StateBase::all_aligned_pairs(StateBase const &w, bool only_one) const {
    vec<PairList> out;
    // Make a list of the {index in *this, position in *this, length of strand} for w's strands
    auto data = reserved<vec<std::array<iseq, 3>>>(len(sys.strands()));
    for_each_aligned_pairs(0, data, sys.strands(), w.sys.strands(), [&] {
        // Make a base to base map of index in w to index in *this
        auto map = reserved<vec<iseq>>(len(sys));
        for (auto const &t : data)
            for (auto i : range(t[1], t[1] + t[2])) map.push_back(i);

        // Copy the data into the right places in a new pair vector
        auto pairs = w.pairs;
        iseq n = 0;
        for (auto const &t : data)
            for (auto i : range(t[1], t[1] + t[2])) pairs[i] = map[w.pairs[n++]];
        
        out.emplace_back(pairs);
        return only_one;
    });
    return out;
}

PairList StateBase::aligned_pairs(StateBase const &w) const {
    auto v = all_aligned_pairs(w, true);
    NUPACK_REQUIRE(len(v), >, 0, "no matches found", *this, w);
    return std::move(v.at(0));
}

/******************************************************************************************/

std::size_t StateBase::symmetry() const {
    std::size_t const sym1 = rotational_symmetry(sys.strands());
    if (sym1 == 1) return 1;

    auto n = sys.data->nicks;
    std::adjacent_difference(n.begin(), n.end(), n.begin());
    std::size_t const sym2 = rotational_symmetry(n);

    if (sym2 == 1) return 1;

    std::size_t const sym3 = pairs.symmetry();
    if (sym3 == 1) return 1;

    return std::gcd(std::gcd(sym1, sym2), sym3);
}

/******************************************************************************************/

string StateBase::dp() const {
    string ret;
    for (auto const &x : complexes.complex_indices) {
        for (auto i : indices<iseq>(x)) {
            auto s = x[i];
            for (auto p = sys.data->nicks[s] + 1; p != sys.data->nicks[s + 1] - 1; ++p) {
                auto q = pairs[p];
                if (q == p) ret.push_back('.');
                else if (q < sys.data->nicks[s] || q >= sys.data->nicks[s + 1]) {
                    auto t = find_if(x, [&](iseq z) {return q < sys.data->nicks[z + 1] && q >= sys.data->nicks[z];}) - begin_of(x);
                    ret.push_back(t < i ? ')' : '(');
                }
                else ret.push_back(q < p ? ')' : '(');
            }
            ret.push_back('+');
        }
        ret.back() = ' ';
    }
    NUPACK_REQUIRE(count(ret, ')'), ==, count(ret, '('), ret);
    ret.pop_back(); return ret;
}

/******************************************************************************************/

/// Returns a vector of (index with null bases) -> (index without null bases)
pair_data_type make_pairs_map(StateBase const &w) {
    pair_data_type map(len(w.pairs), -1);
    auto ppos = begin_of(map);
    auto mpos = 0;
    for (auto s : indices<iseq>(w.sys.strands())) {
        ++ppos;
        std::iota(ppos, ppos + len(w.sys.strands()[s]) - 2, mpos);
        ppos += len(w.sys.strands()[s]) - 2; mpos += len(w.sys.strands()[s]) - 2;
        ++ppos;
    };
    return map;
}

/******************************************************************************************/

string StateBase::sequence() const {
    std::ostringstream os;
    izip(complexes.complex_indices, [&](auto i, auto const &x) {
    	if (i) os << ' ';
    	izip(x, [&](auto i, auto s) {
    		if (i) os << '+';
        	dump_os(os, sys.data->alphabet.to_string(sys.strands()[s].offset(1, -1)));
    	});
    });
    return os.str();
}

/******************************************************************************************/

void initialize_system_data(SystemData &sys, SequenceList const &v) {
    NUPACK_ASSERT(!v.empty());

    sys.total_sequence = Sequence(sum(v, len) + 2 * n_strands(v), [&](auto p, auto n) {
        for (auto const &s : v) {
            if (front(s) != Base::null()) *(p++) = (Base::null());
            p = std::copy(s.begin(), s.end(), p);
            if (back(s) != Base::null()) *(p++) = (Base::null());
        }
    });

    sys.nicks = {0};
    for (auto const &s : v) {
        auto const n = len(s) + iseq(back(s) != Base::null()) + iseq(front(s) != Base::null());
        extend(sys.strand_map, n, len(sys.nicks) - 1);
        sys.nicks.push_back(sys.nicks.back() + n);
    }
    
    for (auto i : range(len(sys.nicks) - 1))
        sys.strands.emplace_back(sys.total_sequence.begin() + sys.nicks[i], 
                                  sys.total_sequence.begin() + sys.nicks[i+1]);
}

System::System(SequenceList const &v) {
    auto ptr = std::make_shared<SystemData>();
    initialize_system_data(*ptr, v);
    data = std::move(ptr);
}

/******************************************************************************************/

}
