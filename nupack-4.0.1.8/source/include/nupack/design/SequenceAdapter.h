/** @brief Contains classes for interfacing with ConstraintHandler and holding and modifying sequence state. */
#pragma once
// #include "../design/constraint_handler.h"
#include "Constraints.h"

#include "../iteration/Transform.h"
#include "../iteration/Range.h"

#include "../types/Matrix.h"
#include "../types/Complex.h"
#include "../types/Named.h"
#include "../reflect/SerializeMatrix.h"

#include <map>

namespace nupack::design { 

// using custom_csp::ConstraintHandler;
// using custom_csp::CompConstraint;
// using custom_csp::NUPACK_CS_STRONG;


/** @brief view of iterators into a larger Sequence definining the domain */
struct DomainView : MemberOrdered {
    span indices;
    NUPACK_REFLECT(DomainView, indices);

    DomainView() = default;

    DomainView(uint beg, uint end) : indices(beg, end) {}

    /** @brief forwarding to underlying span */
    auto start() const {return indices.start();}
    /** @brief forwarding to underlying span */
    auto stop() const {return indices.stop();}
    /** @brief forwarding to underlying span */
    auto size() const {return len(indices);}

    /**
     * @brief return the substring of the Sequence that is the domain
     *
     * @param s a Sequence containing a substring that is the domain
     * @return the domain Sequence
     */
    Sequence to_sequence(Subsequence const &s) const {
        return Sequence(subview(s, indices));
    }

    /**
     * @brief expand the span of indices into a full vector
     */
    vec<uint> to_indices() const {return vec<uint>(indices);}
};


/** @brief collection of DomainViews that can be converted into a Sequence of the strand. */
struct StrandView : MemberOrdered {
    vec<DomainView> domains;
    NUPACK_REFLECT(StrandView, domains);

    StrandView() = default;
    StrandView(vec<DomainView> domains) : domains(std::move(domains)) {}

    auto size() const {return sum(domains, len);}

    StrandView slice(uint beg, uint end) const;

    /**
     * @brief concatenate the domains of the strand from a given Sequence
     *
     * @param s a Sequence containing contiguous regions that are the domains of this strand
     * @return the strand Sequence
     */
    Sequence to_sequence(Subsequence const &s) const {
        auto seq = join(vmap(domains, [&](auto const &d) {return d.to_sequence(s);}));
        return seq;
    }

    /**
     * @brief return the concatenation of the domain indices
     */
    vec<uint> to_indices() const {
        auto vecs = vmap(domains, [](auto const &d) {return d.to_indices();});
        return join(vecs);
    }
};


/** @brief Underlying sequence which other sequence elements in design have
  views to and which forwards update requests to the ConstraintHandler. */
struct DesignSequence : MemberOrdered {
    template <class K, class V> 
    using map_type = std::unordered_map<K, V>;
    
    /** @brief underlying sequence that is mutated to match ConstraintHandler */
    Sequence nucleotides; // need a previous sequence as well to step back to...revert() or something
    /** @brief constraint handler */
    Constraints constraints;
    /** @brief alphabet of bases */
    Alphabet alphabet;
    /** @brief base pair possibilities */
    BasePairing base_pairing;
    /** @brief internal map for returning views of strands */
    map_type<string, StrandView> strands;
    /** @brief internal map for returning views of domains */
    map_type<string, DomainView> domains;
    /** @brief number of times each nucleotide was chosen for mutation */
    vec<uint> times_mutated;
    // real_mat correlation_matrix;
    uint real_variables = 0;
    bool wobble_mutations = false;

    DesignSequence() = default;
    explicit DesignSequence(Alphabet a, BasePairing p, bool wobble) : alphabet(std::move(a)), base_pairing(std::move(p)), wobble_mutations(wobble) {}

    NUPACK_REFLECT(DesignSequence, nucleotides, constraints, strands, domains, times_mutated, alphabet, base_pairing, wobble_mutations);

    /**
     * @brief get a StrandView corresponding to the name if it exists,
     *     otherwise throw an exception
     *
     * @throws std::out_of_range if name is not a strand
     * @param name name of strand to return
     * @return the corresponding view into the nucleotides representing the strand
     */
    StrandView get_strand(string const &name) const {return strands.at(name);}
    /**
     * @brief same as get_strand(), but for domains
     */
    DomainView get_domain(string const &name) const {return domains.at(name);}


    void set_domain(string const &name, Sequence const &in);

    /* constraint-related methods */
    void add_complementarity_constraints();


    string json_domains(Sequence s={}) const;

    void make_sequence(vec<NamedDomain> const &domains, vec<TargetStrand> const &strands);

    void initialize_sequence();

    bool mutate_sequence(vec<uint> const &vars);

    /**
     * @brief validation and setting of nucleotides to new Sequence
     *
     * @throws std::invalid_argument if incoming sequence is not of matching length
     * @param s incoming sequence
     */
    void set_sequence(Sequence s) {
        if (len(s) == len(nucleotides)) nucleotides = std::move(s);
        else throw std::invalid_argument("incoming sequence is incorrect length");
    }

    /**
     * @brief print names and sequences for each domain and the same for each sequence
     */
    void print_components() {
        for (auto const &d: domains) print(d.first, d.second.to_sequence(nucleotides));
        for (auto const &d: strands) print(d.first, d.second.to_sequence(nucleotides));
    }
};

/** @brief extract sequence elements for strands in vector and return implied ::nupack::Complex */
inline ::nupack::Complex to_nick_sequence(vec<StrandView> const &strands, Sequence const &s) {
    return ::nupack::Complex(indirect_view(strands, [&](auto const &d) {return d.to_sequence(s);}));
}

}
