#include <nupack/design/SequenceAdapter.h>
#include <nupack/reflect/Serialize.h>

namespace nupack::design {




/**
 * @brief return a view of the substrand between position begin and end, INCLUSIVE
 */
StrandView StrandView::slice(uint begin, uint end) const {
    auto length = size();
    if (begin > end) throw std::runtime_error("invalid slice. begin must be <= end.");
    if (begin >= length) throw std::out_of_range("begin is not in range.");
    if (end >= length) throw std::out_of_range("end is not in range.");

    auto lengths = prefixes<vec<uint>>(false, indirect_view(domains, len));
    auto prev = prefixes<vec<uint>>(true, indirect_view(domains, len));

    auto first = upper_bound(lengths, begin) - begin_of(lengths);
    auto last = upper_bound(lengths, end) - begin_of(lengths);
    auto b = begin - prev[first]; // begin relative to first domain
    auto e = end - prev[last] + 1; // end relative to last domain


    // single domain edge case
    if (first == last) return decltype(domains) {DomainView(domains[first].start() + b, domains[first].start() + e)};

    decltype(domains) new_domains;
    new_domains.emplace_back(domains[first].start() + b, domains[first].stop());
    for (auto i : range(first+1, last)) new_domains.emplace_back(domains[i]);
    new_domains.emplace_back(domains[last].start(), domains[last].start() + e);
    return new_domains;
}


/**
 * @brief add any missing domain complements. this function is idempotent once
 *     non-complement domains are constant.
 */
// void DesignSequence::add_domain_complements() {
    // auto cur_domains = domain_specs;
    // for (auto const &dom : cur_domains) {
    //     if (back(dom.name) == '*') continue;

    //     auto complement_name = dom.name + "*";
    //     if (none_of(cur_domains, [&](auto const &c) {return complement_name == c.name;})) {
    //         // can eventually add funcitonality to restrict domains for
    //         // possible performance improvement, but it shouldn't matter for
    //         // correctness of applying constraints.
    //         add_domain({complement_name, alphabet.domain(string(len(dom), 'N'))});
    //     }
    // }
// }


/**
 * @brief convert current set of domain and strand specs into underlying
 *     sequence and views on the sequence for the domains and strands.
 */
void DesignSequence::make_sequence(vec<NamedDomain> const &domain_specs, vec<TargetStrand> const &strand_specs) {
    domains.clear(); strands.clear();
    // constraints = ConstraintHandler();

    auto total_domain = join<vec<Wildcard>>(domain_specs);

    nucleotides = indirect_view(total_domain, 
        [](Wildcard b) {NUPACK_ASSERT(+b, "empty wildcard index"); return b.first();});

    int last = 0;
    Alphabet a;
    for (auto const & d: domain_specs) {
        if (a.data) {NUPACK_REQUIRE(d.alphabet, ==, a);}
        else {a = d.alphabet;}
        auto beg = last; last += len(d); auto end = last;
        domains.emplace(d.name, DomainView(beg, end));
    }

    for (auto const & s: strand_specs) {
        strands.emplace(s.name, vmap<vec<DomainView>>(s.domains, 
            [&](NamedDomain const &d) {return get_domain(d.name);}));
    }

    // add variables to ConstraintHandler
    // for_each(nucleotides, [&](auto const &n) {constraints.add_nucleotide_variable(n);});
    NUPACK_ASSERT(a.data);
    constraints = Constraints(a, total_domain);
    add_complementarity_constraints();
}


/**
 * @brief initialize nucleotides by initializing constraints and
 *     converting variables to nucleotides
 * @details nucleotides are resized and blanked out as the addition of any
 *     WordConstraints increases the number of variables in the
 *     ConstraintHandler. This resize only needs to be done once at
 *     initialization assuming no further constraints are added after this
 *     point.
 */
void DesignSequence::initialize_sequence() {
    // auto vars = constraints.init_random();
    // nucleotides = Sequence(len(vars), 'N');
    // transform(vars, nucleotides, [](auto i) {return i < 16 ? Base::from_index(i) : Base::null();});
    auto result = constraints.initial_sequence();
    if (result) {
        nucleotides = result.value();
    } else {
        NUPACK_ERROR("unable to find sequence satisfying all constraints");
    }
}


/**
 * @brief Add complementarity constraints between domain x and x* for all domains x
 */
void DesignSequence::add_complementarity_constraints() {
    for (auto const &dom : domains) {
        auto name = dom.first;
        if (back(name) != '*') {
            auto comp_name = name + "*";
            auto it = domains.find(comp_name);
            if (it != domains.end())
                zip(dom.second.indices, ~(it->second.indices), [&](auto i, auto j) {
                    // constraints.add_constraint(CompConstraint(i, j, NUPACK_CS_STRONG));
                    constraints.complementarity_constraint(alphabet, base_pairing, i, j, wobble_mutations);
                });
        }
    }
}

/**
 * @brief set a given domain to the sequence
 *
 * @param name the domain name
 * @param in the sequence to change the domain to
 */
void DesignSequence::set_domain(string const &name, Sequence const &in) {
    auto domain = get_domain(name);
    NUPACK_REQUIRE(len(domain), ==, len(in), "Input sequence does not match domain length", name);
    
    // Check gecode variables to make sure it's compatible
    auto const &allowed = constraints.initial->nucs;
    izip(domain.indices, in, [&](uint index, uint i, Base b) {
        NUPACK_ASSERT(allowed[i].in(+b), "Nucleotide is not compatible with constraints", name, index);
    });

    nucleotides = Sequence(len(nucleotides), [&](Base *p, Ignore) {
        std::copy(nucleotides.begin(), nucleotides.end(), p);
        izip(domain.indices, [&](auto i, auto j) {p[j] = in[i];});
    });
    NUPACK_QUICK_REQUIRE(domain.to_sequence(nucleotides), ==, in);
}


/**
 * @brief mutate sequence at the given positions but maintain constraint satisfaction.
 *
 * @param vars positions to change
 */
bool DesignSequence::mutate_sequence(vec<uint> const &vars) {
    vec<int> muts{view(vars)};
    // vec<int> current_variables{view(nucleotides)};
    auto old_sequence = nucleotides;
    // transform(constraints.make_mutation(muts, current_variables), nucleotides, [](auto i) {
    //     /** handling WordConstraint non-nucleotide variables with values exceeding Base
    //       * representation
    //       */
    //     return i < len(Base::names) ? Base::from_index(i) : Base::null();
    // });

    if (auto result = constraints.make_mutation(nucleotides, muts)) {
        nucleotides = std::move(*result);
        return true;
    } else {
        return false;
    }
}


string DesignSequence::json_domains(Sequence s) const {
    if (len(s) == 0) s = nucleotides;

    std::map<string, Sequence> temp;
    for (auto const &domain : domains) temp.emplace(domain.first, domain.second.to_sequence(s));

    std::ostringstream ss;
    ss << json(temp);
    return ss.str();
}


}
