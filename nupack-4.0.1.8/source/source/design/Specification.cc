#include <nupack/design/Specification.h>
#include <nupack/design/OutputResult.h>
#include <nupack/types/Named.h>
// #include <spdlog/spdlog.h>
// #include <spdlog/sinks/basic_file_sink.h>

namespace nupack::design {

// using custom_csp::CompConstraint;
// using custom_csp::CompConstraint;
// using custom_csp::IdentConstraint;
// using custom_csp::PatternConstraint;
// using custom_csp::WordConstraint;
// using custom_csp::MatchConstraint;
// using custom_csp::NUPACK_CS_STRONG;
// using custom_csp::trinary;

void CompareConstraint::check_constraint() const {
    NUPACK_REQUIRE(sum(left, len), ==, sum(right, len), *this);
}

void PatternConstraint::check_constraint() const {
    // NUPACK_REQUIRE(sum(domains, len), ==, len(pattern), *this);
}

void DiversityConstraint::check_constraint() const {
    NUPACK_REQUIRE(word_length, >, 0, *this);
    NUPACK_REQUIRE(min_nucleotide_types, >, 0, *this);
}

void SimilarityConstraint::check_constraint() const {
    NUPACK_ASSERT(range.first >= 0 && range.second <= 1 && range.first <= range.second, *this);
    NUPACK_REQUIRE(sum(domains, len), ==, len(reference), *this);
}

void LibraryConstraint::check_constraint() const {
    std::size_t const size = sum(domains, len);
    std::size_t cum = 0;
    for (auto const &d : comparisons) {
        NUPACK_ASSERT(!d.empty());
        cum += d[0].size();
        for (auto const &c : d) NUPACK_REQUIRE(c.size(), ==, d[0].size());
    }
    NUPACK_REQUIRE(cum, ==, size, "Words in the library constraint do not agree with total length of constrained domains", *this);
}

void WindowConstraint::check_constraint() const {
    std::size_t const size = sum(domains, len);
    NUPACK_REQUIRE(len(comparisons), ==, 1, *this);
    for (auto const &c: comparisons[0]) 
        NUPACK_REQUIRE(size, <=, c.size(), *this);
}

template <class T>
T as_objective(T const &t) {return t;}

PatternObjective as_objective(PatternConstraint const &c) {
    return {vmap(c.domains, [](auto const &d) {return d.name;}), c.patterns};
}

SimilarityObjective as_objective(SimilarityConstraint const &c) {
    return {c.domains, c.reference, c.range};
}

void Specification::replace(DomainMap const &map) {
    for (auto &t : tubes) t.replace(map);
    for (auto &t : hard_constraints) fork(t, [&](auto &t) {t.replace(map);});
    for (auto &t : soft_constraints) fork(t, [&](auto &t) {t.replace(map);});
}


Designer Specification::create_designer() const {
    DesignSequence seqs(model.alphabet(), model.pairing(), options.wobble_mutations);

    std::unordered_map<string, NamedDomain> domains;
    std::unordered_map<string, TargetStrand> strands;
    std::map<vec<TargetStrand>, std::pair<TargetComplex, std::size_t>> complexes;
    std::unordered_map<string, vec<TargetStrand>> complex_check;

    auto put_domain = [&](NamedDomain const &d) {
        if (auto [it, p] = domains.try_emplace(d.name, d); !p)
            NUPACK_REQUIRE(it->second, ==, d, "Different domains with the same name");
    };

    auto put_complex = [&](TargetComplex const &x) {
        auto const v = lowest_rotation(x.strands);
        if (auto [it, p] = complex_check.try_emplace(x.name, v); !p)
            NUPACK_REQUIRE(it->second, ==, v, "Different complexes with the same name");

        if (auto [it, p] = complexes.try_emplace(v, x, 0); !p)
            it->second.first.resolve(x);

        for (auto const &s : x.strands) {
            if (auto [it, p] = strands.try_emplace(s.name, s); !p)
                NUPACK_REQUIRE(it->second, ==, s, "Different strands with the same name");

            for (auto const &d : s.domains) put_domain(d);
        }
    };

    for (auto const &t : tubes) for (auto const &x : t.complexes) put_complex(x);

    for (auto const &c : hard_constraints) fork(c, [&](auto const &c) {c.check_constraint(); c.register_definitions(put_domain, put_complex);});
    for (auto const &c : soft_constraints) fork(c, [&](auto const &c) {c.check_constraint(); c.register_definitions(put_domain, put_complex);});

    seqs.make_sequence(item_view(domains), item_view(strands));

    for (auto const &c : hard_constraints) fork(c, [&](auto const &c) {c.implement_hard(seqs, model);});

    /** Design level operations */
    Design design(std::move(seqs));
    
    DecompositionParameters params {options.H_split, options.N_split, options.f_split, options.f_sparse, options.dG_clamp};

    // Set the indices of the complexes
    izip(complexes, [](auto i, auto &p) {p.second.second = i;});

    for (auto const &[k, v] : complexes) {
        design.add_complex(indirect_view(v.first.strands, [](auto const &s) {return s.name;}), 
            model, v.first.name, v.first.structure, params);
    }

    design.add_structure_complementarity();

    /* add tubes */
    for (auto const &t : tubes) {
        vec<uint> indices;
        for (auto const &x : t.complexes) 
            indices.push_back(complexes.at(complex_check.at(x.name)).second);
        vec<real> concs(indices.size(), 0.0);
        std::copy(t.concentrations.begin(), t.concentrations.end(), concs.begin());
        design.add_tube(indices, std::move(concs), t.name);
    };
    
    Weights weighting;
    for (auto &w : defect_weights) weighting.add(std::move(w));

    vec<Objective> objectives;
    objectives.emplace_back(MultitubeObjective());
    weighting.add_objective_weight(objective_weight);

    for (auto const &c : soft_constraints)
        fork(c, [&](auto const &c) {
            objectives.emplace_back(as_objective(c));
            weighting.add_objective_weight(c.weight);    
        });

    Designer ret(std::move(design), std::move(objectives), std::move(weighting), options);
    return ret;
}



// Specification::operator Designer() const {
//     DesignSequence seqs(model.alphabet());

//     /** Sequence-level operations */
//     /* add domains */
//     for_each(domains, [&](auto const &x) {seqs.add_domain(x);});
//     /* add strands */
//     for_each(strands, [&](auto const &x) {seqs.add_strand(x);});

//     seqs.make_sequence();
//     /** Design level operations */
//     Design design(std::move(seqs), std::make_shared<thermo::Storage>(parameters.cache_bytes_of_RAM));
//     /* add complexes */
//     auto comp_name = [](auto const &c) {
//         if (c.name != "") return c.name;
//         std::stringstream s; s << at(c.strands, 0);
//         for (auto i : range(1, len(c.strands))) s << "-" << c.strands[i];
//         return s.str();
//     };

//     DecompositionParameters params {parameters.H_split, parameters.N_split, parameters.f_split, parameters.f_sparse, parameters.dG_clamp};

//     for_each(complexes, [&](auto const &x) {design.add_complex(x.strands, model, comp_name(x), x.structure, params);});


/**
 * @brief returns a pair of the variables for the left and right concatenated sequences
 *
 * @param seqs the DesignSequence to pull the variable mapping from
 * @return a pair with the left indices first and the right indices second
 */
std::pair<vec<int>, vec<int>> CompareConstraint::get_variables(DesignSequence const &seqs) const {
    return {extract_variables(left, seqs), extract_variables(right, seqs)};
}

void MatchConstraint::implement_hard(DesignSequence &seqs, Ignore) const {
    auto vars = get_variables(seqs);
    zip(vars.first, vars.second, [&](int i, int j) {
        // seqs.constraints.add_constraint(IdentConstraint(i, j));
        seqs.constraints.match_constraint(i, j);
    });
}

void ComplementarityConstraint::implement_hard(DesignSequence &seqs, Model<> const &model) const {
    auto vars = get_variables(seqs);
    zip(vars.first, reversed(vars.second), [&](int i, int j) {
        // seqs.constraints.add_constraint(CompConstraint(i, j, NUPACK_CS_STRONG));
        seqs.constraints.complementarity_constraint(model.alphabet(), model.pairing(), i, j, wobble_mutations);
    });
}

void PairingConstraint::implement_hard(DesignSequence &seqs, Model<> const &model) const {
    auto vars = get_variables(seqs);
    zip(vars.first, reversed(vars.second), [&](int i, int j) {
        // seqs.constraints.add_constraint(CompConstraint(i, j, NUPACK_CS_STRONG));
        seqs.constraints.pairing_constraint(model.pairing(), i, j);
    });
}

/**
 * @brief converts the specification into a PatternConstraint and adds it to the
 *     list of constraints in seqs
 * @details if name is the empty string, pattern is prevented in every strand
 *     in seqs. Otherwise, name is looked up as either a strand or a domain
 *     and the pattern is prevented in that one element.
 *
 * @param seqs the DesignSequence to which to add the constraint implied by
 *     the specification
 */
void PatternConstraint::implement_hard(DesignSequence &seqs, Ignore) const {
    /* handle global case */
    // auto poss = seqs.constraints.get_possible_nucleotides();
    for (auto const &pattern : patterns) {
        if (domains.empty()) {
            auto vars = vmap(seqs.strands, [&](auto const &s) {
                return vmap(s.second.to_indices(), [](auto i) {return int(i);});
            });
            for_each(vars, [&](auto const &v) {
                // seqs.constraints.add_constraint(PatternConstraint(v, pattern, poss));
                seqs.constraints.pattern_constraint(v, pattern);
            });
        } else {
            // seqs.constraints.add_constraint(PatternConstraint(extract_element(name, seqs), pattern, poss));
            seqs.constraints.pattern_constraint(extract_variables(domains, seqs), pattern);
        }
    }
}


void DiversityConstraint::implement_hard(DesignSequence &seqs, Ignore) const {
    if (domains.empty()) {
        auto vars = vmap(seqs.strands, [&](auto const &s) {
            return vmap(s.second.to_indices(), [](auto i) {return int(i);});
        });
        for_each(vars, [&](auto const &v) {
            seqs.constraints.diversity_constraint(v, word_length, min_nucleotide_types);
        });
    } else {
        seqs.constraints.diversity_constraint(extract_variables(domains, seqs), word_length, min_nucleotide_types);
    }
}

/**
 * @brief converts the specification into a WordConstraint and adds it to the
 *     list of constraints in seqs
 * @details Adding a WordConstraint must be preceded by adding a variable to
 *     the ConstraintHandler to represent which of the enumerated words is
 *     still accessible during propagation. Both adding the variable and then
 *     adding the constraint are handled in this function.
 *
 * @param seqs the DesignSequence to which to add the constraint implied by
 *     the specification
 */
void LibraryConstraint::implement_hard(DesignSequence &seqs, Ignore) const {
    auto vars = extract_variables(domains, seqs);
    int i = 0;
    for_each(comparisons, [&](DomainList const &c) {
        auto length = std::min<int>(int(vars.size()) - i, len(c.at(0)));
        // auto supp_var = seqs.constraints.add_variable(vec<trinary>(len(c), true));
        auto cur = vec<int>(subview(vars, span{i, i+length}));

        // seqs.constraints.add_constraint(WordConstraint(cur, c, supp_var));
        seqs.constraints.word_constraint(cur, c);
        i += length;
    });
}


/**
 * @brief converts the specification into a MatchConstraint and adds it to the
 *     list of constraints in seqs
 *
 * @param seqs the DesignSequence to which to add the constraint implied by
 *     the specification
 */
void SimilarityConstraint::implement_hard(DesignSequence &seqs, Ignore) const {
    // seqs.constraints.add_constraint(MatchConstraint(extract_element(name, seqs), reference, {range.first}, {range.second}));
    seqs.constraints.similarity_constraint(extract_variables(domains, seqs), reference, range);
}


/**
 * @brief return a concatenation of the indices associated with the elements
 *     in names if they are valid elements of seqs
 *
 * @param names the names of strands or domains
 * @param seqs the DesignSequence from which to pull the elements
 *
 * @return the concatenation of the individual indices from each of the
 *     elements named in names
 */
vec<int> extract_variables(vec<NamedDomain> const &domains, DesignSequence const &seqs) {
    auto temp = vmap(domains, [&](auto const &domain) {
        return extract_element(domain, seqs);
    });
    return join(temp);
}


/**
 * @brief return the variable indices corresponding to the named element
 * @details attempt to find either a domain or strand with the given name,
 *     convert the DomainSpec or StrandSpec into its vector of indices and
 *     return. If no element is found with the given name, an exception is
 *     thrown.
 *
 * @param name the name of a strand or domain
 * @param seqs the DesignSequence from which to pull the DomainSpec or
 *     StrandSpec
 * @throws nupack::Error if name is neither a strand nor a domain in seqs
 *
 * @return a vector of the indices corresponding to the named strand or domain
 */
vec<int> extract_element(NamedDomain const &domain, DesignSequence const &seqs) {
    auto const domain_view = seqs.get_domain(domain.name);
    return vmap(domain_view.to_indices(), [](auto i) {return int(i);});
}




vec<uint> Specification::ensure_compatibility(SingleResult const &res) const {
    // existing result: map from name to index
    std::unordered_map<std::string, uint> res_comp;
    for(uint i = 0; i < res.complexes.size(); i++){
        res_comp[res.complexes[i].name] = i;
    }
    
    std::unordered_map<string, vec<TargetStrand>> complex_check;

    for (auto const &t : tubes) {
        for (auto const &x : t.complexes) {
            auto const v = lowest_rotation(x.strands);
            if (auto [it, p] = complex_check.try_emplace(x.name, v); !p)
                NUPACK_REQUIRE(it->second, ==, v, "Different complexes with the same name");
        }
    }
    return indices(complex_check);
}

}
